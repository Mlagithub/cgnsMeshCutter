#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>

#include <cgns_io.h>
#include <parmetis.h>

#include "cmdLine.h"
#include "metisCutter.h"
#include "mpiAdapter.h"
#include "stringUtil.h"

static const int DataTAG = 1000;

static std::vector<idx_t> distribute(const idx_t len, const int n, idx_t &nCellThisRank)
{
    vector<idx_t> dist(n + 1, 0);
    idx_t base = len / n;
    idx_t remainder = len % n;

    // 计算当前 rank 应该处理的数量
    int rank = MPIAdapter::rank();
    nCellThisRank = base + (rank < remainder ? 1 : 0);

    // 同步所有 rank 的数量
    MPI_Allgather(&nCellThisRank, 1, GetMPIDataType<idx_t>(), dist.data() + 1, 1, GetMPIDataType<idx_t>(), MPI_COMM_WORLD);

    // 计算累积偏移
    for (int iProc = 1; iProc < n + 1; ++iProc) { dist[iProc] += dist[iProc - 1]; }

    return dist;
}

namespace MeshCut
{
MetisCutter::MetisCutter() : bigMesh_(nullptr) {}

MetisCutter::~MetisCutter() {}

void MetisCutter::cut(int argc, char **argv)
{
    CmdLine cl{};
    cl.regist<std::string, CmdLine::MustOffer>("m", "mesh", "mesh file name.", " ");
    cl.regist<size_t, CmdLine::MustOffer>("np", "npart", "part number of sub-mesh.", " ");
    cl.parse(argc, argv);

    MPIAdapter::Initialize(&argc, &argv);
#ifdef DEBUG_MODE
    if (!MPIAdapter::isParallel() || (MPIAdapter::isParallel() && MPIAdapter::isMaster()))
    {
        std::cout << "DEBUG_MODE: type any char to continue\n";
        auto a = std::getchar();
    }
    MPIAdapter::Barrier();
#endif
    MPIAdapter::isParallel() ? this->cut_parmetis(cl.get<std::string>("mesh"), cl.get<size_t>("npart")) : this->cut_metis(cl.get<std::string>("mesh"), cl.get<size_t>("npart"));
}

// collect interface from all other threads
// 优化：使用二进制数据传输，减少字符串操作开销
void MetisCutter::collect_interface()
{
    // 打包：二进制格式 [fileId(int)][count(int)][faceStr1][faceStr2]...
    std::vector<char> sendBuf;
    sendBuf.reserve(outerFace_.size() * 256);

    for (const auto &ifile : outerFace_)
    {
        if (ifile.second.empty()) continue;

        int fileId = ifile.first;
        int count = ifile.second.size();

        size_t pos = sendBuf.size();
        sendBuf.resize(pos + sizeof(int));
        std::memcpy(sendBuf.data() + pos, &fileId, sizeof(int));

        pos = sendBuf.size();
        sendBuf.resize(pos + sizeof(int));
        std::memcpy(sendBuf.data() + pos, &count, sizeof(int));

        for (const auto &face : ifile.second)
        {
            const std::string &faceStr = face.first;
            pos = sendBuf.size();
            sendBuf.resize(pos + faceStr.size() + 1);
            std::memcpy(sendBuf.data() + pos, faceStr.c_str(), faceStr.size() + 1);
        }
    }

    MPIAdapter::Barrier();

    int sendCnt = sendBuf.size();
    vector<int> recvCnt(SIZE, 0);
    MPI_Allgather(&sendCnt, 1, MPI_INT, recvCnt.data(), 1, MPI_INT, MPI_COMM_WORLD);

    vector<int> rdisp(SIZE, 0);
    int totalRecv = recvCnt[0];
    for (int i = 1; i < SIZE; ++i)
    {
        rdisp[i] = rdisp[i - 1] + recvCnt[i - 1];
        totalRecv += recvCnt[i];
    }

    std::vector<char> recvBuf(totalRecv);
    MPI_Allgatherv(sendBuf.data(), sendCnt, MPI_CHAR, recvBuf.data(), recvCnt.data(), rdisp.data(), MPI_CHAR, MPI_COMM_WORLD);

    MPIAdapter::Barrier();

    // 解包 - 只处理不属于当前进程的 fileId
    size_t pos = 0;
    while (pos + 2 * sizeof(int) <= recvBuf.size())
    {
        int fileId;
        std::memcpy(&fileId, recvBuf.data() + pos, sizeof(int));
        pos += sizeof(int);

        int count;
        std::memcpy(&count, recvBuf.data() + pos, sizeof(int));
        pos += sizeof(int);

        // 如果这个文件属于当前进程，跳过（因为已经有了）
        if (ownerFile_.count(fileId) != 0)
        {
            for (int i = 0; i < count; ++i)
            {
                while (pos < recvBuf.size() && recvBuf[pos] != '\0') ++pos;
                ++pos;
            }
            continue;
        }

        check(true, " ", format("  Thread %d: recv %d interface of file %d\n", RANK, count, fileId));

        auto &nbrOuterface = outerFace_[fileId];
        for (int i = 0; i < count; ++i)
        {
            std::string faceStr(recvBuf.data() + pos);
            nbrOuterface.insert({faceStr, {}});
            pos += faceStr.size() + 1;
        }
    }
}

// collect nodeIdG2L_ from all other threads for interface writing
void MetisCutter::collect_nodeIdG2L(const int np)
{
    // 打包：格式 [fileId(int)][count(int)][gId1][lId1][gId2][lId2]...
    std::vector<char> sendBuf;
    sendBuf.reserve(ownerFile_.size() * 1024);

    for (auto fileId : ownerFile_)
    {
        auto &g2l = nodeIdG2L_[fileId];
        if (g2l.empty()) continue;

        size_t pos = sendBuf.size();
        sendBuf.resize(pos + sizeof(int));
        std::memcpy(sendBuf.data() + pos, &fileId, sizeof(int));

        int count = g2l.size();
        pos = sendBuf.size();
        sendBuf.resize(pos + sizeof(int));
        std::memcpy(sendBuf.data() + pos, &count, sizeof(int));

        for (const auto &kv : g2l)
        {
            cgsize_t gId = kv.first;
            cgsize_t lId = kv.second;
            pos = sendBuf.size();
            sendBuf.resize(pos + 2 * sizeof(cgsize_t));
            std::memcpy(sendBuf.data() + pos, &gId, sizeof(cgsize_t));
            std::memcpy(sendBuf.data() + pos + sizeof(cgsize_t), &lId, sizeof(cgsize_t));
        }
    }

    int sendCnt = sendBuf.size();
    vector<int> recvCnt(SIZE, 0);
    MPI_Allgather(&sendCnt, 1, MPI_INT, recvCnt.data(), 1, MPI_INT, MPI_COMM_WORLD);

    vector<int> rdisp(SIZE, 0);
    int totalRecv = recvCnt[0];
    for (int i = 1; i < SIZE; ++i)
    {
        rdisp[i] = rdisp[i - 1] + recvCnt[i - 1];
        totalRecv += recvCnt[i];
    }

    std::vector<char> recvBuf(totalRecv);
    MPI_Allgatherv(sendBuf.data(), sendCnt, MPI_CHAR, recvBuf.data(), recvCnt.data(), rdisp.data(), MPI_CHAR, MPI_COMM_WORLD);

    // 解包 - 跳过当前进程已拥有的
    size_t pos = 0;
    while (pos + 2 * sizeof(int) <= recvBuf.size())
    {
        int fileId;
        std::memcpy(&fileId, recvBuf.data() + pos, sizeof(int));
        pos += sizeof(int);

        int count;
        std::memcpy(&count, recvBuf.data() + pos, sizeof(int));
        pos += sizeof(int);

        // 如果这个文件当前进程已拥有，跳过
        if (ownerFile_.count(fileId) != 0)
        {
            pos += count * 2 * sizeof(cgsize_t);
            continue;
        }

        auto &g2l = nodeIdG2L_[fileId];
        for (int i = 0; i < count; ++i)
        {
            cgsize_t gId, lId;
            std::memcpy(&gId, recvBuf.data() + pos, sizeof(cgsize_t));
            std::memcpy(&lId, recvBuf.data() + pos + sizeof(cgsize_t), sizeof(cgsize_t));
            pos += 2 * sizeof(cgsize_t);
            g2l.insert({gId, lId});
        }
    }
}

// collect maxElemId_ from all other threads for interface writing
void MetisCutter::collect_maxElemId(const int np)
{
    // 打包：格式 [fileId(int)][maxId(cgsize_t)]...
    std::vector<char> sendBuf;
    sendBuf.reserve(ownerFile_.size() * (sizeof(int) + sizeof(cgsize_t)));

    for (auto fileId : ownerFile_)
    {
        if (maxElemId_.count(fileId) == 0) continue;
        cgsize_t maxId = maxElemId_[fileId];

        size_t pos = sendBuf.size();
        sendBuf.resize(pos + sizeof(int));
        std::memcpy(sendBuf.data() + pos, &fileId, sizeof(int));

        pos = sendBuf.size();
        sendBuf.resize(pos + sizeof(cgsize_t));
        std::memcpy(sendBuf.data() + pos, &maxId, sizeof(cgsize_t));
    }

    int sendCnt = sendBuf.size();
    vector<int> recvCnt(SIZE, 0);
    MPI_Allgather(&sendCnt, 1, MPI_INT, recvCnt.data(), 1, MPI_INT, MPI_COMM_WORLD);

    vector<int> rdisp(SIZE, 0);
    int totalRecv = recvCnt[0];
    for (int i = 1; i < SIZE; ++i)
    {
        rdisp[i] = rdisp[i - 1] + recvCnt[i - 1];
        totalRecv += recvCnt[i];
    }

    std::vector<char> recvBuf(totalRecv);
    MPI_Allgatherv(sendBuf.data(), sendCnt, MPI_CHAR, recvBuf.data(), recvCnt.data(), rdisp.data(), MPI_CHAR, MPI_COMM_WORLD);

    // 解包
    size_t pos = 0;
    while (pos + sizeof(int) + sizeof(cgsize_t) <= recvBuf.size())
    {
        int fileId;
        std::memcpy(&fileId, recvBuf.data() + pos, sizeof(int));
        pos += sizeof(int);

        cgsize_t maxId;
        std::memcpy(&maxId, recvBuf.data() + pos, sizeof(cgsize_t));
        pos += sizeof(cgsize_t);

        if (ownerFile_.count(fileId) == 0)
        {
            maxElemId_[fileId] = maxId;
        }
    }
}

// collect cell id from all other threads
// 优化：使用连续内存和 MPI 打包发送，减少通信次数
vector<MetisCutter::Cell> MetisCutter::collect_subBody(const MetisCutter::DecomposeResult &decomposeResult, const vector<int> &ownerThread)
{
    // Step 1: 本地分组 - 按目标线程分组
    vector<vector<Cell>> sendBuffers(SIZE);
    vector<Cell> localCells;

    // 预估每个缓冲区大小
    size_t avgCellsPerThread = (decomposeResult.nCellThisRank + SIZE - 1) / SIZE;
    for (auto &buf : sendBuffers) buf.reserve(avgCellsPerThread);
    localCells.reserve(avgCellsPerThread);

    for (idx_t i = 0; i < decomposeResult.nCellThisRank; ++i)
    {
        int partId = cellPartition_[i];
        auto threadId = ownerThread[partId];
        auto cellId = decomposeResult.start + i;
        Cell c{cellId, partId};

        if (threadId == RANK)
            localCells.push_back(c);
        else
            sendBuffers[threadId].push_back(c);
    }
    vector<idx_t>{}.swap(cellPartition_); // free

    // Step 2: 交换发送/接收计数
    vector<int> sendCounts(SIZE), recvCounts(SIZE);
    for (int i = 0; i < SIZE; ++i) sendCounts[i] = sendBuffers[i].size() * sizeof(Cell);

    MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Step 3: 计算总接收大小并分配缓冲区
    int totalRecvBytes = 0;
    vector<int> recvDispls(SIZE, 0);
    for (int i = 0; i < SIZE; ++i)
    {
        recvDispls[i] = totalRecvBytes;
        totalRecvBytes += recvCounts[i];
    }

    // Step 4: 打包发送数据
    int totalSendBytes = 0;
    vector<int> sendDispls(SIZE, 0);
    for (int i = 0; i < SIZE; ++i)
    {
        sendDispls[i] = totalSendBytes;
        totalSendBytes += sendCounts[i];
    }

    vector<char> sendBuf(totalSendBytes);
    for (int i = 0; i < SIZE; ++i)
    {
        if (sendCounts[i] > 0) { std::memcpy(sendBuf.data() + sendDispls[i], sendBuffers[i].data(), sendCounts[i]); }
        vector<Cell>{}.swap(sendBuffers[i]); // free
    }

    // Step 5: Alltoallv 一次通信
    vector<char> recvBuf(totalRecvBytes);
    MPI_Alltoallv(sendBuf.data(), sendCounts.data(), sendDispls.data(), MPI_CHAR, recvBuf.data(), recvCounts.data(), recvDispls.data(), MPI_CHAR, MPI_COMM_WORLD);

    // Step 6: 解包接收数据
    for (int i = 0; i < SIZE; ++i)
    {
        if (recvCounts[i] > 0)
        {
            int numCells = recvCounts[i] / sizeof(Cell);
            Cell *cells = reinterpret_cast<Cell *>(recvBuf.data() + recvDispls[i]);
            for (int j = 0; j < numCells; ++j) localCells.push_back(cells[j]);
        }
    }

    return localCells;
}

void MetisCutter::cut_metis(string meshFilename, const int np)
{
    // open mesh file to read/write
    bigMesh_ = std::make_shared<CGFile>(meshFilename);
    this->updateOwnerThread(meshFilename, np);

    // load body-data into memory from file
    auto bigBody = bigMesh_->loadSection(bigMesh_->bodySectionIdList());

    check(true, "", format("Call METIS_V3_PartMeshKway decompose: %s\n", bigBody.name));
    // cut
    check(cut_metis(bigBody.data, (idx_t)np, cellPartition_, nodePartition_) == METIS_OK, format("Call METIS_V3_PartMeshKway decompose: %s return error %s", bigBody.name), format("Call METIS_V3_PartMeshKway decompose: %s successfuly\n", bigBody.name));

    // write sub-mesh data to file
    for (auto ifile = 0; ifile < np; ++ifile)
    {
        // open sub-mesh to write
        this->openSubMeshToWrite(meshFilename, ifile, np);

        // collect cell/node id of sub-mesh
        nodeIdG2L_.insert({ifile, {}}); // 1-base id
        nodeIds_.clear();               // 1-base id
        cellIds_.clear();               // 0-base id
        for (auto i = 0; i < cellPartition_.size(); ++i)
        {
            if (cellPartition_[i] != ifile) continue;

            cellIds_.insert(i + 1);
            for (auto it : bigBody.cell(i + 1))
                if (nodeIds_.count(it) == 0) nodeIds_.insert(it);
        }
        for (auto it : nodeIds_) nodeIdG2L_[ifile].insert({it, nodeIdG2L_[ifile].size() + 1});

        // write coordinate
        this->rwNode(ifile);

        // write body-section
        this->rwBody(ifile, bigBody);

        // write bdy-section
        this->rwBoundary(ifile);

        // write interface
        this->rwInterface(ifile);

        // clear
        for (auto i = 0; i <= ifile; ++i)
        {
            if (outerFace_[i].empty())
            {
                nodeIdG2L_[i].clear();
                smallMesh_[i]->close();
            }
        }
    }

    // clear
    bigMesh_->close();
    bigBody.clear();
}

void MetisCutter::cut_parmetis(string meshFilename, const int np)
{
    RANK = MPIAdapter::rank(), SIZE = MPIAdapter::size();

    if (SIZE > np)
    {
        if (MPIAdapter::isMaster()) check(true, "", format("MPI size should not bigger than subMesh number: %d > %d\n", SIZE, np));
        MPIAdapter::Finalize();
        exit(EXIT_SUCCESS);
    }

    // open mesh file to read/write
    bigMesh_ = std::make_shared<CGFile>(meshFilename);
    this->updateOwnerThread(meshFilename, np);

    // handle bodysection
    auto bigBody = bigMesh_->bodySection();

    // ===== 新并行算法：避免数据交换 =====
    // Step 1: 并行分区
    auto decomposeResult = this->decompose_body(bigBody, np);

    // Step 2: 本地构建 subBody，按 partId 分组
    // 每个进程只处理自己负责的分区
    auto comPartId = [](const Cell &lhs, const Cell &rhs) { return lhs.partId < rhs.partId; };
    auto comId = [](const Cell &lhs, const Cell &rhs) { return lhs.id < rhs.id; };

    // 本地分组
    unordered_map<int, vector<Cell>> localCellsByPart;
    for (idx_t i = 0; i < decomposeResult.nCellThisRank; ++i)
    {
        int partId = cellPartition_[i];
        auto cellId = decomposeResult.start + i;
        localCellsByPart[partId].push_back(Cell{cellId, partId});
    }
    vector<idx_t>{}.swap(cellPartition_); // free

    // Step 3: 对于不属于当前进程的分区，发送给目标进程
    // 对于属于当前进程的分区，从其他进程接收
    vector<Cell> myCells;

    // 3.1 统计每个进程需要发送/接收的数据量
    vector<int> sendCounts(SIZE, 0), recvCounts(SIZE, 0);
    for (const auto &kv : localCellsByPart)
    {
        int partId = kv.first;
        int targetRank = ownerThread_[partId];
        if (targetRank != RANK) { sendCounts[targetRank] += kv.second.size() * sizeof(Cell); }
        else
        {
            myCells.insert(myCells.end(), kv.second.begin(), kv.second.end());
        }
    }
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // 3.2 打包发送数据
    vector<int> sendDispls(SIZE, 0), recvDispls(SIZE, 0);
    int totalSend = 0, totalRecv = 0;
    for (int i = 0; i < SIZE; ++i)
    {
        sendDispls[i] = totalSend;
        totalSend += sendCounts[i];
        recvDispls[i] = totalRecv;
        totalRecv += recvCounts[i];
    }

    vector<char> sendBuf(totalSend);
    for (const auto &kv : localCellsByPart)
    {
        int partId = kv.first;
        int targetRank = ownerThread_[partId];
        if (targetRank != RANK)
        {
            char *ptr = sendBuf.data() + sendDispls[targetRank];
            std::memcpy(ptr, kv.second.data(), kv.second.size() * sizeof(Cell));
            sendDispls[targetRank] += kv.second.size() * sizeof(Cell);
        }
    }

    // 重置位移
    for (int i = 0; i < SIZE; ++i) { sendDispls[i] = (i == 0) ? 0 : sendDispls[i - 1] + sendCounts[i - 1]; }

    // 3.3 交换数据
    vector<char> recvBuf(totalRecv);
    MPI_Alltoallv(sendBuf.data(), sendCounts.data(), sendDispls.data(), MPI_CHAR, recvBuf.data(), recvCounts.data(), recvDispls.data(), MPI_CHAR, MPI_COMM_WORLD);

    // 3.4 解包接收数据
    for (int i = 0; i < SIZE; ++i)
    {
        if (recvCounts[i] > 0)
        {
            int numCells = recvCounts[i] / sizeof(Cell);
            Cell *cells = reinterpret_cast<Cell *>(recvBuf.data() + recvDispls[i]);
            for (int j = 0; j < numCells; ++j) myCells.push_back(cells[j]);
        }
    }

    // 清理
    localCellsByPart.clear();
    sendBuf.clear();
    recvBuf.clear();

    // Step 4: 处理本地数据
    if (myCells.empty()) { std::cout << format("Thread %d has no subBody\n", RANK); }

    std::sort(myCells.begin(), myCells.end(), comPartId);
    auto pos = myCells.begin();
    while (pos != myCells.end())
    {
        auto lastPos = pos;
        const auto ifile = lastPos->partId;
        pos = std::upper_bound(lastPos, myCells.end(), *lastPos, comPartId);
        std::sort(lastPos, pos, comId);

        auto indexLB = lastPos->id, indexUB = (pos - 1)->id;
        auto curBody = bigMesh_->loadSection(bigMesh_->bodySectionIdList(), indexLB, indexUB);

        // open sub-mesh file for writing
        this->openSubMeshToWrite(meshFilename, ifile, np);

        // cell and node
        cellIds_.clear(), nodeIds_.clear();
        for (auto it = lastPos; it != pos; ++it)
        {
            cellIds_.insert(it->id);
            for (auto node : curBody.cell(it->id)) nodeIds_.insert(node);
        }

        nodeIdG2L_.insert({ifile, {}});
        if (!nodeIds_.empty())
            for (auto it : nodeIds_) nodeIdG2L_[ifile].insert({it, nodeIdG2L_[ifile].size() + 1});

        // write coordinate
        this->rwNode(ifile);

        // write body
        this->rwBody(ifile, curBody);

        // write bdy
        this->rwBoundary(ifile);
    }
    myCells.clear();
    bigBody.clear();

    // Step 5: 收集并写入界面
    this->collect_interface();

    // 收集所有分区的 nodeIdG2L_，用于写入 interface 到不属于自己的分区
    this->collect_nodeIdG2L(np);

    // 收集所有分区的 maxElemId_
    this->collect_maxElemId(np);

    // 先关闭所有文件，然后再以修改模式打开写入 interface
    bigMesh_->close();
    for (auto it : smallMesh_)
        if (it != nullptr)
        {
            it->close();
            it.reset();
        }

    // 确保所有进程都完成关闭后再写入界面
    MPIAdapter::Barrier();

    // write interface
    for (auto ifile : ownerFile_) this->writeInterface(ifile, np);
}

int MetisCutter::cut_metis(const vector<vector<idx_t>> &cellToplogy, idx_t np, vector<idx_t> &cellPartition, vector<idx_t> &nodePartition)
{
    // metis options
    idx_t options[METIS_NOPTIONS];
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING] = 0;
    // options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
    // options[METIS_OPTION_IPTYPE] = METIS_IPTYPE_EDGE;
    // options[METIS_OPTION_RTYPE] = METIS_RTYPE_GREEDY;
    // options[METIS_OPTION_NCUTS] = 1;
    // options[METIS_OPTION_NITER] = 10;
    // options[METIS_OPTION_UFACTOR] = 30;
    METIS_SetDefaultOptions(options);

    // eind-eptr
    std::vector<idx_t> eind, eptr(cellToplogy.size() + 1, 0);
    cgsize_t offset = 0;
    for (auto i = 0; i < cellToplogy.size(); ++i)
    {
        for (auto jt : cellToplogy[i]) { eind.push_back(jt - 1); }
        offset += cellToplogy[i].size();
        eptr[i + 1] = offset;
    }
    idx_t nCell = cellToplogy.size();
    idx_t nNode = bigMesh_->nNode();
    cellPartition.assign(nCell, 0);
    nodePartition.assign(nNode, 0);
    idx_t ncommon = 4;
    idx_t objval;

    auto rst = METIS_PartMeshDual(&nCell, &nNode, eptr.data(), eind.data(), NULL, NULL, &ncommon, &np, NULL, options, &objval, cellPartition.data(), nodePartition.data());

    nodePartition.clear();
    vector<idx_t>{}.swap(nodePartition);

    return rst;
}

int MetisCutter::cut_parmetis(idx_t np, vector<idx_t> &elmdist, vector<idx_t> &eptr, vector<idx_t> &eind, vector<idx_t> &cellPartition)
{
    // parallel read element connectivity data
    idx_t wgtflag = 0;      // 0: no weights(elmwgt is null), 2: weights on the vertices only
    idx_t numflag = 0;      // 0: C-style numbering, 1: Fortran-style numbering
    idx_t ncon = 1;         // the number of weights that each vertex has.
    idx_t ncommonnodes = 4; // 对于3D四面体/六面体网格，4更合适，减少对偶图边数
    idx_t edgecut = 0;
    idx_t opts[3] = {0, 1, static_cast<idx_t>(std::time(nullptr))};
    real_t tpwgts[ncon * np];
    real_t ubvec[ncon];

    for (auto i = 0; i < ncon * np; ++i) { tpwgts[i] = 1.0 / (real_t)np; }
    for (auto i = 0; i < ncon; ++i) { ubvec[i] = 1.05; }

    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    return ParMETIS_V3_PartMeshKway(elmdist.data(), eptr.data(), eind.data(),
                                    NULL, // weights of the elements
                                    &wgtflag, &numflag, &ncon, &ncommonnodes, &np, tpwgts, ubvec, opts, &edgecut, cellPartition.data(), &mpi_comm);
}

MetisCutter::DecomposeResult MetisCutter::decompose_body(const CGFile::Section &bigBody, const int np)
{
    DecomposeResult result;

    auto elmdist = distribute(bigBody.end - bigBody.start + 1, SIZE, result.nCellThisRank);
    result.start = bigBody.start + elmdist[RANK];

    auto shareBigbody = bigMesh_->loadSection(bigMesh_->bodySectionIdList(), result.start, result.start + result.nCellThisRank - 1);

    // 预计算节点总数并预分配内存，避免多次重分配
    size_t totalNodes = 0;
    idx_t beg = shareBigbody.isMixed() ? 1 : 0;
    for (auto id = result.start; id < result.start + result.nCellThisRank; ++id) { totalNodes += shareBigbody.cell(id).size() - beg; }

    std::vector<idx_t> eind, eptr;
    eind.reserve(totalNodes);
    eptr.reserve(result.nCellThisRank + 1);
    eptr.push_back(0);

    for (auto id = result.start; id < result.start + result.nCellThisRank; ++id)
    {
        auto &tmp = shareBigbody.cell(id);
        for (int j = beg; j < tmp.size(); ++j) { eind.push_back(tmp[j] - 1); }
        eptr.push_back(eptr.back() + tmp.size() - beg);
    }
    check(true, "", format("ParMETIS_V3_PartMeshKway begin divide: %s\n", shareBigbody.name));
    cellPartition_.assign(result.nCellThisRank, 0);
    check(cut_parmetis(np, elmdist, eptr, eind, cellPartition_) == METIS_OK, format("Process %d: ParMETIS_V3_PartMeshKway return error while divide section: %s\n", RANK, bigBody.name), format("ParMETIS_V3_PartMeshKway done\n"));

    return result;
}

void MetisCutter::updateOwnerThread(string bigFileName, const int np)
{
    idx_t dummy;

    // calculate partId collection of each thread
    ownerThread_.assign(np, 0);
    {
        auto tmp = distribute(np, SIZE, dummy);
        for (auto i = 1; i <= SIZE; ++i)
        {
            for (auto j = tmp[i - 1]; j < tmp[i]; ++j) ownerThread_[j] = i - 1;
        }
    }

    smallMesh_.assign(np, nullptr);
    // for(auto i=0; i<np; ++i)
    // {
    //     if(ownerThread[i] == RANK)
    //     {
    //         smallMesh_[i] = std::make_shared<CGFile>(smallMeshName(bigFileName, i, np), CG_MODE_WRITE);
    //         ownerFile_.insert(i);
    //     }
    // }
}

void MetisCutter::openSubMeshToWrite(string filename, const int id, const int np)
{
    if (ownerThread_[id] == RANK)
    {
        smallMesh_[id] = std::make_shared<CGFile>(smallMeshName(filename, id, np), CG_MODE_WRITE);
        ownerFile_.insert(id);
    }
}

void MetisCutter::rwBody(const int ifile, const CGFile::Section &bigBody)
{
    auto &subFile = smallMesh_[ifile];

    // new a sectionk
    auto &curS = subFile->addSection();
    strcpy(curS.name, bigBody.name);
    curS.cellType = bigBody.cellType;
    for (auto id : cellIds_) { curS.addCell(id, bigBody.cell(id), bigBody.flag(id)); }

    // write data
    subFile->writeSection(curS, nodeIdG2L_[ifile]);

    // 更新最大元素 ID
    maxElemId_[ifile] = curS.end;

    // global info
    auto len = curS.end - curS.start + 1;
    subFile->writeGlobalInfo(curS, bigBody.data.size(), globalOffset_, globalOffset_ + len);
    globalOffset_ += len;

    // update outerface
    this->updateOuterFace(curS, ifile);

    // clear
    curS.clear();
}

void MetisCutter::rwBoundary(const int ifile)
{
    auto &subFile = smallMesh_[ifile];
    auto &curOuterFace = outerFace_[ifile];

    for (auto ithSection : bigMesh_->bdySectionIdList())
    {
        auto &subBdy = subFile->addSection();
        auto &bigBdy = bigMesh_->loadSection(ithSection);
        strcpy(subBdy.name, bigBdy.name);
        subBdy.cellType = bigBdy.cellType;

        cgsize_t ID = bigBdy.start;
        for (auto i = bigBdy.start; i <= bigBdy.end; ++i)
        {
            auto face = bigBdy.cell(i);
            auto faceStr = CGFile::stringAFace(face);
            if (curOuterFace.find(faceStr) != curOuterFace.end())
            {
                curOuterFace.erase(faceStr);
                subBdy.addCell(ID++, face, bigBdy.flag(i));
            }
        }

        if (!subBdy.data.empty())
        {
            subFile->writeSection(subBdy, nodeIdG2L_[ifile]);
            maxElemId_[ifile] = subBdy.end;
        }

        // clear
        subBdy.clear();
    }
}

void MetisCutter::rwInterface(const int id)
{
    auto &subFile = smallMesh_[id];
    auto &curOuterFace = outerFace_[id];

    for (auto nbr = 0; nbr < id; ++nbr)
    {
        if (outerFace_[nbr].empty()) continue;

        auto &nbrOuterFace = outerFace_[nbr];

        set<cgsize_t> typeFlags;
        CGFile::Section curS, nbrS;
        for (auto face = nbrOuterFace.begin(); face != nbrOuterFace.end();)
        {
            auto it = curOuterFace.find(face->first);
            if (it != curOuterFace.end())
            {
                curS.data.emplace_back(it->second);
                nbrS.data.push_back(face->second);

                auto t = CGFile::CellType(it->second.size(), 2);
                curS.typeFlag.push_back(t);
                nbrS.typeFlag.push_back(t);

                typeFlags.insert(it->second.size());

                it = curOuterFace.erase(it);
                face = nbrOuterFace.erase(face);
            }
            else
            {
                ++face;
            }
        }
        if (!curS.data.empty())
        {
            strcpy(curS.name, format("Interface_%d", nbr).c_str());
            strcpy(nbrS.name, format("Interface_%d", id).c_str());

            if (typeFlags.size() > 1)
            {
                curS.cellType = ElementType_t::MIXED;
                nbrS.cellType = ElementType_t::MIXED;
                auto curOffset = 0;
                for (auto it : curS.data)
                {
                    curS.offset.push_back(curOffset);
                    nbrS.offset.push_back(curOffset);
                    curOffset += it.size() + 1;
                }
                // 补全最后一个 offset
                curS.offset.push_back(curOffset);
                nbrS.offset.push_back(curOffset);
            }
            else if (typeFlags.size() == 1)
            {
                auto type = (*typeFlags.begin() == 3) ? ElementType_t::TRI_3 : ElementType_t::QUAD_4;
                curS.cellType = type;
                nbrS.cellType = type;
            }
            subFile->writeSection(curS, nodeIdG2L_[id]);
            smallMesh_[nbr]->writeSection(nbrS, nodeIdG2L_[nbr]);
        }

        // clear
        vector<vector<cgsize_t>>{}.swap(curS.data);
        vector<cgsize_t>{}.swap(curS.offset);
        vector<cgsize_t>{}.swap(curS.typeFlag);
    }
}

void MetisCutter::rwNode(const int ifile)
{
    auto &subFile = smallMesh_[ifile];

    // zone info
    subFile->nCell(cellIds_.size());
    subFile->nNode(nodeIds_.size());

    // node data
    auto dataTmp = bigMesh_->loadCoordinate(*nodeIds_.begin(), *nodeIds_.rbegin());
    vector<vector<double>> data(3, vector<double>(nodeIds_.size(), 0.0));
    auto start = *nodeIds_.begin();
    cgsize_t i = 0;
    for (auto id : nodeIds_)
    {
        auto loc = id - start;
        data[0][i] = dataTmp[0][loc];
        data[1][i] = dataTmp[1][loc];
        data[2][i] = dataTmp[2][loc];
        ++i;
    }
    smallMesh_[ifile]->writeCoordinate(data);
}

std::string MetisCutter::smallMeshName(string bigMeshName, const int id, const int np)
{
    char smallFilename[128];
    sprintf(smallFilename, "decomposed_mesh_%01d.cgns", id);
    return std::string(smallFilename);
}

void MetisCutter::updateOuterFace(const CGFile::Section &curS, const int ifile)
{
    auto &curOuterFace = outerFace_[ifile];

    for (size_t i = 0; i < curS.data.size(); ++i)
    {
        auto cellType = curS.isMixed() ? CGFile::CellType(curS.typeFlag[i]) : curS.cellType;
        for (const auto &face : CGFile::allFaceInCell(curS.data[i], cellType))
        {
            auto val = CGFile::stringAFace(face);
            auto it = curOuterFace.find(val);
            if (it != curOuterFace.end()) { curOuterFace.erase(it); }
            else
            {
                curOuterFace.insert({std::move(val), face});
            }
        }
    }
}

void MetisCutter::writeGlobalInfo(const CGFile::Section &curS, const int id)
{
    auto len = curS.end - curS.start + 1;
    smallMesh_[id]->writeGlobalInfo(curS, 10, globalOffset_, globalOffset_ + len);
    globalOffset_ += len;
}

void MetisCutter::writeInterface(const int id, const int np)
{
    auto &curOuterFace = outerFace_[id];

    // 当前分区文件
    std::string curFilename = this->smallMeshName("", id, np);
    CGFile curFile(curFilename, CG_MODE_MODIFY);
    curFile.setIdOffset(maxElemId_[id] + 1);

    // 遍历所有其他分区
    for (auto nbr = 0; nbr < np; ++nbr)
    {
        if (nbr == id) continue;
        if (outerFace_.find(nbr) == outerFace_.end()) continue;
        if (outerFace_[nbr].empty()) continue;

        auto &nbrOuterFace = outerFace_[nbr];

        set<cgsize_t> typeFlags;
        CGFile::Section curS, nbrS;
        for (auto face = nbrOuterFace.begin(); face != nbrOuterFace.end();)
        {
            auto it = curOuterFace.find(face->first);
            if (it != curOuterFace.end())
            {
                curS.data.emplace_back(it->second);
                nbrS.data.push_back(face->second);

                auto t = CGFile::CellType(it->second.size(), 2);
                curS.typeFlag.push_back(t);
                nbrS.typeFlag.push_back(t);

                typeFlags.insert(it->second.size());

                it = curOuterFace.erase(it);
                face = nbrOuterFace.erase(face);
            }
            else
            {
                ++face;
            }
        }

        if (!curS.data.empty())
        {
            // 设置单元类型
            if (typeFlags.size() > 1)
            {
                curS.cellType = ElementType_t::MIXED;
                nbrS.cellType = ElementType_t::MIXED;

                // 设置 offset 数组
                auto curOffset = 0;
                for (auto it : curS.data)
                {
                    curS.offset.push_back(curOffset);
                    nbrS.offset.push_back(curOffset);
                    curOffset += it.size() + 1;
                }
                // 补全最后一个 offset
                curS.offset.push_back(curOffset);
                nbrS.offset.push_back(curOffset);
            }
            else if (typeFlags.size() == 1)
            {
                auto type = (*typeFlags.begin() == 3) ? ElementType_t::TRI_3 : ElementType_t::QUAD_4;
                curS.cellType = type;
                nbrS.cellType = type;
            }

            // 写入当前分区的 interface
            strcpy(curS.name, format("Interface_%d", nbr).c_str());
            curFile.writeSection(curS, nodeIdG2L_[id]);

            // 写入邻居分区的 interface (只有当前进程拥有邻居分区时才写入)
            // 注意: 并行版本中,每个进程只写入自己拥有的分区文件
            // 邻居分区的 interface 由拥有该分区的进程写入
            if (ownerFile_.count(nbr) != 0)
            {
                strcpy(nbrS.name, format("Interface_%d", id).c_str());
                std::string nbrFilename = this->smallMeshName("", nbr, np);
                CGFile nbrFile(nbrFilename, CG_MODE_MODIFY);
                nbrFile.setIdOffset(maxElemId_[nbr] + 1);
                nbrFile.writeSection(nbrS, nodeIdG2L_[nbr]);
                nbrFile.close();
            }
        }

        // clear
        vector<vector<cgsize_t>>{}.swap(curS.data);
        vector<cgsize_t>{}.swap(curS.offset);
        vector<cgsize_t>{}.swap(curS.typeFlag);
        vector<vector<cgsize_t>>{}.swap(nbrS.data);
        vector<cgsize_t>{}.swap(nbrS.offset);
        vector<cgsize_t>{}.swap(nbrS.typeFlag);
    }

    curFile.close();
}

} // namespace MeshCut
