#include <algorithm>
#include <cmath>
#include <numeric>
#include <cstring>
#include <fstream>
#include <ctime>

#include <parmetis.h>
#include <cgns_io.h>

#include "cmdLine.h"
#include "metisCutter.h"
#include "mpiAdapter.h"
#include "stringUtil.h"

static const int DataTAG = 1000;

static std::vector<idx_t> distribute(const idx_t len, const int n, idx_t &nCellThisRank)
{
    vector<idx_t> dist(n + 1, 0);
    auto step = (len + n - 1) / n;
    nCellThisRank = MPIAdapter::isLast() ? len - (n - 1) * step : step;
    MPI_Allgather(&nCellThisRank, 1, GetMPIDataType<idx_t>(), dist.data() + 1, 1, GetMPIDataType<idx_t>(), MPI_COMM_WORLD);
    for (int iProc = 1; iProc < n + 1; ++iProc) { dist[iProc] += dist[iProc - 1]; }
    return dist;
}

namespace MeshCut
{
MetisCutter::MetisCutter() : bigMesh_(nullptr)
{
}

MetisCutter::~MetisCutter()
{
}

void MetisCutter::cut(int argc, char** argv)
{
    CmdLine cl{};
    cl.regist<std::string, CmdLine::MustOffer>("m", "mesh", "mesh file name."," ");
    cl.regist<size_t, CmdLine::MustOffer>("np", "npart", "part number of sub-mesh."," ");
    cl.parse(argc, argv);

    MPIAdapter::Initialize(&argc, &argv);
#ifdef DEBUG_MODE
    if(!MPIAdapter::isParallel() || (MPIAdapter::isParallel() && MPIAdapter::isMaster())){
        std::cout << "DEBUG_MODE: type any char to continue\n"; 
        auto a = std::getchar();
    } 
    MPIAdapter::Barrier();
#endif
    MPIAdapter::isParallel() ? this->cut_parmetis(cl.get<std::string>("mesh"), cl.get<size_t>("npart")) : this->cut_metis(cl.get<std::string>("mesh"), cl.get<size_t>("npart"));
}

// collect interface from all other threads
void MetisCutter::collect_interface()
{
    auto packInterface = [&, this]() -> std::string
    {
        std::string sendStr;
        for(auto ifile : outerFace_)
        {
            if(ifile.second.empty()) continue;

            sendStr += (std::to_string(ifile.first) + ":");
            for(auto face : ifile.second)
            {
                sendStr += (face.first + ";");
            }
            sendStr += "#";
        }

        return sendStr;
    };

    auto unpackInterface = [&, this](const string& recvStr)
    {
        for(auto ifile : stringSplit(recvStr, "#"))
        {
            auto cur = stringSplit(ifile, ":");
            auto nbrFile = std::stoi(cur[0]);
            if(ownerFile_.count(nbrFile)!=0) continue;

            check(true, " ", format("  Thread %d: recv %d interface of file %d\n", RANK, stringSplit(cur[1], ";").size(), nbrFile));

            auto& nbrOuterface = outerFace_[nbrFile];
            for(auto face : stringSplit(cur[1], ";"))
            {
                nbrOuterface.insert({face, {}});
            }
        }
    };

    MPIAdapter::Barrier();

    string sendTmp = packInterface(), recvStr, sendStr;
    int sendCnt = sendTmp.size();
    vector<int> sendCnts(SIZE, sendCnt), sdisp(SIZE, 0), rdisp(SIZE, 0), recvCnt(SIZE, 0);

    MPI_Alltoall(sendCnts.data(), 1, MPI_INT, recvCnt.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for(auto i=1; i<SIZE; ++i) rdisp[i] = rdisp[i-1] + recvCnt[i-1];
    recvStr.assign(std::accumulate(recvCnt.begin(), recvCnt.end(), 0)+1, ' ');
    for(auto i=0; i<SIZE; ++i) sendStr+=sendTmp;
    for(auto i=1; i<SIZE; ++i) sdisp[i] = sdisp[i-1] + sendCnt;
    MPI_Alltoallv(sendStr.data(), sendCnts.data(), sdisp.data(), MPI_CHAR, (char*)(recvStr.data()), recvCnt.data(), rdisp.data(), MPI_CHAR, MPI_COMM_WORLD);
    
    MPIAdapter::Barrier();
    unpackInterface(stringTrim(recvStr));
}

// collect cell id from all other threads
vector<MetisCutter::Cell> MetisCutter::collect_subBody(const MetisCutter::DecomposeResult& decomposeResult, const vector<int>& ownerThread)
{
    map<int, vector<Cell>> cellInfo;
    {
        for (idx_t i = 0; i < decomposeResult.nCellThisRank; ++i)
        {
            int partId = cellPartition_[i];
            auto threadId = ownerThread[partId];
            auto cellId = decomposeResult.start + i;
            cellInfo[threadId].push_back(Cell{cellId, partId});
        }
    }
    vector<idx_t>{}.swap(cellPartition_); // free

    // all process collect send-recv info from all process
    vector<long> sendInfo(SIZE, 0), recvInfo(SIZE, 0);
    {
        for (auto it : cellInfo)
        {
            if (it.first != RANK) sendInfo[it.first] += it.second.size();
        } 
        MPI_Alltoall(sendInfo.data(), 1, MPI_LONG, recvInfo.data(), 1, MPI_LONG, MPI_COMM_WORLD);
    }

    // send-recv data
    MPI_Datatype cellType;
    {
        const int cnt = 2;
        int blockLen[cnt] = {1, 1};
        MPI_Aint disp[cnt] = {0, sizeof(idx_t)};
        MPI_Datatype types[cnt] = {MPI_LONG, MPI_INT};
        MPIAdapter::TypeStruct(cnt, blockLen, disp, types, &cellType);
        MPIAdapter::TypeCommit(&cellType);
    }
    vector<MPI_Request> requests;
    vector<MPI_Status> statuses;
    vector<vector<Cell>> recvData;
    // send
    for(int i=0; i<SIZE; ++i)
    {
        if(sendInfo[i] != 0)
        {
            requests.push_back(MPI_Request{});
            MPIAdapter::ISend(cellInfo[i].data(), cellInfo[i].size(), cellType, i, DataTAG + i, MPI_COMM_WORLD, &requests.back());
        }
    }
    // recv 
    for(int i=0; i<SIZE; ++i)
    {
        auto recvDataSize = recvInfo[i];
        if(recvDataSize != 0)
        {
            recvData.push_back(vector<Cell>(recvDataSize, Cell{}));
            requests.push_back(MPI_Request{});
            MPI_Irecv(recvData.back().data(), recvDataSize, cellType, i, RANK + DataTAG, MPI_COMM_WORLD, &(requests.back()));
        }
    }
    statuses.assign(requests.size(), MPI_Status{});
    MPIAdapter::WaitAll(requests.size(), requests.data(), statuses.data());
    requests.clear(), statuses.clear();

    // update cellInfo
    for (auto &it : recvData)
    {
        for (auto &jt : it) cellInfo[RANK].emplace_back(jt);
    }
    MPIAdapter::TypeFree(&cellType);

    return cellInfo[RANK];
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
    check(
        cut_metis(bigBody.data, (idx_t)np, cellPartition_, nodePartition_) == METIS_OK,
        format("Call METIS_V3_PartMeshKway decompose: %s return error %s", bigBody.name),
        format("Call METIS_V3_PartMeshKway decompose: %s successfuly\n", bigBody.name));

    // write sub-mesh data to file
    for(auto ifile=0; ifile<np; ++ifile)
    {
        // open sub-mesh to write 
        this->openSubMeshToWrite(meshFilename, ifile, np);

        // collect cell/node id of sub-mesh
        nodeIdG2L_.insert({ifile, {}}); // 1-base id
        nodeIds_.clear(); // 1-base id
        cellIds_.clear(); // 0-base id
        for(auto i=0; i<cellPartition_.size(); ++i)
        {
            if (cellPartition_[i] != ifile) continue;

            cellIds_.insert(i+1);
            for(auto it : bigBody.cell(i+1)) if(nodeIds_.count(it)==0) nodeIds_.insert(it);
        }
        for(auto it : nodeIds_) nodeIdG2L_[ifile].insert({it, nodeIdG2L_[ifile].size()+1});

        // write coordinate
        this->rwNode(ifile);

        // write body-section
        this->rwBody(ifile, bigBody);

        // write bdy-section
        this->rwBoundary(ifile);

        // write interface
        this->rwInterface(ifile);

        // clear
        for(auto i = 0; i<=ifile; ++i)
        {
            if(outerFace_[i].empty())
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

    if(SIZE>np)
    {
        if(MPIAdapter::isMaster()) check(true, "", format("MPI size should not bigger than subMesh number: %d > %d\n", SIZE, np));

        MPIAdapter::Finalize();
        exit(EXIT_SUCCESS);
    }

    // open mesh file to read/write
    bigMesh_ = std::make_shared<CGFile>(meshFilename);
    this->updateOwnerThread(meshFilename, np);

    // handle bodysection
    auto bigBody = bigMesh_->bodySection();

    // decompose body-section into np parts, 
    // and collect subBody belong to each process
    auto comPartId = [](const Cell& lhs, const Cell& rhs){ return lhs.partId < rhs.partId; };
    auto comId = [](const Cell& lhs, const Cell& rhs){ return lhs.id < rhs.id; };
    auto subBody = collect_subBody(this->decompose_body(bigBody, np), ownerThread_);
    if(subBody.size()==0) std::cout << format("Thread %d has no subBody\n", RANK);
    std::sort(subBody.begin(), subBody.end(), comPartId);
    auto pos = subBody.begin();
    while (pos != subBody.end())
    {
        // relod body-section
        auto lastPos = pos;
        const auto ifile = lastPos->partId;
        pos = std::upper_bound(lastPos, subBody.end(), *lastPos, comPartId);
        std::sort(lastPos, pos, comId);
        auto indexLB = lastPos->id, indexUB = (pos-1)->id;
        auto curBody = bigMesh_->loadSection(bigMesh_->bodySectionIdList(), indexLB, indexUB);
        // check(true, "", format("  Thread %d: load %d cells [%d-%d] of part %d\n", RANK, indexUB-indexLB+1, indexLB, indexUB, ifile));

        // cell and node
        cellIds_.clear(), nodeIds_.clear();
        for(auto cell : subBody)
        {
            if(cell.partId != ifile) continue;
            cellIds_.insert(cell.id);
            for(auto it : curBody.cell(cell.id)) if(nodeIds_.count(it)==0) nodeIds_.insert(it);
        }
        nodeIdG2L_.insert({ifile, {}});
        if(!nodeIds_.empty()) for(auto it : nodeIds_) nodeIdG2L_[ifile].insert({it, nodeIdG2L_[ifile].size()+1});

        // write coordinate
        this->rwNode(ifile);

        // write body
        this->rwBody(ifile, curBody);

        // write bdy
        this->rwBoundary(ifile);
    }
    bigBody.clear();

    // collect interface
    this->collect_interface();

    // write interface
    for(auto ifile : ownerFile_) this->writeInterface(ifile, np);

    // clear
    bigMesh_->close();
    for(auto it : smallMesh_) if(it != nullptr) it->close();

    
    MPIAdapter::Barrier();
}

int MetisCutter::cut_metis(const vector<vector<idx_t>>& cellToplogy, idx_t np, vector<idx_t>& cellPartition, vector<idx_t>& nodePartition)
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
    std::vector<idx_t> eind, eptr(cellToplogy.size()+1, 0);
    cgsize_t offset = 0;
    for(auto i=0; i<cellToplogy.size(); ++i)
    {
        for(auto jt : cellToplogy[i])
        {
            eind.push_back(jt-1);
        }
        offset += cellToplogy[i].size();
        eptr[i+1] = offset;
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

int MetisCutter::cut_parmetis(idx_t np, vector<idx_t>& elmdist, vector<idx_t>& eptr, vector<idx_t>& eind, vector<idx_t>& cellPartition)
{
    // parallel read element connectivity data
    idx_t wgtflag = 0;      // 0: no weights(elmwgt is null), 2: weights on the vertices only
    idx_t numflag = 0;      // 0: C-style numbering, 1: Fortran-style numbering
    idx_t ncon = 1;         // the number of weights that each vertex has.
    idx_t ncommonnodes = 2; //should be greater than 0, 2 is ok for most meshes.
    idx_t edgecut = 0;
    idx_t opts[3] = {0, 1, static_cast<idx_t>(std::time(nullptr))};
    real_t tpwgts[ncon * np];
    real_t ubvec[ncon];

    for (auto i = 0; i < ncon * np; ++i)
    {
        tpwgts[i] = 1.0 / (real_t)np;
    }
    for (auto i = 0; i < ncon; ++i)
    {
        ubvec[i] = 1.05;
    }

    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    return ParMETIS_V3_PartMeshKway(
        elmdist.data(),
        eptr.data(),
        eind.data(),
        NULL, // weights of the elements
        &wgtflag,
        &numflag,
        &ncon,
        &ncommonnodes,
        &np,
        tpwgts,
        ubvec,
        opts,
        &edgecut,
        cellPartition.data(),
        &mpi_comm
    );
}

MetisCutter::DecomposeResult MetisCutter::decompose_body(const CGFile::Section& bigBody, const int np)
{
    DecomposeResult result;

    auto elmdist = distribute(bigBody.end - bigBody.start + 1, SIZE, result.nCellThisRank);
    result.start = bigBody.start + elmdist[RANK];

    auto shareBigbody = bigMesh_->loadSection(bigMesh_->bodySectionIdList(), result.start, result.start+result.nCellThisRank-1);

    std::vector<idx_t> eind, eptr{0};
    idx_t beg = shareBigbody.isMixed() ? 1 : 0;
    for (auto id = result.start; id < result.start + result.nCellThisRank; ++id)
    {
        auto &tmp = shareBigbody.cell(id);
        for (int j = beg; j < tmp.size(); ++j) { eind.push_back(tmp[j] - 1); }
        eptr.push_back(eptr.back()+tmp.size()-beg);
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
        for(auto i=1; i<=SIZE; ++i)
        {
            for(auto j=tmp[i-1]; j<tmp[i]; ++j) ownerThread_[j] = i-1;
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
    if(ownerThread_[id] == RANK)
    {
        smallMesh_[id] = std::make_shared<CGFile>(smallMeshName(filename, id, np), CG_MODE_WRITE);
        ownerFile_.insert(id);
    }
}

void MetisCutter::rwBody(const int ifile, const CGFile::Section& bigBody)
{
    auto &subFile = smallMesh_[ifile];

    // new a sectionk
    auto &curS = subFile->addSection();
    strcpy(curS.name, bigBody.name);
    curS.cellType = bigBody.cellType;
    for (auto id : cellIds_)
    {
        curS.addCell(id, bigBody.cell(id), bigBody.flag(id));
    }

    // write data
    subFile->writeSection(curS, nodeIdG2L_[ifile]);

    // global info  
    auto len = curS.end - curS.start + 1;
    subFile->writeGlobalInfo(curS, bigBody.data.size(), globalOffset_, globalOffset_+len);
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

    for(auto ithSection : bigMesh_->bdySectionIdList())
    {
        auto &subBdy = subFile->addSection();
        auto &bigBdy = bigMesh_->loadSection(ithSection);
        strcpy(subBdy.name, bigBdy.name);
        subBdy.cellType = bigBdy.cellType;

        cgsize_t ID = bigBdy.start;
        for(auto i = bigBdy.start; i<=bigBdy.end; ++i)
        {
            auto face = bigBdy.cell(i);
            auto faceStr = CGFile::stringAFace(face);
            if(curOuterFace.find(faceStr) != curOuterFace.end())
            {
                curOuterFace.erase(faceStr);
                subBdy.addCell(ID++, face, bigBdy.flag(i));
            }
        }

        if(!subBdy.data.empty()) subFile->writeSection(subBdy, nodeIdG2L_[ifile]);

        // clear
        subBdy.clear();
    }
}

void MetisCutter::rwInterface(const int id)
{
    auto &subFile = smallMesh_[id];
    auto &curOuterFace = outerFace_[id];

    for(auto nbr = 0; nbr < id; ++nbr)
    {
        if (outerFace_[nbr].empty()) continue;

        auto &nbrOuterFace = outerFace_[nbr];

        set<cgsize_t> typeFlags;
        CGFile::Section curS, nbrS;
        for(auto face = nbrOuterFace.begin(); face!=nbrOuterFace.end();)
        {
            auto it = curOuterFace.find(face->first);
            if(it != curOuterFace.end())
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
        if(!curS.data.empty())
        {
            strcpy(curS.name, format("Interface_%d", nbr).c_str());
            strcpy(nbrS.name, format("Interface_%d", id).c_str());

            if(typeFlags.size()>1){
                curS.cellType = ElementType_t::MIXED;
                nbrS.cellType = ElementType_t::MIXED;
                auto curOffset = 0;
                for(auto it : curS.data)
                {
                    curS.offset.push_back(curOffset);
                    nbrS.offset.push_back(curOffset);
                    curOffset += it.size()+1;
                }
            } 
            else if(typeFlags.size()==1){
                auto type = (*typeFlags.begin()==3) ? ElementType_t::TRI_3 : ElementType_t::QUAD_4;
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
    for(auto id : nodeIds_)
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
    char fmt[128], smallFilename[128];
    sprintf(fmt, "%%s.%d.%%0%dd.cgns", np, int(std::log10(np) + 1));
    sprintf(smallFilename, fmt, bigMeshName.substr(0, bigMeshName.size() - 5).c_str(), id);

    return smallFilename;
}

void MetisCutter::updateOuterFace(const CGFile::Section& curS, const int ifile)
{
    auto &curOuterFace = outerFace_[ifile];

    for (auto i = 0; i < curS.data.size(); ++i)
    {
        auto cellType = curS.isMixed() ? CGFile::CellType(curS.typeFlag[i]) : curS.cellType;
        for (auto it : CGFile::allFaceInCell(curS.data[i], cellType))
        {
            auto val = CGFile::stringAFace(it);
            if (curOuterFace.find(val) != curOuterFace.end()) { curOuterFace.erase(val); }
            else
            {
                curOuterFace.insert({val, it});
            }
        }
    }
}

void MetisCutter::writeGlobalInfo(const CGFile::Section& curS, const int id)
{
    auto len = curS.end - curS.start + 1;
    smallMesh_[id]->writeGlobalInfo(curS, 10, globalOffset_, globalOffset_+len);
    globalOffset_ += len;
}

void MetisCutter::writeInterface(const int id, const int np)
{
    auto &subFile = smallMesh_[id];
    auto &curOuterFace = outerFace_[id];

    for(auto nbr = 0; nbr < np; ++nbr)
    {
        if (nbr == id) continue;
        auto nbrIsOwner = ownerFile_.count(nbr) != 0;

        auto &nbrOuterFace = outerFace_[nbr];

        set<cgsize_t> typeFlags;
        CGFile::Section curS, nbrS;
        for(auto face = nbrOuterFace.begin(); face!=nbrOuterFace.end();)
        {
            auto it = curOuterFace.find(face->first);
            if(it != curOuterFace.end())
            {
                curS.data.emplace_back(it->second);
                auto t = CGFile::CellType(it->second.size(), 2);
                curS.typeFlag.push_back(t);
                typeFlags.insert(it->second.size());
                it = curOuterFace.erase(it);

                if(nbrIsOwner)
                {
                    nbrS.data.push_back(face->second);
                    nbrS.typeFlag.push_back(t);
                    face = nbrOuterFace.erase(face);
                }
            }
            else
            {
                ++face;
            }
        }
        if(!curS.data.empty())
        {
            strcpy(curS.name, format("Interface_%d", nbr).c_str());
            if(nbrIsOwner) strcpy(nbrS.name, format("Interface_%d", id).c_str());

            if(typeFlags.size()>1){
                curS.cellType = ElementType_t::MIXED;
                if(nbrIsOwner) nbrS.cellType = ElementType_t::MIXED;
            } 
            else if(typeFlags.size()==1){
                auto type = (*typeFlags.begin()==3) ? ElementType_t::TRI_3 : ElementType_t::QUAD_4;
                curS.cellType = type;
                if(nbrIsOwner) nbrS.cellType = type;
            }
            subFile->writeSection(curS, nodeIdG2L_[id]);
            if(nbrIsOwner) smallMesh_[nbr]->writeSection(nbrS, nodeIdG2L_[nbr]);   
        }

        // clear
        vector<vector<cgsize_t>>{}.swap(curS.data);
        vector<cgsize_t>{}.swap(curS.offset);
        vector<cgsize_t>{}.swap(curS.typeFlag);            
    } 
}

} // namespace MeshCut
