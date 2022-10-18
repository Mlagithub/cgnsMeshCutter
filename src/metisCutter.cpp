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

static const int LenTAG = 100, DataTAG = 1000;

auto distribute = [&](const idx_t len, const int n, idx_t &nCellThisRank) {
    vector<idx_t> dist(n + 1, 0); // same for all threads
    {
        auto step = (len + n - 1) / n;
        nCellThisRank = MPIAdapter::isLast() ? len - (n - 1) * step : step;
        MPI_Allgather(&nCellThisRank, 1, GetMPIDataType<idx_t>(), dist.data() + 1, 1, GetMPIDataType<idx_t>(), MPI_COMM_WORLD);
        for (int iProc = 1; iProc < n + 1; ++iProc) { dist[iProc] += dist[iProc - 1]; }
    }
    return dist;
};

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

// collect cell id from all other threads
vector<MetisCutter::Cell> MetisCutter::collect_subBody(const MetisCutter::DecomposeResult& decomposeResult, const vector<int>& ownerThread)
{
    map<int, vector<Cell>> cellInfo;
    {
        for (idx_t i = 0; i < decomposeResult.nCellThisRank; ++i)
        {
            auto partId = cellPartition_[i];
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
    this->openSubMeshToWrite(meshFilename, np);

    for(auto it : bigMesh_->bodySections())
    {
        // load data into memory from file
        auto &bigBody = bigMesh_->loadSection(it);
        
        check(true, "", format("Call METIS_V3_PartMeshKway decompose: %s\n", bigBody.name));
        // cut
        check(
            cut_metis(bigBody.data, np, cellPartition_, nodePartition_) == METIS_OK,
            format("Call METIS_V3_PartMeshKway decompose: %s return error %s", bigBody.name),
            format("Call METIS_V3_PartMeshKway decompose: %s successfuly\n", bigBody.name));

        // write sub-mesh data to file
        for(auto ifile=0; ifile<np; ++ifile)
        {
            // collect cell/node id of sub-mesh
            nodeIdG2L_.insert({ifile, {}}); // 1-base id
            nodeIds_.clear(); // 1-base id
            cellIds_.clear(); // 0-base id
            for(auto i=0; i<cellPartition_.size(); ++i)
            {
                if (cellPartition_[i] != ifile) continue;

                cellIds_.insert(i);
                for(auto it : bigBody.data[i]) if(nodeIds_.count(it)==0) nodeIds_.insert(it);
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
    }

    // clear
    bigMesh_->close();
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
    auto ownerThread = this->openSubMeshToWrite(meshFilename, np);
    for(auto i=0; i<SIZE; ++i) if(ownerThread[i] == RANK) nodeIdG2L_.insert({i, {}});

    // handle bodysection
    for(auto it : bigMesh_->bodySections())
    {
        auto &bigBody = bigMesh_->section(it);

        // decompose body-section into np parts, 
        // and collect subBody belong to each process
        auto subBody = collect_subBody(this->decompose_body(bigBody, np), ownerThread);
        if(subBody.size()==0) continue;

        // write coordinate
        auto comPartId = [](const Cell& lhs, const Cell& rhs){ return lhs.partId < rhs.partId; };
        auto comId = [](const Cell& lhs, const Cell& rhs){ return lhs.id < rhs.id; };
        std::sort(subBody.begin(), subBody.end(), comPartId);
        auto pos = subBody.begin();
        while (pos != subBody.end())
        {
            auto lastPos = pos;
            pos = std::upper_bound(lastPos, subBody.end(), *lastPos, comPartId);
            std::sort(lastPos, pos, comId);
            auto indexLB = lastPos->id, indexUB = (pos-1)->id;
            bigMesh_->loadSection(it, indexLB, indexUB);
            check(true, "", format("Thread %d: load %d cells [%d-%d] of part %d\n", RANK, indexUB-indexLB+1, indexLB, indexUB, lastPos->partId));

            // cell and node
            cellIds_.clear(), nodeIds_.clear();
            for(auto cell : subBody)
            {
                if(cell.partId != lastPos->partId) continue;
                auto id = cell.id;
                cellIds_.insert(id);
                for(auto it : bigBody.data[id - indexLB]) if(nodeIds_.count(it)==0) nodeIds_.insert(it);
            }
            if(!nodeIds_.empty()) for(auto it : nodeIds_) nodeIdG2L_[lastPos->partId].insert({it, nodeIdG2L_[lastPos->partId].size()+1});

            this->rwNode(lastPos->partId);

            // write body
            this->rwBody(lastPos->partId, bigBody);

            // write bdy
            this->rwBoundary(lastPos->partId);

            // write interface            
        }

        // clear
        for(auto it : smallMesh_) if(it != nullptr) it->close();

    }
    MPIAdapter::Barrier();
}

int MetisCutter::cut_metis(const vector<vector<cgsize_t>>& cellToplogy, idx_t np, vector<idx_t>& cellPartition, vector<idx_t>& nodePartition)
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
    idx_t opts[3] = {0, 1, std::time(0)};
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

    idx_t nCell = bigBody.end - bigBody.start + 1;
    auto elmdist = distribute(nCell, SIZE, result.nCellThisRank);
    result.start = bigBody.start + elmdist[RANK];
    auto &curS = bigMesh_->loadSection(bigBody.id, result.start, result.start + result.nCellThisRank - 1);
    auto dataSize = std::accumulate(curS.data.begin(), curS.data.end(), 0, [](idx_t sum, vector<idx_t> vec) { return sum + vec.size(); });

    idx_t dataOffset = 0, beg = curS.isMixed() ? 1 : 0;
    std::vector<idx_t> eind(dataSize - result.nCellThisRank * beg, 0), eptr(result.nCellThisRank + 1, 0);
    for (int i = 0; i < result.nCellThisRank; ++i)
    {
        auto &tmp = curS.data[i];
        for (int j = beg; j < tmp.size(); ++j) { eind[dataOffset++] = tmp[j] - 1; }
        eptr[i + 1] = eptr[i] + tmp.size() - beg;
    }
    check(true, "", format("ParMETIS_V3_PartMeshKway begin divide: %s\n", bigBody.name));
    cellPartition_.assign(result.nCellThisRank, 0);
    check(cut_parmetis(np, elmdist, eptr, eind, cellPartition_) == METIS_OK, format("Process %d: ParMETIS_V3_PartMeshKway return error while divide section: %s\n", RANK, bigBody.name), format("ParMETIS_V3_PartMeshKway done\n"));

    return result;
}

vector<int> MetisCutter::openSubMeshToWrite(string bigFileName, const int np)
{
    idx_t dummy;

    // calculate partId collection of each thread
    vector<int> ownerThread(np, 0);
    {
        auto tmp = distribute(np, SIZE, dummy);
        for(auto i=1; i<=SIZE; ++i)
        {
            for(auto j=tmp[i-1]; j<tmp[i]; ++j) ownerThread[j] = i-1;
        }
    }

    smallMesh_.assign(np, nullptr);
    for(auto i=0; i<np; ++i)
    {
        if(ownerThread[i] == RANK) smallMesh_[i] = std::make_shared<CGFile>(smallMeshName(bigFileName, i, np), CG_MODE_WRITE);
    }

    return ownerThread;
}

void MetisCutter::rwBody(const int ifile, const CGFile::Section& bigBody)
{
    auto &subFile = smallMesh_[ifile];

    // new a section
    auto &curS = subFile->addSection();
    strcpy(curS.name, bigBody.name);
    curS.cellType = bigBody.cellType;
    int n = 0, curOffset = 0;
    cgsize_t start = *(cellIds_.begin());
    for (auto id : cellIds_)
    {
        curS.data.emplace_back(bigBody.data[id-start]);
        if (curS.isMixed())
        {
            curS.typeFlag.push_back(bigBody.typeFlag[id]);
            curS.offset.push_back(curOffset);
            curOffset += curS.data.back().size()+1;
        }
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
    vector<vector<cgsize_t>>{}.swap(curS.data);
    vector<cgsize_t>{}.swap(curS.offset);
    vector<cgsize_t>{}.swap(curS.typeFlag);
}

void MetisCutter::rwBoundary(const int ifile)
{
    auto &subFile = smallMesh_[ifile]; 
    auto &curOuterFace = outerFace_[ifile];

    for(auto ithSection : bigMesh_->bdySections())
    {
        auto &subBdy = subFile->addSection();
        auto &bigBdy = bigMesh_->loadSection(ithSection);
        strcpy(subBdy.name, bigBdy.name);
        subBdy.cellType = bigBdy.cellType;

        int curOffset = 0, i=0;
        for(auto face : bigBdy.data)
        {
            auto faceStr = CGFile::stringAFace(face);
            if(curOuterFace.find(faceStr) != curOuterFace.end())
            {
                curOuterFace.erase(faceStr);
                subBdy.data.emplace_back(face);
                if(bigBdy.isMixed())
                {
                    subBdy.typeFlag.push_back(bigBdy.typeFlag[i]);
                    subBdy.offset.push_back(curOffset);
                    curOffset += (face.size() + 1);
                }
            }
            ++i;
        }

        if(!subBdy.data.empty()) subFile->writeSection(subBdy, nodeIdG2L_[ifile]);

        // clear
        vector<vector<cgsize_t>>{}.swap(subBdy.data);
        vector<cgsize_t>{}.swap(subBdy.offset);
        vector<cgsize_t>{}.swap(subBdy.typeFlag);
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
            strcpy(curS.name, format("Interface_%d_%d", 0, nbr).c_str());
            strcpy(nbrS.name, format("Interface_%d_%d", 0, id).c_str());

            if(typeFlags.size()>1){
                curS.cellType = ElementType_t::MIXED;
                nbrS.cellType = ElementType_t::MIXED;
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
    sprintf(fmt, "%%s_%%0%dd.cgns", int(std::log10(np) + 1));
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

} // namespace MeshCut
