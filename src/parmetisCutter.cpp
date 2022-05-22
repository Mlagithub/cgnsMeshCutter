#include "parmetisCutter.h"
#include "timeStat.h"
#include "mpiAdapter.h"
#include "timeStat.h"
#include "format.h"
// #include "fmt/format.h"

#include <parmetis.h>
#include <pcgnslib.h>
#include <cgnslib.h>
#include <cmath>
#include <ctime>
#include <array>
#include <algorithm>
#include <cstring>
#include <numeric>
#include <ostream>
#include <iterator>
#include <unistd.h>


#define IS_FIXED_SIZE(type) ((type >= CGNS_ENUMV(NODE) && \
                              type <= CGNS_ENUMV(HEXA_27)) || \
                              type == CGNS_ENUMV(PYRA_13) || \
                             (type >= CGNS_ENUMV(BAR_4) && \
                              type <= CGNS_ENUMV(HEXA_125)))


// update this func with func nodeIdToBuildFaceInCell
static ElementType_t CellType(const int CGNSCellTypeFalg)
{
    auto rst = ElementType_t::ElementTypeNull;
    switch (CGNSCellTypeFalg)
    {
    case 5:
        rst = ElementType_t::TRI_3;
        break;
    case 7:
        rst = ElementType_t::QUAD_4;
        break;
    case 10:
        rst = ElementType_t::TETRA_4;
        break;
    case 12:
        rst = ElementType_t::PYRA_5;
        break;
    case 14:
        rst = ElementType_t::PENTA_6;
        break;
    case 17:
        rst = ElementType_t::HEXA_8;
        break;
    default:
        //fmt::format("unknown element type", static_cast<ElementType_t>(CGNSCellTypeFalg),'\n');
        break;
    }
    return rst;
}

static int NodeCountOfFaceCell(const ElementType_t FaceCellType)
{
    switch (FaceCellType)
    {
    case ElementType_t::TRI_3 :
        return 3;
        break;
    case ElementType_t::QUAD_4:
        return 4;
        break;
    default:
        //fmt::format("Wrong to get node count of face type: ", FaceCellType, '\n');
        break;
    }
}

static ElementType_t CellType(const int nodeCnt, const int dim)
{
    switch (dim)
    {
    case 2:
        switch (nodeCnt)
        {
        case 3:
            return ElementType_t::TRI_3;
            break;
        case 4:
            return ElementType_t::QUAD_4;
            break;
        default:
            //fmt::format("not supported element type with[dim, nodes] [", dim, nodeCnt, "]\n");
            break;
        }
        break;
    case 3:
        switch (nodeCnt)
        {
        case 4:
            return ElementType_t::TETRA_4;
            break;
        case 5:
            return ElementType_t::PYRA_5;
            break;
        case 6:
            return ElementType_t::PENTA_6;
            break;
        case 8:
            return ElementType_t::HEXA_8;
            break;
        default:
            //fmt::format("not supported element type with[dim, nodes] [", dim, nodeCnt, "]\n");
            break;
        }
        break;
    default:
        break;
    }

    return ElementTypeNull;
}

// update this func with func CellType
// 各种单元类型每个面的节点在单元内部的编号。
// 由于面的法向量具有方向性，这些数据在计算相关向量时需要使用。
const std::unordered_map<ElementType_t, std::vector<std::set<int>>> nodeIdToBuildFaceInCell = {
    {ElementType_t::TETRA_4, {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}},
    {ElementType_t::HEXA_8,  {{0, 1, 5, 4},{1, 2, 6, 5},{2, 3, 7, 6},{3, 0, 4, 7},{0, 3, 2, 1},{4, 5, 6, 7}}},
    {ElementType_t::PENTA_6, {{0, 1, 4, 3}, {1,2,5,4}, {2,0,3,5}, {0,2,1}, {3,4,5}}},
    {ElementType_t::PYRA_5,  {{0, 3, 2, 1}, {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}}}
};


namespace MeshCut
{
    bool isSupportedElementType(const ElementType_t& t);
    
    bool isBoundaryElementType(const ElementType_t& t, const int Dim);

    template<typename ...Args>
    void checkPCGNSCall(int errorCode, Args ... args) {
        if (errorCode != CG_OK)
        {
            //fmt::format("Wrong at parmetisDecomposer. ", std::forward<Args>(args)...);
            // cgp_error_exit();
        }
    }
    
    template<typename ...Args>
    void checkCGNSCall(int errorCode, Args ... args) {
        if (errorCode != CG_OK)
        {
            //fmt::format("Wrong at parmetisDecomposer. ", std::forward<Args>(args)...);
            // cg_error_exit();
        }
    }

    struct ParMetisMeshCutter::SectionPartition
    {
        // section id, build in CGNS file
        int id=-1;

        // section name
        std::string name;

        // 
        ElementType_t elementType = ElementType_t::ElementTypeNull;
        
        // global ele id of this section
        // each section's ele id is continutius
        cgsize_t start = 0, end = 0;
        cgsize_t cell_size_this_rank = 0, nodeSize = 0;

        // different for mixed and no-mixed.
        // no-mixed = ele_size * node_size_per_ele
        // mixed = total node number of all ele + ele_size;
        cgsize_t element_data_size;
        
        // ============= mixed element ==============
        // type flag of first ele, int value, eg,17(HEXA_8)
        // node_1_id,
        // node_2_id,
        // ...
        // node_8_id,
        // type flag of second ele, eg, 14(PENTA_6)
        // node_1_id,
        // node_2_id,
        // ...
        // node_6_id,
        // ...
        // ============= un-mixed element ==============
        // node_1_id of ele_1
        // node_2_id of ele_1
        // ...
        // node_n_id of ele_1,
        // ...
        // node_1_id of ele_m,
        // node_2_id of ele_m,
        // ...
        // node_n_id of ele_m,
        // m is the ele number, n is the node number of element,
        // so all element has same type. 
        void* element_connectivity;
        CGIO::Section::NodeDataType connectivityDataType;

        // only exit for mixed ele type
        // array.size == ele size + 1(for parallel, local ele size)
        // 0
        // node_size_of_ele_1+1,
        // node_size_of_ele_1+1 + node_size_of_ele_2+1,
        // node_size_of_ele_1+1 + node_size_of_ele_2+1 + node_size_of_ele_3+1
        // ...
        std::vector<cgsize_t> connect_offset;

        // element partition
        std::vector<idx_t> elempartition;

        cgsize_t* parent_data = nullptr;
        int nbndry, parent_flag;


        void init(const int id, CGIO::Section& cgio_section)
        {
            auto lowerBnd = cgio_section.startID();
            auto upperBnd = cgio_section.endID();   

            this->id = id;                
            this->name = cgio_section.GetName();
            this->elementType = cgio_section.GetElementType();
            this->cell_size_this_rank = this->end - this->start + 1;

            this->connectivityDataType = cgio_section.GetDataType("ElementConnectivity").first;
            this->element_data_size = cgio_section.GetDataSize("ElementConnectivity", this->start, this->end);

            switch (this->connectivityDataType)
            {
            case CGIO::Section::I4:
                this->element_connectivity = new int[this->element_data_size];
                break;
            case CGIO::Section::I8:
                this->element_connectivity = new cglong_t[this->element_data_size];
                break;
            default:
                break;
            }

            this->connect_offset.assign(this->end-this->start+2, 0);
            cgio_section.ReadElementConnectData(this->element_connectivity, this->connect_offset, this->start, this->end);
        }

    };

    ParMetisMeshCutter::ParMetisMeshCutter(int argc, char** argv)
    {
        int flag;
        MPI_Initialized(&flag);
        if (!flag)
        {
            MPI_Init(&argc, &argv);
        }
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    }

    ParMetisMeshCutter::~ParMetisMeshCutter()
    {
        int flag;
        MPI_Finalized(&flag);
        if (!flag)
        {
            MPI_Finalize();
        }
    }

    bool ParMetisMeshCutter::checkCGNSFile(const std::shared_ptr<CGIO::CGFile> fp)
    {
        // check number of Base
        if (fp->nDataBase() != 1)
        {
            //fmt::format("Thread", this->rank_, "Wrong number of Base:", fp->nDataBase(), "should be 1.\n");
            // cg_error_exit();
            return false;
        }

        // check number of zone
        if (fp->GetDatabase().nZone() != 1)
        {
            //fmt::format("Thread", this->rank_, "Wrong number of Zone:", fp->GetDatabase().nZone(), "should be 1.\n");
            // cg_error_exit();
            return false;
        }

        // check zone type (structure/unstructure)
        // TODO: only support unstuucture
        if (fp->GetDatabase().GetZone().type() != CGIO::ZoneType::Unstructured)
        {
            // fmt::format(
            //     "Thread", this->rank_, "Wrong type of Zone:", [&]() {
            //         std::string rst;
            //         switch (fp->GetDatabase().GetZone().type())
            //         {
            //         case ZoneType_t::Structured:
            //             rst = "Structured";
            //             break;
            //         case ZoneType_t::Unstructured:
            //             rst = "Unstructured";
            //             break;
            //         case ZoneType_t::ZoneTypeNull:
            //             rst = "ZoneTypeNull";
            //             break;
            //         case ZoneType_t::ZoneTypeUserDefined:
            //             rst = "ZoneTypeUserDefined";
            //             break;
            //         default:
            //             break;
            //         }
            //         return rst;
            //     }(),
            //     "should be Unstructured\n");
            // cg_error_exit();
            return false;
        }
        return true;
    }

    void ParMetisMeshCutter::cut(const std::string &mesh, const int npart, const int nx, const int ny, const int nz)
    {
        this->isMaster_ = (this->rank_==0);
        if(this->isMaster_) //fmt::format("Thread", this->rank_, "Begin ParmetisDecomposer\n");
        nparts_ = npart;
        if(this->nparts_<2)
        {
            //fmt::format("Do nothing with decompose number", this->nparts_, " Please check your config file.\n");
            return;
        }

        TimeStat start;
        // open original mesh file
        bigFile_ = std::make_shared<CGIO::CGFile>(CGIO::CGFile());
        bigFile_->Open(mesh.c_str());

        // check support
        if(!this->checkCGNSFile(bigFile_)) //fmt::format("check", bigFile_->name(), "field, from", __FILE__, __LINE__, "\n");

        // 
        this->totalNode_=bigFile_->GetDatabase().GetZone().nVertex();
        this->totalCell_=bigFile_->GetDatabase().GetZone().nCell();
        if(isMaster_) //fmt::format(format("Thread %d zone %-20s [node, element]->[%6d, %6d]\n", this->rank_, bigFile_->GetDatabase().GetZone().GetName().c_str(), this->totalNode_, this->totalCell_));
        for(auto ithSection = 0; ithSection < bigFile_->GetDatabase().GetZone().nSections(); ++ithSection)
        {
            auto curSection = bigFile_->GetDatabase().GetZone().GetSection(ithSection);
            // if(isMaster_) //fmt::format(format("Thread %d   section %-20s [%2d, %6d, %6d].  is body-section: %d\n", this->rank_, curSection.GetName().c_str(), curSection.GetElementType(), curSection.startID(), curSection.endID(), curSection.isBodySection()));
        }
        
        this->decompose_body_section();
        this->rebuild_cgnsFile(mesh);


        TimeStat end;
        //fmt::format("Cut time: {}\n", end-start);
        return;
    }

    void ParMetisMeshCutter::decompose_body_section()
    {
        // parallel read element connectivity data
        // MPI_Barrier(MPI_COMM_WORLD);
        idx_t wgtflag = 0;      // 0: no weights(elmwgt is null), 2: weights on the vertices only
        idx_t numflag = 0;      // 0: C-style numbering, 1: Fortran-style numbering
        idx_t ncon = 1;         // the number of weights that each vertex has.
        idx_t ncommonnodes = 2; //should be greater than 0, 2 is ok for most meshes.
        idx_t edgecut = 0;
        idx_t opts[3] = {0, 1, std::time(0)};
        real_t tpwgts[ncon * nparts_];
        real_t ubvec[ncon];

        for (auto i = 0; i < ncon * nparts_; ++i)
        {
            tpwgts[i] = 1.0 / (real_t)nparts_;
        }
        for (auto i = 0; i < ncon; ++i)
        {
            ubvec[i] = 1.05;
        }

        // TODO: cmdline 
        // parmetis use default option
        if (true)
        {
            opts[0] = 1;
        }
        else{
            // opts[1] = (idx_t)cfg.partitionDir().subDir("parmetis").lookup("infoLevel");
            // opts[2] = (idx_t)cfg.partitionDir().subDir("parmetis").lookup("randomSeed");
        }

        // get element connectivity data, support any type of element, include mixed
        auto zones = bigFile_->GetDatabase().GetZone();
        sectionPartition.assign(zones.nSections(), SectionPartition());

        for(auto ithSection = 0; ithSection < zones.nSections(); ++ithSection)
        {
            auto &cgio_section = zones.GetSection(ithSection);
            if(cgio_section.isBodySection())
            {
                bodySections_.emplace_back(ithSection);
            }
            else 
            {
                boundarySections_.push_back(ithSection);
            } 
            
            auto &curSection = sectionPartition[ithSection];

            auto upperBnd = cgio_section.endID(), lowerBnd = cgio_section.startID();
            if(cgio_section.isBodySection())
            {
                auto tmpEleSizeLocal = static_cast<cgsize_t>((upperBnd - lowerBnd + 1) / size_);
                curSection.start = rank_ * tmpEleSizeLocal + lowerBnd;
                curSection.end = curSection.start + tmpEleSizeLocal-1;
                if(rank_ == size_-1) curSection.end += (upperBnd-lowerBnd+1)-size_*tmpEleSizeLocal;
            }
            else{
                curSection.start = lowerBnd;
                curSection.end = upperBnd;
            }

            curSection.init(ithSection, cgio_section);

            //fmt::format(format("[Thread, start, end, cell_size_this_rank] => [%d, %d, %d, %d]\n", rank_, curSection.start, curSection.end, curSection.cell_size_this_rank));

        }
        MPI_Barrier(MPI_COMM_WORLD);
        

        for(auto ithSection : bodySections_)
        {
            auto &curSection = sectionPartition[ithSection];

            // compute a k-way partition
            // distribute elements of the mesh to all processors
            // elmdist, same for all threads
            std::vector<idx_t> elmdist(this->size_+1, 0);
            idx_t sv = static_cast<idx_t>(curSection.cell_size_this_rank);
            MPI_Allgather(&sv, 1, GetMPIDataType<idx_t>(), elmdist.data()+1, 1, GetMPIDataType<idx_t>(), MPI_COMM_WORLD);
            for(int iProc = 1; iProc < this->size_+1; ++iProc)
            {
                elmdist[iProc] += elmdist[iProc-1];
            }
            
            std::vector<idx_t> eptr(curSection.cell_size_this_rank + 1, 0);
            std::vector<idx_t> eind;

            // for mixed element type
            if(curSection.elementType == CGNS_ENUMV(ElementType_t::MIXED))
            {
                eind.assign(curSection.element_data_size-curSection.cell_size_this_rank, 0);
                auto dataOffset=0;
                switch (curSection.connectivityDataType)
                {
                case CGIO::Section::NodeDataType::I4:
                    for(int i=0; i<curSection.cell_size_this_rank; ++i)
                    {
                        auto n=0;
                        for(auto j=curSection.connect_offset[i]+1; j<curSection.connect_offset[i+1]; ++j)
                        {
                            eind[dataOffset++] = ((int*)curSection.element_connectivity)[j];
                            ++n;
                        }
                        eptr[i+1] = eptr[i] + n;
                    }
                    break;
                case CGIO::Section::NodeDataType::I8:
                    for(int i=0; i<curSection.cell_size_this_rank; ++i)
                    {
                        auto n=0;
                        for(auto j=curSection.connect_offset[i]+1; j<curSection.connect_offset[i+1]; ++j)
                        {
                            eind[dataOffset++] = ((cglong_t*)curSection.element_connectivity)[j];
                            ++n;
                        }
                        eptr[i+1] = eptr[i] + n;
                    }
                    break;
                default:
                    break;
                } 
            }
            else
            {
                int nodeSizePerCell = 0;
                cg_npe((ElementType_t)curSection.elementType, &nodeSizePerCell);
                
                for(int i=1; i< curSection.cell_size_this_rank+1;++i)
                {
                    eptr[i] = eptr[i-1] + nodeSizePerCell;
                }

                eind.assign(curSection.cell_size_this_rank*nodeSizePerCell, 0);
                switch (curSection.connectivityDataType)
                {
                case CGIO::Section::NodeDataType::I4:
                    for (int i = 0; i < curSection.cell_size_this_rank * nodeSizePerCell; ++i)
                    {
                        eind[i] = ((int*)curSection.element_connectivity)[i];
                    }
                    break;
                case CGIO::Section::NodeDataType::I8:
                    for (int i = 0; i < curSection.cell_size_this_rank * nodeSizePerCell; ++i)
                    {
                        eind[i] = ((cglong_t*)curSection.element_connectivity)[i];
                    }
                    break;
                default:
                    break;
                }
            }

            TimeStat stime;
            curSection.elempartition.assign(curSection.cell_size_this_rank,0);// result array.
            MPI_Comm mpi_comm = MPI_COMM_WORLD;
            auto stat = ParMETIS_V3_PartMeshKway(
                elmdist.data(),
                eptr.data(),
                eind.data(),
                NULL, // weights of the elements
                &wgtflag,
                &numflag,
                &ncon,
                &ncommonnodes,
                &nparts_,
                tpwgts,
                ubvec,
                opts,
                &edgecut,
                curSection.elempartition.data(),
                &mpi_comm
            );

            if(stat!=METIS_OK)
            {
                //fmt::format("Thread", this->rank_, "ParMETIS_V3_PartMeshKway return error while divide mesh:", curSection.name, '\n');
                exit(EXIT_FAILURE);
            }
            else{
                //fmt::format(format("Thread %4d decompose %-15s successfully within %10.3G seconds\n", this->rank_, curSection.name.c_str(), TimeStat()-stime));
            }
        }   

        // close Mesh file.
        MPI_Barrier(MPI_COMM_WORLD);
        // bigFile_->Close();
    }
    
    void ParMetisMeshCutter::rebuild_cgnsFile(const std::string meshName)
    {
        
        if(bodySections_.empty())
        {
            if(this->isMaster_) //fmt::format("body section is empty.\n");
            return;
        } 
        if(bodySections_.size()>1)
        {
            //fmt::format("Thread",this->rank_, "Multiple body-sections are not supported yet.\n");
            return;
        }

        // map part cell/node id to original-mesh id 
        this->get_zone_cell_node(bodySections_);
        for(auto ithPart = 0; ithPart < this->nparts_; ++ithPart) this->get_boundary_cell(std::forward<std::vector<cgsize_t>&&>(bodySections_), ithPart);

        // read coordinate from original mesh  
        std::vector<std::vector<std::vector<double>>> localCoords(this->nparts_);
        this->load_coordinate(meshName, localCoords);

        // open file
        // file name of every parts
        char formatStr[128];            
        sprintf(formatStr, "%%s_%%0%dd.cgns", (int)std::log10(this->nparts_) + 1);
        std::vector<int> fp(this->nparts_, 0);
        for(int i=0; i<this->nparts_;++i)
        {
            cgp_mpi_comm(MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            auto nthPartFileName = new char[128];
            sprintf(nthPartFileName, formatStr, meshName.c_str(), i);
            checkPCGNSCall(cgp_open(nthPartFileName, CG_MODE_WRITE, &fp[i]));
            delete[] nthPartFileName;
            printf("Thread %d open file %d\n", rank_, i);
        }

        // write body section
        MPI_Barrier(MPI_COMM_WORLD);
        std::vector<cgsize_t> cellStartOffset(this->nparts_, 1);
        std::vector<int> B_(this->nparts_), Z_(this->nparts_);
        for (cgsize_t ithPart = 0; ithPart < this->nparts_; ++ithPart)
        {
            // base
            checkCGNSCall(cg_base_write(fp[ithPart], "CGNSBASE", 3, 3, &B_[ithPart]));

            // zone and coordinate
            for(auto bodySection : bodySections_)
            {
                // zone
                cgsize_t zoneInfo[3]{nodeInfo_[bodySection][ithPart].size_all_rank, cellInfo_[bodySection][ithPart].size_all_rank, 0};
                checkCGNSCall(cg_zone_write(
                        fp[ithPart],
                        B_[ithPart],
                        meshName.c_str(),
                        zoneInfo,
                        CGNS_ENUMV(Unstructured),
                        &Z_[ithPart]));
                
                // coordinate
                constexpr char coordname[3][12]{"CoordinateX", "CoordinateY", "CoordinateZ"};
                for (int icoords = 0; icoords < 3; ++icoords)
                {
                    // write
                    int C_;
                    checkPCGNSCall(cgp_coord_write(fp[ithPart], 1, 1, CGNS_ENUMV(RealDouble), coordname[icoords], &C_));
                    cgsize_t rmin = nodeInfo_[bodySection][ithPart].id_start_this_rank + 1;
                    cgsize_t rmax = nodeInfo_[bodySection][ithPart].id_end_this_rank;
                    checkPCGNSCall(cgp_coord_write_data(fp[ithPart], B_[ithPart], Z_[ithPart], C_, &rmin, &rmax, localCoords[ithPart][icoords].data()));
                    printf("Thread %4d write %-15s [  %-5ld - %-5ld  ] into %3ld'th file \n", this->rank_, coordname[icoords], rmin, rmax, ithPart);                
                }
                localCoords[ithPart].clear();
                localCoords[ithPart].shrink_to_fit();
            }

            // body section
            for(auto ithSection : bodySections_)
            {
                auto &curSection = sectionPartition[ithSection];

                auto sectionStart = cellStartOffset[ithPart];
                auto sectionEnd = sectionStart + cellInfo_[ithSection][ithPart].size_all_rank - 1;
                cellStartOffset[ithPart] += cellInfo_[ithSection][ithPart].size_all_rank;

                auto start = cellInfo_[ithSection][ithPart].id_start_this_rank+1;
                start = (start > sectionEnd ? sectionEnd : start); 
                auto end = start + cellInfo_[ithSection][ithPart].size_this_rank - 1;
                end = (end < start ? start : end);

                // mixed cell
                if (curSection.elementType == CGNS_ENUMV(MIXED))
                {
                    std::vector<cgsize_t> element_connectivity;
                    std::vector<cgsize_t> offsets(cellInfo_[ithSection][ithPart].size_this_rank+1, 0);

                    switch (curSection.connectivityDataType)
                    {
                    case CGIO::Section::I4:
                        for(auto ithCell = 0; ithCell< cellInfo_[ithSection][ithPart].size_this_rank; ++ithCell)
                        {
                            auto cellId = cellPartIdMapToGlobalId_[ithSection][ithPart][ithCell] - curSection.start + 1;
                            auto begin = curSection.connect_offset[cellId];
                            auto end = curSection.connect_offset[cellId+1];
                            offsets[ithCell+1] = offsets[ithCell] + end - begin;

                            element_connectivity.push_back(((int*)curSection.element_connectivity)[begin]);
                            for(auto j = begin+1; j<end; ++j)
                            {
                                element_connectivity.push_back(nodeGlobalIdMapToPartId_[ithSection][ithPart][((int*)curSection.element_connectivity)[j]-1] + 1);
                            }
                        }
                        break;
                    case CGIO::Section::I8:
                        for(auto ithCell = 0; ithCell< cellInfo_[ithSection][ithPart].size_this_rank; ++ithCell)
                        {
                            auto cellId = cellPartIdMapToGlobalId_[ithSection][ithPart][ithCell] - curSection.start + 1;
                            auto begin = curSection.connect_offset[cellId];
                            auto end = curSection.connect_offset[cellId+1];
                            offsets[ithCell+1] = offsets[ithCell] + end - begin;

                            element_connectivity.push_back(((cglong_t*)curSection.element_connectivity)[begin]);
                            for(auto j = begin+1; j<end; ++j)
                            {
                                element_connectivity.push_back(nodeGlobalIdMapToPartId_[ithSection][ithPart][((cglong_t*)curSection.element_connectivity)[j]-1] + 1);
                            }
                        }
                        break;
                    default:
                        break;
                    }

                    int S_;
                    std::vector<cgsize_t> offsetSizeAllRank(this->size_, 0);
                    MPI_Allgather(&offsets[cellInfo_[ithSection][ithPart].size_this_rank], 1, GetMPIDataType<cgsize_t>(), offsetSizeAllRank.data(), 1, GetMPIDataType<cgsize_t>(), MPI_COMM_WORLD);
                    auto offsetsTotalSize = std::accumulate(offsetSizeAllRank.begin(), offsetSizeAllRank.end(), decltype(offsetSizeAllRank)::value_type(0));
                    cgp_poly_section_write(fp[ithPart], B_[ithPart], Z_[ithPart], curSection.name.c_str(), CGNS_ENUMV(MIXED), sectionStart, sectionEnd, offsetsTotalSize, 0, &S_);

                    cgsize_t offsetVal = 0;
                    for(auto i=0; i<this->rank_; ++i)
                    {
                        offsetVal += offsetSizeAllRank[i];
                    }
                    for(auto i=0; i<cellInfo_[ithSection][ithPart].size_this_rank+1; ++i)
                    {
                        offsets[i] += offsetVal;
                    }
                    cgp_poly_elements_write_data(fp[ithPart], B_[ithPart], Z_[ithPart], S_, start, end, element_connectivity.data(), offsets.data());
                    printf("Thread %4d write %3ld'th file of section: %-15s [  %-5ld - %-5ld,  %-5ld - %-5ld  ]\n", this->rank_, ithPart, curSection.name.c_str(), start, end, sectionStart, sectionEnd);
                }
                // fixed size cell
                else
                {
                    int S_;
                    checkPCGNSCall(cgp_section_write(fp[ithPart], B_[ithPart], Z_[ithPart], curSection.name.c_str(), curSection.elementType, sectionStart, sectionEnd, 0, &S_));

                    int nodeSizePerCell = 0;
                    checkCGNSCall(cg_npe((ElementType_t)curSection.elementType, &nodeSizePerCell));

                    std::vector<cgsize_t> element_connectivity(nodeSizePerCell*cellInfo_[ithSection][ithPart].size_this_rank);
                    cgsize_t id = 0;

                    switch (curSection.connectivityDataType)
                    {
                    case CGIO::Section::I4:
                        for(cgsize_t ithCell=0; ithCell<cellInfo_[ithSection][ithPart].size_this_rank;++ithCell)
                        {
                            auto cellId = cellPartIdMapToGlobalId_[ithSection][ithPart][ithCell] - curSection.start + 1;
                            for(auto j = nodeSizePerCell*cellId; j<nodeSizePerCell*(cellId+1); ++j)
                            {
                                element_connectivity[id++] = (nodeGlobalIdMapToPartId_[ithSection][ithPart][((int*)curSection.element_connectivity)[j]-1] + 1);
                            }
                        }
                        break;
                    case CGIO::Section::I8:
                        for(cgsize_t ithCell=0; ithCell<cellInfo_[ithSection][ithPart].size_this_rank;++ithCell)
                        {
                            auto cellId = cellPartIdMapToGlobalId_[ithSection][ithPart][ithCell] - curSection.start + 1;
                            for(auto j = nodeSizePerCell*cellId; j<nodeSizePerCell*(cellId+1); ++j)
                            {
                                element_connectivity[id++] = (nodeGlobalIdMapToPartId_[ithSection][ithPart][((cglong_t*)curSection.element_connectivity)[j]-1] + 1);
                            }
                        }
                        break;
                    default:
                        break;
                    }

                    checkPCGNSCall(cgp_elements_write_data(fp[ithPart], B_[ithPart], Z_[ithPart], S_, start, end, element_connectivity.data()));
                    printf("Thread %4d write %3ld'th file of section: %-15s [  %-5ld - %-5ld,  %-5ld - %-5ld  ]\n", this->rank_, ithPart, curSection.name.c_str(), start, end, sectionStart, sectionEnd);
                }
            }
        }
        
        // write boundary section
        MPI_Barrier(MPI_COMM_WORLD);
        std::vector<std::set<std::set<cgsize_t>>> outerFaces;
        for (auto i = 0; i < nparts_; ++i)
        {
            outerFaces.emplace_back(this->get_outside_faceCell(0, i));
        }
        for (auto ithSection : boundarySections_)
        {
            auto &curSection = sectionPartition[ithSection];

            for(auto ithPart=0; ithPart<this->nparts_; ++ithPart)
            {
                auto &outerFace = outerFaces[ithPart];

                std::vector<cgsize_t> connectivity;
                std::vector<cgsize_t> offset; offset.push_back(0);
                
                auto nodeG2L = nodeGlobalIdMapToPartId_[0][ithPart];
                if(curSection.elementType==ElementType_t::MIXED)
                {
                    // loop all cell of this part 
                    for(auto ithCell = 0; ithCell < curSection.cell_size_this_rank; ++ithCell)
                    {
                        std::set<cgsize_t> curFace;
                        auto begin = curSection.connect_offset[ithCell];
                        auto end = curSection.connect_offset[ithCell + 1];

                        for (auto ithNode = begin + 1; ithNode < end; ++ithNode)
                        {
                            cgsize_t idG = ((int *)curSection.element_connectivity)[ithNode] - 1;
                            curFace.insert(idG);
                        }
                        if(outerFace.find(curFace)!=outerFace.end())
                        {
                            // offset
                            offset.push_back(curFace.size()+offset.back()+1);
                            // cell connectivity
                            connectivity.push_back(((int *)curSection.element_connectivity)[begin]);
                            // for(auto id : curFace) connectivity.push_back(nodeG2L[id]+1);

                            for (auto ithNode = begin+1; ithNode < end; ++ithNode)
                            {
                                cgsize_t idG = ((int *)curSection.element_connectivity)[ithNode]-1;
                                connectivity.push_back(nodeG2L[idG]+1);
                            }
                            outerFace.erase(curFace);
                        }
                    }
                }
                else
                {
                    // loop all cell of this part 
                    for(auto ithCell = 0; ithCell < curSection.cell_size_this_rank; ++ithCell)
                    {
                        std::set<cgsize_t> curFace;
                        auto begin = curSection.connect_offset[ithCell];
                        auto end = curSection.connect_offset[ithCell + 1];

                        for (auto ithNode = begin; ithNode < end; ++ithNode)
                        {
                            cgsize_t idG = ((int *)curSection.element_connectivity)[ithNode] - 1;
                            curFace.insert(idG);
                        }

                        if(outerFace.find(curFace)!=outerFace.end())
                        {
                            // offset
                            offset.push_back(curFace.size()+offset.back());

                            // cell connectivity
                            // for(auto id : curFace) connectivity.push_back(nodeG2L[id]+1);
                            for (auto ithNode = begin; ithNode < end; ++ithNode)
                            {
                                cgsize_t idG = ((int *)curSection.element_connectivity)[ithNode]-1;
                                connectivity.push_back(nodeG2L[idG]+1);
                            }
                            outerFace.erase(curFace);
                        }
                    }
                } // end if cell type

                // write 
                cellInfo_[ithSection][ithPart].size_this_rank = std::max(0UL, offset.size()-1);
                this->get_boundary_cell({ithSection}, ithPart);
                auto cellSizeAllRank = cellInfo_[ithSection][ithPart].size_all_rank;
                if(cellSizeAllRank==0) continue;

                auto sectionStart = cellStartOffset[ithPart];
                auto sectionEnd = sectionStart + cellSizeAllRank - 1;
                cellStartOffset[ithPart] += cellSizeAllRank;

                auto start = cellInfo_[ithSection][ithPart].id_start_this_rank+1;
                start = (start > sectionEnd ? sectionEnd : start); 
                auto end = start + cellInfo_[ithSection][ithPart].size_this_rank - 1;
                end = (end < start ? start : end);

                if(curSection.elementType==ElementType_t::MIXED)
                {
                    int S_;
                    std::vector<cgsize_t> offsetSizeAllRank(this->size_, 0);
                    MPI_Allgather(&offset[cellInfo_[ithSection][ithPart].size_this_rank], 1, GetMPIDataType<cgsize_t>(), offsetSizeAllRank.data(), 1, GetMPIDataType<cgsize_t>(), MPI_COMM_WORLD);
                    auto offsetsTotalSize = std::accumulate(offsetSizeAllRank.begin(), offsetSizeAllRank.end(), decltype(offsetSizeAllRank)::value_type(0));
                    cgp_poly_section_write(fp[ithPart], B_[ithPart], Z_[ithPart], curSection.name.c_str(), CGNS_ENUMV(MIXED), sectionStart, sectionEnd, offsetsTotalSize, 0, &S_);

                    cgsize_t offsetVal = 0;
                    for(auto i=0; i<this->rank_; ++i)
                    {
                        offsetVal += offsetSizeAllRank[i];
                    }
                    for(auto i=0; i<cellInfo_[ithSection][ithPart].size_this_rank+1; ++i)
                    {
                        offset[i] += offsetVal;
                    }
                    cgp_poly_elements_write_data(fp[ithPart], B_[ithPart], Z_[ithPart], S_, start, end, connectivity.data(), offset.data());
                    printf("Thread %4d write %3ld'th file of section: %-15s [  %-5ld - %-5ld,  %-5ld - %-5ld  ]\n", this->rank_, ithPart, curSection.name.c_str(), start, end, sectionStart, sectionEnd);
                }
                else
                {
                    int S_;
                    checkPCGNSCall(cgp_section_write(fp[ithPart], B_[ithPart], Z_[ithPart], curSection.name.c_str(), curSection.elementType, sectionStart, sectionEnd, 0, &S_));
                    checkPCGNSCall(cgp_elements_write_data(fp[ithPart], B_[ithPart], Z_[ithPart], S_, start, end, connectivity.data()));
                    printf("Thread %4d write %3ld'th file of section: %-15s [  %-5ld - %-5ld,  %-5ld - %-5ld  ]\n", this->rank_, ithPart, curSection.name.c_str(), start, end, sectionStart, sectionEnd);
                } // end if cell type
            } // end loop part
        }// end loop section

        // close file
        MPI_Barrier(MPI_COMM_WORLD);
        for(int i=0; i<this->nparts_;++i)
        {
            checkPCGNSCall(cgp_close(fp[i]));
            //fmt::format("{} close file: {}\n", this->rank_, fp[i]);
        }


        // find_outside_faceCell();
        // generate_boundary_section();
        // generate_interface();

        /*
        // write interface-section
        MPI_Barrier(MPI_COMM_WORLD);
        auto curPart = rank_;
        while (curPart<this->nparts_)
        {
            int fp;
            auto nthPartFileName = new char[128];
            sprintf(nthPartFileName, formatStr, meshName.lessExtension().c_str(), curPart);
            checkCGNSCall(cg_open(nthPartFileName, CG_MODE_MODIFY, &fp));

            // map global id to thread-section-part-local id
            for (auto nbrPart = 0; nbrPart < this->nparts_; nbrPart++)
            {
                if(nbrPart==curPart) continue;
                
                auto sectionStart = cellStartOffset[curPart];
                auto sizeThisPart = 0;
                for(auto it : interface_[curPart][nbrPart]) {sizeThisPart += it.size();}
                if(sizeThisPart<1) continue;
                auto sectionEnd = sectionStart + sizeThisPart - 1; 
                cellStartOffset[curPart] += sizeThisPart;

                std::vector<cgsize_t> element_connectivity;
                std::vector<cgsize_t> offsets(sizeThisPart+1, 0);
                
                auto id = 1;
                for(auto &ithThread : interface_[curPart][nbrPart])
                {
                    for(auto &ithCell : ithThread)
                    {
                        offsets[id] = offsets[id-1] + ithCell.nodeList.size() + 1; 
                        element_connectivity.push_back(CellType(ithCell.nodeList.size(), 2)); //TODO 
                        for(auto it : ithCell.nodeList)
                        {
                            element_connectivity.push_back(it+1);
                        }
                        ++id;
                    }
                }

                int S_;
                char name[32];
                sprintf(name, "Interface_%ld", nbrPart);
                checkCGNSCall(cg_poly_section_write(fp, B_[curPart], Z_[curPart], name, CGNS_ENUMV(MIXED), sectionStart, sectionEnd, 0, element_connectivity.data(), offsets.data(), &S_));

                printf("Thread %4d write %3ld'th file of section: %-15s [  %-5ld - %-5ld  ]\n", this->rank_, curPart, name, sectionStart, sectionEnd);
            }
            // next part to process 
            curPart += this->size_;

            checkCGNSCall(cg_close(fp));
        }
        

        //fmt::format(this->rank_, "clear memory and finish buildMesh\n");
        */
       
    }

    void ParMetisMeshCutter::charArrayBitOR(char *in, char *inout, int *len, MPI_Datatype *dptr)
    {
        int n = this->nparts_ / 8 + 1;
        for (auto i = 0; i < (*len); ++i)
        {
            for (auto j = 0; j < n; ++j)
            {
                inout[i * n + j] |= in[i * n + j];
            }
        }
    }
    
    void ParMetisMeshCutter::get_zone_cell_node(const std::vector<cgsize_t>& sectionIDs)
    {
        auto const n = sectionIDs.size();
        cellPartIdMapToGlobalId_.assign(n, std::vector<std::unordered_map<cgsize_t, cgsize_t>>(this->nparts_));
        nodeGlobalIdMapToPartId_.assign(n, std::vector<std::map<cgsize_t, cgsize_t>>(this->nparts_));
        nodePartIdMapToGlobalId_.assign(n, std::vector<std::map<cgsize_t, cgsize_t>>(this->nparts_));

        // std::vector<std::vector<BitArray>> nodeGlobalIDFlag(n, std::vector<BitArray>(this->totalNode_, BitArray(this->nparts_)));
        std::vector<std::vector<std::vector<cgsize_t>>> nodeGlobalIdFlag(n, vvcgt(this->nparts_,vcgt(this->totalNode_,0)));
        std::vector<std::vector<std::vector<cgsize_t>>> nodeGlobalIdMerge(n, vvcgt(this->nparts_, vcgt(this->totalNode_,this->size_))); 

        cellInfo_.assign(sectionPartition.size(), std::vector<CellInfo>(this->nparts_, CellInfo()));

        // auto fp = std::fstream(format("thread-%d-node-part.log", rank_), std::ios_base::out);

        for(auto ithSection : sectionIDs)
        {
            auto& curSection = sectionPartition[ithSection];
            if(curSection.elementType == CGNS_ENUMV(ElementType_t::MIXED))
            {
                switch (curSection.connectivityDataType)
                {
                case CGIO::Section::I4:
                    for (auto i = 0; i < curSection.cell_size_this_rank; ++i)
                    {
                        auto ithPart = curSection.elempartition[i];
                        // fp << "Cell " << i  << "    " << ithPart << "    " << i + curSection.start << "\t";
                        cellPartIdMapToGlobalId_[ithSection][ithPart].insert({cellInfo_[ithSection][ithPart].size_this_rank++, i + curSection.start - 1});
                        for (auto j = curSection.connect_offset[i] + 1; j < curSection.connect_offset[i + 1]; ++j)
                        {
                            // fp << ((int*)curSection.element_connectivity)[j] << "\t";
                            nodeGlobalIdFlag[ithSection][ithPart][((int*)curSection.element_connectivity)[j]-1] = 1;
                        }
                        // fp << '\n';
                    }
                    break;
                case CGIO::Section::I8:
                    for (auto i = 0; i < curSection.cell_size_this_rank; ++i)
                    {
                        auto ithPart = curSection.elempartition[i];
                        cellPartIdMapToGlobalId_[ithSection][ithPart].insert({cellInfo_[ithSection][ithPart].size_this_rank++, i + curSection.start - 1});
                        for (auto j = curSection.connect_offset[i] + 1; j < curSection.connect_offset[i + 1]; ++j)
                        {
                            nodeGlobalIdFlag[ithSection][ithPart][((cglong_t*)curSection.element_connectivity)[j]-1] = 1;
                        }
                    }
                    break;
                default:
                    break;
                }
            }
            else
            {
                int nodeSizePerCell = 0;
                checkCGNSCall(cg_npe((ElementType_t)curSection.elementType, &nodeSizePerCell));
                switch (curSection.connectivityDataType)
                {
                case CGIO::Section::I4:
                    for(auto i=0; i<curSection.cell_size_this_rank; ++i)
                    {
                        auto ithPart = curSection.elempartition[i];
                        cellPartIdMapToGlobalId_[ithSection][ithPart].insert({cellInfo_[ithSection][ithPart].size_this_rank++, i + curSection.start - 1});
                        for(auto j=i*nodeSizePerCell; j<nodeSizePerCell*(i+1); ++j)
                        {
                            nodeGlobalIdFlag[ithSection][ithPart][((int*)curSection.element_connectivity)[j]-1] = 1;
                        }
                    }
                    break;
                case CGIO::Section::I8:
                    for(auto i=0; i<curSection.cell_size_this_rank; ++i)
                    {
                        auto ithPart = curSection.elempartition[i];
                        cellPartIdMapToGlobalId_[ithSection][ithPart].insert({cellInfo_[ithSection][ithPart].size_this_rank++, i + curSection.start - 1});
                        for(auto j=i*nodeSizePerCell; j<nodeSizePerCell*(i+1); ++j)
                        {
                            nodeGlobalIdFlag[ithSection][ithPart][((cglong_t*)curSection.element_connectivity)[j]-1] = 1;
                        }
                    }
                    break;
                default:
                    break;
                }
            }
        }
   
        for (auto ithSection = 0; ithSection < bodySections_.size(); ++ithSection)
        {
            for(auto ithPart = 0; ithPart < this->nparts_; ++ithPart)
            {
                MPI_Allreduce(nodeGlobalIdFlag[ithSection][ithPart].data(), nodeGlobalIdMerge[ithSection][ithPart].data(), this->totalNode_, GetMPIDataType<cgsize_t>(), MPI_BOR, MPI_COMM_WORLD);

                nodeInfo_.emplace_back(std::vector<NodeInfo>(this->nparts_, NodeInfo()));
                for (cgsize_t ithGlobalNodeId = 0; ithGlobalNodeId < this->totalNode_; ++ithGlobalNodeId)
                {
                    if (nodeGlobalIdMerge[ithSection][ithPart][ithGlobalNodeId]==1)
                    {
                        nodeGlobalIdMapToPartId_[ithSection][ithPart].insert({ithGlobalNodeId, nodeInfo_[ithSection][ithPart].size_all_rank});
                        nodePartIdMapToGlobalId_[ithSection][ithPart].insert({nodeInfo_[ithSection][ithPart].size_all_rank++, ithGlobalNodeId});
                    }
                }
            }

            for(auto ithPart=0; ithPart < this->nparts_; ++ithPart)
            {
                // start/end of part id
                auto tmp = nodeInfo_[ithSection][ithPart].size_all_rank / this->size_;
                nodeInfo_[ithSection][ithPart].id_start_this_rank = this->rank_ * tmp;
                nodeInfo_[ithSection][ithPart].id_end_this_rank = (this->rank_ == this->size_ - 1) ? nodeInfo_[ithSection][ithPart].size_all_rank : (this->rank_ + 1) * tmp;
                nodeInfo_[ithSection][ithPart].size_this_rank = nodeInfo_[ithSection][ithPart].id_end_this_rank - nodeInfo_[ithSection][ithPart].id_start_this_rank;

                // start end of original mesh global id
                auto iter = nodeGlobalIdMapToPartId_[ithSection][ithPart].begin();
                for (int i = 0; i < nodeInfo_[ithSection][ithPart].id_start_this_rank; ++i)
                {
                    ++iter;
                }
                nodeInfo_[ithSection][ithPart].id_start_global = iter->first;
                for (int i = 0; i < nodeInfo_[ithSection][ithPart].size_this_rank - 1; ++i)
                {
                    ++iter;
                }
                nodeInfo_[ithSection][ithPart].id_end_global = iter->first;
            }
        }
    }

    void ParMetisMeshCutter::get_boundary_cell(std::vector<cgsize_t>&& sectionIDs, const cgsize_t ithPart)
    {
        for(auto ithSection : sectionIDs)
        {
            // for (int ithPart = 0; ithPart < this->nparts_; ++ithPart)
            // {
                // size all rank
                MPI_Allreduce(&cellInfo_[ithSection][ithPart].size_this_rank, &cellInfo_[ithSection][ithPart].size_all_rank, 1, GetMPIDataType<cgsize_t>(), MPI_SUM, MPI_COMM_WORLD);

                // id start
                std::vector<cgsize_t> cellSizeAllRank(this->size_, 0);
                MPI_Allgather(&cellInfo_[ithSection][ithPart].size_this_rank, 1, GetMPIDataType<cgsize_t>(), cellSizeAllRank.data(), 1, GetMPIDataType<cgsize_t>(), MPI_COMM_WORLD);
                for (int iProc = 0; iProc < this->rank_; iProc++)
                {
                    cellInfo_[ithSection][ithPart].id_start_this_rank += cellSizeAllRank[iProc];
                }
                cgsize_t cellIdOffset = 0;
                for (auto i = 0; i < ithSection; ++i)
                {
                    cellIdOffset += (cellInfo_[i][ithPart].size_all_rank);
                }
                cellInfo_[ithSection][ithPart].id_start_this_rank += cellIdOffset;
                // id end
                cellInfo_[ithSection][ithPart].id_end_this_rank = cellInfo_[ithSection][ithPart].id_start_this_rank + cellInfo_[ithSection][ithPart].size_this_rank;
                printf("Thread %4d has %15s's %3d'th part with cell range [  %-10ld - %-10ld  ] and size(this,all) [  %-10ld, %-10ld  ]\n", this->rank_, sectionPartition[ithSection].name.c_str(), ithPart, cellInfo_[ithSection][ithPart].id_start_this_rank, cellInfo_[ithSection][ithPart].id_end_this_rank, cellInfo_[ithSection][ithPart].size_this_rank, cellInfo_[ithSection][ithPart].size_all_rank);
            // }
        }
    }

    std::set<std::set<cgsize_t>> ParMetisMeshCutter::get_outside_faceCell(const cgsize_t ithSection, const cgsize_t ithPart)
    {
        auto & curSection = sectionPartition[ithSection];

        // record all face of current part
        std::set<std::set<cgsize_t>> hash;

        if (curSection.elementType == CGNS_ENUMV(MIXED))
        {
            for (auto ithCell = 0; ithCell < cellInfo_[ithSection][ithPart].size_this_rank; ++ithCell)
            {
                auto cellId = cellPartIdMapToGlobalId_[ithSection][ithPart][ithCell] - curSection.start + 1;
                auto begin = curSection.connect_offset[cellId];
                auto end = curSection.connect_offset[cellId + 1];
                auto faces = nodeIdToBuildFaceInCell.at(CellType(((int*)curSection.element_connectivity)[begin]));

                for (auto face : faces)
                {
                    std::set<cgsize_t> tmp;
                    std::for_each(face.begin(), face.end(), [&](const int node)
                                  { tmp.insert(((int*)curSection.element_connectivity)[begin + 1 + node] - 1); });
                    if (hash.find(tmp) == hash.end())
                    {
                        hash.insert(tmp);
                    }
                    else
                    {
                        hash.erase(tmp);
                    }
                }
            }
        }
        else // if-elementType
        {
            int nodeSizePerCell = 0;
            checkCGNSCall(cg_npe((ElementType_t)curSection.elementType, &nodeSizePerCell));
            auto faces = nodeIdToBuildFaceInCell.at(curSection.elementType);

            for (auto ithCell = 0; ithCell < cellInfo_[ithSection][ithPart].size_this_rank; ++ithCell)
            {
                auto cellId = cellPartIdMapToGlobalId_[ithSection][ithPart][ithCell] - curSection.start + 1;
                auto begin = nodeSizePerCell * cellId;

                for (auto face : faces)
                {
                    std::set<cgsize_t> tmp;
                    std::for_each(face.begin(), face.end(), [&](const int node)
                                  { tmp.insert(((int*)curSection.element_connectivity)[begin + node] - 1); });
                    if (hash.find(tmp) == hash.end())
                    {
                        hash.insert(tmp);
                    }
                    else
                    {
                        hash.erase(tmp);
                    }
                }
            }
        } // end-if-elementType

        //
        printf("Thread %d has outside face-cell cnt: %d\n", rank_, hash.size());

        return hash;
    }

    // TODO: for multiply zone section, localCoords needs to be 4-dimensiton vector.
    void ParMetisMeshCutter::load_coordinate(const FileName& meshName, std::vector<std::vector<std::vector<double>>>& localCoords)
    {
        int fpMesh; 
        cgp_mpi_comm(MPI_COMM_WORLD);
        checkPCGNSCall(cgp_open(meshName.c_str(), CG_MODE_READ, &fpMesh), "Thread", this->rank_, "Wrong to open CGNS file:", meshName.c_str(),'\n');
        for(auto bodySection : bodySections_)
        {
            for (cgsize_t ithPart = 0; ithPart < this->nparts_; ++ithPart)
            {
                localCoords[ithPart].assign(3, {});
                for (int icoords = 0; icoords < 3; ++icoords)
                {
                    auto rmin = nodeInfo_[bodySection][ithPart].id_start_global + 1;
                    auto rmax = nodeInfo_[bodySection][ithPart].id_end_global + 1;
                    std::vector<double> globalCoords(rmax - rmin + 1, 0.0);
                    localCoords[ithPart][icoords].assign(nodeInfo_[bodySection][ithPart].size_this_rank, 0.0);
                    cgp_coord_read_data(fpMesh, 1, 1, icoords + 1, &rmin, &rmax, globalCoords.data());
                    for (cgsize_t id = 0; id < nodeInfo_[bodySection][ithPart].size_this_rank; ++id)
                    {
                        auto partId = id + nodeInfo_[bodySection][ithPart].id_start_this_rank;
                        auto globalId = nodePartIdMapToGlobalId_[bodySection][ithPart][partId] - nodeInfo_[bodySection][ithPart].id_start_global;
                        localCoords[ithPart][icoords][id] = globalCoords[globalId];
                    }
                }
            }
        }
        checkPCGNSCall(cgp_close(fpMesh));
        MPI_Barrier(MPI_COMM_WORLD);
    }

   // TODO: for multiply zone section, this function is not correct.
    void ParMetisMeshCutter::find_interface_perThread()
    {
        // different part maybe has interface cell 
        for(auto ithPart=0; ithPart < this->nparts_; ++ithPart)
        {
            // record all face of current part
            // std::map<std::set<cgsize_t>, cgsize_t> hash;
            std::set<std::set<cgsize_t>> hash;
            cgsize_t index = 0;

            // find non-internal face, include interface-face and boundary-face 
            for (auto &curSection : sectionPartition)
            {
                // if(!curSection.isBodySection) continue;

                auto cellCnt = cellInfo_[curSection.id-1][ithPart].size_this_rank;

                if (curSection.elementType == CGNS_ENUMV(MIXED))
                {
                    for (auto ithCell = 0; ithCell < cellCnt; ++ithCell)
                    {
                        auto cellId = cellPartIdMapToGlobalId_[curSection.id-1][ithPart][ithCell] - curSection.start + 1;
                        auto begin = curSection.connect_offset[cellId];
                        auto end = curSection.connect_offset[cellId+1];
                        // auto faces = nodeIdToBuildFaceInCell.at(CellType(curSection.element_connectivity[begin]));

                        // for (auto face : faces)
                        // {
                        //     std::set<cgsize_t> tmp;
                        //     std::for_each(face.begin(), face.end(), [&](const int node)
                        //     {
                        //         tmp.insert(curSection.element_connectivity[begin + 1 + node] - 1);
                        //     });
                        //     if(hash.find(tmp) == hash.end())
                        //     {
                        //         hash.insert(tmp);
                        //     }
                        //     // if (hash.find(tmp) == hash.end())
                        //     // {
                        //     //     hash[tmp] = index++;
                        //     // }
                        //     // else
                        //     // {
                        //     //     auto n = hash.erase(tmp);
                        //     // }
                        // }
                    }
                }
                else // if-elementType
                {
                    int nodeSizePerCell = 0;
                    checkCGNSCall(cg_npe((ElementType_t)curSection.elementType, &nodeSizePerCell));
                    auto faces = nodeIdToBuildFaceInCell.at(curSection.elementType);

                    for (auto ithCell = 0; ithCell < cellCnt; ++ithCell)
                    {
                        auto cellId = cellPartIdMapToGlobalId_[curSection.id-1][ithPart][ithCell] - curSection.start + 1;
                        auto begin = nodeSizePerCell * cellId;

                        for (auto face : faces)
                        {
                            // std::set<cgsize_t> tmp;
                            // std::for_each(face.begin(), face.end(), [&](const int node)
                            // {
                            //     tmp.insert(curSection.element_connectivity[begin + node] - 1);
                            // });
                            // if(hash.find(tmp) == hash.end())
                            // {
                            //     hash.insert(tmp);
                            // }
                            // // if (hash.find(tmp) == hash.end())
                            // // {
                            // //     hash[tmp] = index++;
                            // // }
                            // // else
                            // // {
                            // //     auto n = hash.erase(tmp);
                            // // }
                        }
                    }
                } // end-if-elementType
            
            } // end-for-sectionPartition

            // printf("Thread %d find non-interface cell of part %d, cnt: %d\n", rank_, ithPart, hash.size());

            
            // remove boundary face
            for(auto &curSection : sectionPartition)
            {
                // if(curSection.isBodySection) continue;

                auto cellCnt = cellInfo_[curSection.id-1][ithPart].size_this_rank;

                if (curSection.elementType == CGNS_ENUMV(MIXED))
                {
                    for (auto ithCell = 0; ithCell < cellCnt; ++ithCell)
                    {
                        auto cellId = cellPartIdMapToGlobalId_[curSection.id-1][ithPart][ithCell] - curSection.start + 1;
                        auto begin = curSection.connect_offset[cellId];
                        auto end = curSection.connect_offset[cellId + 1];

                        // std::set<cgsize_t> tmp;
                        // for(auto i=0; i< NodeCountOfFaceCell(CellType(curSection.element_connectivity[begin])); ++i)
                        // {
                        //     tmp.insert(curSection.element_connectivity[begin + 1 + i] - 1);
                        // }
                        // if (hash.find(tmp) != hash.end())
                        // {
                        //     auto n = hash.erase(tmp);
                        // }
                    }
                }
                else // if-elementType
                {
                    int nodeSizePerCell = 0;
                    checkCGNSCall(cg_npe((ElementType_t)curSection.elementType, &nodeSizePerCell));

                    for (auto ithCell = 0; ithCell < cellCnt; ++ithCell)
                    {
                        auto cellId = cellPartIdMapToGlobalId_[curSection.id-1][ithPart][ithCell] - curSection.start + 1;
                        auto begin = nodeSizePerCell * cellId;

                        // std::set<cgsize_t> tmp;
                        // for(auto i=0; i< nodeSizePerCell; ++i)
                        // {
                        //     tmp.insert(curSection.element_connectivity[begin + i] - 1);
                        // }
                        // if (hash.find(tmp) != hash.end())
                        // {
                        //     auto n = hash.erase(tmp);
                        // }  
                    }
                } // end-if-elementType
        
            } // end-for-sectionPartition


            // printf("Thread %d remove boundary cell of part %d, cnt: %d\n", rank_, ithPart, hash.size());

            //TODO: only support single body-section.
            auto bodySectinId = 0;
            for(auto it : sectionPartition)
            {
                // if(it.isBodySection) 
                // {
                //     bodySectinId = it.id-1;
                //     break;
                // }
            }
            //
            faceInfo_.push_back(FaceInfo());
            faceInfo_.back().size_this_rank = hash.size();
            // faceInfo_.back().node_id_offset.reserve(hash.size());
            index = 0;
            for (auto &it : hash)
            {
                Face_ face;
                // face.id = index++;
                face.nodeList = it;
                // face.type = CellType(it.size(), 2);
                outsideFaceCell_[ithPart].insert(face);
                faceInfo_.back().node_size_this_rank += it.size();
                for(auto n : it)
                {
                    outsideNode_[ithPart].insert({n, nodeGlobalIdMapToPartId_[bodySectinId][ithPart][n]}); //TODO: only support single body-section.
                }
            }
        }

        return;
    }

    std::ostream &operator<<(std::ostream &os, std::set<cgsize_t> val)
    {
        // os << "[";
        for (auto it : val)
        {
            os << it << ',';
        }
        os << "\n";
        return os;
    }

    void ParMetisMeshCutter::map_interface_allThread()
    {

        auto totalNodeOnInterfaceAllPart = std::accumulate(faceInfo_.begin(), faceInfo_.end(), 0, [](int curSum, const FaceInfo& f){
            return curSum + f.node_size_this_rank;
        });
        auto totalInterFaceAllPart = std::accumulate(faceInfo_.begin(), faceInfo_.end(), 0, [](int curSum, const FaceInfo& f){
            return curSum + f.size_this_rank;
        });

        // header
        //   offset of first part, equal this->nparts_
        //   offset of sectond part, equal length of part 0
        //   ...
        //   offset of this->nparts part
        // data
        //   face_number_of_part_0
        //   cell_0: node_count, node1, node2, ... nodeN;
        //   cell_1: node_count, node1, node2, ... nodeN;
        //   ...
        //   face_number_of_part_n
        //   cell_0: node_count, node1, node2, ... nodeN;
        //   cell_1: node_count, node1, node2, ... nodeN;
        //   ...
        // total len = header + data(part_0 + part_1 + part_n)
        int lenLocal = 2*this->nparts_+totalNodeOnInterfaceAllPart+totalInterFaceAllPart;
        cgsize_t* sendBuf = new cgsize_t[lenLocal];

        // setup send-buffer
        cgsize_t curPart = 0;
        cgsize_t curPos = this->nparts_;
        for(auto faces : outsideFaceCell_) // loop all part
        {
            sendBuf[curPart] = curPos;
            sendBuf[curPos++] = faces.size();
            std::for_each(faces.begin(), faces.end(), [&](const Face_& face)
            {
                sendBuf[curPos++] = face.nodeList.size();
                for(auto node : face.nodeList)
                {
                    sendBuf[curPos++] = node;
                }
            });
            ++curPart;
        }

        // mpi all-gather-v communicate
        // all thread get all interface data
        int *lenAll = new int[this->size_];
        MPI_Allgather(&lenLocal, 1, GetMPIDataType<int>(), lenAll, 1, GetMPIDataType<int>(), MPI_COMM_WORLD);      
        cgsize_t lenTotal = std::accumulate(lenAll, lenAll+this->size_, 0);
        cgsize_t *recvBuf = new cgsize_t[lenTotal];
        int *disp = new int[this->size_];
        disp[0] = 0;
        for(int i=1; i<this->size_; ++i)
        {
            disp[i] = disp[i-1] + lenAll[i-1];
        }
        MPI_Allgatherv(sendBuf, lenLocal, GetMPIDataType<cgsize_t>(), recvBuf, lenAll, disp, GetMPIDataType<cgsize_t>(), MPI_COMM_WORLD);
        
        // free 
        delete[] lenAll;
        delete[] sendBuf;
        sendBuf = nullptr;
        lenAll = nullptr;


        // each thread process specific parts
        // thread 0: part_0, part_(0+1*(rankSize_-1)), part_(0+2*(rankSize_-1)) ... 
        // thread 1: part_1, part_(1+1*(rankSize_-1)), part_(1+2*(rankSize_-1)) ...
        // ...
        curPart = rank_;
        while (curPart<this->nparts_)
        {
            // collect all curPart-interface cell from all thread
            std::unordered_set<Face_> hash;
            for (auto ithThread = 0; ithThread < this->size_; ++ithThread)
            {
                auto threadOffset = disp[ithThread];
                auto partOffset = threadOffset + recvBuf[threadOffset + curPart];
                auto faceCnt = recvBuf[partOffset];
                auto curPos = partOffset + 1;
                for (int ithFace = 0; ithFace < faceCnt; ++ithFace)
                {
                    Face_ face;
                    auto nNode = recvBuf[curPos++];
                    for (int i = 0; i < nNode; ++i)
                    {
                        face.nodeList.insert(recvBuf[curPos++]);
                    }
                    if (hash.find(face) == hash.end())
                    {
                        hash.insert(face);
                    }
                }
            }

            // find common cell
            for (auto nbrPart = 0; nbrPart < this->nparts_; nbrPart++)
            {
                if(nbrPart==curPart) continue;

                cgsize_t index = 0;
                for(auto ithThread = 0; ithThread<this->size_; ++ithThread)
                {
                    auto threadOffset = disp[ithThread];
                    auto partOffset = threadOffset+recvBuf[threadOffset+nbrPart];
                    auto faceCnt = recvBuf[partOffset];
                    auto curPos = partOffset+1;
                    for(int ithFace = 0; ithFace<faceCnt; ++ithFace)
                    {
                        Face_ face;
                        auto nNode = recvBuf[curPos++];
                        for(int i=0; i<nNode; ++i)
                        {
                            face.nodeList.insert(recvBuf[curPos++]);
                        }
                        if(hash.find(face) != hash.end())
                        {
                            // face.nbrPart = nbrPart;
                            // face.id = ithThread;
                            interface_[curPart][nbrPart][ithThread].push_back(face);
                        }
                    }
                }
            }
            // next part to process 
            curPart += this->size_;
        }

        // free 
        delete[] recvBuf;
        delete[] disp; 
        recvBuf = nullptr;
        disp = nullptr;

        this->map_node_id();

        // for(auto ipart = 0; ipart<this->nparts_; ++ipart)
        // {
        //     for(auto nbrPart = 0; nbrPart < this->nparts_; ++nbrPart)
        //     {
        //         for(auto nThread = 0; nThread < this->size_; ++nThread)
        //         {
        //             std::fstream fp;
        //             if(interface_[ipart][nbrPart][nThread].size()>0) fp.open(("thread-" + std::to_string(rank_) + "-part-" + std::to_string(ipart) + "-nbrPart-" + std::to_string(nbrPart)+ "-nThread-" + std::to_string(nThread) + ".log").c_str(), std::ios_base::out);
        //             for(auto ithface : interface_[ipart][nbrPart][nThread])
        //             {
        //                 fp <<  ithface.nbrPart << ",\t" << ithface.id << ",\t" << ithface.nodeList;
        //             }
        //             fp.close();
        //         }
        //     }
        // }
    }

    void ParMetisMeshCutter::map_node_id()
    {
        // map global node id to local id 
        auto nodeCntAllPart = std::accumulate(outsideNode_.begin(), outsideNode_.end(), 0, [](int curSum, const std::unordered_map<cgsize_t, cgsize_t>& m){
            return curSum + m.size();
        });

        // header
        //   offset of first part, equal this->nparts_
        //   offset of sectond part, equal length of part 0
        //   ...
        //   offset of this->nparts part
        // data
        //   node_number_part_0
        //   node_0_global_id, node_0_local_id;
        //   node_1_global_id, node_1_local_id;
        //   ...
        //   node_number_part_n
        //   node_0_global_id, node_n_local_id;
        //   node_1_global_id, node_1_local_id;
        //   ... 
        // total len = header + data(part_0 + part_1 + part_n)
        int lenLocal = 2*this->nparts_+nodeCntAllPart*2;
        cgsize_t* sendBuf = new cgsize_t[lenLocal];

        // setup send-buffer
        cgsize_t curPart = 0;
        cgsize_t curPos = this->nparts_;
        for(auto nodePairs : outsideNode_) // loop all part
        {
            sendBuf[curPart] = curPos;
            sendBuf[curPos++] = nodePairs.size();
            for(auto nodePair : nodePairs)
            {
                sendBuf[curPos++] = nodePair.first;
                sendBuf[curPos++] = nodePair.second;
            }
            ++curPart;
        }


        // mpi all-gather-v communicate
        // all thread get all node data
        int *lenAll = new int[this->size_];
        MPI_Allgather(&lenLocal, 1, GetMPIDataType<int>(), lenAll, 1, GetMPIDataType<int>(), MPI_COMM_WORLD);      
        cgsize_t lenTotal = std::accumulate(lenAll, lenAll+this->size_, 0);
        cgsize_t *recvBuf = new cgsize_t[lenTotal];
        int *disp = new int[this->size_];
        disp[0] = 0;
        for(int i=1; i<this->size_; ++i)
        {
            disp[i] = disp[i-1] + lenAll[i-1];
        }
        MPI_Allgatherv(sendBuf, lenLocal, GetMPIDataType<cgsize_t>(), recvBuf, lenAll, disp, GetMPIDataType<cgsize_t>(), MPI_COMM_WORLD);

        // free 
        delete[] lenAll;
        delete[] sendBuf;
        sendBuf = nullptr;
        lenAll = nullptr;


        // each thread process specific parts
        // thread 0: part_0, part_(0+1*(rankSize_-1)), part_(0+2*(rankSize_-1)) ... 
        // thread 1: part_1, part_(1+1*(rankSize_-1)), part_(1+2*(rankSize_-1)) ...
        // ...
        curPart = rank_;
        while (curPart<this->nparts_)
        {
            std::unordered_map<cgsize_t, cgsize_t> idMap;

            for(auto ithThread = 0; ithThread<this->size_; ++ithThread)
            {
                auto threadOffset = disp[ithThread];
                auto partOffset = threadOffset+recvBuf[threadOffset+curPart];
                auto nodeCnt = recvBuf[partOffset];
                auto curPos = partOffset+1;
                for(int ithNode = 0; ithNode<nodeCnt; ++ithNode)
                {
                    idMap.insert({recvBuf[curPos], recvBuf[curPos+1]});
                    curPos+=2;
                }
            }

            // map global id to thread-section-part-local id
            for (auto nbrPart = 0; nbrPart < this->nparts_; nbrPart++)
            {
                if(nbrPart==curPart) continue;

                for(auto ithThread = 0; ithThread<this->size_; ++ithThread)
                {
                    for(auto &ithFace : interface_[curPart][nbrPart][ithThread])
                    {
                        std::set<cgsize_t> nodes;
                        for(auto it : ithFace.nodeList)
                        {
                            nodes.insert(idMap[it]);
                        }
                        ithFace.nodeList.swap(nodes);
                    }
                }
            }
            // next part to process 
            curPart += this->size_;
        }

        // free 
        delete[] recvBuf;
        delete[] disp; 
        recvBuf = nullptr;
        disp = nullptr;
    }

    bool isBoundaryElementType(const ElementType_t& t, const int Dim)
    {
        bool rst = false;
        switch (Dim)
        {
        case 2:
            switch(t)
            {
            // 1-D
            case ElementType_t::BAR_2 :
            case ElementType_t::BAR_3 :
            case ElementType_t::BAR_4 :
            case ElementType_t::BAR_5 :
                rst = true;
                break;
            default:
                break;
            }
            break;
        case 3:
            switch(t)
            {
            // 2-D
            case ElementType_t::TRI_3:
            case ElementType_t::TRI_6:
            case ElementType_t::TRI_9:
            case ElementType_t::TRI_10:
            case ElementType_t::TRI_12:
            case ElementType_t::TRI_15:
            case ElementType_t::QUAD_4:
            case ElementType_t::QUAD_8:
            case ElementType_t::QUAD_9:
            case ElementType_t::QUAD_12:
            case ElementType_t::QUAD_16:
            case ElementType_t::QUAD_P4_16:
            case ElementType_t::QUAD_25:
                rst = true;
                break;
            default:
                break;
            }
            break;
        }
        return rst;
    };

    bool isSupportedElementType(const ElementType_t& t)
    {
        bool rst = false;
        switch (t)
        {
        case ElementType_t::TRI_3 :
        case ElementType_t::QUAD_4 :
        case ElementType_t::TETRA_4 :
        case ElementType_t::HEXA_8 :
        case ElementType_t::PENTA_6 :
        case ElementType_t::PYRA_5 :
            rst = true;
            break;
        default:
            rst = false;
            break;
        }
        if(!rst)
        {
            //fmt::format("Unsupported element type:", t, '\n');
        }
        return rst;
    };

    ::std::ostream& operator<<(::std::ostream& o, ElementType_t t)
    {
        switch (t)
        {
        case ElementTypeNull:
            return o << "ElementTypeNull";
        case ElementTypeUserDefined:
            return o << "ElementTypeUserDefined";
        case NODE:
            return o << "NODE";
        case BAR_2:
            return o << "BAR_2";
        case BAR_3:
            return o << "BAR_3";
        case TRI_3:
            return o << "TRI_3";
        case TRI_6:
            return o << "TRI_6";
        case QUAD_4:
            return o << "QUAD_4";
        case QUAD_8:
            return o << "QUAD_8";
        case QUAD_9:
            return o << "QUAD_9";
        case TETRA_4:
            return o << "TETRA_4";
        case TETRA_10:
            return o << "TETRA_10";
        case PYRA_5:
            return o << "PYRA_5";
        case PYRA_14:
            return o << "PYRA_14";
        case PENTA_6:
            return o << "PENTA_6";
        case PENTA_15:
            return o << "PENTA_15";
        case PENTA_18:
            return o << "PENTA_18";
        case HEXA_8:
            return o << "HEXA_8";
        case HEXA_20:
            return o << "HEXA_20";
        case HEXA_27:
            return o << "HEXA_27";
        case MIXED:
            return o << "MIXED";
        case PYRA_13:
            return o << "PYRA_13";
        case NGON_n:
            return o << "NGON_n";
        case NFACE_n:
            return o << "NFACE_n";
        case BAR_4:
            return o << "BAR_4";
        case TRI_9:
            return o << "TRI_9";
        case TRI_10:
            return o << "TRI_10";
        case QUAD_12:
            return o << "QUAD_12";
        case QUAD_16:
            return o << "QUAD_16";
        case TETRA_16:
            return o << "TETRA_16";
        case TETRA_20:
            return o << "TETRA_20";
        case PYRA_21:
            return o << "PYRA_21";
        case PYRA_29:
            return o << "PYRA_29";
        case PYRA_30:
            return o << "PYRA_30";
        case PENTA_24:
            return o << "PENTA_24";
        case PENTA_38:
            return o << "PENTA_38";
        case PENTA_40:
            return o << "PENTA_40";
        case HEXA_32:
            return o << "HEXA_32";
        case HEXA_56:
            return o << "HEXA_56";
        case HEXA_64:
            return o << "HEXA_64";
        case BAR_5:
            return o << "BAR_5";
        case TRI_12:
            return o << "TRI_12";
        case TRI_15:
            return o << "TRI_15";
        case QUAD_P4_16:
            return o << "QUAD_P4_16";
        case QUAD_25:
            return o << "QUAD_25";
        case TETRA_22:
            return o << "TETRA_22";
        case TETRA_34:
            return o << "TETRA_34";
        case TETRA_35:
            return o << "TETRA_35";
        case PYRA_P4_29:
            return o << "PYRA_P4_29";
        case PYRA_50:
            return o << "PYRA_50";
        case PYRA_55:
            return o << "PYRA_55";
        case PENTA_33:
            return o << "PENTA_33";
        case PENTA_66:
            return o << "PENTA_66";
        case PENTA_75:
            return o << "PENTA_75";
        case HEXA_44:
            return o << "HEXA_44";
        case HEXA_98:
            return o << "HEXA_98";
        case HEXA_125:
            return o << "HEXA_125";

        default:
            return o << t;
        }
    }

} // namespace MeshCut