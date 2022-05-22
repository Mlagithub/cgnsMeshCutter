#include <algorithm>
#include <numeric>
#include <cstring>
#include <functional>
#include <sstream>

#include "mpi.h"
#include "cgnslib.h"
#include "pcgnslib.h"

#include "cartesianCutter.h"
#include "format.h"
#include "stringUtil.h"


#if defined(__GUNC__) || defined(__GNUG__)
#include <features.h>
#if __GNUC_PREREQ(5,2)
#else
namespace std
{
template<>
struct hash<DataType_t>
{
    typedef DataType_t argument_type;
    typedef size_t result_type;

    result_type operator()(const argument_type& x) const
    {
        using type = typename std::underlying_type<argument_type>::type;
        return std::hash<type>()(static_cast<type>(x));
    }
};
}
#endif
#endif


namespace MeshCut
{

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
            std::cout << format("not supported element type with[dim, nodes] [ %d, %d]\n", dim, nodeCnt);
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
            std::cout << format("not supported element type with[dim, nodes] [ %d, %d]\n", dim, nodeCnt);
            break;
        }
        break;
    default:
        break;
    }

    return ElementTypeNull;
}
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
        std::cout << format("unknown element type %d\n", static_cast<ElementType_t>(CGNSCellTypeFalg));
        break;
    }
    return rst;
}
static unordered_map<DataType_t, std::function<void(void *&, const size_t len)>> New{
    {CGNS_ENUMV(RealSingle), [](void *&ptr, const size_t len) { ptr = new float[len]; }},
    {CGNS_ENUMV(RealDouble), [](void *&ptr, const size_t len) { ptr = new double[len]; }}};

static std::unordered_map<DataType_t, std::function<void(void *&)>> Delete{
    {CGNS_ENUMV(RealSingle), [](void *&ptr) { float* tmp = ((float*)ptr); delete []tmp; ptr=nullptr; }},
    {CGNS_ENUMV(RealDouble), [](void *&ptr) { double* tmp = ((double*)ptr); delete []tmp; ptr=nullptr; }}};

static std::unordered_map<DataType_t, std::function<double(void *&, const size_t id)>> Get{
    {CGNS_ENUMV(RealSingle), [](void *&ptr, const size_t id) { return ((float *)ptr)[id]; }},
    {CGNS_ENUMV(RealDouble), [](void *&ptr, const size_t id) { return ((double *)ptr)[id]; }}};

// 各种单元类型每个面的节点在单元内部的编号。
// 由于面的法向量具有方向性，这些数据在计算相关向量时需要使用。
const unordered_map<ElementType_t, vector<vector<int>>>
    nodeIdToBuildFaceInCell = {
        {ElementType_t::TETRA_4, {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}},
        {ElementType_t::HEXA_8,
         {{0, 3, 2, 1},
          {0, 1, 5, 4},
          {1, 2, 6, 5},
          {2, 3, 7, 6},
          {0, 4, 7, 3},
          {4, 5, 6, 7}}},
        {ElementType_t::PENTA_6,
         {{0, 2, 1}, 
          {3, 4, 1, 0},
          {5, 2, 1, 4},
          {2, 5, 3, 0},
          {5, 4, 3}}},
        {ElementType_t::PYRA_5,
         {{0, 3, 2, 1},
          {0, 1, 4},
          {1, 2, 4},
          {2, 3, 4},
          {3, 0, 4}}}};


CartesianCutter::CellLoader::CellLoader(Section& curS, vector<cgsize_t>& data, vector<cgsize_t>& offset, const int fp)
    : curS(curS), data(data), offset(offset)
{
    cg_section_read(fp, 1, 1, curS.id, curS.name, &curS.cellType, &curS.start, &curS.end, &curS.nBdy, &curS.flag);
    cg_ElementDataSize(fp, 1, 1, curS.id, &curS.dataSize);
    len = curS.end - curS.start + 1;
    id = 0;

    if(curS.cellType == CGNS_ENUMV(MIXED))
    {
        data.assign(curS.dataSize, 0);
        offset.assign(len+1, 0);
        cg_poly_elements_read(fp, 1, 1, curS.id, data.data(), offset.data(), NULL);
    }
    else
    {
        cg_npe(curS.cellType, &npe);
        data.assign(curS.dataSize, 0);
        cg_elements_read(fp, 1, 1, curS.id, data.data(), NULL);
        flag = curS.cellType;
    }
}

bool CartesianCutter::CellLoader::nextCell(vector<cgsize_t>& nodes)
{
    if(id>=len) return false;

    cgsize_t start, end;
    if(curS.cellType == CGNS_ENUMV(MIXED))
    {
        start = offset[id] + 1;
        end = offset[id+1];
        flag = data[start-1];
    }
    else
    {
        start = npe * id;
        end = start + npe;
    }

    for(auto j=start; j<end; ++j) { nodes.push_back(data[j]); }

    ++id;

    return true;
}

CartesianCutter::Face::Face(const vector<cgsize_t>& val)
    : nodeList(val) {}

CartesianCutter::Face::Face(vector<cgsize_t>&& val) noexcept
    : nodeList(std::move(val)){}

CartesianCutter::Face::Face(const Face& other)
    : nodeList(other.nodeList) {}

CartesianCutter::Face::Face(Face&& other) noexcept
    : nodeList(std::move(other.nodeList)) {}

CartesianCutter::Face& CartesianCutter::Face::operator=(const Face& other)
{
    if(other == *this) return *this;
    nodeList = other.nodeList;

    return *this;
}

CartesianCutter::Face& CartesianCutter::Face::operator=(Face&& other) noexcept
{
    nodeList = std::move(other.nodeList);
    
    return *this;
}

bool CartesianCutter::Face::operator==(const Face &other) const
{
    if (this->nodeList.empty() || other.nodeList.empty())
        return false;
    if(nodeList.size() == other.nodeList.size() ) return nodeList == other.nodeList;

    auto nums1 = nodeList, nums2 = other.nodeList;
    std::sort(nums1.begin(), nums1.end());
    std::sort(nums2.begin(), nums2.end());
    auto length1 = nums1.size(), length2 = nums2.size();
    size_t index1 = 0, index2 = 0;
    int common = 0;
    while (index1 < length1 && index2 < length2) {
        auto num1 = nums1[index1], num2 = nums2[index2];
        if (num1 == num2) {
            ++common;
            ++index1;
            ++index2;
        } else if (num1 < num2) {
            ++index1;
        } else {
            ++index2;
        }
    }
    return common>=3 ? true : false;
}

bool CartesianCutter::Face::operator<(const Face& rhs) const
{
    std::set<cgsize_t> lhsList{this->nodeList.begin(), this->nodeList.end()};
    std::set<cgsize_t> rhsList{rhs.nodeList.begin(), rhs.nodeList.end()};\
    return lhsList < rhsList;
}

vector<vector<cgsize_t>> CartesianCutter::allFaceInCell(const vector<cgsize_t>& idList, ElementType_t type)
{
    check(
        (nodeIdToBuildFaceInCell.count(type)!=0),
        format("not supported cell type %d\n", type));

    vector<vector<cgsize_t>> rst;

    for(auto nodes : nodeIdToBuildFaceInCell.at(type))
    {
        vector<cgsize_t> face;
        for(auto it : nodes) face.push_back(idList[it]);
        rst.emplace_back(face);
    }

    return rst;
}

vector<double> CartesianCutter::boundingBoxBigMesh()
{
    char coordname[33];
    int ncoords;
    vector<double> bbox(6,0);
    cgsize_t range_min = 1, range_max = zoneInfo[0];    

    if(!cg_grid_bounding_box_read(bigFile_, 1, 1, 1, coordDataType_, bbox.data()))
    {
        cg_ncoords(bigFile_, ibase, izone, &ncoords);
        check(cg_coord_info(bigFile_, 1, 1, 1, &coordDataType_, coordname)); 

        if(coordDataType_ == CGNS_ENUMV(RealDouble))
        {
            vector<double> coordsDouble(range_max, 0.0);
            for(int icoord = 1; icoord<=ncoords; ++icoord)
            {
                cg_coord_info(bigFile_, ibase, izone, icoord, &coordDataType_, coordname);
                cg_coord_read(bigFile_, ibase, izone, coordname, coordDataType_, &range_min, &range_max, coordsDouble.data());

                auto rst = std::minmax_element(coordsDouble.begin(), coordsDouble.end());
                bbox[(icoord-1) * 2 + 0] = *rst.first;
                bbox[(icoord-1) * 2 + 1] = *rst.second;

                coordname_.push_back(coordname);
            }
        }  
        else
        {
            vector<float> coordsSingle(range_max, 0.0);
            for(int icoord = 1; icoord<=ncoords; ++icoord)
            {
                cg_coord_info(bigFile_, ibase, izone, icoord, &coordDataType_, coordname);
                cg_coord_read(bigFile_, ibase, izone, coordname, coordDataType_, &range_min, &range_max, coordsSingle.data());
                auto rst = std::minmax_element(coordsSingle.begin(), coordsSingle.end());
                bbox[(icoord-1) * 2 + 0] = *rst.first;
                bbox[(icoord-1) * 2 + 1] = *rst.second;

                coordname_.push_back(coordname);
            }
        }
    }

    return bbox;
}

void CartesianCutter::boundingBoxSmallMesh()
{
    auto bbox = this->boundingBoxBigMesh();
    auto stepI = (bbox[1] - bbox[0]) / nx_;
    auto stepJ = (bbox[3] - bbox[2]) / ny_;
    auto stepK = (bbox[5] - bbox[4]) / nz_;

    for (auto rank = 0; rank < n_; ++rank)
    {
        auto &box = bbox_[rank];
        auto loc = id2Location(rank);
        box.min_x = bbox[0] + stepI * loc.i;
        box.max_x = bbox[0] + stepI * (loc.i+1);
        box.min_y = bbox[2] + stepJ * loc.j;
        box.max_y = bbox[2] + stepJ * (loc.j+1);
        box.min_z = bbox[4] + stepK * loc.k;
        box.max_z = bbox[4] + stepK * (loc.k+1);
    }
}

void CartesianCutter::check(bool state, string&& errStr, string&& okStr)
{
    if(!state)
    {
        if(!errStr.empty())
        {
            std::cerr << (errStr.c_str());
        }
        exit(EXIT_FAILURE);
    }
    else
    {
        if(!okStr.empty()) std::cout << (okStr.c_str());
    }
}

void CartesianCutter::check(const int runCode, string&& errStr, string&& okStr)
{
    if(runCode != CG_OK)
    {
        if(!errStr.empty())
        {
            std::cerr << "ERROR: " << (errStr.c_str());
        }
        exit(EXIT_FAILURE);
    }
    else
    {
        if(!okStr.empty()) std::cout << (okStr.c_str());
    }
}

void CartesianCutter::checkFile(string fileName)
{   
    char basename[33], zonename[33];
    int nbases, nzones, ngrids, ncoords, phys_dim;
    int filetype;
    ZoneType_t zoneType;
    int &fp = bigFile_;
    check(
        cg_is_cgns(fileName.c_str(), &filetype),
        format("%s is not a CGNS file\n", fileName.c_str()));

    check(
        cg_open(fileName.c_str(), CG_MODE_READ, &fp), 
        format("can not open file: %s \n", fileName.c_str()));

    check(
        (cg_nbases(fp, &nbases), nbases==1), 
        format("now only support one Base, you have %d\n", nbases));

    check(
        (cg_base_read(fp, ibase, basename, &cell_dim_, &phys_dim), phys_dim==3),
        format("now only support 3-Dim mesh, you have %d\n", phys_dim));

    check(
        (cg_nzones(fp, ibase, &nzones), nzones==1),
        format("now only support one Zone, you have %d\n", nzones));
    
    check(
        (cg_zone_type(fp, ibase, izone, &zoneType), zoneType == ZoneType_t::Unstructured),
        format("now only support UnStructured Mesh\n"));

    check(cg_zone_read(fp, ibase, izone, zonename, zoneInfo));

    check(
        (cg_ngrids(fp, ibase, izone, &ngrids), ngrids==1),
        format("now only support single layer Mesh, you have %d\n", ngrids));

    int nSection;
    Section section;
    check(
        (cg_nsections(fp, ibase, izone, &nSection), nSection>0),
        format("no section found in file %s", fileName.c_str()));
    for(int iSection = 1; iSection<=nSection; ++iSection)
    {
        section.id = iSection;
        cg_section_read(fp, ibase, izone, iSection, section.name, &section.cellType, &section.start, &section.end, &section.nBdy, &section.flag);
        cg_ElementDataSize(fp, ibase, izone, iSection, &section.dataSize);
        if(string{section.name}.find(interiorFacePrefix_) != string::npos) interiorSection_.push_back(iSection);
        else if(this->isBodySection(section)) bodySection_.push_back(iSection);
        else bdySection_.push_back(iSection);

        sections_.insert({iSection, section});
        globalNumber_.insert({iSection, section.end-section.start+1});
    }

    coordDataType_ = (this->precision() == 32) ? CGNS_ENUMV(RealSingle) : CGNS_ENUMV(RealDouble);
}

Vector<double> CartesianCutter::cellCenter(const vector<cgsize_t> &cell, const ElementType_t type)
{
    auto nodeCoordinate = [&](cgsize_t id) {
        --id;
        return Vector<double>{nodesData_[0][id], nodesData_[1][id], nodesData_[2][id]};
    };

    auto faceCenter = [&](const vector<cgsize_t> & nodelist)
    {
        if(nodelist.size() == 3)
        {
            return (nodeCoordinate(nodelist[0]) + nodeCoordinate(nodelist[1]) + nodeCoordinate(nodelist[2])) / 3.0;
        }
        else if(nodelist.size()==4)
        {
            vector<Vector<double>> nodes;
            for(auto it : nodelist) nodes.push_back(nodeCoordinate(it));
            Scalar v[2];
            auto a = nodes[1] - nodes[0];
            auto b = nodes[2] - nodes[0];
            v[0] = (a * b).norm() * 0.5;
            a = nodes[2] - nodes[0];
            b = nodes[3] - nodes[0];
            v[1] = (a * b).norm() * 0.5;
            auto vs = v[0] + v[1];
            auto center = (v[0] / vs) * (nodes[0] + nodes[1] + nodes[2]) / 3.0;
            center = center + (v[1] / vs) * (nodes[0] + nodes[2] + nodes[3]) / 3.0;
            return center;
        }
    };

    auto nodeLists = [&](const int id)
    {
        vector<cgsize_t> nList;
        for(auto it : nodeIdToBuildFaceInCell.at(type)[id]) nList.push_back(cell[it]);
        return nList;
    };

    Vector<double> center;
    vector<Vector<double>> nodeList;
    for(auto it : cell)
    { 
        nodeList.push_back(nodeCoordinate(it));
    }

    switch (type)
    {
    case ElementType_t::TRI_3:
    {
        center = (nodeList[0] + nodeList[1] + nodeList[2]) / 3.;
    }
    break;
    case ElementType_t::QUAD_4:
    {
        // nodeList.push_back(*nodeList.begin());
        // for(int i=0; i<4; ++i)
        // {
        //     auto t1 = nodeList[i] + nodeList[i+1];
        //     auto t2 = nodeList[i] * nodeList[i+1];
        //     center.x() += t1.x() *  t2.x();
        //     center.y() += t1.y() *  t2.y();
        //     center.z() += t1.z() *  t2.z();
        // }
        // auto area = ((nodeList[2]-nodeList[0]) * (nodeList[3]-nodeList[1])).norm();
        // center = center / (6 * area);
        Scalar v[2];
        ScalarVector a = nodeList[1] - nodeList[0];
        ScalarVector b = nodeList[2] - nodeList[0];
        v[0] = (a * b).norm() * 0.5;
        a = nodeList[2] - nodeList[0];
        b = nodeList[3] - nodeList[0];
        v[1] = (a * b).norm() * 0.5;
        Scalar vs = v[0] + v[1];
        center = (v[0] / vs) * (nodeList[0] + nodeList[1] + nodeList[2]) / 3.;
        center = center + (v[1] / vs) * (nodeList[0] + nodeList[2] + nodeList[3]) / 3.;
    }
    break;
    case ElementType_t::TETRA_4 :
    {
        center = (nodeList[0] + nodeList[1] + nodeList[2] + nodeList[3]) / 4;
    }
    break;
    case ElementType_t::PYRA_5:
    {
        auto faceCenter0 = faceCenter(nodeLists(0));
        center = (nodeList[4] - faceCenter0) / 3 + faceCenter0;
    }
    break;
    case ElementType_t::PENTA_6:
    {
        center = 0.5 * (faceCenter(nodeLists(0)) + faceCenter(nodeLists(4)));
    }
    break;
    case ElementType_t::HEXA_8:
    {
        // 将六面体分成五个四面体，分别求出中心和体积；
        // 六面体的中心为各个四面体的中心按体积加权平均。
        ScalarVector a, b, c;
        Scalar v[5];
        a = nodeList[0] - nodeList[4];
        b = nodeList[5] - nodeList[4];
        c = nodeList[7] - nodeList[4];
        v[0] = fabs(a & (b * c)) / 6;
        a = nodeList[2] - nodeList[6];
        b = nodeList[5] - nodeList[6];
        c = nodeList[7] - nodeList[6];
        v[1] = fabs(a & (b * c)) / 6;
        a = nodeList[0] - nodeList[1];
        b = nodeList[2] - nodeList[1];
        c = nodeList[5] - nodeList[1];
        v[2] = fabs(a & (b * c)) / 6;
        a = nodeList[0] - nodeList[3];
        b = nodeList[2] - nodeList[3];
        c = nodeList[7] - nodeList[3];
        v[3] = fabs(a & (b * c)) / 6;
        a = nodeList[2] - nodeList[0];
        b = nodeList[5] - nodeList[0];
        c = nodeList[7] - nodeList[0];
        v[4] = fabs(a & (b * c)) / 6;
        auto volume = v[0] + v[1] + v[2] + v[3] + v[4];
        center = 0.25 * (v[0] / volume) *
                  (nodeList[0] + nodeList[4] +
                   nodeList[5] + nodeList[7]);
        center += 0.25 * (v[1] / volume) *
                   (nodeList[2] + nodeList[5] +
                    nodeList[6] + nodeList[7]);
        center += 0.25 * (v[2] / volume) *
                   (nodeList[0] + nodeList[1] +
                    nodeList[2] + nodeList[5]);
        center += 0.25 * (v[3] / volume) *
                   (nodeList[0] + nodeList[2] +
                    nodeList[3] + nodeList[7]);
        center += 0.25 * (v[4] / volume) *
                   (nodeList[0] + nodeList[2] +
                    nodeList[5] + nodeList[7]);
    }
    break;
    default:
        break;
    }

    return center;
}

void CartesianCutter::cut(const std::string &mesh, const int npart, const int nx, const int ny, const int nz)
{
    string bigFileName = mesh;
    interiorFacePrefix_ = "INTERFACE";

    // check partition parameter
    n_ = npart, nx_ = nx, ny_ = ny, nz_ = nz;
    {
        check(n_>=2, format("partition number should >=2.\n"));
        check((n_ == nx_ * ny_ * nz_), format("partition and nproc mismatch, check your config file\n"));
    }
    fluidDomain_ = "FLUID";
    check(!fluidDomain_.empty(), format("setup fluid domain name at YAML\n"));

    this->checkFile(bigFileName);

    // calculate bbox
    this->boundingBoxSmallMesh();

    this->readNode();
    
    // handle all file
    for(auto i=0; i<sections_.size(); ++i) globalOffset_.insert({sections_[i].id, 0});
    idOffset_.assign(n_, 1);
    for (auto iFile = 0; iFile < n_; ++iFile)
    {
        set<int> sectionInBox;
        auto smallFilename = this->openToWrite(bigFileName, iFile);
        auto &curBox = bbox_[iFile];
        nodeIdInBox_.clear();
        auto &curNodeIdMap = nodeIdMap_[iFile];

        // pick body section
        for (auto iBody : bodySection_)
        {
            auto &curS = sections_[iBody];
            curS.data.clear();
            curS.offset.clear();
            curS.typeFlag.clear();
            
            if (!this->pickBodyCellInBox(curBox, curS)) continue;
            else std::cout << format("  %d: pick section %s \n", iFile, curS.name);
            
            sectionInBox.insert(iBody);
            updateNodeId(curS);
        }

        // write coordiante
        std::cout << format("  %d: write coordinate %d\n", iFile, nodeIdInBox_.size());
        for (auto it : nodeIdInBox_)
        {
            curNodeIdMap.insert({it, curNodeIdMap.size() + 1});
        }
        zoneInfo[0] = curNodeIdMap.size();
        zoneInfo[1] = std::accumulate(sectionInBox.begin(), sectionInBox.end(), 0, [&](cgsize_t sum, const int ithSection){
            return sum + sections_[ithSection].data.size();});
        writeNode(iFile);

        // body section
        for(auto it : sectionInBox)
        {
            auto &curS = sections_[it];
            writeSection(curS, iFile);
            writeGlobalInfo(curS, iFile);
            updateOuterFace(curS, iFile);
            curS.data.clear();
            curS.offset.clear();
        }

        // bdy section
        this->rwBoundary(iFile);

        // interior section
        this->rwInterior(iFile);

        // interface
        this->rwInterface(iFile);


        if(outerFace_[iFile].empty()) nodeIdMap_[iFile].clear();
    }

    check(cg_close(bigFile_), "error close bigFile\n");
    for(auto it : smallFiles_) check(cg_close(it.second), "error close smallFile\n");
}

bool CartesianCutter::isBodySection(const Section& section)
{   
    vector<cgsize_t> bodyConn, offset;
    int cell_dim;
    switch (section.cellType)
    {
    case NODE:
        cell_dim = 0;
        break;
    case BAR_2:
        cell_dim = 1;
        break;
    case TRI_3:
    case QUAD_4:
        cell_dim = 2;
        break;
    case TETRA_4:
    case PYRA_5:
    case PENTA_6:
    case HEXA_8:
        cell_dim = 3;
        break;
    case MIXED:
        cgsize_t flag;
        bodyConn.assign(section.dataSize, 0);
        offset.assign(section.end-section.start+1, 0);
        cg_poly_elements_read(bigFile_, ibase, izone, section.id, bodyConn.data(), offset.data(), &flag);
        switch(bodyConn[0])
        {
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
            cell_dim = 2;
            break;
        default:
            cell_dim = 3;
            break;
        }
        break;
    default:
        std::cerr << ("Unsupported element type\n");
        break;
    }
    return cell_dim == cell_dim_;
}

bool CartesianCutter::inBox(const BoundBox &box, const vector<cgsize_t> &cell, const ElementType_t type) 
{
    auto center = this->cellCenter(cell, type);
    auto rst = (center.x() >= box.min_x && center.x() <= box.max_x &&
            center.y() >= box.min_y && center.y() <= box.max_y &&
            center.z() >= box.min_z && center.z() <= box.max_z );
    return rst;
}

CartesianCutter::Location CartesianCutter::id2Location(const int x) const
{
    Location loc;
    loc.k = x / (nx_ * ny_);
    loc.i = (x - loc.k*(nx_*ny_)) / ny_;
    loc.j = x % ny_;

    return loc;
}

int CartesianCutter::location2Id(const CartesianCutter::Location& loc) const
{
    return loc.k*nx_*ny_ + loc.i*ny_ + loc.j;
}

int CartesianCutter::location2Id(CartesianCutter::Location&& loc) const
{
    return loc.k*nx_*ny_ + loc.i*ny_ + loc.j;
}

vector<int> CartesianCutter::neighbors(const int x)
{
    std::vector<int> rst;

    auto loc = id2Location(x);
    for(auto i = loc.i-1; i<=loc.i+1; ++i)
    {
        if(i<0 || i>=nx_) continue;
        for(auto j=loc.j-1; j<=loc.j+1; ++j)
        {
            if(j<0 || j>=ny_) continue;
            for(auto k=loc.k-1; k<=loc.k+1; ++k)
            {
                if(k<0 || k>=nz_) continue;
                if(i==loc.i && j==loc.j && k==loc.k) continue;
                rst.push_back(location2Id(Location{i, j, k}));
            }
        }
    }

    std::sort(rst.begin(), rst.end());
    return rst;
}

string CartesianCutter::openToWrite(string bigFileName, const int id)
{
    char fmt[128], smallFilename[128];
    sprintf(fmt, "%%s_%%0%dd.cgns", int(std::log10(n_) + 1));
    sprintf(smallFilename, fmt, bigFileName.substr(0, bigFileName.size() - 5).c_str(), id);

    check(
        cg_open(smallFilename, CG_MODE_WRITE, &smallFiles_[id]), 
        format("can not open file %s to write\n", smallFilename),
        format("open file %s to write\n", smallFilename));

    return string{smallFilename};
}

std::ostream &operator<<(std::ostream &os, const CartesianCutter::BoundBox& bbox)
{
    os << "BoundingBox: \n";
    os << "RangeX: [" << bbox.min_x << ", " << bbox.max_x << "]\n";
    os << "RangeY: [" << bbox.min_y << ", " << bbox.max_y << "]\n";
    os << "RangeZ: [" << bbox.min_z << ", " << bbox.max_z << "]\n";
    return os;
}

bool CartesianCutter::pickBodyCellInBox(const BoundBox& box, Section &curSection)
{
    int curOffset = 0;
    vector<cgsize_t> data, offset, tmp;
    auto isMixed = (curSection.cellType == CGNS_ENUMV(MIXED));
    auto loader = CellLoader{curSection, data, offset, bigFile_};
    while (loader.nextCell(tmp))
    {
        if (this->inBox(box, tmp, CellType(loader.flag)))
        {
            curSection.data.push_back(tmp);
            if (isMixed)
            {
                curSection.typeFlag.push_back(loader.flag);
                curSection.offset.push_back(curOffset);
                curOffset += (tmp.size() + 1);
            }
        }
        tmp.clear();
    }
    return !curSection.data.empty();
}

int CartesianCutter::precision()
{
    int rst;
    check(cg_precision(bigFile_, &rst));

    return rst;
}

void CartesianCutter::readNode()
{
    cgsize_t range_min = 1, range_max = zoneInfo[0];
    auto len = range_max - range_min + 1;
    nodesData_.assign(3, vector<double>(len, 0.0));
    void *data;
    New[coordDataType_](data, len);
    for (auto iCoord = 0; iCoord < 3; ++iCoord)
    {
        const char *coordname = coordname_[iCoord].c_str();
        std::cout << format("Read %s [%d - %d]\n", coordname, range_min, range_max);
        cg_coord_read(bigFile_, ibase, izone, coordname, coordDataType_, &range_min, &range_max, data);
        auto &ndata = nodesData_[iCoord];
        for (cgsize_t i = 0; i < len; ++i)
        {
            ndata[i] = Get[coordDataType_](data, i);
        }
    }
    Delete[coordDataType_](data);
}

void CartesianCutter::rwBoundary(const int id)
{
    auto &curOuterFace = outerFace_[id];

    for(auto ithSection : bdySection_)
    {
        auto &curSection = sections_[ithSection];

        int curOffset = 0;
        vector<cgsize_t> data, offset, tmp;
        auto isMixed = (curSection.cellType == CGNS_ENUMV(MIXED));
        auto loader = CellLoader{curSection, data, offset, bigFile_};
        while (loader.nextCell(tmp))
        {
            auto face = Face{tmp};
            if(curOuterFace.find(face) != curOuterFace.end())
            {
                curOuterFace.erase(face);
                curSection.data.emplace_back(tmp);
                if(isMixed)
                {
                    curSection.typeFlag.push_back(loader.flag);
                    curSection.offset.push_back(curOffset);
                    curOffset += (tmp.size() + 1);
                }
            }
            tmp.clear();
        }
        
        if(!curSection.data.empty())
        {
            this->writeSection(curSection, id);
        }

        // clear memory
        vector<vector<cgsize_t>>{}.swap(curSection.data);
    }
}

void CartesianCutter::rwInterior(const int id)
{
    auto curNodeIdMap = nodeIdMap_[id];

    auto isInBox = [&](const vector<cgsize_t>& nodes) -> bool {
        for(auto it : nodes) { if(curNodeIdMap.find(it) == curNodeIdMap.end()) {return false;} }
        return true;
    };

    for(auto ithSection : interiorSection_)
    {
        auto &curSection = sections_[ithSection];

        int curOffset = 0;
        vector<cgsize_t> data, offset, tmp;
        auto isMixed = (curSection.cellType == CGNS_ENUMV(MIXED));
        auto loader = CellLoader{curSection, data, offset, bigFile_};
        while (loader.nextCell(tmp))
        {
            if(isInBox(tmp))
            {
                curSection.data.push_back(tmp);
                if(isMixed)
                {
                    curSection.typeFlag.push_back(loader.flag);
                    curSection.offset.push_back(curOffset);
                    curOffset += (tmp.size() + 1);
                }
            }
            tmp.clear();
        }
        
        if(!curSection.data.empty())
        {
            this->writeSection(curSection, id);
        }

        // clear memory
        vector<vector<cgsize_t>>{}.swap(curSection.data);
    }
}

void CartesianCutter::rwInterface(const int id)
{
    auto &curOuterFace = outerFace_[id];

    for(auto nbr : neighbors(id))
    {
        if (nbr >= id || outerFace_[nbr].empty()) continue;

        auto &nbrOuterFace = outerFace_[nbr];

        set<cgsize_t> typeFlags;
        Section curS, nbrS;
        for(auto face = nbrOuterFace.begin(); face!=nbrOuterFace.end();)
        {
            auto it = curOuterFace.find(face->first);
            if(it != curOuterFace.end())
            {
                curS.data.emplace_back(it->second);
                nbrS.data.push_back(face->second);

                auto t = CellType(it->second.size(), 2);
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
            strcpy(curS.name, format("Interface_%d_%d", fluidDomainId_, nbr).c_str());
            strcpy(nbrS.name, format("Interface_%d_%d", fluidDomainId_, id).c_str());

            if(typeFlags.size()>1){
                curS.cellType = ElementType_t::MIXED;
                nbrS.cellType = ElementType_t::MIXED;
            } 
            else if(typeFlags.size()==1){
                auto type = (*typeFlags.begin()==3) ? ElementType_t::TRI_3 : ElementType_t::QUAD_4;
                curS.cellType = type;
                nbrS.cellType = type;
            }
            writeSection(curS, id);      
            writeSection(nbrS, nbr);   

            // std::fstream fp(format("Interface_%d_%d", rank_, nbr), std::ios_base::out);
            // std::for_each(curS.data.begin(), curS.data.end(), [&](const vector<cgsize_t> &cell) {
            //     auto tmp = this->cellCenter(cell, (cell.size()==3) ? ElementType_t::TRI_3 : ElementType_t::QUAD_4);
            //     fp << tmp << '\n';  
            // });
            // fp.close();
        }
    }
}

void CartesianCutter::updateNodeId(const Section &curS)
{
    for (auto jt : curS.data)
    {
        for (auto it : jt)
        {
            if (nodeIdInBox_.find(it) == nodeIdInBox_.end()) nodeIdInBox_.insert(it);
        }
    }
}

void CartesianCutter::updateOuterFace(const Section& curS, const int id)
{
    auto &curOuterFace = outerFace_[id];

    for (auto i = 0; i < curS.data.size(); ++i)
    {
        auto cellType = curS.cellType == CGNS_ENUMV(MIXED) ? CellType(curS.typeFlag[i]) : curS.cellType;
        for (auto it : allFaceInCell(curS.data[i], cellType))
        {
            auto val = Face{it};
            if (curOuterFace.find(val) != curOuterFace.end()) { curOuterFace.erase(val); }
            else
            {
                curOuterFace.insert({val, it});
            }
        }
    }
}

void CartesianCutter::writeSection(Section &curS, const int id)
{
    auto fp = smallFiles_[id];
    auto &curNodeIdMap = nodeIdMap_[id];

    // update id
    curS.dataSize = std::accumulate(
        curS.data.begin(), 
        curS.data.end(), 
        0, 
        [](cgsize_t sum, vector<cgsize_t> vec){
            return sum + vec.size();
        });
     
    curS.start = idOffset_[id];
    curS.end = curS.start + curS.data.size() - 1;
    idOffset_[id] = curS.end + 1;        

    // write 
    int dummy;
    if(curS.cellType == CGNS_ENUMV(MIXED))
    {
        curS.dataSize += curS.typeFlag.size();
        vector<cgsize_t> data(curS.dataSize, 0);
        int icell=0, i=0;
        for(auto it : curS.data)
        {
            data[i++] = curS.typeFlag[icell++];
            for(auto j : it)
            {
                data[i++] = curNodeIdMap[j];
            }
        }
        check(
            cg_poly_section_write(fp, 1, 1, curS.name, curS.cellType, curS.start, curS.end, 0, data.data(), curS.offset.data(), &dummy),
            format("file %d write section %s\n", id, curS.name));
    }
    else
    {
        vector<cgsize_t> data(curS.dataSize);
        int i=0;
        for(auto it : curS.data)
        {
            for(auto j : it)
            {
                data[i++] = curNodeIdMap[j];
            }
        }
        check(
            cg_section_write(fp, 1, 1, curS.name, curS.cellType, curS.start, curS.end, curS.nBdy, data.data(), &dummy),
            format("file %d write section %s\n", id, curS.name));
    }

    std::cout << format("  %d: write section %s [%d: %d, %d]\n", id, curS.name, curS.end-curS.start+1, curS.start, curS.end);
}

void CartesianCutter::writeNode(const int id)
{
    int dummy;
    auto fp = smallFiles_[id];

    // base
    check(
        cg_base_write(fp, "CGNSBASE", 3, 3, &dummy), 
        "write base\n");

    // zone
    check(
        cg_zone_write(fp, 1, "ZONE", zoneInfo, Unstructured, &dummy),
        "write zone\n");
    
    // grid
    check(
        cg_grid_write(fp, 1, 1, "GridCoordinates", &igrid),
        "write grid\n");

    // coordinate
    std::vector<double> localCoords(zoneInfo[0], 0.0);
    auto &curNodeIdMap = nodeIdMap_[id];
    for(auto iCoord : {0, 1, 2})
    {
        int i=0;
        for(auto it : curNodeIdMap) localCoords[i++] = nodesData_[iCoord][it.first-1];
        const char *coordname = coordname_[iCoord].c_str();
        check(
            cg_coord_write(fp, 1, 1, RealDouble, coordname, localCoords.data(), &dummy),
            format("write coordinate %s\n", coordname));
    }

    // clear memory
    vector<double>{}.swap(localCoords);
}

void CartesianCutter::writeGlobalInfo(const Section& curS, const int id)
{
    auto fp = smallFiles_[id];
    auto len = curS.end - curS.start + 1;
    auto start = globalOffset_[curS.id];
    globalOffset_[curS.id] += len;

    cg_goto(fp, 1, "Zone_t", 1, curS.name, 0, "end");
    std::ostringstream ostr;
    ostr << "GlobalNumber:" << globalNumber_[curS.id] << '\n';
    ostr << "GlobalStart:" << start << '\n';
    ostr << "GlobalEnd:" << start + len << '\n';
    auto str = ostr.str();
    cg_descriptor_write("GlobalDomainInfo", str.c_str());
}


} // namespace MeshCut



    // MPI_Barrier(MPI_COMM_WORLD);
    // {
    //     for(auto it : fluidSection_)
    //     {
    //         auto &curS = sections_[it];
    //         cgsize_t lenSend = curS.end - curS.start + 1, globalOffset = 0;
    //         vector<cgsize_t> cellSizeAllRank(this->size_, 0);
    //         MPI_Allgather(&lenSend, 1, GetMPIDataType<cgsize_t>(), cellSizeAllRank.data(), 1, GetMPIDataType<cgsize_t>(), MPI_COMM_WORLD);
    //         for (int iProc = 0; iProc < this->rank_; iProc++)
    //         {
    //             globalOffset += cellSizeAllRank[iProc];
    //         }

    //         cg_goto(smallFile_, 1, "Zone_t", 1, curS.name, 0, "end");
    //         std::ostringstream ostr;
    //         ostr << "GlobalNumber:" << std::accumulate(cellSizeAllRank.begin(), cellSizeAllRank.end(), 0) << '\n';
    //         ostr << "GlobalStart:" << globalOffset << '\n';
    //         globalOffset += (curS.end-curS.start+1);
    //         ostr << "GlobalEnd:" << globalOffset << '\n';
    //         auto str = ostr.str();
    //         cg_descriptor_write("GlobalDomainInfo", str.c_str());
    //     }

    //     for(auto it : solidSection_)
    //     {
    //         auto &curS = sections_[it];

    //         auto len = curS.end-curS.start+1;
    //         cg_goto(smallFile_, 1, "Zone_t", 1, curS.name, 0, "end");
    //         std::ostringstream ostr;
    //         ostr << "GlobalNumber:" << len << '\n';
    //         ostr << "GlobalStart:" << 0 << '\n';
    //         ostr << "GlobalEnd:" << len << '\n';
    //         auto str = ostr.str();
    //         cg_descriptor_write("GlobalDomainInfo", str.c_str()); 
    //     }
    // }