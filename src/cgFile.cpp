#include <algorithm>
#include <iostream>
#include <numeric>
#include <sstream>
#include <cgns_io.h>
#include <functional>
#include <fstream>
#include <cstring>

#include "format.h"
#include "cgFile.h"
#include "mpiAdapter.h"


static std::unordered_map<DataType_t, std::function<void(void *&, const size_t len)>> New{{DataType_t::Integer, [](void *&ptr, const size_t len) { ptr = new int[len]; }}, {DataType_t::LongInteger, [](void *&ptr, const size_t len) { ptr = new long[len]; }}};

static std::unordered_map<DataType_t, std::function<void(void *&)>> Delete{{DataType_t::Integer,
                                                                          [](void *&ptr) {
                                                                              int *tmp = ((int *)ptr);
                                                                              delete[] tmp;
                                                                              ptr = nullptr;
                                                                          }},
                                                                         {DataType_t::LongInteger, [](void *&ptr) {
                                                                              long *tmp = ((long *)ptr);
                                                                              delete[] tmp;
                                                                              ptr = nullptr;
                                                                          }}};

static std::unordered_map<DataType_t, std::function<long(void *&, const size_t id)>> Get{{DataType_t::Integer, [](void *&ptr, const size_t id) { return ((int *)ptr)[id]; }}, {DataType_t::LongInteger, [](void *&ptr, const size_t id) { return ((long *)ptr)[id]; }}};

static std::unordered_map<DataType_t, std::function<void(void *&, const size_t, void *&, const size_t)>> Copy{{DataType_t::Integer, [](void *&dst, const size_t dstId, void *&src, const size_t srcId) { ((int *)dst)[dstId] = ((int *)src)[srcId]; }},
                                                                                                            {DataType_t::LongInteger, [](void *&dst, const size_t dstId, void *&src, const size_t srcId) { ((long *)dst)[dstId] = ((long *)src)[srcId]; }}};

static std::unordered_map<DataType_t, std::function<long(void *&, const size_t)>> Cast2cgsize_t{{DataType_t::Integer, [](void *&src, const size_t id) -> long { return cgsize_t(((int *)src)[id]); }}, {DataType_t::LongInteger, [](void *&src, const size_t id) -> long { return ((long *)src)[id]; }}};


static cgsize_t BLOCK_LEN_EACH_LOAD = 10000;

namespace MeshCut
{
void check(bool state, string&& errStr, string&& okStr)
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

void check(const int runCode, string&& errStr, string&& okStr)
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

// 各种单元类型每个面的节点在单元内部的编号。
// 由于面的法向量具有方向性，这些数据在计算相关向量时需要使用。
const unordered_map<ElementType_t, vector<vector<int>>> CGFile::nodeIdToBuildFaceInCell
    = {{ElementType_t::TETRA_4, {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}},
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


CGFile::CellLoader::CellLoader(Section& curS, const int fp, const cgsize_t lowerBd, const cgsize_t upperBd)
    : curS(curS), fp(fp)
{
    lowerBound = std::max(lowerBd, curS.start);
    upperBound = std::min(upperBd, curS.end);
    blockSize = std::min(upperBound - lowerBound + 1, BLOCK_LEN_EACH_LOAD);
    start = lowerBound - blockSize;
    loadNextBlock();
}

bool CGFile::CellLoader::nextCell(vector<cgsize_t>& nodes)
{
    if(id>=len && !loadNextBlock()) return false;

    nodes.clear();

    cgsize_t begin, end;
    if(curS.isMixed())
    {
        begin = offset[id] + 1;
        end = offset[id+1];
        flag = data[begin-1];
    }
    else
    {
        begin = npe * id;
        end = begin + npe;
    }

    for(auto j=begin; j<end; ++j) { nodes.push_back(data[j]); }

    ++id;

    return true;
}

bool CGFile::CellLoader::loadNextBlock()
{
    if(start + blockSize>upperBound) return false;

    cgsize_t curBlockDataSize = 0;
    start += blockSize;
    auto end = std::min(start + blockSize - 1, upperBound);
    len = end - start + 1;
    cg_ElementPartialSize(fp, 1, 1, curS.id, start, end, &curBlockDataSize);

    vector<cgsize_t>{}.swap(data);
    if(curS.isMixed())
    {
        vector<cgsize_t>{}.swap(offset);
        data.assign(curBlockDataSize, 0);
        offset.assign(len+1, 0);
        cg_poly_elements_partial_read(fp, 1, 1, curS.id, start, end, data.data(), offset.data(), NULL);
    }
    else
    {
        cg_npe(curS.cellType, &npe);
        data.assign(curBlockDataSize, 0);
        cg_elements_partial_read(fp, 1, 1, curS.id, start, end, data.data(), NULL);
    }

    ++blockId;
    id = 0;

    return true;
}

void CGFile::Section::addCell(const cgsize_t id, const vector<cgsize_t>& idList, const cgsize_t flag)
{
    data.emplace_back(idList);
    IDMap.insert({id, IDMap.size()});
    typeFlag.push_back(flag);
    offset.push_back(offset.back() + idList.size()+1);
}

const vector<cgsize_t>& CGFile::Section::cell(const int id) const
{
    return data[IDMap.at(id)];
}

void CGFile::Section::clear()
{
    vector<vector<cgsize_t>>{}.swap(data);
    vector<cgsize_t>{0}.swap(offset);
    vector<cgsize_t>{}.swap(typeFlag);
    map<cgsize_t, cgsize_t>{}.swap(IDMap);
}

cgsize_t CGFile::Section::flag(const int id) const
{
    return isMixed() ? typeFlag[IDMap.at(id)] : (cgsize_t)cellType;
}

bool CGFile::Section::isMixed() const
{
    return (cellType == CGNS_ENUMV(MIXED));
}

void CGFile::Section::printToFile(string fname)
{
    std::fstream fp(fname, std::ios_base::out);
    if(!fp.is_open()) std::cout << format("Can not open file %s to write\n", fname.c_str());

    fp << format("Cell -> %d:[%d - %d]\n", data.size(), start, end);
    for(auto i=0; i<data.size(); ++i)
    {
        for(auto it : data[i]) fp << it << ',';
        if(isMixed()) fp << format(" => %d,%d\n", offset[i], typeFlag[i]);
    }

    fp.close();
}

CGFile::Section CGFile::Section::subSection(const cgsize_t beg, const cgsize_t end)
{
    Section rst;

    rst.cellType = this->cellType;
    memcpy(rst.name, this->name, 32);
    rst.start = beg, rst.end = end;
    rst.nBdy = this->nBdy;

    for(auto id = beg; id<=end; ++id)
    {
        rst.addCell(id, this->cell(id), this->flag(id));
    }

    return rst;
}

CGFile::CGFile(string filename, int mode) : filename_(filename)
{
    switch (mode)
    {
    case CG_MODE_READ:
        this->checkFile();
        // this->loadCGIOInfo();
        this->loadCoordinateInfo();
        break;
    case CG_MODE_WRITE:
        check(cg_open(filename.c_str(), mode, &fp), format("can not open file %s to write\n", filename.c_str()));
        coordname_.push_back("CoordinateX");
        coordname_.push_back("CoordinateY");
        coordname_.push_back("CoordinateZ");
        break;
    default:
        break;
    }
    isOpen_ = true;
}

vector<vector<cgsize_t>> CGFile::allFaceInCell(const vector<cgsize_t>& idList, ElementType_t type)
{
    check(
        (CGFile::nodeIdToBuildFaceInCell.count(type)!=0),
        format("not supported cell type %d\n", type));

    vector<vector<cgsize_t>> rst;

    for(auto nodes : CGFile::nodeIdToBuildFaceInCell.at(type))
    {
        vector<cgsize_t> face;
        for(auto it : nodes) face.push_back(idList[it]);
        rst.emplace_back(face);
    }

    return rst;
}

CGFile::Section& CGFile::addSection()
{
    int n = sections_.size();
    sections_[n] = Section{};

    return sections_[n];
}

CGFile::Section CGFile::bodySection()
{
    Section rst;
    string name;
    set<ElementType_t> types;
    vector<cgsize_t> startFlag;

    for(auto it : bodySectionIDList_)
    {
        auto &curS = sections_[it];
        name += (curS.name + string("_"));
        startFlag.push_back(curS.start);
        rst.end = std::max(rst.end, curS.end);
        types.insert(curS.cellType);
    }
    rst.start = *(std::min_element(startFlag.begin(), startFlag.end()));
    std::memcpy(rst.name, name.c_str(), std::min((size_t)31, name.rfind('_')));
    rst.cellType = types.size() == 1 ? *types.begin() : ElementType_t::MIXED;

    return rst;
}

const vector<int>& CGFile::bodySectionIdList() const
{
    return bodySectionIDList_;
}

const vector<int>& CGFile::bdySectionIdList() const
{
    return bdySectionIDList_;
}

void CGFile::checkFile()
{
    char basename[33], zonename[33];
    int nbases, nzones, ngrids, ncoords, phys_dim;
    int filetype;
    ZoneType_t zoneType;
    check(
        cg_is_cgns(filename_.c_str(), &filetype),
        format("%s is not a CGNS file\n", filename_.c_str()));

    check(
        (cg_set_file_type(CG_FILE_ADF) == CG_OK && cg_open(filename_.c_str(), CG_MODE_READ, &fp) == CG_OK) || (cg_set_file_type(CG_FILE_HDF5) == CG_OK && cg_open(filename_.c_str(), CG_MODE_READ, &fp) == CG_OK),
        format("can not open file: %s \n", filename_.c_str()));
    isOpen_ = true;

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
        format("no section found in file %s", filename_.c_str()));
    for(int iSection = 1; iSection<=nSection; ++iSection)
    {
        section.id = iSection;
        int flag;
        cg_section_read(fp, ibase, izone, iSection, section.name, &section.cellType, &section.start, &section.end, &section.nBdy, &flag);
        cg_ElementDataSize(fp, ibase, izone, iSection, &section.dataSize);
        if(this->isBodySection(section)) bodySectionIDList_.push_back(iSection);
        else bdySectionIDList_.push_back(iSection);

        sections_.insert({iSection, section});
        globalNumber_.insert({iSection, section.end-section.start+1});
    }

    coordDataType_ = (this->precision() == 32) ? CGNS_ENUMV(RealSingle) : CGNS_ENUMV(RealDouble);
}

ElementType_t CGFile::CellType(const int nodeCnt, const int dim)
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

ElementType_t CGFile::CellType(const int CGNSCellTypeFalg)
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

void CGFile::close()
{
    if(!isOpen_) return;
    check(cg_close(fp), "", format("close file %s\n", filename_.c_str()));
    isOpen_ = false;
}

bool CGFile::isBodySection(const Section& section)
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
        cg_poly_elements_read(fp, ibase, izone, section.id, bodyConn.data(), offset.data(), &flag);
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

vector<vector<double>> CGFile::loadCoordinate(const cgsize_t rangeMin, const cgsize_t rangeMax)
{
    vector<vector<double>> data(coordname_.size(), vector<double>(rangeMax-rangeMin+1, 0.0));
    if(coordDataType_ == CGNS_ENUMV(RealDouble))
    {
        for(auto i = 0; i<3; ++i)
        {
            auto coordname = coordname_[i].c_str();
            check(cg_coord_read(fp, ibase, izone, coordname, coordDataType_, &rangeMin, &rangeMax, data[i].data()), format("read %s\n", coordname));
        }
    }
    else if(coordDataType_ == CGNS_ENUMV(RealSingle))
    {
        for(auto i = 0; i<3; ++i)
        {
            vector<float> tmp;
            auto coordname = coordname_[i].c_str();
            check(cg_coord_read(fp, ibase, izone, coordname, coordDataType_, &rangeMin, &rangeMax, tmp.data()), format("read %s\n", coordname));
            for(auto id = 0; id<tmp.size(); ++id) data[i][id] = tmp[id];
        }
    }

    return data;
}

void CGFile::loadCoordinateInfo()
{
    char coordname[33];
    cgsize_t range_min = 1, range_max = zoneInfo[0];

    cg_ncoords(fp, ibase, izone, &ncoords);
    check(cg_coord_info(fp, 1, 1, 1, &coordDataType_, coordname));

    vector<double> coordsDouble(range_max, 0.0);
    for (int icoord = 1; icoord <= ncoords; ++icoord)
    {
        cg_coord_info(fp, ibase, izone, icoord, &coordDataType_, coordname);
        coordname_.push_back(coordname);
    }
}

void CGFile::loadCGIOInfo()
{
    double rootId;
    int cgioNum;

    auto getChildrenIds = [&](double parentId){
        int numChildren, tmp;
        vector<double> childrenIds;

        cgio_number_children(cgioNum, parentId, &numChildren);
        childrenIds.assign(numChildren, 0.0);
        cgio_children_ids(cgioNum, parentId, 1, numChildren, &tmp, childrenIds.data());

        return childrenIds;
    };

    auto lookupNodeIdByLabel = [&](double nodeId, string label)
    {
        char tmpLabel[CGIO_MAX_LABEL_LENGTH+1];
        vector<double> rst;

        for(auto id : getChildrenIds(nodeId))
        {
            cgio_get_label(cgioNum, id, tmpLabel);
            if(string(tmpLabel) == label)
            {
                rst.push_back(id);
            }
        }

        return rst;
    };

    cg_get_cgio(fp, &cgioNum);
    cgio_get_root_id(cgioNum, &rootId);
    auto baseId = lookupNodeIdByLabel(rootId, "CGNSBase_t");
    auto zoneId = lookupNodeIdByLabel(baseId[0], "Zone_t");
    auto sectionId = lookupNodeIdByLabel(zoneId[0], "Elements_t");

    char name[CGIO_MAX_NAME_LENGTH+1];
    char type[CGIO_MAX_DATATYPE_LENGTH+1];
    for(auto id : sectionId)
    {
        cgio_get_name(cgioNum, id, name);
        // sectionNodeId.insert({name, id});
        cgio_get_data_type(cgioNum, id, type);
        // sectionDataType.insert({name, type});
    }
}

CGFile::Section& CGFile::loadSection(const int id)
{
    auto &curSection = sections_[id];    
    // if(!curSection.data.empty()) return curSection;
    curSection.clear();

    cgsize_t ID = curSection.start;
    vector<cgsize_t> tmp;
    auto loader = CellLoader{curSection, fp, curSection.start, curSection.end};
    while (loader.nextCell(tmp))
    {
        curSection.addCell(ID++, tmp, loader.flag);
    }

    return curSection;
}

CGFile::Section& CGFile::loadSection(const int id, const cgsize_t start, const cgsize_t end)
{
    auto &curSection = sections_[id];
    // if(!curSection.data.empty() && start>=curSection.start && end <= curSection.end) return curSection;
    curSection.clear();
    curSection.start = start;
    curSection.end = end;

    cgsize_t ID = start;
    vector<cgsize_t> tmp;
    auto loader = CellLoader{curSection, fp, start, end};
    while (loader.nextCell(tmp))
    {
        curSection.addCell(ID++, tmp, loader.flag);
    }

    return curSection;
}

CGFile::Section CGFile::loadSection(const vector<int> idList)
{
    Section curS;
    vector<cgsize_t> tmp;
    string name;
    set<ElementType_t> types;

    for(auto bodyId : idList)
    {
        auto &curBody = sections_[bodyId];
        curBody.clear();
        cgsize_t ID = curBody.start;
        name += (curBody.name + string("_"));
        curS.start = std::min(curS.start, curBody.start);
        curS.end = std::max(curS.end, curBody.end);
        types.insert(curBody.cellType);
        auto loader = CellLoader{curBody, fp, curBody.start, curBody.end};
        while (loader.nextCell(tmp))
        {
            curS.addCell(ID++, tmp, loader.flag);
        }
    }
    std::memcpy(curS.name, name.c_str(), std::min((size_t)31, name.rfind('_')));
    curS.cellType = types.size() == 1 ? *types.begin() : ElementType_t::MIXED;

    return curS;
}

CGFile::Section CGFile::loadSection(const vector<int> idList, const cgsize_t start, const cgsize_t end)
{
    Section curS;
    vector<cgsize_t> tmp;
    string name;
    set<ElementType_t> types;
    curS.start = start;
    curS.end = end;
    cgsize_t ID = start;

    for(auto bodyId : idList)
    {
        auto &curBody = sections_[bodyId];
        if(curBody.start>end || curBody.end<start) continue;
        
        curBody.clear();
        name += (curBody.name + string("_"));
        types.insert(curBody.cellType);
        auto loader = CellLoader{curBody, fp, std::max(start, curBody.start), std::min(end, curBody.end)};
        while (loader.nextCell(tmp))
        {
            curS.addCell(ID++, tmp, loader.flag);
        }
    }
    std::memcpy(curS.name, name.c_str(), std::min((size_t)31, name.rfind('_')));
    curS.cellType = types.size() == 1 ? *types.begin() : ElementType_t::MIXED;

    return curS;
}

cgsize_t CGFile::nCell() const
{
    return zoneInfo[1];
}

void CGFile::nCell(const cgsize_t ncell)
{
    zoneInfo[1] = ncell;
}

cgsize_t CGFile::nNode() const
{
    return zoneInfo[0];
}

void CGFile::nNode(const cgsize_t nnode)
{
    zoneInfo[0] = nnode;
}

int CGFile::precision()
{
    int rst;
    check(cg_precision(fp, &rst));

    return rst;
}

CGFile::Section& CGFile::section(const int id)
{
    return sections_[id];
}

string CGFile::stringAFace(const vector<cgsize_t>& face)
{
    if(face.empty()) return "";

    set<cgsize_t> tmp;
    for(auto it : face) tmp.insert(it);

    return std::accumulate(
        std::next(tmp.begin()),
        tmp.end(),
        std::to_string(*tmp.begin()),
        [](string sum, const cgsize_t id){
            return sum + "-" + std::to_string(id);
        });
}

void CGFile::writeCoordinate(const vector<vector<double>>& data)
{
    int dummy;

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
    for(auto iCoord : {0, 1, 2})
    {
        const char *coordname = coordname_[iCoord].c_str();
        check(cg_coord_write(fp, 1, 1, RealDouble, coordname, data[iCoord].data(), &dummy), format("write coordinate %s\n", coordname));
    }
}

void CGFile::writeSection(Section &curS, const map<cgsize_t, cgsize_t>& nodeIdG2L)
{
    // update id
    curS.dataSize = std::accumulate(
        curS.data.begin(),
        curS.data.end(),
        0,
        [](cgsize_t sum, vector<cgsize_t> vec){
            return sum + vec.size();
        });

    curS.start = idOffset_;
    curS.end = curS.start + curS.data.size() - 1;
    idOffset_ = curS.end + 1;

    int dummy;
    if(curS.isMixed())
    {
        curS.dataSize += curS.typeFlag.size();
        vector<cgsize_t> data(curS.dataSize, 0);
        int icell=0, i=0;
        for(auto it : curS.data)
        {
            data[i++] = curS.typeFlag[icell++];
            for(auto j : it)
            {
                data[i++] = nodeIdG2L.at(j);
            }
        }
        check(
            cg_poly_section_write(fp, 1, 1, curS.name, curS.cellType, curS.start, curS.end, 0, data.data(), curS.offset.data(), &dummy),
            format("file %s write section %s\n", filename_.c_str(), curS.name));
    }
    else
    {
        vector<cgsize_t> data(curS.dataSize);
        int i=0;
        for(auto it : curS.data)
        {
            for(auto j : it)
            {
                data[i++] = nodeIdG2L.at(j);
            }
        }
        check(
            cg_section_write(fp, 1, 1, curS.name, curS.cellType, curS.start, curS.end, curS.nBdy, data.data(), &dummy),
            format("file %s write section %s\n", filename_.c_str(), curS.name));
    }

    std::cout << format("  %s: write section %s [%d: %d, %d]\n", filename_.c_str(), curS.name, curS.end-curS.start+1, curS.start, curS.end);
}

void CGFile::writeGlobalInfo(const Section &curS, const cgsize_t n, const cgsize_t start, const cgsize_t end)
{
    cg_goto(fp, 1, "Zone_t", 1, curS.name, 0, "end");
    std::ostringstream ostr;
    ostr << "GlobalNumber:" << n << '\n';
    ostr << "GlobalStart:" << start << '\n';
    ostr << "GlobalEnd:" << end << '\n';
    auto str = ostr.str();
    cg_descriptor_write("GlobalDomainInfo", str.c_str());
}

} // namespace MeshCut
