#include <iostream>
#include <numeric>

#include "format.h"
#include "cgFile.h"

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

CGFile::CGFile(string filename, int mode) : filename_(filename)
{
    switch (mode)
    {
    case CG_MODE_READ:
        this->checkFile();
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

CGFile::CellLoader::CellLoader(Section& curS, const int fp)
    : curS(curS)
{
    cg_section_read(fp, 1, 1, curS.id, curS.name, &curS.cellType, &curS.start, &curS.end, &curS.nBdy, &curS.flag);
    cg_ElementDataSize(fp, 1, 1, curS.id, &curS.dataSize);
    len = curS.end - curS.start + 1;
    id = 0;

    if(curS.isMixed())
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

bool CGFile::CellLoader::nextCell(vector<cgsize_t>& nodes)
{
    if(id>=len) return false;

    cgsize_t start, end;
    if(curS.isMixed())
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

bool CGFile::Section::isMixed() const
{
    return (cellType == CGNS_ENUMV(MIXED));
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

vector<int> CGFile::bodySections() const
{
    return bodySection_;
}

vector<int> CGFile::bdySections() const
{
    return bdySection_;
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
        cg_section_read(fp, ibase, izone, iSection, section.name, &section.cellType, &section.start, &section.end, &section.nBdy, &section.flag);
        cg_ElementDataSize(fp, ibase, izone, iSection, &section.dataSize);
        if(this->isBodySection(section)) bodySection_.push_back(iSection);
        else bdySection_.push_back(iSection); 

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

CGFile::Section& CGFile::loadSection(const int id)
{
    auto &curSection = sections_[id];
    if(!curSection.data.empty()) return curSection;

    int curOffset = 0;
    vector<cgsize_t> tmp;
    auto isMixed = curSection.isMixed();
    auto loader = CellLoader{curSection, fp};
    while (loader.nextCell(tmp))
    {
        curSection.data.push_back(tmp);
        curSection.offset.push_back(curOffset);
        curOffset += (tmp.size() + 1);
        if (isMixed)
        {
            curSection.typeFlag.push_back(loader.flag);
        }
        tmp.clear();
    }

    return curSection;
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

} // namespace MeshCut
