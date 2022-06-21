#include <iostream>

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
}

CGFile::CellLoader::CellLoader(Section& curS, vector<cgsize_t>& data, vector<cgsize_t>& offset, const int fp)
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

bool CGFile::CellLoader::nextCell(vector<cgsize_t>& nodes)
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
            cg_coord_read(fp, ibase, izone, coordname, coordDataType_, &rangeMin, &rangeMax, tmp.data());
            check(cg_coord_read(fp, ibase, izone, coordname, coordDataType_, &rangeMin, &rangeMax, data[i].data()), format("read %s\n", coordname));
            for(auto it : tmp) data[i].push_back(it);
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

void CGFile::loadSection(const int id)
{
    int curOffset = 0;
    vector<cgsize_t> data, offset, tmp;
    auto &curSection = sections_[id];
    auto isMixed = (curSection.cellType == CGNS_ENUMV(MIXED));
    auto loader = CellLoader{curSection, data, offset, fp};
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
}

cgsize_t CGFile::nCell() const
{
    return zoneInfo[1];
}

cgsize_t CGFile::nNode() const
{
    return zoneInfo[0];
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
        check(cg_coord_write(fp, 1, 1, RealDouble, coordname, data.data(), &dummy), format("write coordinate %s\n", coordname));
    }
}

} // namespace MeshCut
