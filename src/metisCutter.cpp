#include <algorithm>
#include <cmath>
#include <numeric>
#include <cstring>

#include "cmdLine.h"
#include "metisCutter.h"

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
#ifndef SUPPORT_METIS
    std::cout << "Not support method: Metis\n";
    return;
#endif

    CmdLine cl{};
    cl.regist<std::string, CmdLine::MustOffer>("m", "mesh", "mesh file name."," ");
    cl.regist<size_t, CmdLine::MustOffer>("np", "npart", "number to be cutting."," ");
    cl.regist<bool, CmdLine::NotMustOffer>("multizone", "multizone_mode", "switch on if mesh is multizone [on, off].", "off");
    cl.regist<string, CmdLine::NotMustOffer>("fluid", "fluid_domain_name", "if multiZone switch on, offer comma separated fluid domain name.","");
    cl.regist<string, CmdLine::NotMustOffer>("interior", "interior_section", "only support multizone mode","");
    cl.regist<string, CmdLine::NotMustOffer>("weight", "weight_filename", "sub-mesh weight at x/y/z direction","");
    cl.parse(argc, argv);

    this->cut(cl.get<std::string>("mesh"), cl.get<size_t>("npart"));
}


void MetisCutter::cut(string meshFilename, const int np)
{
    // open mesh file to read/write
    bigMesh_ = std::make_shared<CGFile>(meshFilename);
    this->openToWrite(meshFilename, np);

    for(auto it : bigMesh_->bodySections())
    {
        // load data into memory from file
        auto &bigBody = bigMesh_->loadSection(it);
        
        std::cout << format("Call METIS_V3_PartMeshKway decompose: %s\n", bigBody.name);
        // cut
        check(
            cut(bigBody.data, np, cellPartition_, nodePartition_) == METIS_OK,
            format("Call METIS_V3_PartMeshKway decompose: %s return error %s", bigBody.name),
            format("Call METIS_V3_PartMeshKway decompose: %s successfuly\n", bigBody.name));

        // write sub-mesh data to file
        for(auto ifile=0; ifile<np; ++ifile)
        {
            // collect cell/node id of sub-mesh
            nodeIdG2L_.clear(); // 1-base id
            nodeIds_.clear(); // 1-base id
            cellIds_.clear(); // 0-base id
            for(auto i=0; i<cellPartition_.size(); ++i)
            {
                if (cellPartition_[i] != ifile) continue;

                cellIds_.insert(i);
                for(auto it : bigBody.data[i]) if(nodeIds_.count(it)==0) nodeIds_.insert(it);
            }
            for(auto it : nodeIds_) nodeIdG2L_.insert({it, nodeIdG2L_.size()+1});

            // write coordinate
            this->rwNode(ifile);

            // write body-section
            this->rwBody(ifile, bigBody);

            // write bdy-section
            this->rwBoundary(ifile);

            // write interface
        }
    }
}

int MetisCutter::cut(const vector<vector<cgsize_t>>& cellToplogy, idx_t np, vector<idx_t>& cellPartition, vector<idx_t>& nodePartition)
{
    // metis options
    idx_t options[METIS_NOPTIONS];
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
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

    return METIS_PartMeshDual(&nCell, &nNode, eptr.data(), eind.data(), NULL, NULL, &ncommon, &np, NULL, options, &objval, cellPartition.data(), nodePartition.data());
}

void MetisCutter::openToWrite(string bigFileName, const int np)
{
    for (auto i = 0; i < np; ++i)
    {
        char fmt[128], smallFilename[128];
        sprintf(fmt, "%%s_%%0%dd.cgns", int(std::log10(np) + 1));
        sprintf(smallFilename, fmt, bigFileName.substr(0, bigFileName.size() - 5).c_str(), i);
        smallMesh_.push_back(std::make_shared<CGFile>(smallFilename, CG_MODE_WRITE));
    }
}

void MetisCutter::rwBody(const int ifile, const CGFile::Section& bigBody)
{
    auto &subFile = smallMesh_[ifile];

    // new a section
    auto &curS = subFile->addSection();
    strcpy(curS.name, bigBody.name);
    curS.cellType = bigBody.cellType;
    int n = 0;
    for (auto id : cellIds_)
    {
        curS.data.emplace_back(bigBody.data[id]);
        if (curS.isMixed()) curS.typeFlag.push_back(bigBody.typeFlag[id]);
    }
    // write data
    subFile->writeSection(curS, nodeIdG2L_);

    // update outerface 
    this->updateOuterFace(curS, ifile);
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

        if(!subBdy.data.empty()) subFile->writeSection(subBdy, nodeIdG2L_);

        // clear memory
        vector<vector<cgsize_t>>{}.swap(subBdy.data);
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



} // namespace MeshCut
