#include <algorithm>
#include <cmath>

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
    for(auto i=0; i<np; ++i) this->openToWrite(meshFilename, np, i);

    for(auto it : bigMesh_->bodySections())
    {
        // load data into memory from file
        bigMesh_->loadSection(it);
        auto &curSection = bigMesh_->section(it);

        // cut
        vector<idx_t> cellPartition, nodePartition;
        check(
            cut(curSection.data, np, cellPartition, nodePartition) == METIS_OK,
            format("METIS_V3_PartMeshKway return error while divide mesh: %s", meshFilename.c_str()),
            format("METIS_V3_PartMeshKway decompose %s successfuly\n", curSection.name));

        // write sub-mesh data to file
        for(auto i=0; i<np; ++i)
        {
            // write coordinate
            this->writeNode(nodePartition, i);

            // write body-section

            // write bdy-section

            // write interface
        }
    }

    auto i = getchar();
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

void MetisCutter::openToWrite(string bigFileName, const int np, const int id)
{
    char fmt[128], smallFilename[128];
    sprintf(fmt, "%%s_%%0%dd.cgns", int(std::log10(np) + 1));
    sprintf(smallFilename, fmt, bigFileName.substr(0, bigFileName.size() - 5).c_str(), id);

    smallMesh_.push_back(std::make_shared<CGFile>(smallFilename, CG_MODE_WRITE));
}

void MetisCutter::writeNode(const vector<idx_t>& nodePartition, const int ifile)
{
    map<cgsize_t, cgsize_t> nodeIdMap;
    vector<cgsize_t> nodeId;
    cgsize_t id = 0;
    for(auto i = 0; i<nodePartition.size(); ++i)
    {
        if(nodePartition[i]==ifile)
        {
            nodeIdMap.insert({i+1, id++});
            nodeId.push_back(i+1);
        }
    }
    smallMesh_[ifile]->writeCoordinate(bigMesh_->loadCoordinate(*nodeId.begin(), *nodeId.rbegin()));
}

} // namespace MeshCut
