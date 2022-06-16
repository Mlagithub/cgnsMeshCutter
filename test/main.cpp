#include <iostream>
#include <string>

#include "cmdLine.h"
#include "cartesianCutter.h"
#include "stringUtil.h"

using MeshCut::CartesianCutter;
using MeshCut::MeshCutter;

int main(int argc, char** argv)
{
    CmdLine cl{};
    cl.regist<std::string, CmdLine::MustOffer>("t", "tool", "cutter tools [cartesian]."," ");
    cl.regist<std::string, CmdLine::MustOffer>("m", "mesh", "mesh file name."," ");
    cl.regist<size_t, CmdLine::MustOffer>("np", "npart", "number to be cutting."," ");
    cl.regist<size_t, CmdLine::MustOffer>("nx", "npart_x", "number to be cutted at x axix.","1");
    cl.regist<size_t, CmdLine::MustOffer>("ny", "npart_y", "number to be cutted at y axix.","1");
    cl.regist<size_t, CmdLine::MustOffer>("nz", "npart_z", "number to be cutted at z axix.","1");
    cl.regist<bool, CmdLine::NotMustOffer>("multizone", "multizone_mode", "switch on if mesh is multizone [on, off].", "off");
    cl.regist<string, CmdLine::NotMustOffer>("fluid", "fluid_domain_name", "if multiZone switch on, offer comma separated fluid domain name.","");
    cl.regist<string, CmdLine::NotMustOffer>("interior", "interior_section", "only support multizone mode","");
    cl.regist<string, CmdLine::NotMustOffer>("weight", "weight_filename", "sub-mesh weight at x/y/z direction","");
    cl.parse(argc, argv);

    MeshCutter *cutter;
    auto tool = cl.get<std::string>("tool");
    if(tool == "cartesian")
    {
        cutter = new CartesianCutter();
    }
    else
    {
        std::cout << "Not support tool: " << tool << "\n";
        return 0;
    }


    cutter->cut(
        cl.get<std::string>("mesh"), 
        cl.get<size_t>("npart"), 
        cl.get<size_t>("nx"), 
        cl.get<size_t>("ny"), 
        cl.get<size_t>("nz"),
        cl.get<string>("weight"),
        cl.get<bool>("multizone"),
        stringSplit(cl.get<string>("fluid"),","),
        cl.get<string>("interior")
        );
    
    delete cutter;
    
    return 0;
}