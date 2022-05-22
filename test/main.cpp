#include <iostream>
#include <string>

#include "cmdLine.h"
#include "parmetisCutter.h"
#include "cartesianCutter.h"

using MeshCut::CartesianCutter;
using MeshCut::MeshCutter;
using MeshCut::ParMetisMeshCutter;

int main(int argc, char** argv)
{
    CmdLine cl{};
    cl.regist<std::string, CmdLine::MustOffer>("t", "tool", "cutter tools [parmetis, cartesian]."," ");
    cl.regist<std::string, CmdLine::MustOffer>("m", "mesh", "mesh file name."," ");
    cl.regist<size_t, CmdLine::MustOffer>("np", "npart", "number to be cutting."," ");
    cl.regist<size_t, CmdLine::MustOffer>("nx", "npart_x", "number to be cutting at x axix.","1");
    cl.regist<size_t, CmdLine::MustOffer>("ny", "npart_y", "number to be cutting at y axix.","1");
    cl.regist<size_t, CmdLine::MustOffer>("nz", "npart_z", "number to be cutting at z axix.","1");
    cl.parse(argc, argv);

    MeshCutter *cutter;
    auto tool = cl.get<std::string>("tool");
    if (tool == "parmetis")
    {
        cutter = new ParMetisMeshCutter(argc, argv);
    }
    else if(tool == "cartesian")
    {
        cutter = new CartesianCutter();
    }
    cutter->cut(cl.get<std::string>("mesh"), cl.get<size_t>("npart"), cl.get<size_t>("nx"), cl.get<size_t>("ny"), cl.get<size_t>("nz"));

    return 0;
}