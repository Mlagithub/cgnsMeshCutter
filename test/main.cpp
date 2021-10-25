#include <iostream>
#include <string>

#include "cmdLine.h"
#include "meshCutter.h"
#include "parmetisCutter.h"


using MeshCut::MeshCutter;
using MeshCut::ParMetisMeshCutter;

int main(int argc, char** argv)
{
    CmdLine cl{};
    cl.regist<std::string, CmdLine::MustOffer>("m", "mesh", "mesh file name."," ");
    cl.regist<size_t, CmdLine::MustOffer>("np", "npart", "number to be cutting."," ");
    cl.parse(argc, argv);

    MeshCutter* cutter = new ParMetisMeshCutter(argc, argv);
    cutter->cut(cl.get<std::string>("mesh"), cl.get<size_t>("npart"));

    return 0;
}