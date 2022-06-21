#include <iostream>
#include <string>
#include <memory>

#include "../src/cmdLine.h"
#include "stringUtil.h"
#include "cartesianCutter.h"
#include "metisCutter.h"

using MeshCut::MeshCutter;
using MeshCut::CartesianCutter;
using MeshCut::MetisCutter;

int main(int argc, char** argv)
{
    CmdLine cl{};
    cl.regist<std::string, CmdLine::MustOffer>("t", "tool", "cutter tools [parmetis, cartesian]."," ");
    cl.parse(argc, argv);

    std::shared_ptr<MeshCutter> cutter;
    auto tool = cl.get<std::string>("tool");
    if(tool == "cartesian")
    {
        cutter = std::make_shared<CartesianCutter>(CartesianCutter{});
    }
    else if(tool == "metis")
    {
        cutter = std::make_shared<MetisCutter>(MetisCutter{});
    }
    else
    {
        std::cout << "Not support tool: " << tool << "\n";
        return 0;
    }

    cutter->cut(argc, argv);
    
    return 0;
}