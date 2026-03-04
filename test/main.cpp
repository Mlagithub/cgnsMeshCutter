#include <iostream>
#include <memory>
#include <string>

#include "../src/cmdLine.h"
#include "cartesianCutter.h"
#include "metisCutter.h"
#include "stringUtil.h"

using MeshCut::CartesianCutter;
using MeshCut::MeshCutter;
using MeshCut::MetisCutter;

int main(int argc, char **argv)
{
    try
    {
        CmdLine cl{};
        cl.regist<std::string, CmdLine::MustOffer>("t", "tool", "cutter tools [metis, cartesian].", "metis");
        cl.regist<bool, CmdLine::NotMustOffer>("u", "usage", "pirnt usage with -t option and exit [on, off].", "off");
        cl.parse(argc, argv);

        std::shared_ptr<MeshCutter> cutter;
        auto tool = cl.get<std::string>("tool");
        if (tool == "cartesian") { cutter = std::make_shared<CartesianCutter>(CartesianCutter{}); }
        else if (tool == "metis") { cutter = std::make_shared<MetisCutter>(MetisCutter{}); }
        else
        {
            std::cout << "Not supported tool: " << tool << "\n";
            return 0;
        }

        if (cl.get<bool>("usage")) argv[1] = const_cast<char *>("-h");
        cutter->cut(argc, argv);

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
