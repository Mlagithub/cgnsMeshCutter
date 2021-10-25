#include <iostream>
#include <string>

#include "meshCutter.h"
#include "parmetisCutter.h"


using MeshCut::MeshCutter;
using MeshCut::ParMetisMeshCutter;

int main(int argc, char** argv)
{
    MeshCutter* cutter = new ParMetisMeshCutter(argc, argv);
    cutter->cut("box.cgns", 2);
    return 0;
}