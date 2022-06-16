#ifndef DECOMPOSER_H
#define DECOMPOSER_H

#include <string>
#include <vector>

using std::vector;
using std::string;

namespace MeshCut
{
class MeshCutter
{
private:
    /* data */
public:
    MeshCutter(/* args */);
    virtual ~MeshCutter();

    virtual void cut(
        const std::string &mesh, 
        const int npart, 
        const int nx = 0, 
        const int ny = 0, 
        const int nz = 0,
        const string weightFilename = "",
        bool multiZone = false, 
        vector<string> fluidDomainNameRule={},
        string interiorSection = "") = 0;
};
} // namespace MeshCut

#endif