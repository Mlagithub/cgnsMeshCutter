#ifndef DECOMPOSER_H
#define DECOMPOSER_H

#include <string>

namespace MeshCut
{
class MeshCutter
{
private:
    /* data */
public:
    MeshCutter(/* args */);
    virtual ~MeshCutter();

    virtual void cut(const std::string &mesh, const int npart, const int nx = 0, const int ny = 0, const int nz = 0) = 0;
};
} // namespace MeshCut

#endif