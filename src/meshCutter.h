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

    virtual void cut(int argc, char** argv) = 0;
};
} // namespace MeshCut

#endif