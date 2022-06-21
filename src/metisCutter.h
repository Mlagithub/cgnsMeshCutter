#pragma once

#include <memory>

#ifdef SUPPORT_METIS
#include <metis.h>
#else 
    typedef int64_t idx_t;
#endif

#include "meshCutter.h"
#include "cgFile.h"

namespace MeshCut
{
class MetisCutter : public MeshCutter
{
public:
    MetisCutter();
    virtual ~MetisCutter();

public:
    void cut(int argc, char** argv) override;

private:
    using cgns_filetype = std::shared_ptr<CGFile>;
    cgns_filetype bigMesh_;
    vector<cgns_filetype> smallMesh_;

private:
    void cut(string mesh, const int np);

    int cut(const vector<vector<idx_t>>& cellToplogy, idx_t np, vector<idx_t>& cellPartition, vector<idx_t>& nodePartition);

    void openToWrite(string filename, const int np, const int id);

    void writeNode(const vector<idx_t>& nodePartition, const int id);
};
} // namespace MeshCut
