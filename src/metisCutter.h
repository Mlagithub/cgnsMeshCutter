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
    map<int, map<cgsize_t, cgsize_t>> nodeIdG2L_;
    set<cgsize_t> nodeIds_, cellIds_;
    vector<idx_t> cellPartition_, nodePartition_;
    map<int, map<string, vector<cgsize_t>>> outerFace_;
    cgsize_t globalOffset_ = 0;

private:
    void cut(string mesh, const int np);

    int cut(const vector<vector<idx_t>>& cellToplogy, idx_t np, vector<idx_t>& cellPartition, vector<idx_t>& nodePartition);

    void openToWrite(string filename, const int np);

    void rwBody(const int ifile, const CGFile::Section& bigBody);

    void rwBoundary(const int ifile);

    void rwInterface(const int ifile);

    void rwNode(const int ifile);

    void updateOuterFace(const CGFile::Section& curS, const int id);

    void writeGlobalInfo(const CGFile::Section& curS, const int id);
};
} // namespace MeshCut
