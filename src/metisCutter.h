#pragma once

#include <memory>
#include <unordered_map>
#include <map>

#ifdef SUPPORT_METIS
#include <metis.h>
#else 
    typedef int64_t idx_t;
#endif

#include "meshCutter.h"
#include "cgFile.h"

using std::map;
using std::unordered_map;


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
    unordered_map<int, map<string, vector<cgsize_t>>> outerFace_;
    cgsize_t globalOffset_ = 0;

private:
    void cut_metis(string mesh, const int np);
    void cut_parmetis(string mesh, const int np);

    int cut_metis(const vector<vector<idx_t>>& cellToplogy, idx_t np, vector<idx_t>& cellPartition, vector<idx_t>& nodePartition);
    int cut_parmetis(idx_t np, vector<idx_t>& elmdist, vector<idx_t>& eptr, vector<idx_t>& eind, vector<idx_t>& cellPartition);

    void openToWrite(string filename, const int np);

    void rwBody(const int ifile, const CGFile::Section& bigBody);

    void rwBoundary(const int ifile);

    void rwInterface(const int ifile);

    void rwNode(const int ifile);

    void updateOuterFace(const CGFile::Section& curS, const int id);

    void writeGlobalInfo(const CGFile::Section& curS, const int id);
};
} // namespace MeshCut
