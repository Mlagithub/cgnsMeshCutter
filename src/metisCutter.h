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
    int RANK = 0, SIZE = 1;
    cgns_filetype bigMesh_;
    vector<cgns_filetype> smallMesh_;
    map<int, map<cgsize_t, cgsize_t>> nodeIdG2L_;
    set<cgsize_t> nodeIds_, cellIds_;
    vector<idx_t> cellPartition_, nodePartition_;
    unordered_map<int, map<string, vector<cgsize_t>>> outerFace_;
    cgsize_t globalOffset_ = 0;
    set<int> ownerFile_;
    struct DecomposeResult
    {
        idx_t nCellThisRank, start;
    };
    struct Cell
    {
        idx_t id;   // begin from 1
        int partId; // begin from 0;
    };

private:
    void collect_interface();
    
    vector<Cell> collect_subBody(const DecomposeResult& decomposeResult, const vector<int>& ownerThread);

    void cut_metis(string mesh, const int np);

    int cut_metis(const vector<vector<idx_t>>& cellToplogy, idx_t np, vector<idx_t>& cellPartition, vector<idx_t>& nodePartition);

    void cut_parmetis(string mesh, const int np);

    int cut_parmetis(idx_t np, vector<idx_t>& elmdist, vector<idx_t>& eptr, vector<idx_t>& eind, vector<idx_t>& cellPartition);

    DecomposeResult decompose_body(CGFile::Section& bigBody, const int np);

    vector<int> openSubMeshToWrite(string filename, const int np);

    void rwBody(const int ifile, const CGFile::Section& bigBody);

    void rwBoundary(const int ifile);

    void rwInterface(const int ifile);

    void rwNode(const int ifile);

    std::string smallMeshName(string bigMeshName, const int id, const int np);

    void updateOuterFace(const CGFile::Section& curS, const int id);

    void writeGlobalInfo(const CGFile::Section& curS, const int id);

    void writeInterface(const int ifile, const int np);
};
} // namespace MeshCut
