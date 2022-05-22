#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <set>
#include <iostream>

#include <cgnslib.h>

#include "vec.h"
#include "meshCutter.h"
#include "format.h"

using std::set;
using std::map;
using std::string;
using std::vector;
using std::unordered_map;

namespace MeshCut
{

class CartesianCutter : public MeshCutter
{
public:
    void cut(const std::string &mesh, const int npart, const int nx = 1, const int ny = 1, const int nz = 1) override;

    struct BoundBox
    {
        double max_x = -std::numeric_limits<double>::max();
        double max_y = -std::numeric_limits<double>::max();
        double max_z = -std::numeric_limits<double>::max();
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        double min_z = std::numeric_limits<double>::max();
        double range_x() const { return max_x - min_x; }
        double range_y() const { return max_y - min_y; }
        double range_z() const { return max_z - min_z; }
        double min(const int i)
        {
            switch (i)
            {
            case 1:
                return min_x;
                break;
            case 2:
                return min_y;
                break;
            case 3:
                return min_z;
                break;
            default:
                std::cout << format("error to get min of bbox.\n");
                break;
            }
        }
        double max(const int i)
        {
            switch (i)
            {
            case 1:
                return max_x;
                break;
            case 2:
                return max_y;
                break;
            case 3:
                return max_z;
                break;
            default:
                std::cout << format("error to get max of bbox.\n");
                break;
            }
        }
    };

private:
    struct Location
    {
        Location(){}
        Location(int i, int j, int k) : i(i), j(j), k(k) {}
        int i, j, k;
    };

    struct Section
    {
        int id = 1;
        char name[33];
        ElementType_t cellType;
        cgsize_t start=1, end=1, dataSize=0;
        int nBdy=0, flag;
        vector<vector<cgsize_t>> data;
        vector<cgsize_t> offset, typeFlag;
    };
    
    struct Face
    {
        Face(const vector<cgsize_t>& val);
        Face(vector<cgsize_t>&& val) noexcept;
        
        Face(const Face& other);
        Face(Face&& other) noexcept;

        Face& operator=(const Face& other);
        Face& operator=(Face&& other) noexcept;

        vector<cgsize_t> nodeList; 
        bool operator==(const Face &other) const;
        bool operator<(const Face& rhs) const;
    };

    struct CellLoader
    {
        CellLoader(Section& curS, vector<cgsize_t>& data, vector<cgsize_t>& offset, const int fp);
        bool nextCell(vector<cgsize_t>& nodes);
        Section& curS;
        vector<cgsize_t>& data;
        vector<cgsize_t>& offset;
        cgsize_t len;
        int npe, flag = -1;
        cgsize_t id = 0;
    };

private:
    int nx_, ny_, nz_, n_;
    int bigFile_;
    map<int, int> smallFiles_;
    map<int, BoundBox> bbox_;
    int fluidDomainId_ = 0;
    string fluidDomain_;
    cgsize_t lenOfFluid_ = 0;

    int cell_dim_, ibase=1, izone=1, igrid=1;
    cgsize_t zoneInfo[3];
    DataType_t coordDataType_;
    vector<string> coordname_;
    string interiorFacePrefix_;

    map<int, cgsize_t> globalNumber_;
    map<int, cgsize_t> globalOffset_;
    map<int, Section> sections_;
    map<int, map<cgsize_t, cgsize_t>> nodeIdMap_;
    vector<vector<double>> nodesData_;
    set<cgsize_t> nodeIdInBox_;
    map<int, map<Face, vector<cgsize_t>>> outerFace_;
    vector<int> bodySection_, bdySection_, interiorSection_;
    set<int> fluidSection_, solidSection_;
    vector<cgsize_t> idOffset_;

private:
    vector<vector<cgsize_t>> allFaceInCell(const vector<cgsize_t>& idList, ElementType_t type);

    vector<double> boundingBoxBigMesh();

    void boundingBoxSmallMesh();

    void check(const int runCode, string&& errStr="", string&& okStr="");

    void check(bool state, string&& errStr="", string&& okStr="");

    void checkFile(string filename);

    Vector<double> cellCenter(const vector<cgsize_t>& cell, const ElementType_t type);

    bool isBodySection(const Section& section);

    bool inBox(const BoundBox& box, const vector<cgsize_t>& cell, const ElementType_t type);

    Location id2Location(const int x) const;

    int location2Id(const Location& loc) const;

    int location2Id(Location&& loc) const;

    vector<int> neighbors(const int x);

    string openToWrite(string filename, const int id);

    bool pickBodyCellInBox(const BoundBox& box, Section& curS);

    int precision();

    void readNode();

    void rwBoundary(const int id);

    void rwInterior(const int id);

    void rwInterface(const int id);

    void updateNodeId(const Section& curS);

    void updateOuterFace(const Section& curS, const int id);

    void writeSection(Section &curS, const int id);

    void writeNode(const int id);

    void writeGlobalInfo(const Section& curS, const int id);
};

std::ostream &operator<<(std::ostream &os, const CartesianCutter::BoundBox &bbox);

} // namespace MeshCut
