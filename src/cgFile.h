#pragma once


#include <string>
#include <limits>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <set>

#include <cgnslib.h>

using std::set;
using std::map;
using std::string;
using std::vector;
using std::unordered_map;

namespace MeshCut
{
void check(const int runCode, string&& errStr="", string&& okStr="");

void check(bool state, string&& errStr="", string&& okStr="");

class CGFile
{
private:
    struct Section;
    struct Face;
    struct CellLoader;

public:
    CGFile(string filename, int mode = CG_MODE_READ);

public:
    vector<int> bodySections() const;

    vector<int> bdySections() const;
    
    void loadSection(const int id);

    vector<vector<double>> loadCoordinate(const cgsize_t rangeMin, const cgsize_t rangeMax);

    cgsize_t nCell() const;

    cgsize_t nNode() const;

    Section& section(const int id);

    void writeCoordinate(const vector<vector<double>>& data);

private:
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
    string filename_;
    int fp;
    map<int, int> smallFiles_;

    int cell_dim_, ibase=1, izone=1, igrid=1, ncoords;
    cgsize_t zoneInfo[3];
    DataType_t coordDataType_;
    vector<string> coordname_;
    
    map<int, cgsize_t> globalNumber_;
    map<int, cgsize_t> globalOffset_;
    map<int, Section> sections_;
    map<int, map<cgsize_t, cgsize_t>> nodeIdMap_;
    vector<vector<double>> nodesData_;
    set<cgsize_t> nodeIdInBox_;
    map<int, map<string, vector<cgsize_t>>> outerFace_;
    vector<int> bodySection_, bdySection_, interiorSection_;
    vector<cgsize_t> idOffset_;

private:

    void checkFile();

    bool isBodySection(const Section& section);

    void loadCoordinateInfo();

    int precision();
};
} // namespace MeshCut
