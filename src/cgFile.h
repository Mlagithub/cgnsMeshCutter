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
public:
    static const unordered_map<ElementType_t, vector<vector<int>>> nodeIdToBuildFaceInCell;
    struct Section;

private:
    struct Face;
    struct CellLoader;

public:
    CGFile(string filename, int mode = CG_MODE_READ);

public:
    static vector<vector<cgsize_t>> allFaceInCell(const vector<cgsize_t>& idList, ElementType_t type);

    Section& addSection();

    vector<int> bodySections() const;

    vector<int> bdySections() const;

    static ElementType_t CellType(const int nodeCnt, const int dim);

    static ElementType_t CellType(const int CGNSCellTypeFalg);

    void close();

    Section& loadSection(const int id);

    vector<vector<double>> loadCoordinate(const cgsize_t rangeMin, const cgsize_t rangeMax);

    cgsize_t nCell() const;

    void nCell(const cgsize_t ncell);

    cgsize_t nNode() const;

    void nNode(const cgsize_t nnode);

    Section& section(const int id);

    static string stringAFace(const vector<cgsize_t>& face);

    void writeCoordinate(const vector<vector<double>>& data);

    void writeSection(Section &curS, const map<cgsize_t, cgsize_t>& nodeIdG2L);

public:
    struct Section
    {
        int id = 1;
        char name[33] = {0};
        ElementType_t cellType;
        cgsize_t start=1, end=1, dataSize=0;
        int nBdy=0, flag;
        vector<vector<cgsize_t>> data;
        vector<cgsize_t> offset, typeFlag;

        bool isMixed() const;
    };

private:
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
        CellLoader(Section& curS, const int fp);
        bool nextCell(vector<cgsize_t>& nodes);
        Section& curS;
        vector<cgsize_t> data;
        vector<cgsize_t> offset;
        cgsize_t len;
        int npe, flag = -1;
        cgsize_t id = 0;
    };

private:
    string filename_;
    int fp;
    map<int, int> smallFiles_;
    bool isOpen_ = false;

    int cell_dim_, ibase=1, izone=1, igrid=1, ncoords;
    cgsize_t zoneInfo[3] = {0, 0, 0};
    DataType_t coordDataType_;
    vector<string> coordname_;
    
    map<int, cgsize_t> globalNumber_;
    map<int, cgsize_t> globalOffset_;
    map<int, Section> sections_;
    vector<int> bodySection_, bdySection_, interiorSection_;
    cgsize_t idOffset_ = 1;

private:

    void checkFile();

    bool isBodySection(const Section& section);

    void loadCoordinateInfo();

    int precision();
};
} // namespace MeshCut
