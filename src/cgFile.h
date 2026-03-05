#pragma once

#include <cgnslib.h>

#include <string>
#include <limits>
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <set>

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

    Section bodySection();

    const vector<int>& bodySectionIdList() const;

    const vector<int>& bdySectionIdList() const;

    static ElementType_t CellType(const int nodeCnt, const int dim);

    static ElementType_t CellType(const int CGNSCellTypeFalg);

    void close();

    vector<vector<double>> loadCoordinate(const cgsize_t rangeMin, const cgsize_t rangeMax);
    
    Section& loadSection(const int id);

    Section& loadSection(const int id, const cgsize_t start, const cgsize_t end);
        
    Section loadSection(const vector<int> idList);

    Section loadSection(const vector<int> idList, const cgsize_t start, const cgsize_t end);

    cgsize_t nCell() const;

    void nCell(const cgsize_t ncell);

    cgsize_t nNode() const;

    void nNode(const cgsize_t nnode);

    Section& section(const int id);

    static string stringAFace(const vector<cgsize_t>& face);

    void writeCoordinate(const vector<vector<double>>& data);

    void setIdOffset(cgsize_t offset);

    void writeSection(Section &curS, const map<cgsize_t, cgsize_t>& nodeIdG2L);

    void writeGlobalInfo(const Section &curS, const cgsize_t n, const cgsize_t start, const cgsize_t end);

public:
    struct Section
    {
        int id = 1;
        char name[33] = {0};
        ElementType_t cellType;
        cgsize_t start=1, end=1, dataSize=0;
        int nBdy=0;
        vector<vector<cgsize_t>> data;
        vector<cgsize_t> offset{0}, typeFlag;
        map<cgsize_t, cgsize_t> IDMap;

        void addCell(const cgsize_t id, const vector<cgsize_t>& idList, const cgsize_t flag = CG_Null);
        void clear();
        const vector<cgsize_t>& cell(const int id) const;
        cgsize_t flag(const int id) const;
        bool isMixed() const;
        void printToFile(string fname);
        Section subSection(const cgsize_t beg, const cgsize_t end);
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
        /**
         * @brief Construct a new Cell Loader object
         * 
         * @param curS CGNS 文件结构 Section 引用
         * @param fp 文件句柄
         * @param lowerBound 当前 Section 的 ID 下界，包含该边界值 
         * @param upperBound 当前 Section 的 ID 上界，包含该边界值
         */
        CellLoader(Section& curS, const int fp, const cgsize_t lowerBound = -1, const cgsize_t upperBound = std::numeric_limits<cgsize_t>::max());

        /**
         * @brief 将入参赋值为当前加载器构造时指定范围内的组成下一个Cell的 Node ID
         * 
         * @param nodes Node 列表引用
         * @return true 赋值成功
         * @return false 已达到末尾，nodes 不可用
         */
        bool nextCell(vector<cgsize_t>& nodes);
        
        int flag = CG_Null;

    private:
        Section& curS;
        vector<cgsize_t> data;
        vector<cgsize_t> offset;
        int fp, npe;
        cgsize_t id = 0, blockSize = 0, start = 0, blockId = 0, len = 0, lowerBound = 0, upperBound = 0;
        bool loadNextBlock();
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
    vector<int> bodySectionIDList_, bdySectionIDList_, interiorSection_;
    cgsize_t idOffset_ = 1;

private:

    void checkFile();

    bool isBodySection(const Section& section);

    void loadCoordinateInfo();

    void loadCGIOInfo();

    int precision();
};
} // namespace MeshCut
