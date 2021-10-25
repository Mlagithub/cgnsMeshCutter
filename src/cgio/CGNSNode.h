#ifndef CGNSNODE_H
#define CGNSNODE_H

#include <iostream>
#include <map>
#include <string>
#include <cgnslib.h>
#include <cgns_io.h>

namespace CGIO
{
/* 
 * @brief 将字符串表示的CGNS节点数据类型转换为内部枚举类型。
 * 
 * @param type 
 * @return constexpr short 
 */
constexpr short cg_typeof(const char *type)
{
    return (short)((type[0] << 8) + type[1]);
}

/* 
 * @brief 表示一个CGNS数据库节点
 * 
 */
class CGNSNode
{
public:
    /* 
     * @brief 节点数据类型
     * 
     */
    enum class DataType : unsigned short
    {
        // 无数据
        MT = cg_typeof("MT"),

        // 32有符号位整数 -- int32_t
        I4 = cg_typeof("I4"),

        // 64有符号位整数 -- int64_t
        I8 = cg_typeof("I8"),

        // 32位无符号整数 -- uint32_t
        U4 = cg_typeof("U4"),

        // 64位无符号整数 -- uint64_t
        U8 = cg_typeof("U8"),

        // 单精度浮点数 -- float
        R4 = cg_typeof("R4"),

        // 双精度浮点数 -- double
        R8 = cg_typeof("R8"),

        // 有符号字符 -- char
        C1 = cg_typeof("C1"),

        // 无符号字符 -- unsigned char
        B1 = cg_typeof("B1"),

        // 链接
        LK = cg_typeof("LK")
    };

private:
    int _cgioNum;
    double _id;
    char _label[33];
    char _name[33];
    void *_data;
    size_t _dataSize;
    int _nDims;
    cgsize_t _Dims[12];
    DataType _type;
    std::map<double, CGNSNode *> _childrenByID;
    std::map<std::string, CGNSNode *> _childrenByName;
    //std::multimap<std::string, CGNSNode *> _childrenByLable;

public:
    CGNSNode(int cgio_num, double id);

    /* 
     * @brief 通过ID访问子节点
     * 
     * @param id 
     * @return CGNSNode& 
     */
    CGNSNode &operator[](double id);

    /* 
     * @brief 通过名字访问子节点
     * 
     * @param name 
     * @return CGNSNode& 
     */
    CGNSNode &operator[](const std::string &name);

    std::string name();

    void ReadData();
    void ReleaseData();

    template<class TR>
    TR* Data()
    {
        return reinterpret_cast<TR *>(_data);
    }

    template <class TR>
    TR Data(size_t index)
    {
        return reinterpret_cast<TR *>(_data)[index];
    }

    ~CGNSNode();
};

} // namespace CGIO

#endif