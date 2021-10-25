#ifndef CGIONODE_H
#define CGIONODE_H

#include "cgns_io.h"

#include <cstdint>
#include <list>
#include <map>
#include <string>

namespace CGIO
{
    namespace
    {
        constexpr short operator "" _cgdt(const char * type, size_t size)
        {
            return (short)((type[1] << 8) + type[0]);
        }
    }
    
    enum DataType : short
    {
        // 无数据
        MT = "MT"_cgdt,

        // 32有符号位整数 -- int32_t
        I4 = "I4"_cgdt,

        // 64有符号位整数 -- int64_t
        I8 = "I8"_cgdt,

        // 32位无符号整数 -- uint32_t
        U4 = "U4"_cgdt,

        // 64位无符号整数 -- uint64_t
        U8 = "U8"_cgdt,

        // 单精度浮点数 -- float
        R4 = "R4"_cgdt,

        // 双精度浮点数 -- double
        R8 = "R8"_cgdt,

        // 有符号字符 -- char
        C1 = "C1"_cgdt,

        // 无符号字符 -- unsigned char
        B1 = "B1"_cgdt,

        // 链接
        LK = "LK"_cgdt,
    };

    struct CGIONode
    {
        CGIONode();
        CGIONode(int,double);
        int cgio_num;
        double id;
        char name[33];
        char label[33];
        union
        {
            char type_in_char[3];
            DataType type;
        };
        int ndims;
        int dims[12];

        // std::list<CGIONode*> children;
        std::map<double, CGIONode*> children_by_id;
        std::map<std::string, CGIONode*> children_by_name;
        std::multimap<std::string, CGIONode*> children_by_label;

        inline CGIONode &operator[](const std::string &name) const
        {
            auto iter = children_by_name.find(name);
            if (iter == children_by_name.end()) throw "Cannot find key";
            return *(iter->second);
        }

        inline int data_size() const
        {
            int size = 1;
            for (int i = 0; i < ndims; i++)
            {
                size *= dims[i];
            }
            return size;
        }

        void ReadNodeInfo();
    };
} // namespace CGIO

#endif