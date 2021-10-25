#ifndef SECTION_H
#define SECTION_H
#include <string>
#include <unordered_map>
#include <vector>

#include "cgnslib.h"

namespace CGIO
{
    class CGIONode;
    class CGNSBase;
    class Zone;

    class Section
    {
    public:
        Section();
        Section(CGNSBase &base, Zone &zone, CGIONode &node);

    public:
        enum NodeDataType{
            C1 = 0,
            I4 = 1,
            I8 = 2,
            R4 = 3,
            R8 = 4
        };

    public:
        std::string GetName() const;

        int nElement() const;
        int GetElemDataSize() const;
        int GetDataSize(const std::string& nodeName) const;
        int GetDataSize(const std::string& nodeName, const int start, const int end);
        inline ElementType_t GetElementType() const { return static_cast<ElementType_t>(element_type); }
        void ReadElementConnectData(void *buffer) const;
        void ReadElementConnectData(void* buffer, std::vector<cgsize_t>& offset, const int start, const int end);
        bool isBodySection() const;

        NodeDataType initConnectData(void** buffer, const int len);

        int startID();
        int endID();

        std::pair<NodeDataType, std::string> GetDataType(const std::string& nodeName) const;

    private:
        bool isCellSection() const;
        
    private:
        CGNSBase *base_ptr;
        // Zone *zone_ptr;
        CGIONode *node_ptr;

        union
        {
            int data[2];
            struct
            {
                int element_type;
                int element_size_boundary;
            };
        };
        union
        {
            int range[2];
            struct
            {
                int elem_start;
                int elem_end;
            };
        };

        bool isBodySection_ = true;

        void* connecBuffer_;
    };

} // namespace CGIO

#endif