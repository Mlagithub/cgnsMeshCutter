#ifndef ZONE_H
#define ZONE_H

#include "CGEnums.h"
#include "Section.h"

#include <vector>
#include <cstdint>
#include <string>
#include <map>

namespace CGIO
{
    class CGIONode;
    class CGNSBase;

    class Zone
    {
        friend class CGNSBase;
    public:
        Zone();
        Zone(CGNSBase& base, CGIONode& node);

        std::string GetName() const;

        int64_t nCell();
        int64_t nVertex();

        inline size_t nSections() const { return sections_by_name.size(); }
        inline Section& GetSection(size_t i) { return sections[i]; }
        Section& GetSection(const std::string& name);

        void ReadCoordinate(int dim, void* data);

        ZoneType type();

    private:
        CGNSBase* base_ptr;
        CGIONode* node_ptr;

        int data[9];

        std::vector<Section> sections;
        std::map<std::string, Section*> sections_by_name;

        ZoneType type_;
    };
} // namespace CGIO

#endif