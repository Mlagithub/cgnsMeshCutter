
#ifndef CGNSBASE_H
#define CGNSBASE_H
#include "Zone.h"

#include <vector>

namespace CGIO
{
    class CGIONode;

    class CGNSBase
    {
        friend class CGFile;
    public:
        CGNSBase();
        CGNSBase(CGIONode &node);

        std::string GetName() const;
        inline int GetCellDim() const { return cell_dim; }
        inline int GetPhysicalDim() const { return phy_dim; }

        inline Zone& GetZone() { return zones[0]; }

        int nZone();
        
    private:

        CGIONode* node_ptr;
        union
        {
            int data[2];
            struct
            {
                int cell_dim;
                int phy_dim;
            };
        };
        
        std::vector<Zone> zones;
    };
} // namespace CGIO

#endif