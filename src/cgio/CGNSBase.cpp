#include "CGNSBase.h"
#include "CGIONode.h"
#include "CGError.h"

#include <cgns_io.h>

namespace CGIO
{
    CGNSBase::CGNSBase() : node_ptr(nullptr)
    {
    }

    CGNSBase::CGNSBase(CGIONode &node) : node_ptr(&node)
    {
        int ierr = CG_OK;
        char data_type[CGIO_MAX_DATATYPE_LENGTH+1];
        cgio_get_data_type(node.cgio_num, node.id, data_type);
        ierr = cgio_read_all_data_type(node.cgio_num, node.id, data_type, data);

        auto zone_begin = node.children_by_label.lower_bound("Zone_t");
        auto zone_end = node.children_by_label.upper_bound("Zone_t");
        for (auto izone = zone_begin; izone != zone_end; izone++)
        {
            zones.emplace_back(*this, *(izone->second));
        }
    }

    std::string CGNSBase::GetName() const
    {
        return std::string(node_ptr->name);
    }

    int CGNSBase::nZone()
    {
        return zones.size();
    }
} // namespace CGIO
