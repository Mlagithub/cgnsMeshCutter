#include "CGIONode.h"
#include "CGError.h"

#include <cgns_io.h>

namespace CGIO
{
    CGIONode::CGIONode() :
        cgio_num(0),
        id(0),
        name{0},
        label{0},
        type_in_char{"MT"},
        ndims(0),
        dims{0},
        // children(),
        children_by_id(),
        children_by_name(),
        children_by_label()
    {}

    CGIONode::CGIONode(int cgio_num, double id) :
        cgio_num(cgio_num),
        id(id)
    {
        this->ReadNodeInfo();

        int ierr = CG_OK;
        int children_num, num_ret;
        ierr = cgio_number_children(cgio_num, id, &children_num);
        double* children_ids = new double[children_num];
        ierr = cgio_children_ids(cgio_num, id, 1, children_num, &num_ret, children_ids);
        for (int i = 0; i < children_num; i++)
        {
            // children.emplace_back(new CGIONode(cgio_num, children_ids[i]));
            // auto new_child = children.back();
            auto new_child = new CGIONode(cgio_num, children_ids[i]);
            children_by_id.insert({new_child->id, new_child});
            children_by_name.insert({new_child->name, new_child});
            children_by_label.insert({new_child->label, new_child});
        }
        
    }

    void CGIONode::ReadNodeInfo()
    {
        int ierr = CG_OK;
        ierr = cgio_get_name(cgio_num, id, name); CheckError(ierr);
        ierr = cgio_get_label(cgio_num, id, label); CheckError(ierr);
        ierr = cgio_get_data_type(cgio_num, id, type_in_char); CheckError(ierr);
        ierr = cgio_get_dimensions(cgio_num, id, &ndims, (cgsize_t*)(&dims)); CheckError(ierr);
    }
} // namespace CGIO
