#include "Zone.h"
#include "CGError.h"
#include "CGIONode.h"
#include "CGNSBase.h"

#include <cgns_io.h>

namespace CGIO
{
    Zone::Zone() :
        base_ptr(nullptr),
        node_ptr(nullptr)
    {}
    
    Zone::Zone(CGNSBase& base, CGIONode& node):
        base_ptr(&base),
        node_ptr(&node)
    {
        int ierr = CG_OK;
        char zone_type[12]{0};
        auto& type_node = *(node.children_by_label.find("ZoneType_t")->second);
        char data_type[CGIO_MAX_DATATYPE_LENGTH+1];
        cgio_get_data_type(type_node.cgio_num, type_node.id, data_type);
        ierr = cgio_read_all_data_type(type_node.cgio_num, type_node.id, data_type, zone_type);
        if (*zone_type == 'U')
        {
            type_ = ZoneType::Unstructured;
        }
        else if (*zone_type == 'S')
        {
            type_ = ZoneType::Structured;
        }
        else
        {
            std::string msg = "Unknown zone type ";
            msg.append(zone_type);
            throw msg.c_str();
        }

        cgio_get_data_type(node.cgio_num, node.id, data_type);
        ierr = cgio_read_all_data_type(node.cgio_num, node.id, data_type, data); CheckError(ierr);

        auto section_begin = node.children_by_label.lower_bound("Elements_t");
        auto section_end = node.children_by_label.upper_bound("Elements_t");

        for (auto isection = section_begin; isection != section_end; isection++)
        {
            sections.emplace_back(base, *this, *(isection->second));
            // sections_by_name.insert({isection->second->name, &sections.back()});

            sections_by_name.insert({isection->second->name, new Section(base, *this, *(isection->second))});
        }
    }

    std::string Zone::GetName() const
    {
        return std::string(node_ptr->name);
    }

    int64_t Zone::nCell()
    {
        switch (type_)
        {
        case ZoneType::Unstructured:
            return data[1];
            break;
        case ZoneType::Structured:
        {
            switch (base_ptr->GetCellDim())
            {
            case 2:
                return data[2] * data[3];
                break;
            case 3:
                return data[3] * data[4] * data[5];
            default:
                break;
            }
        }
            break;
        default:
            break;
        }
        return 0;
    }

    int64_t Zone::nVertex()
    {
        switch (type_)
        {
        case ZoneType::Unstructured:
            return data[0];
            break;
        case ZoneType::Structured:
        {
            switch (base_ptr->GetCellDim())
            {
            case 2:
                return data[0] * data[1];
                break;
            case 3:
                return data[0] * data[1] * data[2];
            default:
                break;
            }
        }
        break;
        default:
            break;
        }
        return 0;
    }

    Section& Zone::GetSection(const std::string& name)
    {
        return *sections_by_name[name];
    }

    void Zone::ReadCoordinate(int dim, void *data)
    {
        if (!data) return;
        const char* coordname;
        switch (dim)
        {
        case 0: coordname = "CoordinateX"; break;
        case 1: coordname = "CoordinateY"; break;
        case 2: coordname = "CoordinateZ"; break;
        default: break;
        }
        auto &coord_node = *(node_ptr->children_by_name["GridCoordinates"]->children_by_name[coordname]);
        char data_type[CGIO_MAX_DATATYPE_LENGTH+1];
        cgio_get_data_type(coord_node.cgio_num, coord_node.id, data_type);
        cgio_read_all_data_type(coord_node.cgio_num, coord_node.id, data_type, data);
        return;
    }

    ZoneType Zone::type()
    {
        return type_;
    }
} // namespace CGIO

