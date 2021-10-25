#include "Section.h"
#include "CGError.h"
#include "CGIONode.h"
#include "CGNSBase.h"
#include "Zone.h"
#include "format.h"

#include "fmt/format.h"


#include <cstring>
#include <cgns_io.h>
#include <cgnslib.h>

using MeshCut::format;

namespace CGIO
{

    std::unordered_map<ElementType_t, int> NodeCountOfElement{
        {ElementType_t::HEXA_8,8},
        {ElementType_t::TETRA_4,4},
        {ElementType_t::TRI_3,3},
        {ElementType_t::PENTA_6,6},
        {ElementType_t::QUAD_4,4},
        {ElementType_t::PYRA_5,5},
    };

    std::unordered_map<int, ElementType_t> ElementTypeFlag{
        {17, ElementType_t::HEXA_8},
        {10, ElementType_t::TETRA_4},
        {5,  ElementType_t::TRI_3},
        {14, ElementType_t::PENTA_6},
        {7,  ElementType_t::QUAD_4},
        {12, ElementType_t::PYRA_5},
        {20, ElementType_t::MIXED},
    };


    Section::Section() :
        base_ptr(nullptr),
        // zone_ptr(nullptr),
        node_ptr(nullptr)
    {
    }

    Section::Section(CGNSBase &base, Zone &zone, CGIONode &node) :
        base_ptr(&base),
        // zone_ptr(&zone),
        node_ptr(&(node))
    {
        int ierr = CG_OK;
        char data_type[CGIO_MAX_DATATYPE_LENGTH+1];
        cgio_get_data_type(node.cgio_num, node.id, data_type);
        ierr = cgio_read_all_data_type(node.cgio_num, node.id, data_type, data);

        auto &range_node = *(node.children_by_name.find("ElementRange")->second);
        cgio_get_data_type(range_node.cgio_num, range_node.id, data_type);
        ierr = cgio_read_all_data_type(range_node.cgio_num, range_node.id, data_type, range);

        isBodySection_ = this->isCellSection();
    }

    std::string Section::GetName() const
    {
        return std::string(node_ptr->name);
    }

    int Section::nElement() const
    {
        return elem_end - elem_start + 1;
    }

    int Section::GetElemDataSize() const
    {
        return node_ptr->children_by_name["ElementConnectivity"]->data_size();
    }

    int Section::GetDataSize(const std::string& nodeName) const
    {
        return node_ptr->children_by_name[nodeName.c_str()]->data_size();
    }

    int Section::GetDataSize(const std::string& nodeName, const int start, const int end)
    {
        // constexpr char* nodeName = "ElementConnectivity";
        auto& conn_node = *(node_ptr->children_by_name[nodeName]);

        // check range 
        if(start < this->elem_start || end > this->elem_end)
        {
            CheckError(CG_ERROR, format("block range [%d,%d] out of range [%d,%d]\n", start, end, elem_start, elem_end));
        }

        auto len = this->GetDataSize(nodeName);
        auto type = this->GetDataType(nodeName);
        switch (type.first)
        {
        case I4:
            connecBuffer_ = new int[len];
            break;
        case I8:
            connecBuffer_ = new cglong_t[len];
            break;
        default:
            CheckError(CG_ERROR, format("Data type %s not supported\n", type.second));
            break;
        }

        CheckError(cgio_read_all_data_type(conn_node.cgio_num, conn_node.id, type.second.c_str(), connecBuffer_), format("cgio_read_all_data_type() of node: %s\n", conn_node.name));

        int dataSize = 0;
        if(this->element_type == static_cast<int>(ElementType_t::MIXED))
        {
            int i=0;
            for(auto n=this->elem_start; n <=this->elem_end; ++n)
            {
                if(n>end) break;
                auto nodeCnt = NodeCountOfElement.at(ElementTypeFlag.at(((int*)connecBuffer_)[i]));
                i += (nodeCnt+1);
                if(n>=start)
                {
                    dataSize += (nodeCnt+1);
                }
            }
        }
        else
        {
            dataSize = NodeCountOfElement.at(ElementTypeFlag.at(element_type)) * (end - start + 1);
        }

        return dataSize;
    }

    Section::NodeDataType Section::initConnectData(void** buffer, const int len)
    {
        auto type = this->GetDataType("ElementConnectivity");
        switch (type.first)
        {
        case I4:
            *buffer = new int[len];
            break;
        case I8:
            *buffer = new cglong_t[len];
            break;
        default:
            fmt::format(format("Data type %s not supported for ElementRange\n", type.second));
            break;
        }
        return type.first;
    }

    std::pair<CGIO::Section::NodeDataType, std::string> Section::GetDataType(const std::string& nodeName) const
    {
        auto& node = *(node_ptr->children_by_name[nodeName.c_str()]);
        char data_type[CGIO_MAX_DATATYPE_LENGTH+1];
        cgio_get_data_type(node.cgio_num, node.id, data_type);

        if (0 == std::strcmp(data_type, "I8"))
        {
            return {NodeDataType::I8, data_type};
        }
        else if (0 == std::strcmp(data_type, "I4"))
        {
            return {NodeDataType::I4, data_type};
        }
        else if (0 == std::strcmp(data_type, "R8"))
        {
            return {NodeDataType::R8, data_type};
        }
        else if (0 == std::strcmp(data_type, "R4"))
        {
            return {NodeDataType::R4, data_type};
        }
        else if (0 == std::strcmp(data_type, "C1"))
        {
            return {NodeDataType::C1, data_type};
        }
        else
        {
            fmt::format(format("Data type %s not supported for %s\n", data_type, nodeName));
        }
    }


    void Section::ReadElementConnectData(void* buffer) const
    {
        auto& conn_node = *(node_ptr->children_by_name["ElementConnectivity"]);
        int ierr = CG_OK;
        char data_type[CGIO_MAX_DATATYPE_LENGTH+1];
        cgio_get_data_type(conn_node.cgio_num, conn_node.id, data_type);
        ierr = cgio_read_all_data_type(conn_node.cgio_num, conn_node.id, data_type, buffer); CheckError(ierr);
        return;
    }

    // TODO: cgns-develop branch(time at 2021-08-01) has bug with function cgio_read_block_data_type()
    // so get all data with cgio_read_all_data_type() and manually slice range of [start,end].
    void Section::ReadElementConnectData(void* buffer, std::vector<cgsize_t>& offset, const int start, const int end)
    {
        if(connecBuffer_ == nullptr){
            this->GetDataSize("ElementConnectivity", start, end);
        }

        auto type = this->GetDataType("ElementConnectivity");

        cgsize_t id_blk = 0;
        cgsize_t id_buf = 0;
        cgsize_t cnt = 0;

        if(this->element_type == static_cast<int>(ElementType_t::MIXED))
        {
            for(auto n=this->elem_start; n <=this->elem_end; ++n)
            {
                if(n>end) break;

                auto nodeCnt = NodeCountOfElement.at(ElementTypeFlag.at(((int*)connecBuffer_)[id_buf]));
                if(n<start)
                {
                    id_buf += (nodeCnt+1);
                    continue;
                }
                else
                {   
                    offset[cnt+1] = offset[cnt] + nodeCnt + 1;
                    cnt++;
                    switch (type.first)
                    {
                    case I4:
                        ((int*)buffer)[id_blk++] = (((int*)connecBuffer_)[id_buf++]);
                        for(int j=0; j<nodeCnt; ++j)
                        {
                            ((int*)buffer)[id_blk++] = (((int*)connecBuffer_)[id_buf++]);
                        }
                        break;
                    case I8:
                        ((int*)buffer)[id_blk++] = (((cglong_t*)connecBuffer_)[id_buf++]);

                        for(int j=0; j<nodeCnt; ++j)
                        {
                            ((int*)buffer)[id_blk++] = (((cglong_t*)connecBuffer_)[id_buf++]);
                        }
                        break;
                    default:
                        break;
                    }
                }
            }
        }
        else
        {
            auto nodeCnt = NodeCountOfElement.at(ElementTypeFlag.at(element_type));
            for(auto n=this->elem_start; n <=this->elem_end; ++n)
            {
                if(n>end) break;

                if(n<start)
                {
                    id_buf += nodeCnt;
                    continue;
                }
                else
                {   
                    offset[cnt+1] = offset[cnt] + nodeCnt;
                    cnt++;
                    switch (type.first)
                    {
                    case I4:
                        for(int j=0; j<nodeCnt; ++j)
                        {
                            ((int*)buffer)[id_blk++] = (((int*)connecBuffer_)[id_buf++]);
                        }
                        break;
                    case I8:
                        for(int j=0; j<nodeCnt; ++j)
                        {
                            ((int*)buffer)[id_blk++] = (((cglong_t*)connecBuffer_)[id_buf++]);
                        }
                        break;
                    default:
                        break;
                    }
                }
            }
        }

        return;
    }

    bool Section::isBodySection() const
    {
        return isBodySection_;
    }

    // TODO: only support 3-dimension mesh 
    bool Section::isCellSection() const
    {
        std::vector<int> bodyConn;
        int cell_dim;

        switch (element_type)
        {
        case NODE:
            cell_dim = 0;
            break;
        case BAR_2:
            cell_dim = 1;
        case TRI_3:
        case QUAD_4:
            cell_dim = 2;
            break;
        case TETRA_4:
        case PYRA_5:
        case PENTA_6:
        case HEXA_8:
            cell_dim = 3;
            break;
        case MIXED:
            bodyConn.assign(this->GetElemDataSize(), 0);
            this->ReadElementConnectData(bodyConn.data());

            switch(bodyConn[0])
            {
            // 2-D
            case ElementType_t::TRI_3:
            case ElementType_t::TRI_6:
            case ElementType_t::TRI_9:
            case ElementType_t::TRI_10:
            case ElementType_t::TRI_12:
            case ElementType_t::TRI_15:
            case ElementType_t::QUAD_4:
            case ElementType_t::QUAD_8:
            case ElementType_t::QUAD_9:
            case ElementType_t::QUAD_12:
            case ElementType_t::QUAD_16:
            case ElementType_t::QUAD_P4_16:
            case ElementType_t::QUAD_25:
                cell_dim = 2;
                break;
            default:
                cell_dim = 3;
                break;
            }
            break;
        default:
            // throw "Unsupported element type";
            fmt::format("Unsupported element type\n");
            break;
        }
        return cell_dim == base_ptr->GetCellDim();
    }

    int Section::startID()
    {
        return elem_start;
    }
    int Section::endID()
    {
        return elem_end;
    }

} // namespace CGIO
