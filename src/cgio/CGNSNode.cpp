#include "CGNSNode.h"

namespace CGIO
{
CGNSNode::CGNSNode(
    int cgio_num,
    double id) : _cgioNum(cgio_num),
                 _id(id),
                 _label{0},
                 _name{0},
                 _data(nullptr),
                 _type(DataType::MT),
                 _nDims(0),
                 _Dims{0},
                 _childrenByID()
{
    cgio_get_name(cgio_num, id, _name);
    cgio_get_label(cgio_num, id, _label);

    char type[3] = "MT";
    cgio_get_data_type(cgio_num, id, type);
    _type = (DataType)(cg_typeof(type));
    if (_type != DataType::MT || _type != DataType::LK)
    {
        cgio_get_dimensions(cgio_num, id, &_nDims, _Dims);
        _dataSize = 1;
        for (int i = 0; i < _nDims; i++)
        {
            _dataSize *= _Dims[i];
        }
    }

    int childrenNum, numReturned;
    cgio_number_children(cgio_num, id, &childrenNum);
    if (childrenNum > 0)
    {
        double *childrenIDs = new double[childrenNum];
        cgio_children_ids(cgio_num, id, 1, childrenNum, &numReturned, childrenIDs);

        for (int iChild = 0; iChild < childrenNum; iChild++)
        {
            auto result = _childrenByID.insert({childrenIDs[iChild], new CGNSNode(cgio_num, childrenIDs[iChild])});
            _childrenByName.insert({result.first->second->_name, result.first->second});
            //_childrenByLable.insert({result.first->second->_label, result.first->second});
        }
        delete[] childrenIDs;
    }
}

CGNSNode &CGNSNode::operator[](double id)
{
    return *_childrenByID[id];
}

CGNSNode &CGNSNode::operator[](const std::string &name)
{
    return *_childrenByName[name];
}

std::string CGNSNode::name()
{
    return _name;
}

void CGNSNode::ReadData()
{
    if (_type != DataType::MT && _type != DataType::LK && !_data)
    {
        switch (_type)
        {
        case DataType::I4:
            _data = new int32_t[_dataSize];
            break;
        case DataType::I8:
            _data = new int64_t[_dataSize];
            break;
        case DataType::U4:
            _data = new uint32_t[_dataSize];
            break;
        case DataType::U8:
            _data = new uint64_t[_dataSize];
            break;
        case DataType::R4:
            _data = new float[_dataSize];
            break;
        case DataType::R8:
            _data = new double[_dataSize];
            break;
        case DataType::C1:
            _data = new char[_dataSize];
            break;
        case DataType::B1:
            _data = new unsigned char[_dataSize];
            break;
        default:
            break;
        }

#if (CGNS_VERSION <= 3400)
        std::cout << "cgio_read_all_data" << std::endl;
        int tmp = cgio_read_all_data(_cgioNum, _id, _data);
#else
        std::cout << "cgio_read_all_data_type" << std::endl;
        char type[3];
        type[0] = (unsigned short)_type >> 8;
        type[1] = (unsigned short)_type & 0x00ff;
        int tmp = cgio_read_all_data_type(_cgioNum, _id, type, _data);
#endif
        if (tmp != 0)
        {
            std::cout << tmp << std::endl;
            char msg[CGIO_MAX_ERROR_LENGTH + 1];
            cgio_error_message(msg);
            std::cout << msg << std::endl;
        }
    }
}

void CGNSNode::ReleaseData()
{
    if (_data)
    {
        switch (_type)
        {
        case DataType::I4:
            delete[](int32_t *) _data;
            break;
        case DataType::I8:
            delete[](int64_t *) _data;
            break;
        case DataType::U4:
            delete[](uint32_t *) _data;
            break;
        case DataType::U8:
            delete[](uint64_t *) _data;
            break;
        case DataType::R4:
            delete[](float *) _data;
            break;
        case DataType::R8:
            delete[](double *) _data;
            break;
        case DataType::C1:
            delete[](char *) _data;
            break;
        case DataType::B1:
            delete[](unsigned char *) _data;
            break;
        default:
            break;
        }
    }
    _data = nullptr;
}

CGNSNode::~CGNSNode()
{
    ReleaseData();
    for (auto &&i : _childrenByID)
    {
        delete i.second;
    }
}

std::ostream &operator<<(std::ostream &o, CGNSNode::DataType t)
{
    char type[3]{0};
    type[0] = (unsigned short)t >> 8;
    type[1] = (unsigned short)t & 0x00ff;
    return o << type;
}

} // namespace CGIO