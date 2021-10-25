#ifndef CGENUMS_H
#define CGENUMS_H

#include "cgnslib.h"

namespace CGIO
{
    enum FileMode : int
    {
        Read = 0,
        Write = 1,
        Modify = 2,
        Closed = 3,
    };

    enum FileType : int
    {
        UNKNOWN = 0,
        ADF = 1,
        HDF5 = 2,
        ADF2 = 3,
    };

    using ZoneType = ZoneType_t;
    // enum ZoneType
    // {
    //     Structured,
    //     Unstructured
    // };
} // namespace Act::CGIO

#endif