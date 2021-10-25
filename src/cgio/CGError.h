#ifndef CGERROR_H
#define CGERROR_H

#include "format.h"
#include "cgnslib.h"

namespace CGIO
{
    // constexpr int CG_OK = 0;
    void CheckError(int error_code);

    template<typename ...Args>
    void CheckError(int errorCode, Args ... args)
    {
        if (errorCode != CG_OK)
        {
            MeshCut::format("Wrong at CGIO:", std::forward<Args>(args)...);
            cg_error_exit();
        }
    }

} // namespace Act::CGIO

#endif