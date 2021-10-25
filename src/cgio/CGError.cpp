#include "CGError.h"

#include <iostream>
#include <cgns_io.h>


namespace CGIO
{
    void CheckError(int error_code)
    {
        if (error_code != CG_OK)
        {
            cgio_error_exit(nullptr);
        }
    }
} // namespace Act::CGIO
