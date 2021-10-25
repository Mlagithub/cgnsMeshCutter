#ifndef CGFILE_H
#define CGFILE_H
#include "CGIONode.h"
#include "CGNSBase.h"
#include "CGEnums.h"

#include <vector>

namespace CGIO
{
    class CGFile
    {
    public:
        CGFile();

        void Open(const char *filename, FileMode mode = FileMode::Read);

        inline CGNSBase& GetDatabase() { return databases[0]; }

        int nDataBase();

        void Close();

        const std::string& name();

    private:
        bool is_opend;
        std::string name_;

        CGIONode root;
        std::vector<CGNSBase> databases;
    };
} // namespace CGIO

#endif