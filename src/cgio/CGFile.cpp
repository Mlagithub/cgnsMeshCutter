#include "CGFile.h"
#include "CGError.h"

#include <cgns_io.h>
#include <cstring>

namespace CGIO
{
    CGFile::CGFile() :
        is_opend(false)
    {}

    void CGFile::Open(const char* filename, FileMode mode)
    {
        if (is_opend) return;

        int ierr = CG_OK;

        int cgio_num, file_type;
        double root_id;
        ierr = cgio_check_file(filename, &file_type); CheckError(ierr, "wrong file type of:", filename);
        ierr = cgio_open_file(filename, mode, file_type, &cgio_num); CheckError(ierr, "can not open file:", filename);
        is_opend = true;
        this->name_ = filename;

        ierr = cgio_get_root_id(cgio_num, &root_id); CheckError(ierr);

        root = CGIONode(cgio_num, root_id);
        auto base_begin = root.children_by_label.lower_bound("CGNSBase_t");
        auto base_end = root.children_by_label.upper_bound("CGNSBase_t");
        for (auto base = base_begin; base != base_end; base++)
        {
            databases.emplace_back(*(base->second));
        }
    }

    void CGFile::Close()
    {
        if (!is_opend) return;
        int ierr = CG_OK;
        ierr = cgio_close_file(root.cgio_num); CheckError(ierr);
        is_opend = false;
        databases.clear();
    }

    const std::string& CGFile::name()
    {
        return this->name_;
    }

    int CGFile::nDataBase()
    {
        return databases.size();
    }

} // namespace CGIO
