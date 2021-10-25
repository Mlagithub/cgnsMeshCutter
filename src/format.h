#ifndef FORMAT_H
#define FORMAT_H

#include <memory>
#include <stdexcept>
#include <string>

namespace MeshCut
{

template <typename... Args>
std::string format(const std::string &fmtStr, Args... args)
{
    size_t len = snprintf(nullptr, 0, fmtStr.c_str(), args...) + 1;
    if (len <= 0) { throw std::runtime_error("Error during formatting."); }
    std::unique_ptr<char[]> buf(new char[len]);
    snprintf(buf.get(), len, fmtStr.c_str(), args...);
    return std::string(buf.get(), buf.get() + len - 1);
}
} // namespace MeshCut

#endif