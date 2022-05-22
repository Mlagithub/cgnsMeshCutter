#include <algorithm>
#include <regex>

#include "stringUtil.h"


#include <regex>

#if __cplusplus >= 201103L &&                             \
    (!defined(__GLIBCXX__) || (__cplusplus >= 201402L) || \
        (defined(_GLIBCXX_REGEX_DFS_QUANTIFIERS_LIMIT) || \
         defined(_GLIBCXX_REGEX_STATE_LIMIT)           || \
             (defined(_GLIBCXX_RELEASE)                && \
             _GLIBCXX_RELEASE > 4)))
#define HAVE_WORKING_REGEX 1
#else
#define HAVE_WORKING_REGEX 0
#endif

const std::string WHITESPACE=" \n\r\t\f\v";

std::string stringJoin(const std::vector<std::string>& strs, std::string&& connector)
{
    std::string rst;
    if(strs.empty()) return rst;
    rst += strs[0];
    for(auto i=1; i<strs.size(); ++i)
    {
        rst += (connector+strs[i]);
    }
    return rst;
};

std::vector<std::string> stringSplit(const std::string& str, const std::string& delimiter)
{
    std::vector<std::string> rst;
    std::string::size_type beg = 0, end = 0;

    while ((end = str.find(delimiter, beg)) != std::string::npos)
    {
        // avoid empty element (caused by doubled slashes)
        if (beg < end) { rst.push_back(str.substr(beg, end - beg)); }
        beg = end + 1;
    }

    // avoid empty trailing element
    if (beg < str.size()) { rst.push_back(str.substr(beg, std::string::npos)); }

    return rst;
}

std::vector<std::string> stringSplit(std::string&& str, const std::string& delimiter)
{
    std::vector<std::string> rst;
    std::string::size_type beg = 0, end = 0;

    while ((end = str.find(delimiter, beg)) != std::string::npos)
    {
        // avoid empty element (caused by doubled slashes)
        if (beg < end) { rst.push_back(str.substr(beg, end - beg)); }
        beg = end + 1;
    }

    // avoid empty trailing element
    if (beg < str.size()) { rst.push_back(str.substr(beg, std::string::npos)); }

    return rst;
}


std::string stringToLower(const string& str)
{
    string rst;
    std::transform(str.begin(), str.end(), std::back_inserter(rst), [](char c) -> char{ return std::tolower(c); });
    
    return rst;
}

std::string stringToLower(string&& str)
{
    string rst;
    std::transform(str.begin(), str.end(), std::back_inserter(rst), [](char c) -> char{ return std::tolower(c); });
    
    return rst;
}

std::string stringTrimLeft(string str)
{
    auto pos = str.find_first_not_of(WHITESPACE);
    return (pos==std::string::npos) ? "" : str.substr(pos);
}

std::string stringTrimRight(string str)
{
    auto pos = str.find_last_not_of(WHITESPACE);
    return (pos==std::string::npos) ? "" : str.substr(0, pos+1);
}

std::string stringTrim(string str)
{
    return stringTrimLeft(stringTrimRight(str));
}

bool stringRegexFind(string rule, string str)
{
#if HAVE_WORKING_REGEX
    std::regex tmp(rule, std::regex_constants::ECMAScript | std::regex_constants::icase);
    return std::regex_search(str, tmp);
#else
    return str.find(rule) != std::string::npos;
#endif
}