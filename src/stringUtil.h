#ifndef STRING_UTIL_H
#define STRING_UTIL_H

#include <string>
#include <vector>

using std::string;
using std::vector;

std::string stringJoin(const std::vector<std::string> &strs, std::string &&connector);
std::vector<std::string> stringSplit(const std::string &str, const std::string &delimiter);
std::vector<std::string> stringSplit(std::string &&str, const std::string &delimiter);
std::string stringToLower(const string &str);
std::string stringToLower(string &&str);
std::string stringTrimLeft(string str);
std::string stringTrimRight(string str);
std::string stringTrim(string str);
bool stringRegexFind(string rule, string str);

#endif