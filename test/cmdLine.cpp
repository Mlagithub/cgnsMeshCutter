
#include "cmdLine.h"


#include <iomanip>
#include <sstream>

void CmdLine::RegistHelp(void) {

    Items tmp;
    tmp.argvFlag_ = "-h";
    tmp.argvName_ = "--help";
    tmp.argvDesc_ = "Print this help message.";
    tmp.isMustOffer_ = false;
    tmp.isRadioOpt_  = false;
    args_[tmp.argvFlag_] = tmp;
}

bool CmdLine::isOption(const string &val) {
    for(auto it : args_){
        if (it.first.compare(val) == 0 || it.second.argvName_.compare(val)==0) return true;
    }
    return false;
}

bool CmdLine::isAtChoiceList(const string &optStr, const string &valStr) {
    auto item = args_[optStr];
    if (item.isRadioOpt_) {
        if (item.argvValueChoice_.find(valStr) != item.argvValueChoice_.end()) {
            return true;
        } else
            return false;
    } else
        return true;
}

void CmdLine::parse(int argc, char **argv) {

    this->RegistHelp();

    programName_ = argv[0];

    if (argc > 1 && (string(argv[1]).compare("-h") == 0 || string(argv[1]).compare("--help") == 0)) {
        std::cout << this->usage();
        exit(EXIT_SUCCESS);
    }

    for (int i = 1; i < argc; i += 2) {

        string optStr{argv[i]};
        string valStr;
        // find a registed option
        if (this->isOption(optStr)) {
            // not find value of current opton
            if (i + 1 >= argc ||
                ((valStr = argv[i + 1]), this->isOption(valStr))) {
                throw std::runtime_error(MeshCut::format("Missing value of option: " + optStr +
                                            this->usage(),
                                        __FILE__, __LINE__));
            }
            // find value of current option
            else {
                if (!this->isAtChoiceList(optStr, valStr)) {
                    throw std::runtime_error(MeshCut::format("Wrong value of Option: " + optStr +
                                                this->usage(),
                                            __FILE__, __LINE__));
                } else {
                    args_[optStr].argvValue_ = valStr;
                }
            }
        }
        // not find a registed valid option
        else {
            // throw std::exception("Wrong option: " + optStr +
            //                             this->usage(),
            //                         __FILE__, __LINE__);
        }
    }

    for(auto it : args_){
        if (it.second.isMustOffer_ && it.second.argvValue_.empty()) {
            std::cout << "Missing Option: " + it.first + this->usage();
            exit(EXIT_FAILURE);
        }
    }
}

string CmdLine::usage(void) {

    string rst = "\nUsage: " + programName_ + " [-" + [&](){
        string rst;
        for (auto it : args_) rst += it.first.substr(1);
        return rst;
    }() + "]\n";

    std::stringstream ss;
    for (auto it : args_) {
        auto item = it.second;
        ss << "    " << std::setw(25) << std::left
           << item.argvFlag_ + "," + item.argvName_ << item.argvDesc_
           << ((item.argvValueChoice_.size() > 0) ? ([&]() {
                  string rst = " The accepted options: ";
                  for (auto it : item.argvValueChoice_) rst+= it + ", ";
                  return rst + '\n';
              }())
                                                  : "\n");
    }
    rst += ss.str();

    rst += "\nProgram " + appName_ + "\n\n";
    return rst;
}

namespace Details{
template <>
float help<float>(const string &val) {
    return std::strtof(val.c_str(), nullptr);
}

template <>
double help<double>(const string &val) {
    return std::strtod(val.c_str(), nullptr);
}

template <>
int help<int>(const string &val) {
    return int(help<double>(val));
}

template <>
size_t help<size_t>(const string &val) {
    return size_t(help<double>(val));
}

template <>
string help<string>(const string &val) {
    return val;
}
}