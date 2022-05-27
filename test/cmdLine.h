

#ifndef CMDLINE_H
#define CMDLINE_H

#include "format.h"

#include <exception>
#include <initializer_list>
#include <iostream>
#include <map>
#include <set>
#include <string>
using std::cout;
using std::endl;
using std::initializer_list;
using std::map;
using std::set;
using std::string;

namespace Details {
template <typename T>
T help(const string &val);
} // namespace Details

class CmdLine {
public:
    CmdLine(): appName_("Program") {}
    CmdLine(const string &appName) : appName_(appName){}
    CmdLine(const string &appName, const string &appMsg,
            const string &appVersion)
        : appName_(appName), appMsg_(appMsg), appVersion_(appVersion) {}

public:
    using IsMustOffer = bool;
    static const CmdLine::IsMustOffer MustOffer    = true;
    static const CmdLine::IsMustOffer NotMustOffer = false;

    string usage(void);

    template <typename T, IsMustOffer isMustOffer, typename... Choice>
    void regist(const string &flag, const string &name, const string &desc,
                const string &defaultVal, Choice &&... choice);

    template <typename T>
    T get(const string &optName);

    void parse(int argc, char **argv);

private:
    bool isOption(const string &val);
    bool isAtChoiceList(const string &optStr, const string &valStr);
    void RegistHelp(void);

private:
    struct Items {
        inline bool operator < (const Items& rhs) const noexcept {
            return this->argvFlag_ < rhs.argvFlag_;
        }
        string argvFlag_;
        string argvName_;
        string argvDesc_;
        string argvValue_;
        IsMustOffer isMustOffer_ = false;
        set<string> argvValueChoice_;
        bool isRadioOpt_ = false;
    };

private:
    string programName_;
    string appName_;
    string appMsg_;
    string appVersion_;
    map<string, Items> args_;
    map<string, string> flagNamePair_;
};


template <typename T>
T CmdLine::get(const string &optName) {
    string name;
    if (!this->isOption("-" + optName)) {
        if (!this->isOption("--" + optName)) {
            throw std::runtime_error(MeshCut::format("Can not distinguish option name: " + optName));
        } else {
            name = flagNamePair_["--" + optName];
        }
    } else {
        name = "-" + optName;
    }
    return Details::help<T>(args_[name].argvValue_);
}

template <typename T, CmdLine::IsMustOffer isMustOffer, typename... Choice>
void CmdLine::regist(const string &flag, const string &name, const string &desc,
            const string &defaultVal, Choice &&... choice) {
    Items tmp;
    tmp.argvFlag_    = "-" + flag;
    tmp.argvName_    = "--" + name;
    tmp.argvDesc_    = desc;
    tmp.isMustOffer_ = isMustOffer;
    tmp.argvValue_   = defaultVal;
    (void)std::initializer_list<int>{
        (tmp.argvValueChoice_.insert(choice), 0)...};
    tmp.isRadioOpt_ = ((tmp.argvValueChoice_.empty()) ? false : true);
    args_.insert(std::make_pair(tmp.argvFlag_, tmp));

    flagNamePair_[tmp.argvFlag_] = tmp.argvName_;
    flagNamePair_[tmp.argvName_] = tmp.argvFlag_;
}

#endif