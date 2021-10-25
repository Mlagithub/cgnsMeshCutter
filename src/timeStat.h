#ifndef TIMESTAT_H
#define TIMESTAT_H
#include <iostream>
#include <chrono>
using namespace std;


namespace MeshCut
{

class TimeStat
{
public:
    class Time
    {
    private:
        std::chrono::steady_clock::time_point timePoint_;
    public:
        Time();
        virtual ~Time(){};
        void update();
        double operator-(const Time& t);
    };



public:
    TimeStat();

    double elapsedTime();

    double timeStep();

    //- 返回时间步索引
    int timeIndex() const;

    bool running();

    int operator-(const TimeStat& rhs) const;

    int operator++();

    int operator++(int);

    virtual ~TimeStat(){};

private:
    // System clock timepoint for duration calculation
    Time *startTimePoint_;
    Time *curTimePoint_;

    /**
     * @brief Run Time Control
     * 
     */
    double startTime_;
    double endTime_;
    double timeStep_;
    int timeIndex_;
    int totalTimeIter_;

};



} // namespace MeshCut

#endif
