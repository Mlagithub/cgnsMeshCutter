#include "timeStat.h"

using namespace std;
MeshCut::TimeStat::Time::Time()
{
    timePoint_ = std::chrono::steady_clock::now();
}

void MeshCut::TimeStat::Time::update()
{
    timePoint_ = std::chrono::steady_clock::now();
}

double MeshCut::TimeStat::Time::operator-(const Time& t)
{
    std::chrono::duration<double> time_run_cost;
    time_run_cost = timePoint_ - t.timePoint_;
    return time_run_cost.count();
}


MeshCut::TimeStat::TimeStat()
    :startTimePoint_(new Time()),
     curTimePoint_(new Time()),
     startTime_(0)
{
    // Directory controlDir = config_.controlDir();
    // timeStep_ = controlDir.lookup("timeStep");
    // endTime_ = controlDir.lookup("totalTime");
    // timeIndex_= startTime_ / timeStep_;
    // totalTimeIter_ = (endTime_ - startTime_) / timeStep_;
}

double MeshCut::TimeStat::elapsedTime()
{
    curTimePoint_->update();
    return *(curTimePoint_)-*(startTimePoint_);
}

double MeshCut::TimeStat::timeStep()
{
    return timeStep_;
}

int MeshCut::TimeStat::timeIndex() const
{
    return timeIndex_;
}

bool MeshCut::TimeStat::running()
{
    return timeIndex_<totalTimeIter_;
}

int MeshCut::TimeStat::operator++()
{
    ++timeIndex_;
    return timeIndex_;
}

int MeshCut::TimeStat::operator++(int)
{
    int t=timeIndex_;
    timeIndex_++;
    return t;
}

int MeshCut::TimeStat::operator-(const TimeStat& rhs) const
{
    return timeIndex_-rhs.timeIndex();
}
