#include "graphclass.h"
#include <cmath>

// -------------------------------------------------
// --- publics -------------------------------------

GraphClass::GraphClass(QString FName) : Manager(FName.toStdString())
{
    TimeHistory.resize(Manager.getNumberOfPopulations());
    TraitHistory.resize(Manager.getNumberOfPopulations());
    maxMembers = 0;
    maxTime = 0;
    jumpedSteps = 1;
    for(size_t i = 0; i < Manager.getNumberOfPopulations(); ++i){
        TimeHistory[i].push_back(0);
        TraitHistory[i].push_back(Manager.getKMembers(i));
    }
}

int GraphClass::generateEvolution(int maxIt)
{
    calcJumpedSteps(maxIt);
    reserveSize(maxIt);
    return generatePointHistory(maxIt);
}

size_t GraphClass::getNumberOfGraphs() const
{
    return Manager.getNumberOfPopulations();
}

double GraphClass::getMaxMembers() const
{
    return maxMembers;
}

double GraphClass::getMaxTime() const
{
    double max = 0;
    for(size_t i = 0; i < Manager.getNumberOfPopulations(); ++i)
        max = std::max(max, TimeHistory[i].back());
    return max;
}

QVector<double> GraphClass::getTimesOf(const int i) const
{
    return TimeHistory.at(i);
//    return TimeHistory[i];
}

QVector<double> GraphClass::getTraitHistOf(const int i) const
{
    return TraitHistory.at(i);
//    return TimeHistory[i];
}

// -------------------------------------------------
// --- privates ------------------------------------

int GraphClass::generatePointHistory(const int maxIt)
{
    int chosen;
    double time = 0;
    int i = 0;
    for(i = 0; i < maxIt; ++i){
        iterateGraphPoint(time, chosen);
        if(Manager.totalEventRate <= 0) break;
    }
    return i*jumpedSteps;
}

void GraphClass::reserveSize(const int maxIt)
{
    for(size_t i = 0; i < Manager.getNumberOfPopulations(); ++i){
        TimeHistory[i].reserve(maxIt);
        TraitHistory[i].reserve(maxIt);
    }
}

void GraphClass::calcJumpedSteps(int & maxIt)
{
    while(true){
        if(maxIt*Manager.getNumberOfPopulations() > pow(10,7)){
            maxIt /= 10;
            jumpedSteps *= 10;
        }
        else
            break;
    }
}

void GraphClass::makeJumpedEvSteps(int & Chosen, double & Time)
{
    for(int i = 0; i < jumpedSteps; ++i){
        Manager.EvolutionStep();
        Time += Manager.event.eventTime;
    }
    Chosen = Manager.event.chosenTrait;
}

void GraphClass::iterateGraphPoint(double & Time, int & Chosen)
{
    makeJumpedEvSteps(Chosen, Time);
    storeCurrentPoint(Time);
    maxMembers = std::max(maxMembers, Manager.getKMembers(Chosen));
}

void GraphClass::storeCurrentPoint(double &Time)
{
    if(Manager.totalEventRate > 0){
        for(size_t i = 0; i < Manager.getNumberOfPopulations(); ++i){
            TimeHistory[i].push_back(Time);
            TraitHistory[i].push_back(Manager.getKMembers(i));
        }
    }
}

