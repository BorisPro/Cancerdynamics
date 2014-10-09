#ifndef GRAPHCLASS_H
#define GRAPHCLASS_H
#include "traiteventmanager.h"
#include <QVector>
#include <QString>

class GraphClass
{
public:

    GraphClass(QString FName);

    int generateEvolution(int max_It);  // generates Time and TraitHistory, so that plotWindow can access it.

    size_t getNumberOfGraphs() const;
    double getMaxMembers() const;
    double getMaxTime() const;
    QVector<double> getTimesOf(const int i) const;
    QVector<double> getTraitHistOf(const int i) const;

    PopulationManager Manager;
private:

    int generatePointHistory(const int maxIt);
    void calcJumpedSteps(int &maxIt);
    void reserveSize(const int maxIt);

private:
    QVector<QVector<double> > TimeHistory;
    QVector<QVector<double> > TraitHistory;
    double maxMembers, maxTime;
    int jumpedSteps;

    void iterateGraphPoint(double &Time, int &Chosen);
        void makeJumpedEvSteps(int &Chosen, double &Time);
        void storeCurrentPoint(double &Time);

};

#endif // GRAPHCLASS_H
