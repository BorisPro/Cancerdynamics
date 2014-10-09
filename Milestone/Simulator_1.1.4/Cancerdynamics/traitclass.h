#ifndef TRAIT_H
#define TRAIT_H

#include <vector>
#include <array>

using std::vector;

class TraitClass
{
public:
    TraitClass();
    ~TraitClass();

public:
    double Members;
    double BirthRate;
    double DeathRate;
//    double NaturalSwitchRate;
    double ProductionRate;
    double ProductionAmount;
    double TraitDeathRate;
    double TraitBirthRate;
    double TraitProductionRate;
    double TraitSwitchRate;
    double TraitKillRate;
    double TraitRate;
    double Mutation;

private:
    /// FIXME: Make getters and setters after moving the statics
};

#endif // TRAIT_H
