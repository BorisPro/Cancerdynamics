#ifndef TRAITEVENTMANAGER_H
#define TRAITEVENTMANAGER_H

#include <string>
#include <vector>
#include "traitclass.h"
#include "eventclass.h"
#include "randomdice.h"
#include "FileStreamer/ifilestreamer.h"
#include "FileStreamer/cfilestreamer.h"

using std::string;
using std::vector;

class PopulationManager
{
public:
//    vector<size_t> population;  // stores N as N[0], and N1,N2,N3 as the rest.
    size_t populations, tCells, genotypes, phenotypes, TNFsAfterKill;
    double totalEventRate, K;
    vector<vector<double>> competition;
    vector<vector<double>> birthReduction;
    vector<vector<double>> tCellKillRate;
    vector<vector<double>> mutationKernel;
    vector<vector<vector<double>>> naturalSwitch;
    vector<vector<vector<double>>> competitionSwitch;
    vector<TraitClass> traits;
    EventClass event;
    RandomDice dice;

public:
    PopulationManager();
    PopulationManager(string FileName);

    void makeEvolutionStep();
    double getKMembers(int TraitIndex) const;
    size_t getNumberOfPopulations() const;
    void initWithFile(const std::string FName);

private:
    // -----------------------------------------
    // --- Initialisiation ---------------------
    void readPopulations(IFileStreamer& stream);
    void createTraits(IFileStreamer& stream);
    void readDeathBirthProductionRates(IFileStreamer& stream);
    void readCompetition(IFileStreamer& stream);
    void readBirthReduction(IFileStreamer& stream);
    void readNaturalCompetitionSwitchRates(IFileStreamer& stream);
    void readTCellKillRateAndAmount(IFileStreamer& stream);
    void readMutationWithKernel(IFileStreamer& stream);


public: //    All Methods below should be private after testing is over.
    // -----------------------------------------
    // --- Utilities ---------------------------
    size_t getGenotypeOfTraitIndex(size_t index);
    size_t getPhenotypeOfTraitIndex(size_t index);
    size_t getMelanomIndex(size_t genotype, size_t phenotype);

    // -----------------------------------------
    // --- Kernel ------------------------------
    void clearCalculatedRates();
    void addNaturalDeathRates();
    void addCompetitionDeathRates();
    void calculateTraitDeathRates();
    void addNaturalBirthRates();
    void addMutationalBirthRates();
    void subtractBirthReduction();
    void calculateTraitBirthRates();
    void calculateTraitProductionRates();
    void addNaturalSwitchRates();
    void addCompetitionSwitchRates();
    void calculateTraitSwitchRates();
    void calculateTraitKillRates();
    void calculateTraitRates();
    void calculateTotalTraitRates();
    void calculateTotalEventRates();

    void sampleEventTime();

    void choseTraitToChange();
    void choseEventType();
    bool decideIfEventTimeWasHitten(double& HittenEvent, double& ChosenEventRate, EventType type);
    void executeEventTypeOnTrait();
        bool isMutation();
        void choseMutantGenotype();
    void executeBirth();
        double getSpecificSwitchRate(size_t chosenPhenotype, size_t chosenGenotype, size_t switchPhenotype);
        size_t getSwitchedTrait();
    void executeSwitch();
    void executeProduction();
        size_t getKilledTrait();
    void executeKill();
    void changePopulation();
    void applyKToParameters();
};

#endif // TRAITEVENTMANAGER_H
