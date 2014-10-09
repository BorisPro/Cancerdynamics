#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "traiteventmanager.h"

using std::string;
using std::vector;

PopulationManager::PopulationManager():
    traits(),
    event(),
    dice()
{

}

PopulationManager::PopulationManager(string FileName):
    traits(),
    event(),
    dice()
{
    initWithFile(FileName);
}

void PopulationManager::applyKToParameters()
{
    for(TraitClass& trait : traits){ // Members and ProductionRates
        trait.Members *= K;
        trait.ProductionRate /= K;
    }
    for(size_t i = 0; i < populations; ++i){    // competition death
        for(size_t j = 0; j < populations; ++j)
            competition[i][j] /= K;
    }
    for(size_t i = 0; i < populations; ++i){    // birth reduction
        for(size_t j = 0; j < populations; ++j)
            birthReduction[i][j] /= K;
    }
    for(size_t k = 0; k < genotypes; ++k){      // competition switch
        for(size_t i = 0; i < phenotypes; ++i){
            for(size_t j = 0; j < phenotypes; ++j)
                competitionSwitch[k][i][j] /= K;
        }
    }
    for(size_t i = 0; i < tCells; ++i){         // killRate
        for(size_t j = 0; j < phenotypes * genotypes; ++j)
            tCellKillRate[i][j] /= K;
    }
}

void PopulationManager::initWithFile(const string FName)
{
    CFileStreamer object(FName);
    IFileStreamer& stream(object);
    readPopulations(stream);
    createTraits(stream);
    readDeathBirthProductionRates(stream);
    readCompetition(stream);
    readBirthReduction(stream);
    readNaturalCompetitionSwitchRates(stream);
    readTCellKillRateAndAmount(stream);
    readMutationWithKernel(stream);
    K = std::stod(stream.getNextWord());
    applyKToParameters();
}

void PopulationManager::EvolutionStep()
{
    calculateTotalEventRates();
    sampleEventTime();
    changePopulation();
}

double PopulationManager::getKMembers(int TraitIndex) const
{
    return traits[TraitIndex].Members / K;
}

size_t PopulationManager::getNumberOfPopulations() const
{
    return populations;
}

void PopulationManager::readPopulations(IFileStreamer &stream)
{
    tCells = stoi(stream.getNextWord());
    genotypes = stoi(stream.getNextWord());
    phenotypes = stoi(stream.getNextWord());
    populations = 1 + tCells + genotypes * phenotypes;
}

void PopulationManager::createTraits(IFileStreamer& stream)
{
    traits.resize(populations);
    for(size_t i = 0; i < populations; ++i){
        traits[i].Members = std::stod(stream.getNextWord());
    }
}

void PopulationManager::readDeathBirthProductionRates(IFileStreamer& stream)
{
    for(TraitClass& cells : traits)
        cells.DeathRate = std::stod(stream.getNextWord());
    for(TraitClass& cells : traits)
        cells.BirthRate = std::stod(stream.getNextWord());
    for(size_t i = 1; i <= tCells; ++i) // andere Produktionsraten wurden bereits mit 0 initialisiert.
        traits[i].ProductionRate = std::stod(stream.getNextWord());
    for(size_t i = 1; i <= tCells; ++i) // andere Produktionsmengen wurden bereits mit 0 initialisiert.
        traits[i].ProductionAmount = std::stod(stream.getNextWord());
}

void PopulationManager::readCompetition(IFileStreamer& stream)
{
    competition.resize(populations);
    for(size_t i = 0; i < populations; ++i){
        competition[i].resize(populations);
        for(size_t j = 0; j < populations; ++j)
            competition[i][j] = std::stod(stream.getNextWord());
    }
}

void PopulationManager::readBirthReduction(IFileStreamer& stream)
{
    birthReduction.resize(populations);
    for(size_t i = 0; i < populations; ++i){
        birthReduction[i].resize(populations);
        for(size_t j = 0; j < populations; ++j)
            birthReduction[i][j] = std::stod(stream.getNextWord());
    }
}

void PopulationManager::readNaturalCompetitionSwitchRates(IFileStreamer& stream)
{
    naturalSwitch.resize(genotypes);
    for(size_t k = 0; k < genotypes; ++k){
        naturalSwitch[k].resize(phenotypes);
        for(size_t i = 0; i < phenotypes; ++i){
            naturalSwitch[k][i].resize(phenotypes);
            for(size_t j = 0; j < phenotypes; ++j)
                naturalSwitch[k][i][j] = std::stod(stream.getNextWord());
        }
    }
    competitionSwitch.resize(genotypes);
    for(size_t k = 0; k < genotypes; ++k){
        competitionSwitch[k].resize(phenotypes);
        for(size_t i = 0; i < phenotypes; ++i){
            competitionSwitch[k][i].resize(phenotypes);
            for(size_t j = 0; j < phenotypes; ++j)
                competitionSwitch[k][i][j] = std::stod(stream.getNextWord());
        }
    }
}

void PopulationManager::readTCellKillRateAndAmount(IFileStreamer& stream)
{
    tCellKillRate.resize(tCells);
    for(size_t i = 0; i < tCells; ++i){
        tCellKillRate[i].resize(phenotypes * genotypes);
        for(size_t j = 0; j < phenotypes * genotypes; ++j)
            tCellKillRate[i][j] = std::stod(stream.getNextWord());
    }
    TNFsAfterKill = std::stoi(stream.getNextWord());
}

void PopulationManager::readMutationWithKernel(IFileStreamer& stream)
{
    for(size_t i = 0; i < genotypes; ++i){
        double mu = std::stod(stream.getNextWord());
        for(size_t j = 0; j < phenotypes; ++j)
            traits[getMelanomIndex(i,j)].Mutation = mu;
    }
    mutationKernel.resize(genotypes);
    for(size_t i = 0; i < genotypes; ++i){
        mutationKernel[i].resize(genotypes);
        for(size_t j = 0; j < genotypes; ++j)
            mutationKernel[i][j] = std::stod(stream.getNextWord());
    }
}

// -----------------------------------------------
// --- Utilities ---------------------------------

size_t PopulationManager::getGenotypeOfTraitIndex(size_t index)
{
    int genotype = index - 1 - tCells;
    if(genotype < 0) throw (string) "getGenotypeOfTraitIndex hat zu kleinen Index erhalten!";
    genotype = genotype / phenotypes;
    return genotype;
}

size_t PopulationManager::getMelanomIndex(size_t genotype, size_t phenotype)
{
    return 1 + tCells + genotype * phenotypes + phenotype;
}

size_t PopulationManager::getPhenotypeOfTraitIndex(size_t index)
{
    int genotype = index - 1 - tCells;
    if(genotype < 0) throw (string) "getPhenotypeOfTraitIndex hat zu kleinen Index erhalten!";
    genotype = genotype % phenotypes;
    return genotype;
}

// -----------------------------------------------
// --- Kernel ------------------------------------

void PopulationManager::clearCalculatedRates()
{
    for(TraitClass& trait : traits){
        trait.TraitDeathRate = 0;
        trait.TraitBirthRate = 0;
        trait.TraitProductionRate = 0;
//        trait.TraitSwitchRate.assign();
        trait.TraitRate = 0;
    }
}

void PopulationManager::addNaturalDeathRates()
{
    for(TraitClass& trait : traits)
        trait.TraitDeathRate += trait.Members * trait.DeathRate;
}

void PopulationManager::addCompetitionDeathRates()
{
    double sum;
    for(size_t i = 0; i < populations; ++i){
        sum = 0;
        for(size_t j = 0; j < populations; ++j){
            sum += competition[i][j] * traits[j].Members;
        }
        traits[i].TraitDeathRate += sum * traits[i].Members;
    }
}

void PopulationManager::calculateTraitDeathRates()
{
    for(TraitClass& trait : traits)
        trait.TraitDeathRate = 0;
    addNaturalDeathRates();
    addCompetitionDeathRates();
}

void PopulationManager::addNaturalBirthRates()
{
    for(TraitClass& trait : traits)
        trait.TraitBirthRate += trait.Members * trait.BirthRate /* * (1-trait.Mutation)*/;
}

void PopulationManager::addMutationalBirthRates()
{
    for(size_t i = 0; i < genotypes; ++i){
        for(size_t j = 0; j < genotypes; ++j){  // iterating through genotype mutation partners.
            for(size_t k = 0; k < phenotypes; ++k){ // considering each phenotype for mutation.
                traits[getMelanomIndex(i,k)].TraitBirthRate +=
                                            traits[getMelanomIndex(j,k)].Members * traits[getMelanomIndex(j,k)].BirthRate
                                            * traits[getMelanomIndex(j,k)].Mutation * mutationKernel[j][i];
            }
        }
    }
}

void PopulationManager::subtractBirthReduction()
{
    double sum;
    for(size_t i = 0; i < populations; ++i){
        sum = 0;
        for(size_t j = 0; j < populations; ++j){
            sum += birthReduction[i][j] * traits[j].Members;
        }
        traits[i].TraitBirthRate -= sum * traits[i].Members;
    }
}

void PopulationManager::calculateTraitBirthRates()
{
    for(TraitClass& trait : traits)
        trait.TraitBirthRate = 0;
    addNaturalBirthRates();
//    addMutationalBirthRates();
    subtractBirthReduction();
}

void PopulationManager::calculateTraitProductionRates()
{
    vector<double> melCells(tCells,0);
    for(size_t i = 0; i < tCells; ++i){
        for(size_t j = 0; j < genotypes * phenotypes; ++j)
            if(tCellKillRate[i][j] > 0)
                melCells[i] += traits[1 + tCells + j].Members;
    }
    for(size_t i = 1; i <= tCells; ++i){
        traits[i].TraitProductionRate = traits[i].Members * melCells[i-1] * traits[i].ProductionRate;
    }

//    double melanomCells = 0;
//    for(size_t i = 1 + tCells; i <= populations; ++i)
//        melanomCells += traits[i].Members;
//    for(size_t i = 1; i <= tCells; ++i){
//        traits[i].TraitProductionRate = traits[i].Members * melanomCells * traits[i].ProductionRate;
//    }
}

void PopulationManager::addNaturalSwitchRates()
{
    double switchrate;
    for(size_t k = 0; k < genotypes; ++k){
        for(size_t i = 0; i < phenotypes; ++i){
            switchrate = 0;
            for(size_t j = 0; j < phenotypes; ++j)  // Improve: calculate once!
                switchrate += naturalSwitch[k][i][j];
            switchrate *= traits[getMelanomIndex(k,i)].Members;
            traits[getMelanomIndex(k,i)].TraitSwitchRate += switchrate;
        }
    }
}

void PopulationManager::addCompetitionSwitchRates()
{
    size_t traitNumber;
    double switchrate;
    for(size_t k = 0; k < genotypes; ++k){
        traitNumber = k * phenotypes;
        for(size_t i = 0; i < phenotypes; ++i){
            traitNumber += i;
            switchrate = 0;
            for(size_t j = 0; j < phenotypes; ++j)
                switchrate += competitionSwitch[k][i][j];
            switchrate *= traits[getMelanomIndex(k,i)].Members * traits[0].Members;
            traits[getMelanomIndex(k,i)].TraitSwitchRate += switchrate;
        }
    }
}

void PopulationManager::calculateTraitSwitchRates()
{
    for(TraitClass& trait : traits)
        trait.TraitSwitchRate = 0;
    addNaturalSwitchRates();
    addCompetitionSwitchRates();
}

void PopulationManager::calculateTraitKillRates()
{
    for(size_t i = 0; i < tCells; ++i){
        traits[1 + i].TraitKillRate = 0;
        for(size_t j = 0; j < genotypes * phenotypes; ++j){
            traits[1 + i].TraitKillRate += tCellKillRate[i][j] * traits[1 + tCells + j].Members;
        }
        traits[1 + i].TraitKillRate *= traits[1 + i].Members;
    }
}

void PopulationManager::calculateTraitRates()
{
    for(TraitClass& trait: traits){
        trait.TraitRate = 0;
        trait.TraitRate += trait.TraitDeathRate;
        trait.TraitRate += trait.TraitBirthRate;
        trait.TraitRate += trait.TraitProductionRate;
        trait.TraitRate += trait.TraitSwitchRate;
        trait.TraitRate += trait.TraitKillRate;
    }
}

void PopulationManager::calculateTotalTraitRates()
{
    calculateTraitDeathRates();
    calculateTraitBirthRates();
    calculateTraitProductionRates();
    calculateTraitSwitchRates();
    calculateTraitKillRates();
    calculateTraitRates();
}

void PopulationManager::calculateTotalEventRates()
{
    calculateTotalTraitRates();
    totalEventRate = 0;
    for(TraitClass trait: traits)
        totalEventRate += trait.TraitRate;
}

//void TraitEventManager::addTotalCompDeathRateOf(int i)
//{
//    double extDeath = 0;
//    for(int j = 0; j < TraitClass::Size; j++)
//        extDeath += TraitClass::CompDeathRate.at(i).at(j) * Trait[j].Members;
//    extDeath *= Trait[i].Members;
//    Trait[i].TotalDeathRate += extDeath;
//}

//void TraitEventManager::addTotalIntrisicDeathRateOf(int i)
//{
//    Trait[i].TotalDeathRate += (Trait[i].DeathRate)*(Trait[i].Members);
//}

//void TraitEventManager::calculateTotalDeathRateOf(int TraitIndex)
//{
//    Trait[TraitIndex].TotalDeathRate = 0;
//    addTotalCompDeathRateOf(TraitIndex);
//    addTotalIntrisicDeathRateOf(TraitIndex);
//}

//void TraitEventManager::calculateTotalDeathRates()
//{
//    for(int k = 0; k < TraitClass::Size; ++k){
//        Trait[k].TotalDeathRate = 0;
//        addTotalCompDeathRateOf(k);
//        addTotalIntrisicDeathRateOf(k);
//    }
//}

//void TraitEventManager::calculateTotalBirthRates()
//{
//    int n = TraitClass::Size;
//    for(int k = 0; k < n; ++k)
//        Trait[k].TotalBirthRate = (Trait[k].Members)*(Trait[k].BirthRate)*(1-TraitClass::Mutation);
//    for(int k = 1; k < n-1; ++k){
//        Trait[k].TotalBirthRate += 0.5*(Trait[k-1].Members)*(Trait[k-1].BirthRate)*(TraitClass::Mutation);
//        Trait[k].TotalBirthRate += 0.5*(Trait[k+1].Members)*(Trait[k+1].BirthRate)*(TraitClass::Mutation);
//    }
//    Trait[n-1].TotalBirthRate += 0.5*(Trait[n-2].Members)*(Trait[n-2].BirthRate)*(TraitClass::Mutation);
//    Trait[0].TotalBirthRate += 0.5*(Trait[1].Members)*(Trait[1].BirthRate)*(TraitClass::Mutation);
//}


//void TraitEventManager::calculateTotalEventRate()
//{
//    TraitClass::TotalEventRate = 0;
//    for(int i = 0; i < TraitClass::Size; i++){
//        Trait[i].TotalTraitRate = Trait[i].TotalBirthRate
//                                    + Trait[i].TotalDeathRate;
//        TraitClass::TotalEventRate += Trait[i].TotalTraitRate;
//    }
//}

//void TraitEventManager::calculateEventRates()
//{
//    calculateTotalDeathRates();
//    calculateTotalBirthRates();
//    calculateTotalEventRate();
//}

void PopulationManager::sampleEventTime()
{
    event.eventTime = dice.rollExpDist(totalEventRate);
}

void PopulationManager::choseTraitToChange()
{
    if(totalEventRate == 0){event.chosenTrait = 0; return;}
    double HittenTrait = dice.rollContUnifDist(totalEventRate);
    for(size_t i = 0; i < populations; i++){
        if(HittenTrait <= traits[i].TraitRate){
            event.chosenTrait = i;
//            if(i == 0 && traits[i].Members == 0)
//                return;
            return;
        }
        HittenTrait -= traits[i].TraitRate;
    }
}

void PopulationManager::choseEventType()
{
    size_t i = event.chosenTrait;
    double HittenEvent = dice.rollContUnifDist(traits[i].TraitRate);

    if(decideIfEventTimeWasHitten(HittenEvent,traits[i].TraitDeathRate,DEATH)){}
    else if(decideIfEventTimeWasHitten(HittenEvent,traits[i].TraitBirthRate,BIRTH)){}
    else if(decideIfEventTimeWasHitten(HittenEvent,traits[i].TraitSwitchRate,SWITCH)){}
    else if(decideIfEventTimeWasHitten(HittenEvent,traits[i].TraitProductionRate,PRODUCTION)){}
    else if(decideIfEventTimeWasHitten(HittenEvent,traits[i].TraitKillRate,KILL)){}
    else if(totalEventRate == 0){event.type = NONE;}
    else throw (string) "No Event was hit! Trait rate roll was bigger than the trait rate!";
}

bool PopulationManager::decideIfEventTimeWasHitten(double &HittenEvent, double &ChosenEventRate, EventType type)
{
    if(HittenEvent < ChosenEventRate){
        event.type = type;
        return true;
    }
    HittenEvent -= ChosenEventRate;
    return false;
}

void PopulationManager::executeEventTypeOnTrait()
{
    if(event.type == DEATH){traits[event.chosenTrait].Members -= 1.;}
    else if(event.type == BIRTH){executeBirth();}  // mutation.
    else if(event.type == SWITCH){executeSwitch();} // switchmatrix (look what swiches can be made).
    else if(event.type == PRODUCTION){executeProduction();} // produce rx.
    else if(event.type == KILL){executeKill();}   // kill and produce tnf after some time.
    else if(event.type == NONE){}
    else throw (string) "No event type has been chosen before: executeEventTypeOnTrait()";
}

bool PopulationManager::isMutation()
{
    double roll = dice.rollContUnifDist(1);
    if(roll < traits[event.chosenTrait].Mutation)
        return true;
    else
        return false;
//    return roll < traits[event.chosenTrait].Mutation ? true : false;
}

void PopulationManager::choseMutantGenotype()
{
    size_t chosenGenotype = getGenotypeOfTraitIndex(event.chosenTrait);
    size_t chosenPhenotype = getPhenotypeOfTraitIndex(event.chosenTrait);
    size_t mutantGenotype = 0;
    double HittenTrait = dice.rollContUnifDist(1);

    for(mutantGenotype = 0; mutantGenotype < genotypes; ++mutantGenotype){  // determine the new genotype
        if(HittenTrait <= mutationKernel[chosenGenotype][mutantGenotype])
            break;
        HittenTrait -= mutationKernel[chosenGenotype][mutantGenotype];
    }
    event.chosenTrait = getMelanomIndex(mutantGenotype, chosenPhenotype); // adding the usual shift and Phenotype
}

void PopulationManager::executeBirth()
{
    if(event.chosenTrait > 0 && event.chosenTrait < 1 + tCells)
        executeProduction();
    else if(isMutation()){
        choseMutantGenotype();
        traits[event.chosenTrait].Members++;
    }
    else{
        traits[event.chosenTrait].Members++;
    }
}

double PopulationManager::getSpecificSwitchRate(size_t chosenPhenotype, size_t chosenGenotype, size_t switchPhenotype)
{
    double currentSwitchRate = 0;
    currentSwitchRate += naturalSwitch[chosenGenotype][chosenPhenotype][switchPhenotype];
    currentSwitchRate += competitionSwitch[chosenGenotype][chosenPhenotype][switchPhenotype]*traits[0].Members;
    return currentSwitchRate;
}

size_t PopulationManager::getSwitchedTrait()
{
    double HittenSwitch = dice.rollContUnifDist(traits[event.chosenTrait].TraitSwitchRate/traits[event.chosenTrait].Members);
    size_t chosenGenotype = getGenotypeOfTraitIndex(event.chosenTrait);
    size_t chosenPhenotype = getPhenotypeOfTraitIndex(event.chosenTrait);
    for(size_t switchPhenotype = 0; switchPhenotype < phenotypes; ++switchPhenotype){
        double currentSwitchRate = getSpecificSwitchRate(chosenPhenotype, chosenGenotype, switchPhenotype);
        if(HittenSwitch <= currentSwitchRate)
            return getMelanomIndex(chosenGenotype, switchPhenotype);
        HittenSwitch -= currentSwitchRate;
    }
    throw (string) "getSwitchedTrait() could not find a matching switch!";
    return 0;
}

void PopulationManager::executeSwitch()
{
    size_t switchedTrait = getSwitchedTrait();
    traits[event.chosenTrait].Members--;
    traits[switchedTrait].Members++;
}

void PopulationManager::executeProduction()
{
    traits[0].Members += traits[event.chosenTrait].ProductionAmount;
    traits[event.chosenTrait].Members++;
}

size_t PopulationManager::getKilledTrait()
{
    double HittenCell = dice.rollContUnifDist(traits[event.chosenTrait].TraitKillRate/traits[event.chosenTrait].Members);
    size_t tcell = event.chosenTrait -1;
    for(size_t i = 0; i < genotypes*phenotypes; ++i){
        if(HittenCell < tCellKillRate[tcell][i] * traits[1+ tCells + i].Members)
            return 1+tCells+i;
        HittenCell -= tCellKillRate[tcell][i] * traits[1+ tCells + i].Members;
    }
    throw (string) "no killed trait found!";
    return 0;
}

void PopulationManager::executeKill()
{
    traits[getKilledTrait()].Members--;
    traits[0].Members += TNFsAfterKill;
}

void PopulationManager::changePopulation()
{
    if(totalEventRate > 0){
        choseTraitToChange();
        choseEventType();
        executeEventTypeOnTrait();
    }
}


//void TraitEventManager::choseEventType()
//{
//    int i = Events.ChosenTrait;
//    double EventType = dice.rollContUnifDist(Trait[i].TotalTraitRate);
//    if(EventType < Trait[i].TotalBirthRate)
//        Events.isBirth = true;
//    else
//        Events.isBirth = false;
//}

//void TraitEventManager::executeEventTypeOnTrait()
//{
//    int ChosenTrait = Events.ChosenTrait;
//    if(Events.isBirth){
//        Trait[ChosenTrait].Members += 1.;
//    }
//    else if(Trait[ChosenTrait].Members > 0)
//        Trait[ChosenTrait].Members -= 1.;
//}

//void TraitEventManager::changeATrait()
//{
//    choseTraitToChange();
//    choseEventType();
//    executeEventTypeOnTrait();
//}

// --- Kernel ------------------------------------
// -----------------------------------------------


//void TraitEventManager::ImprovedEvolutionStep()
//{
//    sampleEventTime();
//    changeATrait();
//    adjustNewEventRates();
//}

// ------------------ Performence Improvement start -----------------

//void TraitEventManager::adjustNewEventRates()
//{
//    int currentTrait = Events.ChosenTrait;
//    int sgn;
//    if(Events.isBirth)
//        sgn = 1;
//    else
//        sgn = -1;

//    // Intrinsic Deathrates
//    Trait[currentTrait].TotalDeathRate += sgn*Trait[currentTrait].DeathRate;
//    // Extrinsic Deathrate
//    for(int i = 0; i < TraitClass::Size; ++i)
//        Trait[i].TotalDeathRate += sgn*TraitClass::CompDeathRate.at(i).at(currentTrait);
//    // Birthrates and TraitRates
//    Trait[currentTrait].TotalBirthRate += sgn*Trait[currentTrait].BirthRate;
//    Trait[currentTrait].TotalTraitRate = Trait[currentTrait].TotalBirthRate + Trait[currentTrait].TotalDeathRate;
//    if(currentTrait > 0){
//        Trait[currentTrait-1].TotalBirthRate += sgn*TraitClass::Mutation*Trait[currentTrait].BirthRate*0.5;
//        Trait[currentTrait-1].TotalTraitRate = Trait[currentTrait-1].TotalBirthRate + Trait[currentTrait-1].TotalDeathRate;
//    }
//    if(currentTrait < TraitClass::Size -1){
//        Trait[currentTrait+1].TotalBirthRate += sgn*TraitClass::Mutation*Trait[currentTrait].BirthRate*0.5;
//        Trait[currentTrait+1].TotalTraitRate = Trait[currentTrait+1].TotalBirthRate + Trait[currentTrait+1].TotalDeathRate;
//    }
//    TraitClass::TotalEventRate = 0;
//    for(int i = 0; i < TraitClass::Size; ++i)
//        TraitClass::TotalEventRate += Trait[i].TotalTraitRate;

//}

// ------------------ Performence Improvement end -----------------

//void TraitEventManager::clearData()
//{
//    TraitClass::setTraitSize(0);
//    Events.clearEvents();
//    Trait.clear();
//}

//double TraitEventManager::retStableDimorphKOf(int i) const
//{
//    int j;
//    if(i > 0)
//        j = 0;
//    else
//        j = 1;
//    double expected;

//    expected = (Trait[i].BirthRate - Trait[i].DeathRate)
//                * TraitClass::CompDeathRate[j][j];
//    expected -= (Trait[j].BirthRate - Trait[j].DeathRate)
//                * TraitClass::CompDeathRate[i][j];

//    double devisor = TraitClass::CompDeathRate[j][j]*TraitClass::CompDeathRate[i][i]
//            - TraitClass::CompDeathRate[j][i]*TraitClass::CompDeathRate[i][j];
//    if(devisor != 0)
//        expected /= devisor;
//    else
//        expected = -1.;

//    return expected/TraitClass::K;
//}

//QVector<double> TraitEventManager::retStableDimorphKVector() const
//{
//    QVector<double> expVals(TraitClass::Size);
//    for(int i = 0; i < TraitClass::Size; ++i)
//        expVals[i] = retStableDimorphKOf(i);
//    return expVals;
//}

//QVector<double> TraitEventManager::retStableMonoKVector() const
//{
//    QVector<double> expVals(TraitClass::Size);
//    for(int i = 0; i < TraitClass::Size; ++i)
//        expVals[i] = retStableKMonoOf(i);
//    return expVals;
//}

//double TraitEventManager::retStableKMonoOf(int i) const
//{
//    double expected = Trait[i].BirthRate /** (1-TraitClass::Mutation)*/ - Trait[i].DeathRate;
//    double divisor = TraitClass::CompDeathRate[i][i];
//    if(divisor > 0){
//        expected /= TraitClass::CompDeathRate[i][i];
//        return expected/TraitClass::K;
//    }
//    return -1.;
//}

//bool TraitEventManager::isNear(QVector<double> &Expected)
//{
//    bool close = true;
//    double diff;
//    for(int i = 0; i < TraitClass::Size; i++){
//        diff = fabs(getKMembers(i) - Expected[i]);
//        close = close && diff < 15./sqrt(TraitClass::K) && diff > -15./sqrt(TraitClass::K);
//    }
//    return close;
//}

//void TraitEventManager::calcFitness()
//{
//    for(int i = 0; i < TraitClass::Size; ++i){
//        for(int j = 0; j < TraitClass::Size; ++j){
//            double tmp = Trait[i].BirthRate  - Trait[i].DeathRate - TraitClass::CompDeathRate[i][j] * TraitClass::K * retStableKMonoOf(j);
//            TraitClass::Fitness[i][j] = round(tmp*pow(10,14))/pow(10,14);
//        }
//    }
//}





//void TraitEventManager::addTotalCompDeathRateOf(int i)
//{
//    double extDeath = 0;
//    for(int j = 0; j < TraitClass::Size; j++)
//        extDeath += TraitClass::CompDeathRate.at(i).at(j) * Trait[j].Members;
//    extDeath *= Trait[i].Members;
//    Trait[i].TotalDeathRate += extDeath;
//}

//void TraitEventManager::setTotalIntrisicDeathRateOf(int i)
//{
//    Trait[i].TotalDeathRate = (Trait[i].DeathRate)*(Trait[i].Members);
//}

//void TraitEventManager::calculateTotalBirthRates()
//{
//    int n = TraitClass::Size;
//    for(int k = 0; k < n; ++k)
//        Trait[k].TotalBirthRate = (Trait[k].Members)*(Trait[k].BirthRate);
//    for(int k = 1; k < n-1; ++k){
//        Trait[k].TotalBirthRate += 0.5*(Trait[k-1].Members)*(Trait[k-1].BirthRate)*(TraitClass::Mutation);
//        Trait[k].TotalBirthRate += 0.5*(Trait[k+1].Members)*(Trait[k+1].BirthRate)*(TraitClass::Mutation);
//    }
//    Trait[n-1].TotalBirthRate += 0.5*(Trait[n-2].Members)*(Trait[n-2].BirthRate)*(TraitClass::Mutation);
//    Trait[0].TotalBirthRate += 0.5*(Trait[1].Members)*(Trait[1].BirthRate)*(TraitClass::Mutation);
//}

//void TraitEventManager::calculateTotalEventRate()
//{
//    for(int i = 0; i < TraitClass::Size; i++){
//        Trait[i].TotalTraitRate = Trait[i].TotalBirthRate + Trait[i].TotalDeathRate;
//        TraitClass::TotalEventRate += Trait[i].TotalTraitRate;
//    }
//}

//void TraitEventManager::calculateEventRates()
//{
//    calculateTotalDeathRates();
//    calculateTotalBirthRates();
//    calculateTotalEventRate();
//}
