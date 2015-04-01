#include <QString>
#include <QtTest>
#include <vector>
#include <string>
#include <windows.h>
#include "../../Simulator/Cancerdynamics/graphclass.h"
#include "Testing_Utilities/Utilities.h"

using std::vector;
using std::string;
using std::size_t;
using namespace Utilities;

class TDD_CancerdynamicsTest : public QObject
{
    Q_OBJECT

public:
    TDD_CancerdynamicsTest();

    vector<double> getMembers(PopulationManager manager);
    vector<double> getProductionRates(PopulationManager manager);
    vector<double> getMutationProbabilities(PopulationManager manager);
    vector<double> getProductionAmount(PopulationManager manager);
private Q_SLOTS:
    void fileStream_Test();
    void readTestInstance_Test();
    void Utility_Tester();
    // -----------------------------------------
    // --- Kernel ------------------------------
    void addNaturalDeathRates_Test();
    void addCompetitionDeathRates_Test();
    void calculateTraitDeathRates_Test();
    void addNaturalBirthRates_Test();
    void addMutationalBirthRates_Test();
    void subtractBirthReduction_Test();
    void calculateTraitBirthRates_Test();
    void calculateProductionRates_Test();
    void addNaturalSwitchRates_Test();
    void addCompetitionSwitchRates_Test();
    void calculateSwitchRates_Test();
    void calculateKillRates_Test();
    void calculateTotalAndTraitRates_Test();
    void applyKToParameters_Test();

    void choseTraitToChange_Test();
    void choseEventType_Test();
    void isMutation_Test();
    void choseMutantGenotype_Test();
    void executeBirth_Test();
    void getSwitchedTrait_Test();
    void executeSwitch_Test();
    void executeProduction_Test();
    void getKilledTrait_Test();

    void EvolutionStep_Test();

    void testWritingToFile();
};

TDD_CancerdynamicsTest::TDD_CancerdynamicsTest()
{
}

void TDD_CancerdynamicsTest::fileStream_Test()
{
    CFileStreamer test("../Parameters/Testinstanz");
    IFileStreamer& streamer(test);
    size_t firstEntry = std::stoi(streamer.getNextWord());
    QCOMPARE(firstEntry, (size_t) 2);
    qDebug()<<"Entries:" << QString::number(firstEntry);
    Utilities::printMassage(streamer.getNextWords(firstEntry));
}

vector<double> TDD_CancerdynamicsTest::getMembers(PopulationManager manager)
{
    vector<double> members;
    for(TraitClass trait : manager.traits)
        members.push_back(trait.Members);
    return members;
}

vector<double> TDD_CancerdynamicsTest::getProductionRates(PopulationManager manager)
{
    vector<double> rate;
    for(TraitClass trait : manager.traits)
        rate.push_back(trait.ProductionRate);
    return rate;
}

vector<double> TDD_CancerdynamicsTest::getMutationProbabilities(PopulationManager manager)
{
    vector<double> rate;
    for(TraitClass trait : manager.traits)
        rate.push_back(trait.Mutation);
    return rate;
}

vector<double> TDD_CancerdynamicsTest::getProductionAmount(PopulationManager manager)
{
    vector<double> rate;
    for(TraitClass trait : manager.traits)
        rate.push_back(trait.ProductionAmount);
    return rate;
}

void TDD_CancerdynamicsTest::readTestInstance_Test()
{
    PopulationManager manager("Testinstanz");
    printMassage("T-cells: " + std::to_string(manager.tCells));
    printMassage("Genotypes: " + std::to_string(manager.genotypes));
    printMassage("Phenotypes: " + std::to_string(manager.phenotypes));
    printMassage(getMembers(manager));
    printMassage("Production Rates and Amounts:");
    printMassage(getProductionRates(manager));
    printMassage(getProductionAmount(manager));
//    printMassage("Competitionmatrix:");
//    printMassage(manager.competition);
//    printMassage("Birthreduction matrix:");
//    printMassage(manager.birthReduction);
//    for(auto matrix : manager.naturalSwitch)
//        printMassage(matrix);
//    for(auto matrix : manager.competitionSwitch)
//        printMassage(matrix);
//    printMassage(manager.tCellKillRate);
//    printMassage(getMutationProbabilities(manager));
//    printMassage(manager.mutationKernel);
//    printMassage(manager.K);
    QCOMPARE(manager.K, 1.);
    QCOMPARE(manager.mutationKernel[2][2], 0.);
}

void TDD_CancerdynamicsTest::Utility_Tester()
{
//    PopulationManager manager("Testinstanz");
//    qDebug()<<manager.getSpecificSwitchRate(0,0,1);

//    PopulationManager manager("Testinstanz");
//    vector<double> tmp;
//    for(size_t i = 0; i < manager.genotypes; ++i)
//        for(size_t j = 0; j < manager.phenotypes; ++j)
//            tmp.push_back(manager.getMelanomIndex(i,j));
//    printMassage(tmp);

//    PopulationManager manager("Testinstanz");
//    vector<double> output;
//    for(size_t i = 1 + manager.tCells; i < manager.populations; ++i)
//        output.push_back(manager.getPhenotypeOfTraitIndex(i));
//    printMassage(output);
}

// -----------------------------------------
// --- Kernel ------------------------------

void TDD_CancerdynamicsTest::addNaturalDeathRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.addNaturalDeathRates();
    vector<double> output(1,0);
    for(TraitClass trait : manager.traits){
        output[0] += trait.TraitDeathRate; // 10 * 1 * 3 + 20 * 2 * 6
        output.push_back(trait.TraitDeathRate);
    }
    printMassage(output);
    QCOMPARE(output[0], 270.);
}

void TDD_CancerdynamicsTest::addCompetitionDeathRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.addCompetitionDeathRates();
    vector<double> output(1,0);
    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitDeathRate);
        output[0] += output.back();
    }

    printMassage(output);   // 1 * 10 * 10 * 3 + 1 * 20 * 20 * 5 + (1 * 20*20 + 1 * 20*20)
    QCOMPARE(output[0], 3100.);
}

void TDD_CancerdynamicsTest::calculateTraitDeathRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.addNaturalDeathRates();
    manager.addCompetitionDeathRates();
    vector<double> output(1,0);
    for(TraitClass trait : manager.traits){
        output[0] += trait.TraitDeathRate;
        output.push_back(trait.TraitDeathRate);
    }
    printMassage(output);
    QCOMPARE(output[0], 3370.);
}

void TDD_CancerdynamicsTest::addNaturalBirthRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.addNaturalBirthRates();
    vector<double> output(1,0);
    vector<double> expected = {830, 50, 50, 50, 100, 100, 120, 120, 120, 120};
    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitBirthRate);
        output[0] += output.back();
    }   // 50 + 2*50 + 100*2 + 120*4
    printMassage(output);
    QCOMPARE(output, expected);
}

void TDD_CancerdynamicsTest::addMutationalBirthRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.addMutationalBirthRates();
    vector<double> output(1,0);
    vector<double> expected = {136, 0, 0, 0, 24, 24, 16, 16, 28, 28};
    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitBirthRate);
        output[0] += output.back();
    }
    printMassage(output);
    QCOMPARE(expected, output);
}

void TDD_CancerdynamicsTest::subtractBirthReduction_Test()
{
    PopulationManager manager("Testinstanz");
    manager.subtractBirthReduction();   // tut nichts weil negative Raten abgefangen werden.
    manager.addNaturalBirthRates();     // erzeugt positive Raten
    manager.subtractBirthReduction();
    vector<double> output(1,0);
    vector<double> expected = {830-270, 50-10, 50-10, 50-10, 100-40, 100-40, 120-40, 120-40, 120-40, 120-40};
    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitBirthRate);
        output[0] += output.back();
    }
    printMassage(output);   // 0.1 * 10 * 10 * 3 + 0.1 * 20 * 20 * 6
    QCOMPARE(output, expected);
}

void TDD_CancerdynamicsTest::calculateTraitBirthRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.calculateTraitBirthRates();
    vector<double> output(1,0);
    vector<double> expected = {560, 40, 40, 40, 60, 60, 80, 80, 80, 80};
    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitBirthRate);
        output[0] += output.back();
    }
    printMassage(output);   // 830 - 270
    QCOMPARE(output, expected);
}

void TDD_CancerdynamicsTest::calculateProductionRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.calculateTraitProductionRates();
    vector<double> output(1,0);
    vector<double> expected = {3600, 0, 1200, 2400, 0, 0, 0, 0, 0, 0 };
    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitProductionRate);
        output[0] += output.back();
    }
    printMassage(output);
    QCOMPARE(output, expected);
}

void TDD_CancerdynamicsTest::addNaturalSwitchRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.addNaturalSwitchRates();
    vector<double> output(1,0);
    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitSwitchRate);
        output[0] += output.back();
    }
    printMassage(output);
    QCOMPARE(output[0], 120.);
}

void TDD_CancerdynamicsTest::addCompetitionSwitchRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.addCompetitionSwitchRates();
    vector<double> output(1,0);
    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitSwitchRate);
        output[0] += output.back();
    }
    printMassage(output);
    QCOMPARE(output[0], 700.);
}

void TDD_CancerdynamicsTest::calculateSwitchRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.calculateTraitSwitchRates();
    vector<double> output(1,0);
    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitSwitchRate);
        output[0] += output.back();
    }
    printMassage(output);
    QCOMPARE(output[0], 820.);
}

void TDD_CancerdynamicsTest::calculateKillRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.calculateTraitKillRates();
    vector<double> output(1,0);
    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitKillRate);
        output[0] += output.back();
    }
    printMassage(output);
    QCOMPARE(output[0], 3800.);
}

void TDD_CancerdynamicsTest::calculateTotalAndTraitRates_Test()
{
    PopulationManager manager("Testinstanz");
    manager.calculateTotalEventRates();
    vector<double> expected = {12150, 150, 3350, 4350, 744, 624, 632, 712, 544, 1044};
    vector<double> output(1,0);

    for(TraitClass trait : manager.traits){
        output.push_back(trait.TraitRate);
        output[0] += output.back();
    }
    printMassage(output);
    QCOMPARE(output, expected);
    QCOMPARE(manager.totalEventRate, output[0]);
}

void TDD_CancerdynamicsTest::applyKToParameters_Test()
{
    PopulationManager manager("Testinstanz");
    manager.K = 10;
    manager.applyKToParameters();

    printMassage("Members for K = 10: ");
    printMassage(getMembers(manager));
    printMassage("Production Rates:");
    printMassage(getProductionRates(manager));
    printMassage("Competitionmatrix:");
    printMassage(manager.competition);
    printMassage("Birthreduction matrix:");
    printMassage(manager.birthReduction);
    printMassage("competition switch matrices:");
    for(auto matrix : manager.competitionSwitch)
        printMassage(matrix);
    printMassage("tCellKillRate:");
    printMassage(manager.tCellKillRate);
    printMassage("K:");
    printMassage(manager.K);
}

void TDD_CancerdynamicsTest::choseTraitToChange_Test()
{
    PopulationManager manager("Testinstanz");
    manager.totalEventRate = 100;
    for(TraitClass trait: manager.traits)
        trait.TraitRate = 0;
    manager.traits[0].TraitRate = 25;
    manager.traits[1].TraitRate = 25;
    manager.traits[2].TraitRate = 50;
    vector<double> history(3,0);
    for(int i = 0; i < 1000; ++i){
        manager.choseTraitToChange();
        history[manager.event.chosenTrait]+= 1./1000;
    }
    history[0] = fabs(history[0] - 0.25)/0.25;
    history[1] = fabs(history[0] - 0.25)/0.25;
    history[2] = fabs(history[0] - 0.5)/0.5;
    printMassage(history);
    for(double hist : history){
        QVERIFY(hist < 2.);
    }
}

void TDD_CancerdynamicsTest::choseEventType_Test()
{
    PopulationManager manager("Testinstanz");
    manager.totalEventRate = 100;
    for(TraitClass trait: manager.traits)
        trait.TraitRate = 0;
    manager.event.chosenTrait = 0;
    manager.traits[0].TraitRate = 125;
    manager.traits[0].TraitBirthRate = 35;
    manager.traits[0].TraitDeathRate = 30;
    manager.traits[0].TraitSwitchRate = 25;
    manager.traits[0].TraitProductionRate = 20;
    manager.traits[0].TraitKillRate = 15;
    vector<double> history(5,0);
    for(int i = 0; i < 100000; ++i){
        manager.choseEventType();
        if(manager.event.type == BIRTH)
             history[0]+= 1./100000;
        if(manager.event.type == DEATH)
             history[1]+= 1./100000;
        if(manager.event.type == SWITCH)
             history[2]+= 1./100000;
        if(manager.event.type == PRODUCTION)
             history[3]+= 1./100000;
        if(manager.event.type == KILL)
             history[4]+= 1./100000;
    }
    printMassage(history);
    history[0] = fabs(history[0] - 0.35)/0.35;
    history[1] = fabs(history[1] - 0.30)/0.30;
    history[0] = fabs(history[2] - 0.25)/0.25;
    history[1] = fabs(history[3] - 0.20)/0.20;
    history[0] = fabs(history[4] - 0.15)/0.15;
}

void TDD_CancerdynamicsTest::isMutation_Test()
{
    PopulationManager manager("Testinstanz");
    manager.calculateTotalEventRates();
    manager.sampleEventTime();
    manager.event.chosenTrait = 5;
    double hits = 0;
    for(int i = 0; i < 1000; ++i)
            hits+= 1./1000 * manager.isMutation();
    hits = fabs((hits-0.2)/0.2);
    printMassage(hits);
    QVERIFY( hits < 0.1);
}

void TDD_CancerdynamicsTest::choseMutantGenotype_Test()
{
    PopulationManager manager("Testinstanz");
    vector<double> hits(2,0);
    for(int i = 0; i < 1000; ++i){
        manager.event.chosenTrait = 3;
        manager.choseMutantGenotype();
        if(manager.event.chosenTrait == 3 + manager.phenotypes)
            hits[0] += 1./1000;
        else if(manager.event.chosenTrait == 3 + 2*manager.phenotypes)
            hits[1] += 1./1000;
        else
            QVERIFY(false);
    }
    printMassage(hits);
    QVERIFY(fabs(hits[0]-0.2)/0.2 < 0.1);
    QCOMPARE(hits[0] + hits[1], 1.);
}

void TDD_CancerdynamicsTest::executeBirth_Test()
{
    PopulationManager manager("Testinstanz");
    for(TraitClass& trait : manager.traits)
        trait.Members = 0;
    int itterations = 1000;
    vector<double> expected = {0, 0, 0, 0.8, 0, 0.04, 0, 0.160, 0};
    for(int i = 0; i < itterations; ++i){
        manager.event.chosenTrait = 3;
        manager.executeBirth();
        if((int)manager.event.chosenTrait%2 == 0)
            QVERIFY(false);
    }
    vector<double> hits;
    for(size_t i = 0; i < manager.populations; ++i){
        hits.push_back(fabs(manager.traits[i].Members/itterations - expected[i]));
        if(hits.back() != 0.)
            hits.back() /= expected[i];
    }
    printMassage(hits);
    for(double comp : hits)
        QVERIFY(comp < 0.5);
    QCOMPARE(manager.traits[3].Members + manager.traits[5].Members + manager.traits[7].Members, (double)itterations);
}

void TDD_CancerdynamicsTest::getSwitchedTrait_Test()
{
    PopulationManager manager("Testinstanz");
    manager.calculateTotalEventRates();
    manager.event.chosenTrait = manager.getMelanomIndex(0,0);
    qDebug()<< manager.event.chosenTrait << manager.getSwitchedTrait();
    QSKIP("Can not be tested properly with less than 3 phenotypes!");
}

void TDD_CancerdynamicsTest::executeSwitch_Test()
{
    QSKIP("Can not be tested properly with less than 3 phenotypes!");
}

void TDD_CancerdynamicsTest::executeProduction_Test()
{
    PopulationManager manager("Testinstanz");
    size_t chosen = 2;
    manager.event.chosenTrait = chosen;
    manager.traits[chosen].Members=0;
    manager.traits[0].Members=0;
    manager.executeProduction();
    qDebug()<< manager.traits[manager.event.chosenTrait].ProductionAmount;
    QCOMPARE(manager.traits[chosen].Members,1.);
    QCOMPARE(manager.traits[0].Members,3.);
//    qDebug()<<manager.traits[chosen].Members;
}

void TDD_CancerdynamicsTest::getKilledTrait_Test()
{
    PopulationManager manager("Testinstanz");
    manager.calculateTotalEventRates();
    size_t chosen = 1;
    vector<double> hist(manager.populations,0.);
    for(int i = 0; i < 10000; ++i){
        manager.event.chosenTrait = chosen;
        hist[manager.getKilledTrait()]+=1./1000;
    }
    printMassage(hist);
}

void TDD_CancerdynamicsTest::EvolutionStep_Test()
{
    string fName = "Testinstanz";
    GraphClass graph(QString::fromStdString(fName));
    try{qDebug()<<"iterations:" << graph.generateEvolution(100000);}
    catch(string error){
        qDebug()<<QString::fromStdString(error);
    }
}

void TDD_CancerdynamicsTest::testWritingToFile()
{
    QString filename="Data.txt";
    QFile file( filename );
    if ( file.open(QIODevice::ReadWrite) )
    {
        QTextStream stream( &file );
        stream << "something else" << endl;
    }
    else{
        QVERIFY(false);
    }

//    CFileStreamer object;
//    if(false != object.accessFileToWrite("Boris.txt"))
//        QVERIFY(false);
//    object.wirteToAccessedFile("test");
//    object.closeFileToWrite();
}




QTEST_APPLESS_MAIN(TDD_CancerdynamicsTest)

#include "tst_tdd_cancerdynamicstest.moc"
