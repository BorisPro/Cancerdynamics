#ifndef FILESTREAMING_H
#define FILESTREAMING_H
#include <vector>
#include <QString>
#include <QFile>
#include <QTextStream>
#include "traiteventmanager.h"

class PopulationManager;    // improve me with Include Guards please!
class FileStreaming
{
public:
    FileStreaming();
    FileStreaming(QString FName);
    ~FileStreaming();

    bool initializeWithFile(QString FName, std::vector<TraitClass> &Traits);

private:
    void readSize_Mutation_K(QTextStream&in, std::vector<TraitClass> & Traits);
    void readMembersBirthsDeaths(QTextStream&in, std::vector<TraitClass> & Traits);
    void readCompetitionKernel(QTextStream&in);
    void clearTraitsAndRewriteWith(QFile& mFile, std::vector<TraitClass>& _Traits);
    std::vector<TraitClass> Traits;
    QFile File;

public:
    /// This is an alterntive Way of storing Data. It streams Data into a Storage File and gives
    /// the opportunity to read the vector from this source.
    /// FIXME: It is not clean yet!!! (but works properly)
    bool writeLimitedData(PopulationManager Manager);
    void copyEvolutionIntoData(PopulationManager &Manager, QTextStream &in, double &ElapsedTime);

//    int storeEvolution(TraitEventManager Manager, int max_it);
    /// This AppendToFile function is not in use,
    /// but it will be usefull when Update functionality comes up
    static void appendToFileEndl(QString FileName, QString Value);
};

#endif // FILESTREAMING_H
