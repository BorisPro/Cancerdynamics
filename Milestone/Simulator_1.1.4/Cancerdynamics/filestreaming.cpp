#include "filestreaming.h"
#include <iostream>
#include <QTextStream>
#include <QDebug>
#include <cmath>


FileStreaming::FileStreaming()
{
}

FileStreaming::FileStreaming(QString FName):File(FName)
{
    if(!File.open(QFile::WriteOnly | QFile::Text | QIODevice::Truncate))
        qDebug() << "Datei konnte nicht geoffnet werden :(";
}

FileStreaming::~FileStreaming()
{
    if(File.isOpen())
        File.close();
}

void FileStreaming::clearTraitsAndRewriteWith(QFile &mFile, std::vector<TraitClass>& _Traits)
{
    Traits.clear();
    QTextStream in(&mFile);
    readSize_Mutation_K(in, _Traits);
    readCompetitionKernel(in);
    readMembersBirthsDeaths(in, _Traits);
    mFile.flush();
}

bool FileStreaming::initializeWithFile(QString FName, std::vector<TraitClass> & _Traits)
{
    QFile mFile(FName);

    if(!mFile.open(QFile::ReadOnly)){
        qDebug() << "Datei konnte nicht geoffnet werden :(";
        return false;
    }
    clearTraitsAndRewriteWith(mFile, _Traits);
    mFile.close();
    return true;
}

void FileStreaming::readSize_Mutation_K(QTextStream & in, std::vector<TraitClass> & Traits)
{
    int size;
    in>>size;
    Traits.assign(size,TraitClass());
    TraitClass::setTraitSize(size);
    in>>TraitClass::K;
    if(TraitClass::K == 0) TraitClass::K = 1;
    in>>TraitClass::Mutation;
    TraitClass::Mutation = TraitClass::Mutation/(TraitClass::K * sqrt(TraitClass::K));
}

void FileStreaming::readMembersBirthsDeaths(QTextStream &in, std::vector<TraitClass> & Traits)
{

    for(int i = 0; i < TraitClass::Size; i++){
        in>>Traits.at(i).Members;
        Traits.at(i).Members = Traits.at(i).Members*TraitClass::K;
    }
    for(int i = 0; i < TraitClass::Size; i++)
        in>>Traits.at(i).BirthRate;
    for(int i = 0; i < TraitClass::Size; i++)
        in>>Traits.at(i).DeathRate;
}

void FileStreaming::readCompetitionKernel(QTextStream &in)
{
    double p;
    for(int i = 0; i < TraitClass::Size; i++){
        for(int j = 0; j < TraitClass::Size; j++){
            in>>p;
            TraitClass::CompDeathRate.at(i).at(j) = p;
            p = TraitClass::K;
            TraitClass::CompDeathRate.at(i).at(j) /= TraitClass::K;
            p = TraitClass::CompDeathRate.at(i).at(j);
        }
    }
}

void FileStreaming::appendToFileEndl(QString FileName, QString Value)
{
    QFile mFile(FileName);
    mFile.open(QFile::Append | QFile::Text);
    QTextStream out(&mFile);
    out<<Value<<endl;
    mFile.close();
}

bool FileStreaming::writeLimitedData(TraitEventManager Manager)
{
    if(File.isOpen())
    {
        QTextStream out(&File);
        out<<Manager.Events.EventTimes<<" ";
        for(int i = 0; i < TraitClass::Size; i++)
            out<<Manager.getKMembers(i)<<" ";
        out<<endl;
        return true;
    }
    return false;
}

void FileStreaming::copyEvolutionIntoData(TraitEventManager& Manager, QTextStream& in, double& ElapsedTime)
{
    ElapsedTime +=Manager.Events.EventTimes;
    in<<ElapsedTime<<" ";
    for(int i = 0; i < TraitClass::Size; i++)
        in<<Manager.getKMembers(i)<<" ";
    in<<endl;
}

//int FileStreaming::storeEvolution(TraitEventManager Manager, int max_it)
//{
//    QFile Storage("Storage.txt");
//    Storage.open(QFile::WriteOnly | QFile::Text | QIODevice::Truncate);
//    QTextStream in(&Storage);
//    QVector<double> expected = Manager.retStableDimorphKVector();
//    double ElapsedTime = 0;
//    for (int iterator=0; iterator< max_it; ++iterator)
//    {
//        Manager.EvolutionStep();
//        copyEvolutionIntoData(Manager, in, ElapsedTime);
//        if(Manager.isNear(expected) && max_it > iterator + TraitClass::K * 10)
//            max_it = iterator + TraitClass::K * 10;
//    }
//    Storage.close();
//    return max_it;
//}
