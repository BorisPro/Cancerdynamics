#include "Utilities.h"

// -------------------------------------
// --- public --------------------------

using namespace Utilities;
using namespace UtilityHelper;

using std::vector;
using std::string;

void Utilities::printMassage(const vector<vector<double>> words)
{
    for(auto line : words){
        QString output;
        for(auto word : line)
            output.append(QString::number(word) + " ");
        qDebug()<<output;
    }
}

void Utilities::printMassage(const vector<string> words){
    QString line;
    for(string word : words)
        line.append(QString::fromStdString(word) + " ");
    qDebug()<<line;
}

void Utilities::printMassage(const vector<double> words){
    QString line;
    for(auto word : words)
        line.append(QString::number(word) + " ");
    qDebug()<<line;
}

void Utilities::printMassage(const string word){
    qDebug()<<QString::fromStdString(word);
}

void Utilities::printMassage(const double word){
    qDebug()<<QString::number(word);
}

void Utilities::printInstanceParameters(IFileStreamer& stream){
    vector<size_t> N(4,0);
    printAndReadNi(stream, N);
    Utilities::printMassage(stream.getNextWord());      // TNF-alpha
    Utilities::printMassage(stream.getNextWords(N[1])); // T-cells
    for(uint i = 0; i < N[3]; ++i)                      // Genotypes for i. Phenotype
        printMassage(stream.getNextWords(N[2]));
    printMassage(stream.getNextWords(N[0]));            // Deathrates for all populations
    printMassage(stream.getNextWords(N[0]));            // Birthrates for all populations
    printMassage(stream.getNextWords(N[1]));            // Productionrates for T-cells
    for(uint i = 0; i < N[0]; ++i)                      // Competitionmatrix
        printMassage(stream.getNextWords(N[0]));
    for(uint i = 0; i < N[0]; ++i)                      // Birthreducing matrix
        printMassage(stream.getNextWords(N[0]));
    for(uint i = 0; i < N[2]; ++i)                      // Switchmatrix (natural) for i. Genotype
        for(uint j = 0; j < N[3]; ++j)                  // Switchmatrix between Phenotypes
            printMassage(stream.getNextWords(N[3]));
    for(uint i = 0; i < N[2]; ++i)                      // Switchmatrix (TNF-alpha) for i. Genotype
        for(uint j = 0; j < N[3]; ++j)                  // Switchmatrix between Phenotypes
            printMassage(stream.getNextWords(N[3]));
    for(uint i = 0; i < N[1]; ++i)                      // T-cell Death matrix
        printMassage(stream.getNextWords(N[2]*N[3]));
    for(uint i = 0; i < N[2]; ++i)                      // Mutation matrix
        printMassage(stream.getNextWords(N[2]));
    Utilities::printMassage("K: " + stream.getNextWord());      // TNF-alpha
}

// -------------------------------------
// --- private -------------------------

void UtilityHelper::printAndReadNi(IFileStreamer& stream, vector<size_t>& N)
{
    for(size_t i = 1; i < N.size(); ++i){
        N[i] = stoi(stream.getNextWord());
        Utilities::printMassage(N[i]);
    }
    N[0] = N[1] + N[2] * N[3] + 1;
    printMassage(N[0]);
}


