#ifndef UTILITIES_H
#define UTILITIES_H

#include <QString>
#include <QtTest>
#include <vector>
#include <string>
#include "../../../Simulator/Cancerdynamics/FileStreamer/ifilestreamer.h"

using std::vector;
using std::string;
using std::size_t;

namespace Utilities{
void printMassage(const vector<vector<double>> words);
void printMassage(const vector<string> words);
void printMassage(const vector<double> words);
void printMassage(const string word);
void printMassage(const double word);
void printInstanceParameters(IFileStreamer& stream);
}

namespace UtilityHelper{
void printAndReadNi(IFileStreamer& stream, vector<size_t>& N);
}

//class Utilities {
//public:
//    static void printMassage(const vector<string> words);
//    static void printMassage(const vector<double> words);
//    static void printMassage(const string word);
//    static void printMassage(const double word);
//    static void printInstanceParameters(IFileStreamer& stream);

//private:
//    void printAndReadNi(IFileStreamer& stream, vector<size_t>& N);
//};


#endif // UTILITIES_H


