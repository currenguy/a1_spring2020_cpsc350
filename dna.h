#ifndef DNA_H
#define DNA_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;

class DNA
{
  public:
    //Default constructor
    DNA();

    //Destructor
    ~DNA();

    void readFile(string filename);

  private:
    string m_lineLengths;
    const unsigned int NUCLEOTIDE_POSSIBLITIES = 4;
    const unsigned int BIGRAM_POSSIBILITIES = 16;
    unsigned int m_sum;
    unsigned int m_numLines;
    unsigned int m_A;
    unsigned int m_C;
    unsigned int m_T;
    unsigned int m_G;
    unsigned int m_AA;
    unsigned int m_AC;
    unsigned int m_AT;
    unsigned int m_AG;
    unsigned int m_CA;
    unsigned int m_CC;
    unsigned int m_CT;
    unsigned int m_CG;
    unsigned int m_TA;
    unsigned int m_TC;
    unsigned int m_TT;
    unsigned int m_TG;
    unsigned int m_GA;
    unsigned int m_GC;
    unsigned int m_GT;
    unsigned int m_GG;

    //Increments number of lines
    void incNumLines();

    //Accessor for sum
    int getSum();

    //Accessor for number of lines
    int getNumLines();

    //Returns the average length of strings from the file
    double getMean();

    //Converts input to an uppercase string
    void makeUpper(string& s);

    //Checks if a line in the file is DNA
    bool isValid(string s);

    //Adds length of a line in the file to the sum
    void addToSum(string s);

    void appendToLineLengths(int d);

    //Computes the variance of the string lengths
    double getVariance();

    string getLineLengths();

    //Computes the standard deviation of the string lengths
    double getStdDeviation();

    //Computes the relative probability of each nucleotide
    double getProbA();

    double getProbC();

    double getProbT();

    double getProbG();

    //Computes the relative probability of each nucleiotide bigram
    double getRelProbBigram();

    //Outputs the labeled results to currentaber.out
    string printResults();

    //Generates 1000 strings whose lengths follow the Gaussian distribution
    // - Strings have the same mean and variance as calculated above
    // - Relative frequency of nucleotides will follow as calculated above
    // - Append results to the output results
    void generateDNA();


};

#endif
