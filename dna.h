#ifndef DNA_H
#define DNA_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cstdlib>

using namespace std;

class DNA
{
  public:
    //Default constructor
    DNA();

    //Destructor
    ~DNA();

    //Reads a file provided in the command line
    //Asks the user if they want to enter another file or quit
    void readFile(string filename);

  private:
    string m_lineLengths;
    const unsigned int m_BIGRAM_POSSIBILITIES = 16;
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

    //Returns the sum of letters in the file
    int getSum();

    //Returns the number of DNA lines in the file
    int getNumLines();

    //Returns the average length of strings from the file
    double getMean();

    //Converts an input string to an uppercase string
    void makeUpper(string& s);

    //Checks if a line only contains valid nucleotide letters
    bool isValid(string s);

    //Adds length of a line in the file to the sum
    //Counts the number of nucleotides in a line
    //Counts the number of nucleotide bigrams in a line (non-overlapping)
    //In an odd string, the last nucleotide doens't have a pair
    void addToSum(string s);

    //Appends a line's sum to the string  m_lineLengths
    void appendToLineLengths(int d);

    //Returns the variance of string lengths from the file
    double getVariance();

    //Returns a string that stores each line's length
    //Each integer is delimited by a '_'
    string getLineLengths();

    //Returns the standard deviation of the string lengths
    double getStdDeviation();

    //The following functions return a relative probability of each nucleotide
    double getProbA();
    double getProbC();
    double getProbT();
    double getProbG();

    //The following functions return a relative probability of each bigram
    double getProbAA();
    double getProbAC();
    double getProbAT();
    double getProbAG();
    double getProbCA();
    double getProbCC();
    double getProbCT();
    double getProbCG();
    double getProbTA();
    double getProbTC();
    double getProbTT();
    double getProbTG();
    double getProbGA();
    double getProbGC();
    double getProbGT();
    double getProbGG();

    //Generates 1000 strings whose lengths follow the Gaussian Distribution
    void generateDNA();
};

#endif
