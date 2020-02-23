#include "dna.h"

//Default constructor
DNA::DNA()
{
  m_lineLengths = "";
  m_sum = 0;
  m_numLines = 0;
  m_A = 0;
  m_C = 0;
  m_T = 0;
  m_G = 0;
  m_AA = 0;
  m_AC = 0;
  m_AT = 0;
  m_AG = 0;
  m_CA = 0;
  m_CC = 0;
  m_CT = 0;
  m_CG = 0;
  m_TA = 0;
  m_TC = 0;
  m_TT = 0;
  m_TG = 0;
  m_GA = 0;
  m_GC = 0;
  m_GT = 0;
  m_GG = 0;
}

//Destructor
DNA::~DNA()
{

}


//Accessor for sum
int DNA::getSum()
{
  return m_sum;
}

//Accessor for number of lines
int DNA::getNumLines()
{
  return m_numLines;
}

//Returns the average length of strings from the file
double DNA::getMean()
{
  return (double)m_sum / m_numLines;
}

string DNA::getLineLengths()
{
  return m_lineLengths;
}

//Converts input to an uppercase string
void DNA::makeUpper(string& s)
{
  for (int i = 0; i < s.size(); ++i)
  {
    s[i] = toupper(s[i]);
  }
}

//Checks if a line in the file is DNA
bool DNA::isValid(string s)
{
  for (int i = 0; i < s.size(); ++i)
  {
    if (!((s[i] == 'A') || (s[i] == 'T') || (s[i] == 'C') || (s[i] == 'G')))
    {
      return false;
    }
  }
  return true;
}

//Adds length of a line in the file to the sum
void DNA::addToSum(string s)
{
  bool morePairs = true;
  bool evenIteration = true;
  int sum = 0;
  for (int i = 0; i < s.size(); ++i)
  {
    if (i == s.size() - 1)
    {
      morePairs = false;
    }
    if (i % 2 == 1)
    {
      evenIteration = false;
    }
    else
    {
      evenIteration = true;
    }

    if (s[i] == 'A')
    {
      ++m_A;
      if (evenIteration == true && morePairs == true)
      {
        if (s[i + 1] == 'A')
        {
          ++m_AA;
        }
        else if (s[i + 1] == 'C')
        {
          ++m_AC;
        }
        else if (s[i + 1] == 'T')
        {
          ++m_AT;
        }
        else if (s[i + 1] == 'G')
        {
          ++m_AG;
        }
      }
    }
    else if (s[i] == 'C')
    {
      ++m_C;
      if (evenIteration == true && morePairs == true)
      {
        if (s[i + 1] == 'A')
        {
          ++m_CA;
        }
        else if (s[i + 1] == 'C')
        {
          ++m_CC;
        }
        else if (s[i + 1] == 'T')
        {
          ++m_CT;
        }
        else if (s[i + 1] == 'G')
        {
          ++m_CG;
        }
      }
    }
    else if (s[i] == 'T')
    {
      ++m_T;
      if (evenIteration == true && morePairs == true)
      {
        if (s[i + 1] == 'A')
        {
          ++m_TA;
        }
        else if (s[i + 1] == 'C')
        {
          ++m_TC;
        }
        else if (s[i + 1] == 'T')
        {
          ++m_TT;
        }
        else if (s[i + 1] == 'G')
        {
          ++m_TG;
        }
      }
    }
    else if (s[i] == 'G')
    {
      ++m_G;
      if (evenIteration == true && morePairs == true)
      {
        if (s[i + 1] == 'A')
        {
          ++m_GA;
        }
        else if (s[i + 1] == 'C')
        {
          ++m_GC;
        }
        else if (s[i + 1] == 'T')
        {
          ++m_GT;
        }
        else if (s[i + 1] == 'G')
        {
          ++m_GG;
        }
      }
    }
    ++sum;
  }
  m_sum += sum;
  this->appendToLineLengths(sum);
}

//Appends a line's sum a string
void DNA::appendToLineLengths(int i)
{
  if (m_lineLengths.size() == 0)
  {
    m_lineLengths.append(to_string(i));
  }
  else
  {
    m_lineLengths.append("_");
    m_lineLengths.append(to_string(i));
  }
}

//Returns the variance
double DNA::getVariance()
{
  double numerator = 0;
  double var = 0;
  int tempNum = 0;
  string localLineLengths = this->getLineLengths();
  while (localLineLengths.size() > 0)
  {
    if (localLineLengths.find('_') != string::npos)
    {
      int delimiter = localLineLengths.find('_');
      int stringSize = localLineLengths.size() - 1;
      tempNum = stoi(localLineLengths.substr(0, delimiter));
      localLineLengths = localLineLengths.substr(delimiter + 1, stringSize);
    }
    else
    {
      tempNum = stoi(localLineLengths);
      localLineLengths = "";
    }
    numerator += ((tempNum - this->getMean()) * (tempNum - this->getMean()));
  }
  return numerator / (this->getNumLines());
}

//Returns the standard deviation of the string lengths
double DNA::getStdDeviation()
{
  return sqrt(this->getVariance());
}

//Computes the relative probability of each nucleotide
double DNA::getProbA()
{
  return (double)m_A / this->getSum();
}

//Computes the relative probability of each nucleotide
double DNA::getProbC()
{
  return (double)m_C / this->getSum();
}

//Computes the relative probability of each nucleotide
double DNA::getProbT()
{
  return (double)m_T / this->getSum();
}

//Computes the relative probability of each nucleotide
double DNA::getProbG()
{
  return (double)m_G / this->getSum();
}

double DNA::getProbAA()
{
  return (double)m_AA / (this->getSum() / 2);
}

double DNA::getProbAC()
{
  return (double)m_AC / (this->getSum() / 2);
}

double DNA::getProbAT()
{
  return (double)m_AT / (this->getSum() / 2);
}

double DNA::getProbAG()
{
  return (double)m_AG / (this->getSum() / 2);
}

double DNA::getProbCA()
{
  return (double)m_CA / (this->getSum() / 2);
}

double DNA::getProbCC()
{
  return (double)m_CC / (this->getSum() / 2);
}

double DNA::getProbCT()
{
  return (double)m_CT / (this->getSum() / 2);
}

double DNA::getProbCG()
{
  return (double)m_CG / (this->getSum() / 2);
}

double DNA::getProbTA()
{
  return (double)m_TA / (this->getSum() / 2);
}

double DNA::getProbTC()
{
  return (double)m_TC / (this->getSum() / 2);
}

double DNA::getProbTT()
{
  return (double)m_TT / (this->getSum() / 2);
}

double DNA::getProbTG()
{
  return (double)m_TG / (this->getSum() / 2);
}

double DNA::getProbGA()
{
  return (double)m_GA / (this->getSum() / 2);
}

double DNA::getProbGC()
{
  return (double)m_GC / (this->getSum() / 2);
}

double DNA::getProbGT()
{
  return (double)m_GT / (this->getSum() / 2);
}

double DNA::getProbGG()
{
  return (double)m_AG / (this->getSum() / 2);
}


//Computes the relative probability of each nucleiotide bigram


//Reads a file provided in the command line
// - Doesn't matter if strings are capitalized  or not in file
// - Keep asking until the user doesn't want to
void DNA::readFile(string filename)
{
  //Will repeat process as long as the user doesn't type "quit"
  while (filename != "exit")
  {
    ifstream inFS;
    ofstream outFS;

    string line = "";
    this->m_sum = 0;
    this->m_numLines = 0;
    this->m_lineLengths = "";

    this->m_A = 0;
    this->m_C = 0;
    this->m_T = 0;
    this->m_G = 0;
    this->m_AA = 0;
    this->m_AC = 0;
    this->m_AT = 0;
    this->m_AG = 0;
    this->m_CA = 0;
    this->m_CC = 0;
    this->m_CT = 0;
    this->m_CG = 0;
    this->m_TA = 0;
    this->m_TC = 0;
    this->m_TT = 0;
    this->m_TG = 0;
    this->m_GA = 0;
    this->m_GC = 0;
    this->m_GT = 0;
    this->m_GG = 0;

    //Opens the file
    cout << endl << endl << endl << "Reading input file." << endl << endl;
    inFS.open(filename);
    outFS.open("currentaber.out");

    //Tells the user if it couldn't open the file
    if (!inFS.is_open())
    {
      cout << "Couldn't open the file." << endl;
      cout << "Enter another file name (type \"exit\" to exit):" << endl;
      cin >> filename;
      continue;
    }

    //Reads the contents of the file
    while (!inFS.eof())
    {
      inFS >> line;
      if (!inFS.fail())
      {
        this->makeUpper(line);
        if (this->isValid(line))
        {
          this->m_numLines++;
          this->addToSum(line);
        }
      }
    }

    outFS << "Full Name: Curren Taber" << endl;
    outFS << "Student ID: 002325149" << endl;
    outFS << "Chapman Email: ctaber@chapman.edu" << endl;
    outFS << "Course Number and Section: CPSC 350-01" << endl;
    outFS << "Assignment 1: C++ Review" << endl << endl;

    outFS << setw(25) << setfill('-') << " " << endl << endl;
    outFS << setfill(' ');

    outFS << setw(20) << left << "SUM:";
    outFS << setw(5) << right << this->getSum() << endl;

    outFS << setw(20) << left << "MEAN:";
    outFS << setw(5) << right << this->getMean() << endl;

    outFS << setw(20) << left << "VARIANCE:";
    outFS << setw(5) << right << this->getVariance() << endl;

    outFS << setw(20) << left << "STD DEVIATION:";
    outFS << setw(5) << right << this->getStdDeviation() << endl << endl;

    outFS << setw(20) << left << "PROB A:";
    outFS << setw(5) << right << this->getProbA() << endl;
    outFS << setw(20) << left << "PROB C:";
    outFS << setw(5) << right << this->getProbC() << endl;
    outFS << setw(20) << left << "PROB T:";
    outFS << setw(5) << right << this->getProbT() << endl;
    outFS << setw(20) << left << "PROB G:";
    outFS << setw(5) << right << this->getProbG() << endl << endl;

    outFS << setw(20) << left << "PROB AA:";
    outFS << setw(5) << right << this->getProbAA() << endl;
    outFS << setw(20) << left << "PROB AC:";
    outFS << setw(5) << right << this->getProbAC() << endl;
    outFS << setw(20) << left << "PROB AT:";
    outFS << setw(5) << right << this->getProbAT() << endl;
    outFS << setw(20) << left << "PROB AG:";
    outFS << setw(5) << right << this->getProbAG() << endl;

    outFS << setw(20) << left << "PROB CA:";
    outFS << setw(5) << right << this->getProbCA() << endl;
    outFS << setw(20) << left << "PROB CC:";
    outFS << setw(5) << right << this->getProbCC() << endl;
    outFS << setw(20) << left << "PROB CT:";
    outFS << setw(5) << right << this->getProbCT() << endl;
    outFS << setw(20) << left << "PROB CG:";
    outFS << setw(5) << right << this->getProbCG() << endl;

    outFS << setw(20) << left << "PROB TA:";
    outFS << setw(5) << right << this->getProbTA() << endl;
    outFS << setw(20) << left << "PROB TC:";
    outFS << setw(5) << right << this->getProbTC() << endl;
    outFS << setw(20) << left << "PROB TT:";
    outFS << setw(5) << right << this->getProbTT() << endl;
    outFS << setw(20) << left << "PROB TG:";
    outFS << setw(5) << right << this->getProbTG() << endl;

    outFS << setw(20) << left << "PROB GA:";
    outFS << setw(5) << right << this->getProbGA() << endl;
    outFS << setw(20) << left << "PROB GC:";
    outFS << setw(5) << right << this->getProbGC() << endl;
    outFS << setw(20) << left << "PROB GT:";
    outFS << setw(5) << right << this->getProbGT() << endl;
    outFS << setw(20) << left << "PROB GG:";
    outFS << setw(5) << right << this->getProbGG() << endl;

    cout << "Results from \"" << filename;
    cout << "\" printed to \"currentaber.out\"" << endl << endl;

    //Closes the file
    inFS.close();
    outFS.close();

    cout << "Enter another file name (type \"exit\" to exit):" << endl;
    cin >> filename;
  }
}

//Generates 1000 strings whose lengths follow the Gaussian distribution
// - Strings have the same mean and variance as calculated above
// - Relative frequency of nucleotides will follow as calculated above
// - Append results to the output results
void DNA::generateDNA()
{

}
