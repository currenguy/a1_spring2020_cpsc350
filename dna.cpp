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
  int sum = 0;
  for (int i = 0; i < s.size(); ++i)
  {
    if (s[i] == 'A')
    {
      ++m_A;
    }
    else if (s[i] == 'C')
    {
      ++m_C;
    }
    else if (s[i] == 'T')
    {
      ++m_T;
    }
    else if (s[i] == 'G')
    {
      ++m_G;
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
  double numerator = 0;
  double dev = 0;
  int tempNum = 0;
  string localLineLengths = this->getLineLengths();
  while (localLineLengths.size() > 0)
  {
    if (localLineLengths.find('_') != string::npos)
    {
      tempNum = stoi(localLineLengths.substr(0, localLineLengths.find('_')));
      localLineLengths = localLineLengths.substr(localLineLengths.find('_') + 1, localLineLengths.size() - 1);
    }
    else
    {
      tempNum = stoi(localLineLengths);
      localLineLengths = "";
    }
    numerator += (((tempNum - this->getMean()) * (tempNum - this->getMean())) / this->getNumLines());
  }
  return sqrt(numerator);
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

//Computes the relative probability of each nucleiotide bigram
double DNA::getRelProbBigram()
{

}

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

    outFS << setw(20) << left << "NUMBER OF LINES:";
    outFS << setw(5) << right << this->getNumLines() << endl;

    outFS << setw(20) << left << "MEAN:";
    outFS << setw(5) << right << this->getMean() << endl;

    outFS << setw(20) << left << "LINE LENGTHS STRING:";
    outFS << setw(5) << right << this->getLineLengths() << endl;

    outFS << setw(20) << left << "VARIANCE:";
    outFS << setw(5) << right << this->getVariance() << endl;

    outFS << setw(20) << left << "STANDARD DEVIATION:";
    outFS << setw(5) << right << this->getStdDeviation() << endl;

    outFS << setw(20) << left << "PROB A:";
    outFS << setw(5) << right << this->getProbA() << endl;
    outFS << setw(20) << left << "PROB C:";
    outFS << setw(5) << right << this->getProbC() << endl;
    outFS << setw(20) << left << "PROB T:";
    outFS << setw(5) << right << this->getProbT() << endl;
    outFS << setw(20) << left << "PROB G:";
    outFS << setw(5) << right << this->getProbG() << endl;

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
