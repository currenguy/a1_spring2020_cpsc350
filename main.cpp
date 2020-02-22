#include "dna.h"

int main(int argc, char** argv)
{
  ifstream inFS;
  ofstream outFS;
  string line;
  string file = argv[1];
  DNA* dna = new DNA();

  //Checks if the user provides a file in the command line
  if (argc < 2)
  {
    cout << "No file provided." << endl;
    return 1;
  }

  dna->readFile(file);

  delete dna;

  return 0;
}
