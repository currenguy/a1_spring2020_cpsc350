# Assignment 1: C++ Review

## 1. IDENTIFYING INFO
- Full Name: Curren Taber
- Student ID: 002325149
- Chapman Email: ctaber@chapman.edu
- Course Number and Section: CPSC 350-01
- Assignment 1: C++ Review

## 2. SOURCE FILES
- dna.h
- dna.cpp
- main.cpp

## 3. DESCRIPTION
- First file name is passed via command line argument, and following
  files are passed through user input. Typing "exit" exits the program.
- All results are written to "currentaber.out"
- This program will count every nucleotide in strings that only contain
  valid nucleotides. It will not count any letters of invalid strings.
- When counting nucleotide pairs, the pairs do not overlap. If the
  string has an odd length, then the last nucleotide does not have a
  pair.
- To generate the strings, the program calculates a length based on the
  Gaussian distribution of the input file. It determines the next letter
  by randomly choosing a letter/pair based on the probabilities
  previously calculated.

## 4. REFERENCES
- Data Structures & Algorithms
- Rene German
- https://www.techiedelight.com/convert-string-uppercase-cpp/
- https://stackoverflow.com/questions/7571326/why-does-dividing-two-int-not-yield-the-right-value-when-assigned-to-double
- http://www.cplusplus.com/forum/beginner/181119/
- https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
- https://www.geeksforgeeks.org/round-in-cpp/
- https://stackoverflow.com/questions/48716109/generating-a-random-number-between-0-1-c
