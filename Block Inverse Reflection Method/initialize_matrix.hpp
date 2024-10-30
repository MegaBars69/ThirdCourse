#ifndef HEADER2  
#define HEADER2
#include <iostream>

int ReadMatrixFromFile(const std::string& filename, double* A,  int n, int m);
void BuildE(double* B, int n, int m);
void FormulaMatrixInitialization(double* A, int n, int m, int s);

#endif 