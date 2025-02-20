#ifndef HEADER2  
#define HEADER2
#include <iostream>
#include "algorithm.hpp"


int ReadMatrixFromFile(const std::string& filename, double* A,  int n, int m);
void BuildE(double* B, int n, int m, int p, int K);
void FormulaMatrixInitialization(double* A, int n, int m, int s, int p, int K);
double f(int i, int j, int n, int s);

#endif 