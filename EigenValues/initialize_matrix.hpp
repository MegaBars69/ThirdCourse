#ifndef HEADER2  
#define HEADER2

double f(int i, int j, int n, int s);
int ReadMatrixFromFile(char* filename, double* A,  int n);
void FormulaMatrixInitialization(double* A, int n, int s);
void PrintMatrix(double* matrix, int n, int m = 25);
int CheckMatrix(double* A, int n);

#endif