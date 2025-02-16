#ifndef HEADER
#define HEADER
#include <string>
using namespace std;
void MatMul1(double *a, double* b, double* c, int n);
void Print(double* a, int n);
void MatMul2(double *a, double* b, double* c, int n);
void MatMul3(double *a, double* b, double* c, int n);
void MatMul4(double *a, double* b, double* c, int n);
void Test(double*A, double*B, double* C, int n, void (*func)(double*, double*, double*, int));
#endif