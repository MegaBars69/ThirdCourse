#ifndef HEADER1  
#define HEADER1
#include <iostream>

void PrintMatrix(double* A,  int n, int m, int r = 25, bool exp_format = false, bool okruglenie = false);
void ApplyVector(double* X, double* A, int row_num, int col_num, int k, bool inside = true);
int InverseTriungleBlock(double* A, double* B, int n, double norm);
void ReplaceWith(double*A, double*B, int row_size, int col_size);
void MatrixMinusEqual(double* A, double* B, int row_size, int col_size);
void MatrixPlusEqual(double* A, double* B, int row_size, int col_size);
void DifferenceOfBlockMatrix(double* A, double* B, int n, int m);
void MinusEqualBlockMulEqualBlockMul(double* Main, double *a, double* b, int n1, int m12, int n2);

int Triungulize(double* A, double* U, int row_num, int col_num, double norm);
int InverseMatrix(double* A, double* B, double* U, double* ProductResult, double norm, int n, int m, int s);
void ApplyMatrix(double* U, double* A, int row_num, int col_num, int amount_of_vectors);
void ZeroOut(double* Diag, double* Down, double* U, int m, int row_size, double norm);
void ApplyMatrixToPair(double* U, double* Up, double* Down, int col_size, int row_size, int amount_of_vectors, bool down_is_zero = false);

double Norm(double* A, double* results, int n, int m);
void AddToMatrix(double* A, double* B, int n, int m);
void BlockMul(double *a, double* b, double* c, int n1, int m12, int n2);
void BlockMulOptimized(double *a, double* b, double* c, int n1, int m12, int n2);
double CalcError(double * A, double * InversedA, int n, int m);
double Discrepancy(double * A, double * InversedA, double* Column, double* ProductResult, double* Sum, int n, int m);
#endif 