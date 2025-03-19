#ifndef S_APROX_H
#define S_APROX_H

void CalculateDiferences(double *dF, double *X, double *F, int n);
void CalculateCoeficients(double* c, double* X, double* d, double* dF, double* F, int n);
void CalculateParametrs(double* A, double*d, double*dF, double* X, double left_derivative, double right_derivative, int n);
double derivativeF(double (*f)(double), double x);
double SplineValue(double x, double a, double b, int n, double* c, double* X, int i = -1);
double derivative(double x, int func_id = -1, double (*f)(double) = nullptr);

#endif