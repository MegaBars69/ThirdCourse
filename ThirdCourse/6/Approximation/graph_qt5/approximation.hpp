#ifndef APROX_H
#define APROX_H

void ChebyshovAproximation(int n, double* f, double* a, double*g, double* z);
double ChebyshovValue(double X, double A, double B, int n, double* a);


#endif