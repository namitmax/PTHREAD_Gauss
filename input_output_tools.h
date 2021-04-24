#ifndef INCLUDE_INP_OUP_TOOLS
#define INCLUDE_INP_OUP_TOOLS

#include <cstdio>
#include <cstring>
#include <string>
#include <cstdlib>
#include <math.h>

void PrintBlock(const double *matrix, const int n, const int l,  const int r);
void Print(const double *matrix, const int l, const int n, const int m, const int r);
bool CheckVars(const int n, const int m, const int p, const int r, const int s, const char* filename);
int FileInput(double *matrix, const int n, const int m, const char* filename);
void InputMatrixWithAlg(double *matrix, double* inverse,
                        const int n, const int m, const int s,
                        const int threadNum, const int p);

#endif
