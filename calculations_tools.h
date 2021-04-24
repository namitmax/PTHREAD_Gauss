#ifndef INCLUDE_CALCULATIONS_TOOLS
#define INCLUDE_CALCULATIONS_TOOLS

#include <unistd.h>
#include <sched.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <cstdlib>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>

class Args {
public:
  double *matrix = 0;
  double *inverse = 0;
  int n = 0;
  int m = 0;
  int p = 0;
  int s = 0;
  int i = 0;
  int r = 0;
  int pt = 0;
  double *temp = 0;
  double *x = 0;
  double *result = 0;
  double *newTemp = 0;
  double *bigTemp = 0;
  int *numCols = 0;
  int *index = 0;
  int error = 0;
  int minI = -1;
  int minJ = -1;
  double resultAllNorm = -1.0;
  double minNorm = -1.0;
  double norm = 0.0;
  char *filename = 0;
  double threadTime = 0.;
  double totalTime = 0.;

  Args() = default;
};

void* Inverse(void* ptr);
double get_full_time();
#endif
