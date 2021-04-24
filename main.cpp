#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include "calculations_tools.h"
#include <unistd.h>
#include <sched.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include "input_output_tools.h"

//#include <fstream>

int main(int argc, char* argv[]) {
  int    n;             // размерность матрицы
  int    m;             // размер блока
  int    r;             // кол-во выводимых значений в матрице
  int    s;             // номер применяемого алгоритма
  int    p;             // число потоков
  double norma;         // норма невязки
  char   *filename;     // имя файла
  int p_true;
  double *matrix;       // исходная матрица
  double *solution;     // решение
  //std::ofstream output("tests/log.txt", std::ios::app);
  if (!(argc == 6 || argc == 7)) {
    printf("Invalid number of parametrs,  usage: n, m, p, r, s, filename \n");
    printf("\n");
    return -1;
  }
  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &m);
  sscanf(argv[3], "%d", &p);
  sscanf(argv[4], "%d", &r);
  sscanf(argv[5], "%d", &s);
  if (argv[6])
    filename = argv[6];
  else 
    filename = 0;
  if (!CheckVars(n, m, p, r, s, filename)) {
    printf("Variables error, usage: n, m, p, r, s, filename \n");
    printf("\n");
    return -1;
  }
  matrix = new double[n * n];
  solution = new double [n * n];
  p_true = p;
  Args *args = new Args[p];
  pthread_t *tids = new pthread_t[p];
  int div = n / m;
  cpu_set_t cpu;
  int tp = (div * m != n ? (div + 1): div);
  int* numCols = new int[tp];
  int err = 0;
  for (int i = 0; i < tp; i++) {
    numCols[i] = i;
  }
  p = (tp >= p ? p : tp);
  double* bigTemp = new double[n * m];
  for (int i = 0; i < (p != 0 ? p : 1); i++) {
    args[i].n = n;
    args[i].m = m;
    args[i].p = p;
    args[i].pt = p_true;
    args[i].s = s;
    args[i].i = i;
    args[i].r = r;
    args[i].matrix = matrix;
    args[i].inverse = solution;
    args[i].error = 0;
    args[i].filename = filename;
    args[i].minI = i;
    args[i].minJ = 0;
    args[i].minNorm = -1;
    args[i].numCols = numCols;
    args[i].result = new double[m * m];
    args[i].temp = new double[m * m];
    args[i].x = new double[m * m];
    args[i].index = new int[m];
    args[i].newTemp = new double[m * m];
    args[i].bigTemp = bigTemp; 
  }
  for (int i = 1; i < p; i++) {
    if (pthread_create(tids + i, 0, &Inverse, args + i)) {
      printf("Can not create thread %d\n", i);
    }
    CPU_ZERO(&cpu);
    CPU_SET(i, &cpu);
    pthread_setaffinity_np(tids[i], sizeof(cpu), &cpu);
  }
  CPU_ZERO(&cpu);
  CPU_SET(0, &cpu);
  tids[0] = pthread_self();
  pthread_setaffinity_np(tids[0], sizeof(cpu), &cpu);
  Inverse(args + 0);
  for (int i = 1; i < p; i++) {
    pthread_join(tids[i], 0);
  }
  for (int i = 0; i < p; i++) {
    delete [] args[i].result;
    delete [] args[i].temp;
    delete [] args[i].x;
    delete [] args[i].index;
    delete [] args[i].newTemp;
  }
  delete [] numCols;
  for (int i = 0; i < p; i++) {
    if (args[i].error != 0) {
      err = args[i].error;
      printf("\n");
      if (err == -2) {
        printf("\n");
        printf("Matrix is not invertible !!!\n");
      } else {
        printf("\n");
        printf("some error was occured, may be init file is bad\n");
        printf("\n");
      }
      delete [] matrix;
      delete [] solution;
      delete [] args;
      delete [] tids;
      delete [] bigTemp;
      return err;
    }
  }
  for (int i = 1; i < p + 1; i++) {
    printf("time for thread %d = %.2f \n", i, args[i - 1].threadTime);
  }
  norma = args[0].resultAllNorm;
  double time = args[0].totalTime;
  delete [] args;
  delete [] tids;
  /////////////////////////////////
  printf ("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d p = %d\n", argv[0], norma, time, s, n, m, p_true);
  //output << norma;
  //output << "\n";
  printf("\n");
  delete [] matrix;
  delete [] bigTemp;
  delete [] solution;
  return 0;
}
