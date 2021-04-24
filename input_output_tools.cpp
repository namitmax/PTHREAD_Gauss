#include "input_output_tools.h"
#include <math.h>

#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))

void Print(const double *matrix, const int l, const int n, const int m, const int r) {
  int k = n / m;
  int l1 = (n != 1 ? n - (m * k) : 1);
  int cols = MIN(n, r);
  int rows = MIN(l, r);
  int count = 0;
  int tmp1 = 0, tmp2 = 0;
  for (int i = 0; i < rows * n; i += n * m) { //
    tmp1 = (i != n * m * k ? m : l1); //
    for (int p = 0; (p < tmp1) && (rows > p + (i / n)); p++) {//
      for (int j = 0; j < tmp1 * cols; j += m * tmp1) { //
        tmp2 = (j != tmp1 * k * m ? m : l1); 
        for (int q = 0; (q < tmp2) && ((cols > q + j / tmp1)); q++) {
	  count++;
	  printf(" %10.3e", matrix[p * tmp2 + i + j + q]);
	  if (count == cols) {
	    count = 0;
	    printf("\n");
	  }
	}
      }
    }
  }
}

void PrintBlock(const double *matrix, const int n, const int l,  const int r) {
  int cols = MIN(l, r);
  int rows = MIN(n, r);
  for (int i = 0; i < rows; i++) {
    for(int j = 0; j < cols; j++) {
      printf(" %10.3e", matrix[i * l + j]);
    }
    printf("\n");
  }
}

int FileInput(double *matrix, const int n, const int m, const char* filename) {
  FILE *IN;
  if ((IN = fopen(filename , "rt")) == NULL){
    return -1;
  }
  int k = n / m;
  int l = n - (m * k);
  int tmp1 = 0;
  int t = 0;
  for (int i = 0; i < n * n; i += n * m) {
    tmp1 = (i != k * m * n ? m : l);
    for (int p = 0; p < tmp1; p++) {
      for (int j = 0; j < n * tmp1; j += m * tmp1) {
        for (int q = 0; q < ((j == m * k * tmp1) ? l : m); q++) {
          if (fscanf(IN, "%lf", &matrix[p * ((j == k * tmp1 * m) ? l : m) + i + j + q]) != 1) {
            fclose(IN);
            return -2;
          }
          t++;
        }
      }
    }
  }
  fclose(IN);
  return 0;
}

void SetInverse(double *solution, const int n, const int m,
               const int threadNum, const int p_) {
  int k = n / m;
  int l = n - (m * k);
  int tmp1 = 0;
  int t = threadNum * m * tmp1;
  for (int i = 0; i < n * n; i += n * m) {
    tmp1 = (i != k * m * n ? m : l);
    for (int p = 0; p < tmp1; p++) { // берем строку блока
      for (int j = threadNum * m * tmp1; j < n * tmp1; j += m * tmp1 * p_) { // берем блок из строки блоков
        t = i + p * n + ( j / (m * tmp1 * p_)) * m * p_ + threadNum * m;
        for (int q = 0; q < ((j == m * k * tmp1) ? l : m); q++) { // заполняем строку блока
          solution[p * ((j == k * tmp1 * m) ? l : m) + i + j + q] = (t / n != t % n ? 0. : 1.);
          t++;
        }
      }
    }
  }
}


void InputAlg0(double *matrix, const int n, const int m,
               const int threadNum, const int p_) {
  int k = n / m;
  int l = n - (m * k);
  int tmp1 = 0;
  int t = threadNum * m * tmp1;
  for (int i = 0; i < n * n; i += n * m) {
    tmp1 = (i != k * m * n ? m : l);
    for (int p = 0; p < tmp1; p++) { // берем строку блока
      for (int j = threadNum * m * tmp1; j < n * tmp1; j += m * tmp1 * p_) { // берем блок из строки блоков
        t = i + p * n + ( j / (m * tmp1 * p_)) * m * p_ + threadNum * m;
        for (int q = 0; q < ((j == m * k * tmp1) ? l : m); q++) { // заполняем строку блока
          matrix[p * ((j == k * tmp1 * m) ? l : m) + i + j + q] = 0.;
          t++;
        }
      }
    }
  }
}

int InputAlg1(double *matrix, const int n, const int m,
              const int threadNum, const int p_) {
  int k = n / m;
  int l = n - (m * k);
  int tmp1 = 0;
  int t = threadNum * m * tmp1;
  for (int i = 0; i < n * n; i += n * m) {
    tmp1 = (i != k * m * n ? m : l);
    for (int p = 0; p < tmp1; p++) { // берем строку блока
      for (int j = threadNum * m * tmp1; j < n * tmp1; j += m * tmp1 * p_) { // берем блок из строки блоков
        t = i + p * n + ( j / (m * tmp1 * p_)) * m * p_ + threadNum * m;
	for (int q = 0; q < ((j == m * k * tmp1) ? l : m); q++) { // заполняем строку блока
	  matrix[p * ((j == k * tmp1 * m) ? l : m) + i + j + q] = n - (t / n < t % n ? t % n : t / n);
	  t++;
        }
      }
    }
  }
  return 0;
}

int InputAlg2(double *matrix, const int n, const int m,
              const int threadNum, const int p_) {
  int k = n / m;
  int l = n - (m * k);
  int tmp1 = 0;
  int t = threadNum * m * tmp1;
  for (int i = 0; i < n * n; i += n * m) {
    tmp1 = (i != k * m * n ? m : l);
    for (int p = 0; p < tmp1; p++) { // берем строку блока
      for (int j = threadNum * m * tmp1; j < n * tmp1; j += m * tmp1 * p_) { // берем блок из строки блоков
        t = i + p * n + ( j / (m * tmp1 * p_)) * m * p_ + threadNum * m;
        for (int q = 0; q < ((j == m * k * tmp1) ? l : m); q++) { // заполняем строку блока
          matrix[p * ((j == k * tmp1 * m) ? l : m) + i + j + q] = (t / n < t % n ? t % n : t / n) + 1;
          t++;
        }
      }
    }
  }
  return 0;
}

int InputAlg3(double *matrix, const int n, const int m,
              const int threadNum, const int p_) {
  int k = n / m;
  int l = n - (m * k);
  int tmp1 = 0;
  int t = threadNum * m * tmp1;
  for (int i = 0; i < n * n; i += n * m) {
    tmp1 = (i != k * m * n ? m : l);
    for (int p = 0; p < tmp1; p++) { // берем строку блока
      for (int j = threadNum * m * tmp1; j < n * tmp1; j += m * tmp1 * p_) { // берем блок из строки блоков
        t = i + p * n + ( j / (m * tmp1 * p_)) * m * p_ + threadNum * m;
        for (int q = 0; q < ((j == m * k * tmp1) ? l : m); q++) { // заполняем строку блока
          matrix[p * ((j == k * tmp1 * m) ? l : m) + i + j + q] = fabs(t / n - t % n);
          t++;
        }
      }
    }
  }
  return 0;
}

int InputAlg4(double *matrix, const int n, const int m,
              const int threadNum, const int p_) {
  int k = n / m;
  int l = n - (m * k);
  int tmp1 = 0;
  int t = threadNum * m * tmp1;
  for (int i = 0; i < n * n; i += n * m) {
    tmp1 = (i != k * m * n ? m : l);
    for (int p = 0; p < tmp1; p++) { // берем строку блока
      for (int j = threadNum * m * tmp1; j < n * tmp1; j += m * tmp1 * p_) { // берем блок из строки блоков
        t = i + p * n + ( j / (m * tmp1 * p_)) * m * p_ + threadNum * m;
        for (int q = 0; q < ((j == m * k * tmp1) ? l : m); q++) { // заполняем строку блока
          matrix[p * ((j == k * tmp1 * m) ? l : m) + i + j + q] = 1.0 / (t / n + t % n + 1);
          t++;
        }
      }
    }
  }
  return 0;
}

bool CheckVars(const int n, const int m, const int p, const int r, const int s, const char* filename) {
  return (p > 0) && (n > 0) && (m > 0) && (r > 0) && (n >= m) && (s > -1) && (s < 5) && (!(s > 0) || (filename == 0)) && (m % 3 == 0);
}

void InputMatrixWithAlg(double *matrix, double* inverse,
                        const int n, const int m, const int s, 
                        const int threadNum, const int p) {
  SetInverse(inverse, n, m, threadNum, p);
  switch (s) {
    case 0:
      InputAlg0(matrix, n, m, threadNum, p);
      break;
    case 1:
      InputAlg1(matrix, n, m, threadNum, p);
      break;
    case 2:
      InputAlg2(matrix, n, m, threadNum, p);
      break;
    case 3:
      InputAlg3(matrix, n, m, threadNum, p);
      break;
    case 4:
      InputAlg4(matrix, n, m, threadNum, p);
      break;
    default:
      return;
      break;
  }
  return;
}
