#include "input_output_tools.h"
#include "calculations_tools.h"

#define MIN(a, b) ((a) < (b)?(a):(b))

static inline double Norma(const double* matrix, const int n, const int m) {
  double sum = 0;
  double max = 0;
  for (int i = 0; i < n; i++) {
    sum = 0;
    for (int j = 0; j < m; j++)
      sum += fabs(matrix[i * n + j]);
    if (sum > max)
      max = sum;
  }
  return max;
}

static inline void GetBlock(const double* matrix, double* block,
              const int i, const int j,
              const int d, const int l,
              const int m, const int n) {
  int k = n / m;
  int l1 = n - m * k;
  int start = (i != k ? i * n * m + j * m * m : i * n * m + j * l1 * m);
  if (d == l && d % 3 == 0) {
    for (int t = 0;  t < d; t += 3)
      for (int p = 0; p < l; p += 3) {
        block[t * l + p]           = matrix[start + t * l + p];
	block[t * l + p + 1]       = matrix[start + t * l + p + 1];
	block[t * l + p + 2]       = matrix[start + t * l + p + 2];
	block[(t + 1) * l + p]     = matrix[start + (t + 1) * l + p];
	block[(t + 1) * l + p + 1] = matrix[start + (t + 1) * l + p + 1];
	block[(t + 1) * l + p + 2] = matrix[start + (t + 1) * l + p + 2];
	block[(t + 2) * l + p]     = matrix[start + (t + 2) * l + p];
	block[(t + 2) * l + p + 1] = matrix[start + (t + 2) * l + p + 1];
	block[(t + 2) * l + p + 2] = matrix[start + (t + 2) * l + p + 2];
      }
  } else {
    for (int t = 0;  t < d; t++)
      for (int p = 0; p < l; p++) {
	block[t * l + p] = matrix[start + t * l + p];
      }
  }
}

static inline int InverseBlock(double *a, double *x, const int n, const double tempNorm, int* index) {
  int indMax1;
  int indMax2;
  int i, j, k;
  double tmp;
  int tmp1;
  double max;
  for (i = 0; i < n; i++)
    index[i] = i;
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      x[i * n + j] = (double)(i == j);
  for (i = 0; i < n; ++i) {
    max = -1;
    indMax1 = 0;
    indMax2 = 0;
    tmp1 = (n > 40 ? i + 1 : n);
    for (j = i; j < tmp1; ++j)
      for (k = i; k < n; ++k)
        if (max < fabs(a[j * n + k])) {
          max = fabs(a[j * n + k]);
          indMax1 = j;
          indMax2 = k;
        }
    if (max <= fabs(1e-16 * tempNorm))
      return -1;
    if (indMax1 != i) {
      for (j = 0; j < n; ++j) {
          tmp = a[i * n + j];
          a[i * n + j] = a[indMax1 * n + j];
          a[indMax1 * n + j] = tmp;
      }
      for (j = 0; j < n; ++j) {
        tmp = x[i * n + j];
        x[i * n + j] = x[indMax1 * n + j];
        x[indMax1 * n + j] = tmp;
      }
    }
    if (indMax2 != i) {
      k = index[i];
      index[i] = index[indMax2];
      index[indMax2] = k;
      for (j = 0; j < n; ++j) {
        tmp = a[j * n + i];
        a[j * n + i] = a[j * n + indMax2];
        a[j * n + indMax2] = tmp;
      }
    }
    tmp = 1.0 / a[i * n + i];
    for (j = i; j < n; ++j)
      a[i * n + j] *= tmp;
    for (j = 0; j < n; ++j)
      x[i * n + j] *= tmp;
    for (j = i + 1; j < n; ++j) {
      tmp = a[j * n + i];
      for (k = i; k < n; ++k)
        a[j * n + k] -= a[i * n + k] * tmp;
      for (k = 0; k < n; ++k)
        x[j * n + k] -= x[i * n + k] * tmp;
    }
  }
  for (k = 0; k < n; ++k)
    for (i = n - 1; i >= 0; --i) {
      tmp = x[i * n + k];
      for (j = i + 1; j < n; ++j)
        tmp -= a[i * n + j] * x[j * n + k];
      x[i * n + k] = tmp;
    }
  for (i = 0; i < n; ++i) {
    k = index[i];
    for (j = 0; j < n; ++j) 
      a[k * n + j] = x[i * n + j];
  }
   return 0;
}

static inline void BlockMultiplicationForMult(const double* block1, double* block2, double* result,
                         const int n, const int l, const int m) {
  int r, t, q;
  double sum = 0;
  double c[9];
  if (n % 3 == 0 && n == m && m == l) {
    for (r = 0; r < n; r += 3) {
      for (t = 0; t < m; t += 3) {
        c[0] = 0.;
        c[1] = 0.;
        c[2] = 0.;
        c[3] = 0.;
        c[4] = 0.;
        c[5] = 0.;
        c[6] = 0.;
        c[7] = 0.;
        c[8] = 0.;
        for (q = 0; q < l; ++q) {
          c[0] += block1[r * l + q] * block2[q * m + t];
          c[1] += block1[r * l + q] * block2[q * m + t + 1];
          c[2] += block1[r * l + q] * block2[q * m + t + 2];
          c[3] += block1[(r + 1) * l + q] * block2[q * m + t];
          c[4] += block1[(r + 1) * l + q] * block2[q * m + t + 1];
          c[5] += block1[(r + 1) * l + q] * block2[q * m + t + 2];
          c[6] += block1[(r + 2) * l + q] * block2[q * m + t];
          c[7] += block1[(r + 2) * l + q] * block2[q * m + t + 1];
          c[8] += block1[(r + 2) * l + q] * block2[q * m + t + 2];
        }
        result[r * n + t]             += c[0];
        result[r * n + t + 1]         += c[1];
        result[r * n + t + 2]         += c[2];
        result[(r + 1) * n + t]       += c[3];
        result[(r + 1) * n + (t + 1)] += c[4];
        result[(r + 1) * n + (t + 2)] += c[5];
        result[(r + 2) * n + t]       += c[6];
        result[(r + 2) * n + (t + 1)] += c[7];
        result[(r + 2) * n + (t + 2)] += c[8];
      }
    }
  }
  else {
    sum = 0;
    for(int i = 0; i < n; i++)
      for(int j = 0; j < m; j++) {
        sum = 0;
        for(int k = 0; k < l; k++)
          sum += block1[i * l + k] * block2[j + k * m];
        result[i * m + j] += sum;
      }
   }
}

static inline void BlockMultiplication(const double* block1, const double* block2, double* result,
                         const int n, const int l, const int m) {
  int r, t, q;
  double sum = 0;
  double c[9];
  if (n == l &&  l == m && n % 3 == 0) {
    for (r = 0; r < n; r += 3) {
      for (t = 0; t < m;  t += 3) {
        c[0] = 0.;
	c[1] = 0.;
	c[2] = 0.;
	c[3] = 0.;
	c[4] = 0.;
	c[5] = 0.;
	c[6] = 0.;
	c[7] = 0.;
	c[8] = 0.;
        for (q = 0; q < l; ++q) {
          c[0] += block1[r * l + q] * block2[q * m + t];
          c[1] += block1[r * l + q] * block2[q * m + t + 1];
	  c[2] += block1[r * l + q] * block2[q * m + t + 2];
          c[3] += block1[(r + 1) * l + q] * block2[q * m + t];
          c[4] += block1[(r + 1) * l + q] * block2[q * m + t + 1];
	  c[5] += block1[(r + 1) * l + q] * block2[q * m + t + 2];
          c[6] += block1[(r + 2) * l + q] * block2[q * m + t];
	  c[7] += block1[(r + 2) * l + q] * block2[q * m + t + 1];
	  c[8] += block1[(r + 2) * l + q] * block2[q * m + t + 2];
        }
        result[r * n + t]             = c[0];
        result[r * n + t + 1]         = c[1];
	result[r * n + t + 2]         = c[2];
        result[(r + 1) * n + t]       = c[3];
        result[(r + 1) * n + (t + 1)] = c[4];
	result[(r + 1) * n + (t + 2)] = c[5];
	result[(r + 2) * n + t]       = c[6];
        result[(r + 2) * n + (t + 1)] = c[7];
        result[(r + 2) * n + (t + 2)] = c[8];
      }
    }
  }
  else {
    for(int i = 0; i < n; i++)
      for(int j = 0; j < m; j++) {
        sum = 0;
        for(int k = 0; k < l; k++)
          sum += block1[i * l + k] * block2[j + k * m];
        result[i * m + j] = sum;
      }
   }
}

static inline void PutBlock(double* matrix, const double* block,
              const int i, const int j,
              const int d, const int l,
              const int m, const int n) {
  int k = n / m;
  int l1 = n - m * k;
  int start = (i != k ? i * n * m + j * m * m : i * n * m + j * l1 * m);
  if (d == l && d % 3 == 0) {
    for (int t = 0;  t < d; t += 3)
      for (int p = 0; p < l; p += 3) {
        matrix[start + t * l + p]           = block[t * l + p];
        matrix[start + t * l + p + 1]       = block[t * l + p + 1];
        matrix[start + t * l + p + 2]       = block[t * l + p + 2];
        matrix[start + (t + 1) * l + p]     = block[(t + 1) * l + p];
        matrix[start + (t + 1) * l + p + 1] = block[(t + 1) * l + p + 1];
        matrix[start + (t + 1) * l + p + 2] = block[(t + 1) * l + p + 2];
        matrix[start + (t + 2) * l + p]     = block[(t + 2) * l + p];
        matrix[start + (t + 2) * l + p + 1] = block[(t + 2) * l + p + 1];
        matrix[start + (t + 2) * l + p + 2] = block[(t + 2) * l + p + 2];
      }
  } else {
    for (int t = 0;  t < d; t++)
      for (int p = 0; p < l; p++)
        matrix[start + t * l + p] = block[t * l + p];
  }
}

static inline void ChangeBlock(double* matrix,
                 const int i1, const int j1,
                 const int i2, const int j2,
                 const int n, const int m,
                 const int l, const int d) {
  double buff;
  int k = n / m;
  int l1 = n - m * k;
  int start1 = (i1 != k ? i1 * n * m + j1 * m * m : i1 * n * m + j1 * l1 * m);
  int start2 = (i2 != k ? i2 * n * m + j2 * m * m : i2 * n * m + j2 * l1 * m);
  for(int t = 0; t < d; t++)
    for(int p = 0; p < l; p++) {
      buff = matrix[start1 + t * l + p];
      matrix[start1 + t * l + p] = matrix[start2 + t * l + p];
      matrix[start2 + t * l + p] = buff;
    }
}

static inline void ChangeRow(double* matrix, 
               const int i, const int j, 
               const int m, const int n,
	       const int threadNum, const int p) {
  int k = n / m;
  int last = (k != 0 ? k + 1 : k);
  for(int t = threadNum; t < last; t += p) {
    if (t != k) {
      ChangeBlock(matrix, i, t, j, t, n, m, m, m);
    }
    else {
      ChangeBlock(matrix, i, k, j, k, n, m, n - k * m, m);
    }
  }
}

static inline void DiffBlocks(double* matrix1, const double* matrix2,
               const int d, const int l) {
 if (d == l && d % 3 == 0) {
  for (int i = 0; i < d; i += 3)
    for (int j = 0; j < l; j += 3) {
      matrix1[i * l + j] -= matrix2[i * l + j];
      matrix1[i * l + j + 1] -= matrix2[i * l + j + 1];
      matrix1[i * l + j + 2] -= matrix2[i * l + j + 2];
      matrix1[(i + 1) * l + j] -= matrix2[(i + 1) * l + j];
      matrix1[(i + 1) * l + j + 1] -= matrix2[(i + 1) * l + j + 1];
      matrix1[(i + 1) * l + j + 2] -= matrix2[(i + 1) * l + j + 2];
      matrix1[(i + 2) * l + j ] -= matrix2[(i + 2) * l + j];
      matrix1[(i + 2) * l + j + 1] -= matrix2[(i + 2) * l + j + 1];
      matrix1[(i + 2) * l + j + 2] -= matrix2[(i + 2) * l + j + 2];
    }
 } else {
   for (int i = 0; i < d; i++)
    for (int j = 0; j < l; j++) {
      matrix1[i * l + j] -= matrix2[i * l + j];
    }
 }
}

static inline void ChangeCol(double* matrix, 
               const int i, const int j, 
               const int m, const int n) {
  int k = n / m;
  for(int t = 0; t < k; t++) {
    ChangeBlock(matrix, t, i, t, j, n, m, m, m);
  }
  if (n % m != 0)
    ChangeBlock(matrix, k, i, k, j, n, m, m, n- k* m);
}

static inline void RowMultiplication(double* matrix, const double* x,
                       const int t, const int p, const int d,
                       const int n, const int m,
                       const int threadNum, const int ptr,
                       double* temp, double* result) {
  int div = n / m;
  int t1;
  t1 = (div * m != n ? (div + 1): div);
  int tmp1 = (t != div ? d : n - div * m);
  int l1 = n - div * m;
  int tmp2;
  int jStart = ((threadNum + ptr * (p / ptr) >= p) ? (p / ptr) * ptr + threadNum : (p / ptr + 1) * ptr + threadNum);
  for (int j = jStart; j < t1; j += ptr) { // умножаем обращенный блок на остальные в строке
    tmp2 = (j != div ? m : n - div * m);
    GetBlock(matrix, temp, t, j, tmp1, tmp2, m, n);
    result = matrix + (t != div ? t * n * m + j * m * m : t * n * m + j * l1 * m);
    BlockMultiplication(x, temp, result, tmp1, tmp1, tmp2);
  }
}

static inline void MiniSecondPass(double* matrix, const double* values, const int num, 
		const int i, const int d, const int n, const int m, 
		double* temp1, double* temp2, double* newTemp) {
  int div = n / m;
  int l = n - m * div;
  int j = num;
  int tempDim;
  tempDim = (j != div ? m : l);
  temp1 = matrix + (num != div ? num * n * m + j * m * m : num * n * m + j * l * m);
  BlockMultiplication(values, temp1, newTemp, m, d, tempDim);
  temp2 = matrix + (i != div ? i * n * m + j * m * m : i * n * m + j * l * m);
  DiffBlocks(temp2, newTemp, m, tempDim);
}

static inline void SecondPass(double* matrix, const double* values, const int num, 
		const int p, const int i, const int d, 
		const int n, const int m, const int threadNum, 
		const int ptr, double* temp1, double* temp2, double *newTemp,
		bool f) {
  int div = n / m;
  int l = n - m * div;
  int tp = (l != 0 ? (div + 1): div);
  int jStart = ((threadNum + ptr * (p / ptr) >= p) ? (p / ptr) * ptr + threadNum : (p / ptr + 1) * ptr + threadNum);
  if (jStart == p && (f)) { // чтобы не стереть values
    jStart += ptr;
  }
  int tempDim;
  for (int j = jStart; j < tp; j += ptr) {
    tempDim = (j != div ? m : l);
    temp1 = matrix + (num != div ? num * n * m + j * m * m : num * n * m + j * l * m);
    temp2 = matrix + (i != div ? i * n * m + j * m * m : i * n * m + j * l * m);
    BlockMultiplication(values, temp1, newTemp, m, d, tempDim);
    DiffBlocks(temp2, newTemp, m, tempDim);
  }
}


static inline void ReverseGauss(double* matrix, double *inverse, 
                  const int num, const int m, const int n, const int threadNum, 
		  const int p, double* values, double* temp1, double* temp2, double* newTemp) {
  int div = n / m;
  int l = n - div * m;
  int tempDim;
  for (int i = num - 1; i > -1; i--) {
    tempDim = (num != div ? m : l);
    values = matrix + (i != div ? i * n * m + num * m * m : i * n * m + num * l * m);
    SecondPass(matrix, values, num, num, i, tempDim, n, m, threadNum, p, temp1, temp2, newTemp, true);
    SecondPass(inverse, values, num, 0, i, tempDim, n, m,  threadNum, p, temp1, temp2, newTemp, false);
  }
}

static inline void MiniFirstPass(double* matrix, const double* values, const int num, const int i, const int d, const int n, const int m, 
                   double* temp1, double* temp2, double* newTemp) {
  int div = n / m;
  int l = n - m * div;
  int j = num;
  int tempDim = (j != div ? m : l);
  temp1 = matrix + (num != div ? num * n * m + j * m * m : num * n * m + j * l * m);
  temp2 = matrix + (i != div ? i * n * m + j * m * m : i * n * m + j * l * m);
  BlockMultiplication(values, temp1, newTemp, d, m, tempDim);
  DiffBlocks(temp2, newTemp, d, tempDim);
}

static inline void FirstPass(double* matrix, const double* values, const int num, const int p, const int i, const int d, const int n, const int m, const int threadNum, const int ptr, 
		double* temp1, double* temp2, double* newTemp, const bool f) {
  int div = n / m;
  int l = n - m * div;
  int tp = (l != 0 ? (div + 1): div);
  int jStart = ((threadNum + ptr * (p / ptr) >= p) ? (p / ptr) * ptr + threadNum : (p / ptr + 1) * ptr + threadNum);
  if ((jStart == p && f)) { // чтобы не стереть values
    jStart += ptr;
  }
  int tempDim;
  for (int j = jStart; j < tp; j += ptr) {
    tempDim = (j != div ? m : l);
    temp1 = matrix + (num != div ? num * n * m + j * m * m : num * n * m + j * l * m);
    temp2 = matrix + (i != div ? i * n * m + j * m * m : i * n * m + j * l * m);
    BlockMultiplication(values, temp1, newTemp, d, m, tempDim);
    DiffBlocks(temp2, newTemp, d, tempDim);
  }
}

static inline void ForwardGauss(double* matrix, double *inverse, const int num, const int m, const int n, const int threadNum, const int p, 
                  double* values, double* temp1, double* temp2, double* newTemp) {
  int div = n / m;
  int l = n - div * m;
  int tp = (l != 0 ? (div + 1): div);
  int tempDim;
  for (int i = num + 1; i < tp; i++) {
    tempDim = (i != div ? m : l);
    values = matrix + (i != div ? i * n * m + num * m * m : i * n * m + num * l * m);
    FirstPass(matrix, values, num, num, i, tempDim, n, m, threadNum, p, temp1, temp2, newTemp, true);
    FirstPass(inverse, values, num, 0, i, tempDim, n, m, threadNum, p, temp1, temp2, newTemp, false);
  }
}

static inline void synchronizeMult(double* a, double* c,
                                   const int threadNum, const int p,
                                   const int m, const int n, const int i) {
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
  static int threadsIn = 0;
  static int threadsOut = 0;
  int totalThreads = p;
  int k = n / m;
  int l = n - (m * k);
  int last = (l > 0 ? k + 1 : k);
  int a_hor;
  int a_vert = (i != k ? m : l);
  double *pa, *pc;
  pthread_mutex_lock(&mutex);
  ++threadsIn;
  if (threadsIn >= totalThreads) {
    threadsOut = 0;
    pthread_cond_broadcast(&condvar_in);
  } else {
    while (threadsIn < totalThreads) {
      pthread_cond_wait(&condvar_in, &mutex);
    }
  }
  ++threadsOut;
  if (threadsOut >= totalThreads) {
    threadsIn = 0;
    pthread_cond_broadcast(&condvar_out);
  } else {
    while (threadsOut < totalThreads) {
      pthread_cond_wait(&condvar_out, &mutex);
    }
  }
  for (int q = threadNum; q < last; q += p) {
    a_hor = (q < k ? m : l);
    pa = a + q * m * a_vert;
    pc = c + q * m * m;
    for (int s = 0; s < a_hor * a_vert; s++) {
      pa[s] = pc[s];
    }
  }
  pthread_mutex_unlock(&mutex);
}

static inline void mult(double *a, double *b, 
                        int n, int m, int threadNum, int p,
                        double* bigTemp) {
  int s, k, l, last, a_vert, a_hor, b_hor;
  double *pa, *pb, *pc, *c;
  c = bigTemp;
  k = n / m;
  l = n - (m * k);
  last = (l > 0 ? k + 1 : k);
  for (int i = 0; i < last; i++) {
    a_vert = (i < k ? m : l);
    pc = c;
    for (int j = threadNum; j < last; j += p) {
      b_hor = (j < k ? m : l);
      pc = c + m * m * j;
      for (int s = 0; s < b_hor * a_vert; s++) {
        pc[s] = 0;
      }
      for (s = 0; s < last; s ++) {
        a_hor = (s < k ? m : l);
	pa = a + (i * n * m) + (s * a_vert * m);
        pb = b + (s * n * m) + (j * a_hor * m);
	BlockMultiplicationForMult(pa, pb, pc, a_vert, a_hor, b_hor);
      }
    }
    pa = a + (i * n * m);
    pc = c;
    synchronizeMult(pa, pc, threadNum, p, m, n, i);
  }
}

static inline void SumBlocks(double* matrix1, const double* matrix2, 
               const int d, const int l) {
  for (int i = 0; i < d; i++)
    for (int j = 0; j < l; j++)
      matrix1[i] += fabs(matrix2[i * d + j]);
}

static inline double ResultNorma(const double* matrix, const int n, const int m, const int x) {
  double* temp = (double*)malloc (m * sizeof(double)) ;
  double* temp1 = (double*)malloc (m * m * sizeof(double));
  double max = 0.;
  int k = n / m;
  int l = n - k * m;
  int last = (l != 0 ? k + 1 : k);
  int t;
  for (int i = 0; i < last; i++) {
    for (int j = 0; j < m; j++)
      temp[j] = 0;
    t = (i != k ? m : l);
    for (int j = 0; j < last; j++) {
      if (j != k) {
        GetBlock(matrix, temp1, i, j, t, m, m, n);
        SumBlocks(temp, temp1, t, m);
      } 
      else {
       GetBlock(matrix, temp1, i, j, t, l, m, n);
       SumBlocks(temp, temp1, t, l);
      }
    }
    for (int j = 0; j < t; j++) {
      if ((temp[j] - x) > max) {
        max = (temp[j] - x);
      }
    }
  }
  free(temp);
  free(temp1);
  return max;
}

static inline double get_time() {
  struct rusage buf;
  getrusage(RUSAGE_THREAD,&buf);
  return (double)buf.ru_utime.tv_sec + (double)buf.ru_utime.tv_usec/1000000.;
}

double get_full_time() {
  struct timeval buf;
  gettimeofday(&buf, 0);
  return (double)buf.tv_sec + (double)buf.tv_usec/1000000.;
}

static inline void synchronize(Args* a) {
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
  static int threadsIn = 0;
  static int threadsOut = 0;
  int totalThreads = a->p;
  pthread_mutex_lock(&mutex);
  ++threadsIn;
  if (threadsIn >= totalThreads) {
    threadsOut = 0;
    pthread_cond_broadcast(&condvar_in);
  } else {
    while (threadsIn < totalThreads) {
      pthread_cond_wait(&condvar_in, &mutex);
    }
  }
  ++threadsOut;
  if (threadsOut >= totalThreads) {
    threadsIn = 0;
    pthread_cond_broadcast(&condvar_out);
  } else {
    /*
    printf("MATRIX AND INVERSE AFTER MULT \n");
    Print(a->matrix, a->n, a->n, a->m, 20);
    printf("\n");
    Print(a->inverse, a->n, a->n, a->m, 20);
    */
    while (threadsOut < totalThreads) {
      pthread_cond_wait(&condvar_out, &mutex);
    }
  }
  pthread_mutex_unlock(&mutex);
}

static inline void synchronizeRead(Args* a, 
                     double* matrix, const int n, const int m,
                     const int s, const char* filename, const bool flag) {
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
  static int threadsIn = 0;
  static int threadsOut = 0;
  static int error = 0;
  static double matrixNorm = 0.;
  static double norma = -1.0;
  static bool outputFirst = true;
  static bool outputLast = true;
  static double time;
  //static int error;
  int totalThreads = a->p;
  pthread_mutex_lock(&mutex);
  ++threadsIn;
  if (threadsIn >= totalThreads) {
    threadsOut = 0;
    pthread_cond_broadcast(&condvar_in);
  } else {
    while (threadsIn < totalThreads) {
      pthread_cond_wait(&condvar_in, &mutex);
    }
  }
  if (!flag) {
    if (outputLast) {
      outputLast = false;
      time = get_full_time() - a->totalTime;
      printf("\n");
      printf("-----------SOLUTION----------\n");
      Print(a->matrix, n, n, a->m, a->r);
      printf("\n");
    }
  }
  ++threadsOut;
  if (threadsOut >= totalThreads) {
    threadsIn = 0;
    if (s == 0) {
      error = FileInput(a->matrix, n, m, filename);
    }
    if (error == 0) {
      if (flag) {
        if (outputFirst) {
          outputFirst = false;
	  matrixNorm = ResultNorma(matrix, n, m, 0);
          printf("\n");
          printf("-----------ORIGINAL MATRIX----------\n");
          Print(matrix, n, n, a->m, a->r);
        }
	time = get_full_time();
      }
    }
    pthread_cond_broadcast(&condvar_out);
  } else {
    while (threadsOut < totalThreads) {
      pthread_cond_wait(&condvar_out, &mutex);
    }
  }
  a->error = error;
  a->norm = matrixNorm;
  if (norma > -1) {
    a->resultAllNorm = norma;
  }
  a->totalTime = time;
  pthread_mutex_unlock(&mutex);
}

static inline void synchronizeGetBlock(Args* a,
                         double* matrix,
                         const int t) {
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
  static int threadsIn = 0;
  static int threadsOut = 0;
  int totalThreads = a->p;
  int div = a->n / a->m;
  int tmp1;
  pthread_mutex_lock(&mutex);
  ++threadsIn;
  if (threadsIn >= totalThreads) {
    threadsOut = 0;
    pthread_cond_broadcast(&condvar_in);
  } else {
    while (threadsIn < totalThreads) {
      pthread_cond_wait(&condvar_in, &mutex);
    }
  }
  tmp1 = (t != div ? a->m : a->n - div * a->m);
  GetBlock(matrix, a->temp, t, t, tmp1, tmp1, a->m, a->n);
  ++threadsOut;
  if (threadsOut >= totalThreads) {
    threadsIn = 0;
    pthread_cond_broadcast(&condvar_out);
  } else {
    while (threadsOut < totalThreads) {
      pthread_cond_wait(&condvar_out, &mutex);
    }
  }
  pthread_mutex_unlock(&mutex);
}

static inline void synchronizeForward(Args* a,
                        const int num) {
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
  static int threadsIn = 0;
  static int threadsOut = 0;
  int totalThreads = a->p;
  int div = a->n / a->m;
  int l = a->n - div * a->m;
  int tp = (l != 0 ? (div + 1): div);
  int tmp1;
  pthread_mutex_lock(&mutex);
  ++threadsIn;
  if (threadsIn >= totalThreads) {
    threadsOut = 0;
    pthread_cond_broadcast(&condvar_in);
  } else {
    while (threadsIn < totalThreads) {
      pthread_cond_wait(&condvar_in, &mutex);
    }
  }
  ++threadsOut;
  if (threadsOut >= totalThreads) {
    threadsIn = 0;
    if (a->i == num || ((num / a->p) * a->p + a->i) == num) {
      for (int i = num + 1; i < tp; i++) {
        tmp1 =  (i != div ? a->m : l);
        GetBlock(a->matrix, a->temp, i, num, tmp1, a->m, a->m, a->n);
        MiniFirstPass(a->matrix, a->temp, num, i, tmp1, a->n, a->m, a->x, a->result, a->newTemp);
      }
    }
    /*printf("MATRIX AND INVERSE AFTER FORWARD \n");
    Print(a->matrix, a->n, a->n, a->m, 20);
    printf("\n");
    Print(a->inverse, a->n, a->n, a->m, 20);*/
    pthread_cond_broadcast(&condvar_out);
  } else {
    while (threadsOut < totalThreads) {
      pthread_cond_wait(&condvar_out, &mutex);
    }
  }
  pthread_mutex_unlock(&mutex);
}

static inline void synchronizeReverse(Args* a,
                        const int num) {
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
  static int threadsIn = 0;
  static int threadsOut = 0;
  static bool f = true;
  int totalThreads = a->p;
  int div = a->n / a->m;
  int l = a->n - div * a->m;
  int tmp1;
  pthread_mutex_lock(&mutex);
  ++threadsIn;
  if (threadsIn >= totalThreads) {
    threadsOut = 0;
    pthread_cond_broadcast(&condvar_in);
  } else {
    while (threadsIn < totalThreads) {
      pthread_cond_wait(&condvar_in, &mutex);
    }
  }
  ++threadsOut;
  if (threadsOut >= totalThreads) {
    threadsIn = 0;
    if (f && (a->i == num || ((num / a->p) * a->p + a->i) == num)) {
      f = false;
      for (int i = num - 1; i > -1; i--) {
        tmp1 = (num != div ? a->m : l);
        GetBlock(a->matrix, a->temp, i, num, tmp1, a->m, a->m, a->n);
        //MiniSecondPass(a->matrix, a->temp, num, i, tmp1, a->n, a->m, a->x, a->result, a->newTemp);
      }
    }
    /*printf("MATRIX AND INVERSE AFTER REVERSE \n");
    Print(a->matrix, a->n, a->n, a->m, 20);
    printf("\n");
    Print(a->inverse, a->n, a->n, a->m, 20);*/
    pthread_cond_broadcast(&condvar_out);
  } else {
    while (threadsOut < totalThreads) {
      pthread_cond_wait(&condvar_out, &mutex);
    }
  }
  pthread_mutex_unlock(&mutex);
}

static inline void synchronizeMinBlock(Args* a, 
                         double* matrix, const int t, const int n, const int m) {
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
  static int threadsIn = 0;
  static int threadsOut = 0;
  static int minI = -1;
  static int minJ = -1;
  static double min = -1;
  int totalThreads = a->p;
  double temp;
  pthread_mutex_lock(&mutex);
  if (fabs(min + 1) < 1e-16) {
    minI = -1;
    minJ= -1;
  }
  if ((a->minNorm > 0) && ((min > a->minNorm) || ((fabs(min + 1) < 1e-16)))) {
    min = a->minNorm;
    minI = a->minI;
    minJ = a->minJ;
  }
  ++threadsIn;
  if (threadsIn >= totalThreads) {
    threadsOut = 0;
    pthread_cond_broadcast(&condvar_in);
  } else {
    while (threadsIn < totalThreads) {
      pthread_cond_wait(&condvar_in, &mutex);
    }
  }
  ++threadsOut;
  if (threadsOut >= totalThreads) {
    threadsIn = 0;
    if ((minJ != t) && (minJ != -1)) {
      ChangeCol(matrix, minJ, t, m, n); 
      temp = a->numCols[t];
      a->numCols[t] = a->numCols[minJ];
      a->numCols[minJ] = temp;
    }
    pthread_cond_broadcast(&condvar_out);
  } else {
    while (threadsOut < totalThreads) {
      pthread_cond_wait(&condvar_out, &mutex);
    }
  }
  if (minI != -1 && minJ != -1) {
      a->minI = minI;
      a->minJ = minJ;
  } else {
    a->error = -2;
  }
  min = -1;
  pthread_mutex_unlock(&mutex);
}

static inline void synchronizeResidual(Args* a) {
  static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
  static int threadsIn = 0;
  static int threadsOut = 0;
  int totalThreads = a->p;
  pthread_mutex_lock(&mutex);
  ++threadsIn;
  if (threadsIn >= totalThreads) {
    threadsOut = 0;
    pthread_cond_broadcast(&condvar_in);
  } else {
    while (threadsIn < totalThreads) {
      pthread_cond_wait(&condvar_in, &mutex);
    }
  }
  ++threadsOut;
  if (threadsOut >= totalThreads) {
    threadsIn = 0;
    pthread_cond_broadcast(&condvar_out);
  } else {
    /*printf("MATRIX AND INVERSE AFTER MULT \n");
    Print(a->matrix, a->n, a->n, a->m, 20);
    printf("\n");
    Print(a->inverse, a->n, a->n, a->m, 20);*/
    while (threadsOut < totalThreads) {
      pthread_cond_wait(&condvar_out, &mutex);
    }
  }
  a->resultAllNorm = ResultNorma(a->inverse, a->n, a->m, 1); 
  pthread_mutex_unlock(&mutex);
}

static inline void minBlock(Args& args,
             const double* matrix, 
	     const int t, const int threadNum, const int p,
             const int m, const int n,
	     const double tempNorm,
             double* x, double* temp, int* ind) {
  int minI = -1;
  int minJ = -1;
  double min = -1;
  double newMin;
  int div = n / m;
  int jStart = ((threadNum + p * (t / p) >= t) ? (t / p) * p + threadNum : (t / p + 1) * p + threadNum);
  minI = -1;
  minJ = -1;
  if (jStart != div) {
    for (int i = t; i < div; i++) {
      for (int j = jStart; j < div; j += p) {
        GetBlock(matrix, temp, i, j, m, m, m, n);
        if (InverseBlock(temp, x, m, tempNorm, ind) == 0) {
          newMin = Norma(temp, m, m);
          if ((newMin < min) || ((fabs(min + 1) < 1e-16))) {
            minI = i;
            minJ = j;
            min = newMin;
          }
	}
      }
    }
  }
  else { 
    if ((jStart == div) && (t == div)) {
      GetBlock(matrix, temp, div, div, n - div * m, n - div * m, m, n);
      if (InverseBlock(temp, x, n - div * m, tempNorm, ind) == 0)  {
        min = Norma(temp, n - div * m, n - div * m);
        minI = t;
        minJ = t;
      } else {
        minI = -1;
	minJ = -1;
	min = -1.0;
      }
    }
  }
  args.minNorm = min;
  args.minI = minI;
  args.minJ = minJ;
  return;
}

void* Inverse(void* ptr) {
  Args *args = (Args*)ptr;
  double *matrix = args->matrix;
  double *inverse = args->inverse;
  int threadNum = args->i;
  int n = args->n;
  int m = args->m;
  int s = args->s;
  int p = args->p;
  double *result = args->result;
  double *temp = args->temp;
  double *x = args->x;
  double *newTemp = args->newTemp;
  int *index = args->index;
  int div = n / m;
  int tmp1;
  double tempNorm = 0.;
  InputMatrixWithAlg(matrix, inverse, n, m, s, threadNum, p);
  synchronizeRead(args, matrix, n, m, s, args->filename, true);
  if (args->error != 0) {
    args->error = -1;
    return 0;
  }
  tempNorm = args->norm;
  int t1 = (div * m != n ? (div + 1): div);
  int tk, dk;
  args->threadTime = get_time();
  for (int t = 0; t < t1; t++) {
    minBlock(*args, matrix, t, threadNum, p,
                 m, n, 
                 tempNorm, x, temp, index);
    synchronizeMinBlock(args, matrix, t, n, m);
    if (args->error != 0) {
      //InputMatrixWithAlg(inverse, inverse, n, m, s, threadNum, p);
      //synchronizeRead(args, inverse, n, m, s, args->filename, true);
      args->error = -2;
      return 0;
    }
    if (args->minI != t) {
      ChangeRow(matrix, args->minI, t, m, n, threadNum, p); // поменять строки для смены блока
      ChangeRow(inverse, args->minI, t, m, n, threadNum, p);
    }
    synchronizeGetBlock(args, matrix, t);
    tmp1 = (t != div ? m : n - div * m);
    InverseBlock(temp, x, tmp1, tempNorm, index);
    RowMultiplication(matrix, temp, t, t, tmp1, n, m, threadNum, p, x, result);
    RowMultiplication(inverse, temp, t, 0, tmp1, n, m, threadNum, p, x, result);
    synchronize(args);
    ForwardGauss(matrix, inverse, t, m, n, threadNum, p, x, temp, result, newTemp);
    synchronizeForward(args, t);
  }
  for(int i = t1 - 1; i > 0; i--) {
    ReverseGauss(matrix, inverse, i, m, n, threadNum, p, x, temp, result, newTemp);
    synchronizeReverse(args, i);
  }
  for (int i = 0; i < t1; i++) {
    tk = args->numCols[i];
    dk = (tk != div ? m : n - div * m);
    for (int j = threadNum; j < t1; j += p) {
      tmp1 = (j != div ? m : n - div * m);
      temp = inverse + (i != div ? i * n * m + j * m * m : i * n * m + j * (n - div * m) * m);
      PutBlock(matrix, temp, tk, j, dk, tmp1, m, n);
    }
  }
  synchronize(args);
  args->threadTime = fabs(get_time() - args->threadTime);
  InputMatrixWithAlg(inverse, inverse, n, m, s, threadNum, p);
  synchronizeRead(args, inverse, n, m, s, args->filename, false);
  if (args->n <= 50 || args->p > 1) {
    mult(inverse, matrix, n, m, threadNum, p, args->bigTemp);
    synchronizeResidual(args);
  } else {
    args->resultAllNorm = -1.0;
  }
  return 0;
}
