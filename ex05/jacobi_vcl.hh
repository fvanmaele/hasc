#ifndef JACOBI_VCL_HH
#define JACOBI_VCL_HH
#include <utility>
#include <algorithm>
#include <vectorclass.h>

void jacobi_vectorized_kernel(int n, int iterations, double *__restrict uold,
                              double *__restrict unew) {
  int blocks4 = ((n - 2) / 4) * 4;

  Vec4d QUARTER = Vec4d(0.25);
  Vec4d A, B, C, D, E;

  // do iterations
  for (int i = 0; i < iterations; i++) {
    for (int i1 = 1; i1 < n - 1; i1++) {
      int istart = i1 * n + 1;
      int iend = istart + blocks4;

      for (int index = istart; index < iend; index += 4) {
        A = Vec4d(0.0);
        B.load(&uold[index - n]);
        C.load(&uold[index - 1]);
        D.load(&uold[index + 1]);
        E.load(&uold[index + n]);
        A = mul_add(QUARTER, B, A);
        A = mul_add(QUARTER, C, A);
        A = mul_add(QUARTER, D, A);
        A = mul_add(QUARTER, E, A);
        A.store(&unew[index]);
      }
      for (int index = iend; index < n - 1; index++)
        unew[index] = 0.25 * (uold[index - n] + uold[index - 1] +
                              uold[index + 1] + uold[index + n]);
    }

    using std::swap;
    swap(uold, unew);
  }
}

template <int K>
void jacobi_vectorized_wave_kernel(int n, int iterations, double *__restrict uold,
                                   double *__restrict unew) {
  double *u[2];
  u[0] = uold;
  u[1] = unew;

  int blocks4 = ((n - 2) / 4) * 4;

  Vec4d QUARTER = Vec4d(0.25);
  Vec4d A, B, C, D, E;

  // do iterations
  for (int kk = 0; kk < iterations; kk += K)
    for (int m = 2; m <= n + K - 2; m++)
      for (int k = std::max(1, m - n + 2); k <= std::min(K, m - 1); k++) {
        int i1 = m - k;
        int dst = k % 2;
        int src = 1 - dst;
        int istart = i1 * n + 1;
        int iend = istart + blocks4;

        for (int index = istart; index < iend; index += 4) {
          A = Vec4d(0.0);
          B.load(&u[src][index - n]);
          C.load(&u[src][index - 1]);
          D.load(&u[src][index + 1]);
          E.load(&u[src][index + n]);
          A = mul_add(QUARTER, B, A);
          A = mul_add(QUARTER, C, A);
          A = mul_add(QUARTER, D, A);
          A = mul_add(QUARTER, E, A);
          A.store(&u[dst][index]);
        }
        for (int index = iend; index < n - 1; index++)
          u[dst][index] = 0.25 * (u[src][index - n] + u[src][index - 1] +
                                  u[src][index + 1] + u[src][index + n]);
      }
}

#endif // JACOBI_VCL_HH
