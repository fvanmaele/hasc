#ifndef JACOBI_SEQ_HH
#define JACOBI_SEQ_HH
#include <utility>
#include <algorithm>

int jacobi_vanilla_kernel(int n, int iterations, double *__restrict__ uold,
                          double *__restrict__ unew) {
  // do iterations
  for (int i = 0; i < iterations; i++) {
    for (int i1 = 1; i1 < n - 1; i1++)
      for (int i0 = 1; i0 < n - 1; i0++)
        unew[i1 * n + i0] =
            0.25 * (uold[i1 * n + i0 - n] + uold[i1 * n + i0 - 1] +
                    uold[i1 * n + i0 + 1] + uold[i1 * n + i0 + n]);

    using std::swap;
    swap(uold, unew);
  }
  return iterations;
}

template <int B>
void jacobi_blocked_kernel(int n, int iterations, double *__restrict__ uold,
                           double *__restrict__ unew) {
  int blocksB = ((n - 2) / B) * B;
  int remainderB = (n - 2) % B;
  // std::cout << "blocksB=" << blocksB << " remainderB=" << remainderB <<
  // std::endl;

  // do iterations
  for (int i = 0; i < iterations; i++) {
    for (int I1 = 1; I1 < 1 + blocksB; I1 += B)
      for (int I0 = 1; I0 < 1 + blocksB; I0 += B)
        for (int i1 = I1; i1 < I1 + B; i1++)
          for (int i0 = I0; i0 < I0 + B; i0++)
            unew[i1 * n + i0] =
                0.25 * (uold[i1 * n + i0 - n] + uold[i1 * n + i0 - 1] +
                        uold[i1 * n + i0 + 1] + uold[i1 * n + i0 + n]);
    for (int I0 = 1; I0 < 1 + blocksB; I0 += B)
      for (int i1 = 1 + blocksB; i1 < n - 1; i1++)
        for (int i0 = I0; i0 < I0 + B; i0++)
          unew[i1 * n + i0] =
              0.25 * (uold[i1 * n + i0 - n] + uold[i1 * n + i0 - 1] +
                      uold[i1 * n + i0 + 1] + uold[i1 * n + i0 + n]);
    for (int i1 = 1 + blocksB; i1 < n - 1; i1++)
      for (int i0 = 1 + blocksB; i0 < n - 1; i0++)
        unew[i1 * n + i0] =
            0.25 * (uold[i1 * n + i0 - n] + uold[i1 * n + i0 - 1] +
                    uold[i1 * n + i0 + 1] + uold[i1 * n + i0 + n]);

    using std::swap;
    swap(uold, unew);
  }
}

template <int K>
void jacobi_wave_kernel(int n, int iterations, double *__restrict__ uold,
                        double *__restrict__ unew) {
  double *u[2];
  u[0] = uold;
  u[1] = unew;

  // do iterations
  for (int kk = 0; kk < iterations; kk += K)
    for (int m = 2; m <= n + K - 2; m++)
      for (int k = std::max(1, m - n + 2); k <= std::min(K, m - 1); k++) {
        int i1 = m - k;
        int dst = k % 2;
        int src = 1 - dst;

        for (int i0 = 1; i0 < n - 1; i0++) {
          u[dst][i1 * n + i0] =
              0.25 * (u[src][i1 * n + i0 - n] + u[src][i1 * n + i0 - 1] +
                      u[src][i1 * n + i0 + 1] + u[src][i1 * n + i0 + n]);
        }
      }
}

#endif // JACOBI_SEQ_HH
