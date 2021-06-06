#ifndef JACOBI_SEQ_HH
#define JACOBI_SEQ_HH
#include <algorithm>
#include <utility>
#include <cmath>

double jacobi_residual(int n, const double *__restrict u)
{
  double sqsum = 0;

  for (int i1 = 1; i1 < n - 1; i1++)
    for (int i0 = 1; i0 < n - 1; i0++)
    {
      double Azk = u[i1 * n + i0];
      Azk -= 0.25 * (u[i1 * n + i0 - n] + u[i1 * n + i0 - 1] +
                     u[i1 * n + i0 + 1] + u[i1 * n + i0 + n]);
      sqsum += Azk*Azk;
    }
  return std::sqrt(sqsum);
}

std::pair<int, double>
jacobi_vanilla_kernel(int n, int iterations, double *__restrict uold,
                      double *__restrict unew, double tol, int tol_check)
{
  double residual = 0;
  int i;
  // do iterations
  for (i = 0; i < iterations; i++)
  {
    for (int i1 = 1; i1 < n - 1; i1++)
      for (int i0 = 1; i0 < n - 1; i0++)
        unew[i1 * n + i0] = 0.25 * (uold[i1 * n + i0 - n] + uold[i1 * n + i0 - 1] +
                                    uold[i1 * n + i0 + 1] + uold[i1 * n + i0 + n]);

    using std::swap;
    swap(uold, unew);

    // perform convergence check (b == 0)
    if (i > 0 && i % tol_check == 0) {
      residual = jacobi_residual(n, uold);

      // check tolerance
      if (residual <= tol)
        break;
    }
  }
  return std::make_pair(i, residual);
}

template <int B>
std::pair<int, double>
jacobi_blocked_kernel(int n, int iterations, double *__restrict uold,
                      double *__restrict unew, double tol, int tol_check)
{
  int blocksB = ((n - 2) / B) * B;
  int remainderB = (n - 2) % B; // unused?
  // std::cout << "blocksB=" << blocksB << " remainderB=" << remainderB <<
  // std::endl;

  double residual = 0;
  int i;
  // do iterations
  for (i = 0; i < iterations; i++)
  {
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

    // perform convergence check (b == 0)
    if (i > 0 && i % tol_check == 0) {
      residual = jacobi_residual(n, uold);

      // check tolerance
      if (residual <= tol)
        break;
    }
  }
  return std::make_pair(i, residual);
}

template <int K>
std::pair<int, double>
jacobi_wave_kernel(int n, int iterations, double *__restrict uold,
                   double *__restrict unew, double tol, int tol_check)
{
  double *u[2];
  u[0] = uold;
  u[1] = unew;

  if (tol_check != -1 && tol_check % K != 0)
    throw std::invalid_argument("tol_check must be a multiple of K");

  // do iterations
  double residual = 0;
  int kk;
  for (kk = 0; kk < iterations; kk += K)
  {
    for (int m = 2; m <= n + K - 2; m++)
      for (int k = std::max(1, m - n + 2); k <= std::min(K, m - 1); k++)
      {
        int i1 = m - k;
        int dst = k % 2;
        int src = 1 - dst;

        for (int i0 = 1; i0 < n - 1; i0++) {
          u[dst][i1 * n + i0] = 0.25 * (u[src][i1 * n + i0 - n] + u[src][i1 * n + i0 - 1] +
                                        u[src][i1 * n + i0 + 1] + u[src][i1 * n + i0 + n]);
        }
      }

    if (kk > 0 && kk % tol_check == 0)
    {
      int m = n + K - 2;
      int k = std::min(K, m - 1);
      residual = jacobi_residual(n, u[k % 2]); // xxx: which was last updated?

      // check tolerance
      if (residual <= tol)
        break;
    }
  }
  return std::make_pair(std::min(kk, iterations), residual);
}

#endif // JACOBI_SEQ_HH
