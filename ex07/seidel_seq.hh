#ifndef SEIDEL_SEQ_HH
#define SEIDEL_SEQ_HH
#include "residual.hh"
#include <algorithm>
#include <utility>

std::pair<int, double> seidel_vanilla_kernel(int n, int iterations, double *__restrict u, double tol,
                                             int tol_check)
{
  double res = 0;
  int i;
  // do iterations
  for (i = 0; i < iterations; i++)
  {
    for (int i1 = 1; i1 < n - 1; i1++)
      for (int i0 = 1; i0 < n - 1; i0++)
        u[i1 * n + i0] = 0.25 * (u[i1 * n + i0 - n] + u[i1 * n + i0 - 1] +
                                 u[i1 * n + i0 + 1] + u[i1 * n + i0 + n]);

    // perform convergence check (b == 0)
    if (i > 0 && i % tol_check == 0)
    {
      res = residual(n, u);

      // check tolerance
      if (res <= tol)
        break;
    }
  }
  return std::make_pair(i, res);
}

template <int B>
std::pair<int, double> seidel_blocked_kernel(int n, int iterations, double *__restrict u, double tol,
                                             int tol_check)
{
  int blocksB = ((n - 2) / B) * B;
  //int remainderB = (n - 2) % B;
  // std::cout << "blocksB=" << blocksB << " remainderB=" << remainderB <<
  // std::endl;

  double res = 0;
  int i;
  // do iterations
  for (i = 0; i < iterations; i++)
  {
    // Interior (blocks)
    for (int I1 = 1; I1 < 1 + blocksB; I1 += B)
      for (int I0 = 1; I0 < 1 + blocksB; I0 += B)
        for (int i1 = I1; i1 < I1 + B; i1++)
          for (int i0 = I0; i0 < I0 + B; i0++)
            u[i1 * n + i0] = 0.25 * (u[i1 * n + i0 - n] + u[i1 * n + i0 - 1] +
                                     u[i1 * n + i0 + 1] + u[i1 * n + i0 + n]);

    // Boundary (remainder)
    for (int I1 = 1; I1 < 1 + blocksB; I1 += B)
      for (int i1 = I1; i1 < I1 + B; i1++)
        for (int i0 = 1 + blocksB; i0 < n - 1; i0++)
          u[i1 * n + i0] = 0.25 * (u[i1 * n + i0 - n] + u[i1 * n + i0 - 1] +
                                   u[i1 * n + i0 + 1] + u[i1 * n + i0 + n]);

    for (int I0 = 1; I0 < 1 + blocksB; I0 += B)
      for (int i1 = 1 + blocksB; i1 < n - 1; i1++)
        for (int i0 = I0; i0 < I0 + B; i0++)
          u[i1 * n + i0] = 0.25 * (u[i1 * n + i0 - n] + u[i1 * n + i0 - 1] +
                                   u[i1 * n + i0 + 1] + u[i1 * n + i0 + n]);

    for (int i1 = 1 + blocksB; i1 < n - 1; i1++)
      for (int i0 = 1 + blocksB; i0 < n - 1; i0++)
        u[i1 * n + i0] = 0.25 * (u[i1 * n + i0 - n] + u[i1 * n + i0 - 1] +
                                 u[i1 * n + i0 + 1] + u[i1 * n + i0 + n]);

    // perform convergence check (b == 0)
    if (i > 0 && i % tol_check == 0)
    {
      res = residual(n, u);

      // check tolerance
      if (res <= tol)
        break;
    }
  }
  return std::make_pair(i, res);
}

#endif // SEIDEL_SEQ_HH
