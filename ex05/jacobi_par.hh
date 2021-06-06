#ifndef JACOBI_PAR_HH
#define JACOBI_PAR_HH
#include <algorithm>
#include <utility>
#include <cmath>

double jacobi_parallel_residual(int n, const double *__restrict u)
{
  double sqsum = 0;

  // no dependencies apart from reduction variable (parallel reads)
#pragma omp parallel for collapse(2) reduction(+:sqsum)
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
jacobi_parallel_kernel(int n, int iterations, double *__restrict uold,
                       double *__restrict unew, double tol, int tol_check)
{
  double residual = 0;
  int i;
  // do iterations
  for (i = 0; i < iterations; i++)
  {
  // parallellize inner loops (parallel reads on uold)
#pragma omp parallel for collapse(2)
    for (int i1 = 1; i1 < n - 1; i1++)
      for (int i0 = 1; i0 < n - 1; i0++)
      {
        unew[i1 * n + i0] = 0.25 * (uold[i1 * n + i0 - n] + uold[i1 * n + i0 - 1] +
                                    uold[i1 * n + i0 + 1] + uold[i1 * n + i0 + n]);
      }
    using std::swap;
    swap(uold, unew);

    // perform convergence check (b == 0)
    if (i > 0 && i % tol_check == 0) {
      residual = jacobi_parallel_residual(n, uold);

      // check tolerance
      if (residual <= tol)
        break;
    }
  }
  return std::make_pair(i, residual);
}
#endif // JACOBI_PAR_HH
