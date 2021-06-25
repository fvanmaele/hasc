#ifndef RESIDUAL_HH
#define RESIDUAL_HH
#include <cmath>

double residual(int n, const double *__restrict u,
                       int i1begin, int i1end,
                       int i0begin, int i0end)
{
  double sqsum = 0;

  for (int i1 = i1begin; i1 < i1end; i1++)
    for (int i0 = i0begin; i0 < i0end; i0++)
    {
      double Azk = 4*u[i1*n + i0];
      Azk -= u[i1*n + i0 - n] + u[i1*n + i0 - 1] + u[i1*n + i0 + 1] +
             u[i1*n + i0 + n];
      Azk *= (double)1/((n-1)*(n-1));
      sqsum += Azk*Azk;
    }
  return std::sqrt(sqsum);
}

double residual(int n, const double *__restrict u)
{
  return residual(n, u, 1, n-1, 1, n-1);
}

#endif // RESIDUAL_HH
