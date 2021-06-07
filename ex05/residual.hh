#ifndef RESIDUAL_HH
#define RESIDUAL_HH

double jacobi_residual(int n, const double *__restrict u)
{
  double sqsum = 0;

  for (int i1 = 1; i1 < n - 1; i1++)
    for (int i0 = 1; i0 < n - 1; i0++)
    {
      double Azk = 4*u[i1 * n + i0];
      Azk -= u[i1 * n + i0 - n] + u[i1 * n + i0 - 1] + u[i1 * n + i0 + 1] +
             u[i1 * n + i0 + n];
      Azk *= (double)1/(n-1);
      sqsum += Azk*Azk;
    }
  return std::sqrt(sqsum);
}

#endif // RESIDUAL_HH
