#ifndef JACOBI_VCL_HH
#define JACOBI_VCL_HH
#include <algorithm>
#include <utility>
#include <cmath>
#include <vectorclass.h>

double jacobi_vectorized_residual(int n, const double *__restrict u)
{
//  double sqsum = 0;

//  for (int i1 = 1; i1 < n - 1; i1++) {
//    for (int i0 = 1; i0 < n - 1; i0++) {
//      double Azk = u[i1 * n + i0];
//      Azk -= 0.25 * (u[i1 * n + i0 - n] + u[i1 * n + i0 - 1] +
//                     u[i1 * n + i0 + 1] + u[i1 * n + i0 + n]);
//      sqsum += Azk*Azk;
//    }
//  }
//  return std::sqrt(sqsum);

  int blocks4 = ((n - 2) / 4) * 4;

  Vec4d QUARTER = Vec4d(0.25);
  Vec4d A, B, C, D, E;
  Vec4d sqsum4 = Vec4d(0.0);
  double sqsum = 0.0;

  for (int i1 = 1; i1 < n - 1; i1++)
  {
    int istart = i1 * n + 1;
    int iend = istart + blocks4;

    for (int index = istart; index < iend; index += 4)
    {
      A.load(&u[index]);
      B.load(&u[index - n]);
      C.load(&u[index - 1]);
      D.load(&u[index + 1]);
      E.load(&u[index + n]);

      Vec4d Azk4 = A;
      Azk4 -= mul_add(QUARTER, B, A);
      Azk4 -= mul_add(QUARTER, C, A);
      Azk4 -= mul_add(QUARTER, D, A);
      Azk4 -= mul_add(QUARTER, E, A);
      sqsum4 += Azk4*Azk4;
    }

    sqsum = horizontal_add(sqsum4);
    for (int index = iend; index < n - 1; index++)
    {
      double Azk = u[index];
      Azk -= 0.25 * (u[index - n] + u[index - 1] + u[index + 1] + u[index + n]);
      sqsum += Azk * Azk;
    }
  }
  return std::sqrt(sqsum);
}

std::pair<int, double>
jacobi_vectorized_kernel(int n, int iterations, double *__restrict uold,
                         double *__restrict unew, double tol, int tol_check)
{
  int blocks4 = ((n - 2) / 4) * 4;

  Vec4d QUARTER = Vec4d(0.25);
  Vec4d A, B, C, D, E;

  double residual = 0;
  int i;
  // do iterations
  for (i = 0; i < iterations; i++)
  {
    for (int i1 = 1; i1 < n - 1; i1++)
    {
      int istart = i1 * n + 1;
      int iend = istart + blocks4;

      for (int index = istart; index < iend; index += 4)
      {
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
      {
        unew[index] = 0.25 * (uold[index - n] + uold[index - 1] +
                              uold[index + 1] + uold[index + n]);
      }
    }
    using std::swap;
    swap(uold, unew);

    // perform convergence check (b == 0)
    if (i > 0 && i % tol_check == 0) {
      residual = jacobi_vectorized_residual(n, uold);

      // check tolerance
      if (residual <= tol)
        break;
    }
  }
  return std::make_pair(i, residual);
}

template <int K>
std::pair<int, double>
jacobi_vectorized_wave_kernel(int n, int iterations, double *__restrict uold,
                              double *__restrict unew, double tol, int tol_check)
{
  double *u[2];
  u[0] = uold;
  u[1] = unew;

  if (tol_check != -1 && tol_check % K != 0)
    throw std::invalid_argument("tol_check must be a multiple of K");

  int blocks4 = ((n - 2) / 4) * 4;

  Vec4d QUARTER = Vec4d(0.25);
  Vec4d A, B, C, D, E;

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
        int istart = i1 * n + 1;
        int iend = istart + blocks4;

        for (int index = istart; index < iend; index += 4)
        {
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

    if (kk > 0 && kk % tol_check == 0)
    {
      int m = n + K - 2;
      int k = std::min(K, m - 1);
      residual = jacobi_vectorized_residual(n, u[k%2]);

      // check tolerance
      if (residual <= tol)
        break;
    }
  }
  return std::make_pair(std::min(kk, iterations), residual);
}

#endif // JACOBI_VCL_HH
