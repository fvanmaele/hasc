#ifndef HASC_EX02_TRANSPOSE_HH
#define HASC_EX02_TRANSPOSE_HH
#include <iosfwd>

// B = A^T,
// A,B are nxn matrices stored row-major in a 1d array,
// assume A and B are NOT the same matrix
// B accessed consecutively, A accessed strided
inline void transpose1(int n, const double * __restrict A, double * __restrict B)
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      B[i * n + j] = A[j * n + i];
}

// B = A^T,
// A,B are nxn matrices stored row-major in a 1d array,
// assume A and B are NOT the same matrix
// A accessed consecutively, B accessed strided
inline void transpose2(int n, const double * __restrict A, double * __restrict B)
{
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      B[j * n + i] = A[i * n + j];
}

// B = A^T,
// A,B are nxn matrices stored row-major in a 1d array,
// assume A and B are NOT the same matrix
// B accessed consecutively, A accessed strided
template <int M, int N> // MxN blocks
inline void transpose3(int n, const double * __restrict A, double * __restrict B)
{
  if (n % M != 0)
  {
    std::cout << n << " is not a multiple of " << M << std::endl;
    exit(0);
  }
  if (n % N != 0)
  {
    std::cout << n << " is not a multiple of " << N << std::endl;
    exit(0);
  }
  for (int i = 0; i < n; i += M)
    for (int j = 0; j < n; j += N)
    {
      int Bstart = i * n + j;
      int Astart = j * n + i;
      for (int ii = 0; ii < M; ++ii)
        for (int jj = 0; jj < N; ++jj)
          B[ii * n + jj + Bstart] = A[jj * n + ii + Astart];
    }
}

// B = A^T,
// A,B are nxn matrices stored row-major in a 1d array,
// assume A and B are NOT the same matrix
// B accessed strided, A accessed consecutively
template <int M, int N> // MxN blocks
inline void transpose4(int n, const double * __restrict A, double * __restrict B)
{
  if (n % M != 0)
  {
    std::cout << n << " is not a multiple of " << M << std::endl;
    exit(0);
  }
  if (n % N != 0)
  {
    std::cout << n << " is not a multiple of " << N << std::endl;
    exit(0);
  }
  for (int i = 0; i < n; i += M)
    for (int j = 0; j < n; j += N)
    {
      int Astart = i * n + j;
      int Bstart = j * n + i;
      for (int ii = 0; ii < M; ++ii)
        for (int jj = 0; jj < N; ++jj)
          B[jj * n + ii + Bstart] = A[ii * n + jj + Astart];
    }
}

// B = A^T,
// A,B are nxn matrices stored row-major in a 1d array,
// assume A and B are NOT the same matrix
template <int M> // bxb blocks with M blocking in the strided direction
inline void transpose5(int n, int b, const double * __restrict A, double * __restrict B)
{
  if (b % M != 0)
  {
    std::cout << b << " is not a multiple of " << M << std::endl;
    exit(0);
  }
  if (n % b != 0)
  {
    std::cout << n << " is not a multiple of " << b << std::endl;
    exit(0);
  }
  // first two loops over the blocks
  for (int I = 0; I < n; I += b)
    for (int J = 0; J < n; J += b)
      for (int j = J; j < J + b; j += M)
        for (int i = I; i < I + b; ++i)
          for (int jj = 0; jj < M; ++jj)
          {
            int jjj = j + jj;
            B[i * n + jjj] = A[jjj * n + i]; // consecutive writes
          }
}

// B = A^T,
// A,B are nxn matrices stored row-major in a 1d array,
// assume A and B are NOT the same matrix
template <int M> // bxb blocks with M blocking in the strided direction
inline void transpose6(int n, int b, const double * __restrict A, double * __restrict B)
{
  if (b % M != 0)
  {
    std::cout << b << " is not a multiple of " << M << std::endl;
    exit(0);
  }
  if (n % b != 0)
  {
    std::cout << n << " is not a multiple of " << b << std::endl;
    exit(0);
  }
  // first two loops over the blocks
  for (int I = 0; I < n; I += b)
    for (int J = 0; J < n; J += b)
      for (int j = J; j < J + b; j += M)
        for (int i = I; i < I + b; ++i)
          for (int jj = 0; jj < M; ++jj)
          {
            int jjj = j + jj;
            B[jjj * n + i] = A[i * n + jjj]; // strided writes
          }
}

#endif // HASC_EX02_TRANSPOSE_HH