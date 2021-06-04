#ifndef GLOBAL_CONTEXT_HH
#define GLOBAL_CONTEXT_HH
#include "jacobi_seq.hh"
#include "jacobi_vcl.hh"
#include "jacobi_par.hh"
#include <memory>

#ifdef _WIN32
#define __restrict__ __restrict
#endif

struct GlobalContext {
  // input data
  int n;          // nxn lattice of points
  int iterations; // number of iterations to do
  double *u0;     // the initial guess
  double *u1;     // temporary vector

  // output data

  GlobalContext(int n_) : n(n_) {}
};

void jacobi_vanilla(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  jacobi_vanilla_kernel(context->n, context->iterations, uold, unew);

  // result should be in u1
  if (context->u1 != uold)
    std::swap(context->u0, context->u1);
}

void jacobi_vectorized(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  jacobi_vectorized_kernel(context->n, context->iterations, uold, unew);

  // result should be in u1
  if (context->u1 != uold)
    std::swap(context->u0, context->u1);
}

template <int B>
void jacobi_blocked(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  jacobi_blocked_kernel<B>(context->n, context->iterations, uold, unew);

  // result should be in u1
  if (context->u1 != uold)
    std::swap(context->u0, context->u1);
}

template <int K>
void jacobi_wave(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  jacobi_wave_kernel<K>(context->n, context->iterations, uold, unew);

  // result should be in u1
  if (context->u1 != uold)
    std::swap(context->u0, context->u1);
}

template <int K>
void jacobi_vectorized_wave(std::shared_ptr<GlobalContext> context) {
  double* uold = context->u0;
  double* unew = context->u1;
  jacobi_vectorized_wave_kernel<K>(context->n, context->iterations, uold, unew);

  // result should be in u1
  if (context->u1 != uold)
    std::swap(context->u0, context->u1);
}

#endif // GLOBAL_CONTEXT_HH
