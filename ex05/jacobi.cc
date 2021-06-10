#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "jacobi.hh"
#include "time_experiment.hh"

const int B = 128;
const int K = 40;

// main function runs the experiments and outputs results as csv
int main(int argc, char **argv)
{
  // read parameters
  int n = 1024;
  int iterations = 1000;
  if (argc != 3) {
    std::cout << "usage: ./jacobi_vanilla <size> <iterations>" << std::endl;
    exit(1);
  }
  sscanf(argv[1], "%d", &n);
  sscanf(argv[2], "%d", &iterations);
  // std::cout << "jacobi_vanilla: n=" << n
  //        << " iterations=" << iterations
  //        << " memory (mbytes)=" << (n*n)*8.0*2.0/1024.0/1024.0
  //        << std::endl;

  // check sizes
  if (K % 2 == 1) {
    std::cout << "K must be even" << std::endl;
    exit(1);
  }
  if (iterations % K != 0) {
    std::cout << "iterations must be a multiple of K" << std::endl;
    exit(1);
  }

  // get global context shared by aall threads
  auto context = std::make_shared<GlobalContext>(n);
  context->iterations = iterations;

  // allocate aligned arrays
  context->u0 = new (std::align_val_t(64)) double[n * n];
  context->u1 = new (std::align_val_t(64)) double[n * n];

  // fill boundary values and initial values
  auto g = [&](int i0, int i1) {
    return (i0 > 0 && i0 < n - 1 && i1 > 0 && i1 < n - 1)
               ? 0.0
               : ((double)(i0 + i1)) / n;
  };

  std::cout << "N,";
  std::cout << "vanilla,";
  std::cout << "blocked,";
  std::cout << "wave,";
  std::cout << "vectorized,";
  std::cout << "vectorized_wave,";
  std::cout << "parallel";
  std::cout << std::endl;
  std::cout << n * n;

  std::vector<std::pair<int, double> > convergence(6);

  // warmup
  for (int i1 = 0; i1 < n; i1++)
    for (int i0 = 0; i0 < n; i0++)
      context->u0[i1 * n + i0] = context->u1[i1 * n + i0] = g(i0, i1);
  auto start = get_time_stamp();
  auto [it, res] = jacobi_blocked<B>(context);
  auto stop = get_time_stamp();
  double elapsed = get_duration_seconds(start, stop);
  double updates = context->iterations;
  updates *= (n - 2) * (n - 2);
  // std::cout << "," << updates/elapsed/1e9;

  // vanilla
  for (int i1 = 0; i1 < n; i1++)
    for (int i0 = 0; i0 < n; i0++)
      context->u0[i1 * n + i0] = context->u1[i1 * n + i0] = g(i0, i1);
  start = get_time_stamp();
  std::tie(it, res) = jacobi_vanilla(context);
  stop = get_time_stamp();
  convergence.at(0) = {it, res};
  elapsed = get_duration_seconds(start, stop);
  updates = context->iterations;
  updates *= (n - 2) * (n - 2);
  std::cout << "," << updates / elapsed / 1e9;

  // blocked
  for (int i1 = 0; i1 < n; i1++)
    for (int i0 = 0; i0 < n; i0++)
      context->u0[i1 * n + i0] = context->u1[i1 * n + i0] = g(i0, i1);
  start = get_time_stamp();
  std::tie(it, res) = jacobi_blocked<B>(context);
  stop = get_time_stamp();
  convergence.at(1) = {it, res};
  elapsed = get_duration_seconds(start, stop);
  std::cout << "," << updates / elapsed / 1e9;

  // wave
  for (int i1 = 0; i1 < n; i1++)
    for (int i0 = 0; i0 < n; i0++)
      context->u0[i1 * n + i0] = context->u1[i1 * n + i0] = g(i0, i1);
  start = get_time_stamp();
  std::tie(it, res) = jacobi_wave<K>(context);
  stop = get_time_stamp();
  convergence.at(2) = {it, res};
  elapsed = get_duration_seconds(start, stop);
  std::cout << "," << updates / elapsed / 1e9;

  // vectorized
  for (int i1 = 0; i1 < n; i1++)
    for (int i0 = 0; i0 < n; i0++)
      context->u0[i1 * n + i0] = context->u1[i1 * n + i0] = g(i0, i1);
  start = get_time_stamp();
  std::tie(it, res) = jacobi_vectorized(context);
  stop = get_time_stamp();
  convergence.at(3) = {it, res};
  elapsed = get_duration_seconds(start, stop);
  std::cout << "," << updates / elapsed / 1e9;

  // vectorized wave
  for (int i1 = 0; i1 < n; i1++)
    for (int i0 = 0; i0 < n; i0++)
      context->u0[i1 * n + i0] = context->u1[i1 * n + i0] = g(i0, i1);
  start = get_time_stamp();
  std::tie(it, res) = jacobi_vectorized_wave<K>(context);
  stop = get_time_stamp();
  convergence.at(4) = {it, res};
  elapsed = get_duration_seconds(start, stop);
  std::cout << "," << updates / elapsed / 1e9;

  // parallel
  for (int i1 = 0; i1 < n; i1++)
    for (int i0 = 0; i0 < n; i0++)
      context->u0[i1 * n + i0] = context->u1[i1 * n + i0] = g(i0, i1);
  start = get_time_stamp();
  std::tie(it, res) = jacobi_parallel(context);
  stop = get_time_stamp();
  convergence.at(5) = {it, res};
  elapsed = get_duration_seconds(start, stop);
  std::cout << "," << updates / elapsed / 1e9;

  std::cout << std::endl;

//  std::cout << "vanilla, iterations: " << convergence[0].first
//            << ", residual: " << convergence[0].second << std::endl;
//  std::cout << "blocked, iterations: " << convergence[1].first
//            << ", residual: " << convergence[1].second << std::endl;
//  std::cout << "wave, iterations: " << convergence[2].first
//            << ", residual: " << convergence[2].second << std::endl;
//  std::cout << "vectorized, iterations: " << convergence[3].first
//            << ", residual: " << convergence[3].second << std::endl;
//  std::cout << "vectorized_wave, iterations: " << convergence[4].first
//            << ", residual: " << convergence[4].second << std::endl;
//  std::cout << "parallel, iterations: " << convergence[5].first
//            << ", residual: " << convergence[5].second << std::endl;

  // deallocate arrays
  delete[] context->u1;
  delete[] context->u0;

  return 0;
}
