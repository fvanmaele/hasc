#include <benchmark/benchmark.h>
#include <vector>
#include "jacobi.hh"

double g(int n, int i0, int i1) {
  return (i0 > 0 && i0 < n - 1 && i1 > 0 && i1 < n - 1)
             ? 0.0
             : ((double)(i0 + i1)) / n;
}

static void BM_jacobi_vanilla(benchmark::State& state) {
  // Perform setup here
  size_t n = state.range(0);
  size_t k = state.range(1); // iterations

  auto context = std::make_shared<GlobalContext>(n);
  auto u0 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
  auto u1 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
  context->u0 = u0.get();
  context->u1 = u1.get();
  context->iterations = k;

  for (auto _ : state) {
    // This code gets timed
    jacobi_vanilla(context);
  }
}

template <int B>
static void BM_jacobi_blocked(benchmark::State& state) {
  // Perform setup here
  size_t n = state.range(0);
  size_t k = state.range(1); // iterations

  auto context = std::make_shared<GlobalContext>(n);
  auto u0 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
  auto u1 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
  context->u0 = u0.get();
  context->u1 = u1.get();
  context->iterations = k;

  for (auto _ : state) {
    // This code gets timed
    jacobi_blocked<B>(context);
  }
}

template <int K>
static void BM_jacobi_wave(benchmark::State& state) {
  // Perform setup here
  size_t n = state.range(0);
  size_t k = state.range(1); // iterations

  auto context = std::make_shared<GlobalContext>(n);
  auto u0 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
  auto u1 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
  context->u0 = u0.get();
  context->u1 = u1.get();
  context->iterations = k;

  for (auto _ : state) {
    // This code gets timed
    jacobi_wave<K>(context);
  }
}

static void BM_jacobi_vectorized(benchmark::State& state) {
  // Perform setup here
  size_t n = state.range(0);
  size_t k = state.range(1); // iterations

  auto context = std::make_shared<GlobalContext>(n);
  auto u0 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
  auto u1 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
  context->u0 = u0.get();
  context->u1 = u1.get();
  context->iterations = k;

  for (auto _ : state) {
    // This code gets timed
    jacobi_vectorized(context);
  }
}

template <int K>
static void BM_jacobi_vectorized_wave(benchmark::State& state) {
  // Perform setup here
  size_t n = state.range(0);
  size_t k = state.range(1); // iterations

  auto context = std::make_shared<GlobalContext>(n);
  auto u0 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
  auto u1 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
  context->u0 = u0.get();
  context->u1 = u1.get();
  context->iterations = k;

  for (auto _ : state) {
    // This code gets timed
    jacobi_vectorized_wave<K>(context);
  }
}

//static void BM_jacobi_vanilla_parallel(benchmark::State& state) {
//  // Perform setup here
//  size_t n = state.range(0);
//  size_t k = state.range(1); // iterations

//  auto context = std::make_shared<GlobalContext>(n);
//  auto u0 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
//  auto u1 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
//  context->u0 = u0.get();
//  context->u1 = u1.get();
//  context->iterations = k;

//  for (auto _ : state) {
//    // This code gets timed
//    jacobi_vanilla_parallel(context);
//  }
//}

//static void BM_jacobi_vectorized_parallel(benchmark::State& state) {
//  // Perform setup here
//  size_t n = state.range(0);
//  size_t k = state.range(1); // iterations

//  auto context = std::make_shared<GlobalContext>(n);
//  auto u0 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
//  auto u1 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
//  context->u0 = u0.get();
//  context->u1 = u1.get();
//  context->iterations = k;

//  for (auto _ : state) {
//    // This code gets timed
//    jacobi_vectorized_parallel(context);
//  }
//}

//template <int K>
//static void BM_jacobi_wave_parallel(benchmark::State& state) {
//  // Perform setup here
//  size_t n = state.range(0);
//  size_t k = state.range(1); // iterations

//  auto context = std::make_shared<GlobalContext>(n);
//  auto u0 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
//  auto u1 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
//  context->u0 = u0.get();
//  context->u1 = u1.get();
//  context->iterations = k;

//  for (auto _ : state) {
//    // This code gets timed
//    jacobi_wave_parallel<K>(context);
//  }
//}

//template <int K>
//static void BM_jacobi_vectorized_wave_parallel(benchmark::State& state) {
//  // Perform setup here
//  size_t n = state.range(0);
//  size_t k = state.range(1); // iterations

//  auto context = std::make_shared<GlobalContext>(n);
//  auto u0 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
//  auto u1 = std::unique_ptr<double>(new (std::align_val_t(64)) double[n * n]);
//  context->u0 = u0.get();
//  context->u1 = u1.get();
//  context->iterations = k;

//  for (auto _ : state) {
//    // This code gets timed
//    jacobi_vectorized_wave_parallel<K>(context);
//  }
//}

// Register the function as a benchmark
// TODO: include billion updates per second
BENCHMARK(BM_jacobi_vanilla)
    ->Ranges({{128, 2048}, {128, 1024}})
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_jacobi_blocked, 128)
    ->Ranges({{128, 2048}, {128, 1024}})
    ->Unit(benchmark::kMillisecond); // B = 128
BENCHMARK_TEMPLATE(BM_jacobi_wave, 40)
    ->Ranges({{128, 2048}, {128, 1024}})
    ->Unit(benchmark::kMillisecond); // K = 40
BENCHMARK(BM_jacobi_vectorized)
    ->Ranges({{128, 2048}, {128, 1024}})
    ->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(BM_jacobi_vectorized_wave, 40)
    ->Ranges({{128, 2048}, {128, 1024}})
    ->Unit(benchmark::kMillisecond); // K = 40

//BENCHMARK(BM_jacobi_vanilla_parallel)->Ranges({{16, 2048}, {128, 1024}});
//BENCHMARK(BM_jacobi_vectorized_parallel)->Ranges({{16, 2048}, {128, 1024}});
//BENCHMARK(BM_jacobi_wave_parallel)->Ranges({{16, 2048}, {128, 1024}});
//BENCHMARK(BM_jacobi_vectorized_wave_parallel)->Ranges({{16, 2048}, {128, 1024}});

BENCHMARK_MAIN();
