#include <benchmark/benchmark.h>
#include <tbb/task_scheduler_init.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "scalar_product_async.hh"

static void BM_scalar_product_serial(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_seq(x, y));
    }
}
static void BM_scalar_product_policy_unseq(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_unseq(x, y));
    }
}
static void BM_scalar_product_policy_par(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_par(x, y));
    }
}
static void BM_scalar_product_policy_par_unseq(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_par_unseq(x, y));
    }
}

template <int P>
static void BM_scalar_product_async(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_async<P>(x, y));
    }
}

template <int P>
static void BM_scalar_product_async_unseq(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_async_unseq<P>(x, y));
    }
}

template <int P>
static void BM_scalar_product_packaged_task(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_packaged_task<P>(x, y));
    }
}

template <int P = 0>
static void BM_scalar_product_openmp(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10); // note: not NUMA-aware
    if (P >= 1)
        omp_set_num_threads(P);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_openmp(x, y));
    }
}

template <int P>
static void BM_scalar_product_tbb(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);
    tbb::task_scheduler_init init(P); // avoid measuring overhead from initializing TBB

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_tbb(x, y));
    }
}

// Register the function as a benchmark
// XXX: this always gets the wrong results for functions using std::async and std::packaged_task.
BENCHMARK(BM_scalar_product_serial)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_scalar_product_policy_unseq)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_policy_par)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_policy_par_unseq)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_async, 2)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_async, 4)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_async, 8)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_async_unseq, 2)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_async_unseq, 4)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_async_unseq, 8)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_packaged_task, 2)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_packaged_task, 4)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_packaged_task, 8)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_openmp)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_openmp, 2)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_openmp, 4)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_openmp, 8)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_tbb, 2)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_tbb, 4)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK_TEMPLATE(BM_scalar_product_tbb, 8)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();

BENCHMARK_MAIN();
