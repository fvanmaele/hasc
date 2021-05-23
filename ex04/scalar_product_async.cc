#include <benchmark/benchmark.h>
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
static void BM_scalar_product_async_2(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_async<2>(x, y));
    }
}
static void BM_scalar_product_async_4(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_async<4>(x, y));
    }
}
static void BM_scalar_product_async_8(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_async<8>(x, y));
    }
}
static void BM_scalar_product_async_unseq_2(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_async_unseq<2>(x, y));
    }
}
static void BM_scalar_product_async_unseq_4(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_async_unseq<4>(x, y));
    }
}
static void BM_scalar_product_async_unseq_8(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_async_unseq<8>(x, y));
    }
}
static void BM_scalar_product_packaged_task_2(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_packaged_task<2>(x, y));
    }
}
static void BM_scalar_product_packaged_task_4(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_packaged_task<4>(x, y));
    }
}
static void BM_scalar_product_packaged_task_8(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_packaged_task<8>(x, y));
    }
}
static void BM_scalar_product_openmp(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_openmp(x, y));
    }
}
static void BM_scalar_product_tbb(benchmark::State& state) {
    // Perform setup here
    ptrdiff_t N = state.range(0);
    std::vector<NumberType> x(N);
    std::vector<NumberType> y(N);
    sp_init(x, y, -10, 10);

    for (auto _ : state) {
        benchmark::DoNotOptimize(sp_tbb(x, y));
    }
}

// Register the function as a benchmark
// XXX: this always gets the wrong results for functions using std::async and std::packaged_task.
BENCHMARK(BM_scalar_product_serial)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_policy_unseq)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_policy_par)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_policy_par_unseq)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_async_2)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_async_4)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_async_8)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_async_unseq_2)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_async_unseq_4)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_async_unseq_8)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_packaged_task_2)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_packaged_task_4)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_packaged_task_8)->Range(8<<3, 8<<22)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_openmp)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();
BENCHMARK(BM_scalar_product_tbb)->Range(8<<3, 8<<25)->Unit(benchmark::kMillisecond)->UseRealTime();

BENCHMARK_MAIN();
