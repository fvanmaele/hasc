#include <benchmark/benchmark.h>
#include <vector>
#include <iostream>
#include <random>

static double sum = 0;

// 1. An empty loop that iterates from 0 to n-1 without doing anything.
static void BM_empty_loop(benchmark::State &state)
{
    // Setup
    int n = state.range(0);

    // Timed code
    for (auto _ : state) {
        for (int i = 0; i < n; ++i)
            ;
    }
}

// 2. A loop that accumulates something in a global variable.
static void BM_reduce_global(benchmark::State& state)
{
    // Setup
    std::mt19937_64 gen(42);
    std::normal_distribution<float> d{ 10, 2 };
    int n = state.range(0);
    std::vector<float> v(n);
    for (float& c: v)
        c = d(gen);

    // Timed code
    for (auto _ : state) {
        sum = 0; // reset reduction value
        for (auto &&c : v)
            sum += c;
    }
}

// 3. A loop that fills a vector of size n with the integers 0, ..., n-1;
static void BM_fill_vector(benchmark::State &state)
{
    // Setup
    int n = state.range(0);
    std::vector<size_t> v(n); // 0-initialized

    // Timed code
    for (auto _ : state) {
        for (size_t i = 0; i < v.size(); ++i)
            //benchmark::DoNotOptimize(v[i] = i);
            v[i] = i;
    }
}

// Register the function as a benchmark
BENCHMARK(BM_empty_loop)->Range(2, 2 << 20);
BENCHMARK(BM_reduce_global)->Range(2, 2 << 20);
BENCHMARK(BM_fill_vector)->Range(2, 2 << 20);

// Run the benchmark
BENCHMARK_MAIN();