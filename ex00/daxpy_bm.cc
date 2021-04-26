#include <benchmark/benchmark.h>
#include <vector>
#include <thread>
#include <numeric> // for iota
#include <algorithm> // for min
#include <utility> // for ref
#include <iostream>

// We choose the amount of threads depending on problem size to reduce
// overhead on small vectors.
const int min_size = 32;
const int hardware_threads = std::thread::hardware_concurrency();
const int min_threads = hardware_threads != 0 ? hardware_threads : 2;

// As the elements of the result vector are written directly to disjoint blocks,
// we do not need locking or condition variables as with a scalar product
// and its global sum.
//
// Note: y is passed by non-const reference, requiring std::ref when used
// with std::thread. An alternative is to take a pointer instead of a reference
// (e.g. throughs std::vector::data)
void daxpy(const std::vector<double> &x, std::vector<double> &y, double a,
           int rank, int nproc)
{
    assert(rank < nproc);
    assert(x.size() == y.size());
    const size_t N = x.size();
    
    // Disjoint partition for input and result vectors
    for (size_t i = (N*rank)/nproc; i < (N*(rank+1))/nproc; ++i) {
        y[i] += a*x[i];
    }
}

// Sequential version for comparison purposes
void daxpy_seq(const std::vector<double> &x, std::vector<double> &y, double a){
    for (std::size_t i=0; i<x.size(); ++i){
        y[i] += a*x[i];
    }
}

static void BM_daxpy(benchmark::State &state)
{
    // Setup
    int N = state.range(0);
    int nproc = std::min(min_threads, (N+min_size-1) / min_size);

    // Fill values of x (1, 2, 3, ..., N)
    std::vector<double> x(N);
    std::iota(x.begin(), x.end(), 1);  
    
    // Fill values of y (N, N-1, ..., 1)
    // Note: y will have different values over multiple runs
    std::vector<double> y(N);
    std::iota(y.rbegin(), y.rend(), 1);
    double a = 1.5;

    // Timed code
    for (auto _ : state) {
        std::vector<std::thread> threads;
        threads.reserve(nproc);

        for (int rank = 0; rank < nproc; ++rank) {
            threads.push_back(std::thread{daxpy, std::cref(x), std::ref(y), a, rank, nproc});
        }
        // daxpy(x, y, a, nproc-1, nproc);
        for (int rank = 0; rank < nproc; ++rank) {
            threads[rank].join();
        }
    }
}

static void BM_daxpy_seq(benchmark::State &state)
{
    // Setup
    int N = state.range(0);

    // Fill values of x (1, 2, 3, ..., N)
    std::vector<double> x(N);
    std::iota(x.begin(), x.end(), 1);  
    
    // Fill values of y (N, N-1, ..., 1)
    // Note: y will have different values over multiple runs
    std::vector<double> y(N);
    std::iota(y.rbegin(), y.rend(), 1);
    double a = 1.5;

    // Timed code
    for (auto _ : state) {
        daxpy_seq(x, y, a);
    }
}

// Register the function as a benchmark
BENCHMARK(BM_daxpy)->RangeMultiplier(4)->Range(256, 33554432)->UseRealTime();
BENCHMARK(BM_daxpy_seq)->RangeMultiplier(4)->Range(256, 33554432)->UseRealTime();

// Run the benchmark
BENCHMARK_MAIN();
