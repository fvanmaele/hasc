#include <benchmark/benchmark.h>
#include <vector>
#include <thread>
#include <numeric> // for iota
#include <algorithm> // for min
#include <utility> // for ref
#include <random>
#include <algorithm>
#include <iterator>
#include <functional>

template <typename T>
void initialize_vector(std::vector<T>& x){
    // First create an instance of an engine.
    std::random_device rnd_device;
    // Specify the engine and distribution.
    std::mt19937 mersenne_engine {rnd_device()};  // Generates random integers
    std::uniform_real_distribution<T> dist(-10,10);

    auto gen = [&dist, &mersenne_engine](){
        return dist(mersenne_engine);
    };

    generate(x.begin(), x.end(), gen);
}


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
    // assert(rank < nproc);
    // assert(x.size() == y.size());
    const size_t N = x.size();
    
    // Disjoint partition for input and result vectors
    for (size_t i = (N*rank)/nproc; i < (N*(rank+1))/nproc; ++i) {
        y[i] += a*x[i];
    }
}

static void BM_daxpy(benchmark::State &state)
{
    // Setup
    int N = state.range(0);
    int nproc = std::thread::hardware_concurrency();
    double a = 1.5;
    std::vector<double> x(N);
    std::vector<double> y(N);
    initialize_vector(x);
    initialize_vector(y);

    // Timed code
    for (auto _ : state) {
        std::vector<std::thread> threads;
        for (int rank = 0; rank < nproc; ++rank) {
            threads.push_back(std::thread{daxpy, std::cref(x), std::ref(y), a, rank, nproc});
        }
        for (int rank = 0; rank < nproc; ++rank) {
            threads[rank].join();
        }
        // daxpy_seq(x,y,a);
    }
}

// Register the function as a benchmark
BENCHMARK(BM_daxpy) -> RangeMultiplier(4) -> Range(256, 33554432)->UseRealTime();


void daxpy_seq(const std::vector<double> &x, std::vector<double> &y, double a){
    for (std::size_t i=0; i<x.size(); ++i){
        y[i] += a*x[i];
    }
}

static void BM_daxpy_seq(benchmark::State &state)
{
    // Setup
    int N = state.range(0);
    int nproc = std::thread::hardware_concurrency();
    double a = 1.5;
    std::vector<double> x(N);
    std::vector<double> y(N);
    initialize_vector(x);
    initialize_vector(y);

    // Timed code
    for (auto _ : state) {
        daxpy_seq(x,y,a);
    }
}


BENCHMARK(BM_daxpy_seq) -> RangeMultiplier(4) -> Range(256, 33554432)->UseRealTime();

// Run the benchmark
BENCHMARK_MAIN();
