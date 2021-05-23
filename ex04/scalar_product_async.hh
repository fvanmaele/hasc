#ifndef SCALAR_PRODUCT_ASYNC_HH
#define SCALAR_PRODUCT_ASYNC_HH
#include <future>
#include <thread>
#include <numeric>
#include <span>
#include <execution>
#include <functional>
#include <cassert>
#include <vector>
#include <random>

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

// Use higher precision reduction type for improved numerical stability (large vectors)
using NumberType = float;
using ReduceType = double;


ReduceType sp_seq(std::span<NumberType> x, std::span<NumberType> y) // or std::inner_product
{
    ReduceType sum = 0;
    for (size_t i = 0, n = x.size(); i < n; ++i) {
        sum += x[i] * y[i];
    }
    return sum;
}

ReduceType sp_unseq(std::span<NumberType> x, std::span<NumberType> y)
{
    return std::transform_reduce(std::execution::unseq, x.begin(), x.end(), y.begin(), ReduceType(0),
                                 std::plus<ReduceType>(), std::multiplies<ReduceType>());
}

ReduceType sp_par(std::span<NumberType> x, std::span<NumberType> y)
{
    return std::transform_reduce(std::execution::par, x.begin(), x.end(), y.begin(), ReduceType(0),
                                 std::plus<ReduceType>(), std::multiplies<ReduceType>());
}

ReduceType sp_par_unseq(std::span<NumberType> x, std::span<NumberType> y)
{
    return std::transform_reduce(std::execution::par_unseq, x.begin(), x.end(), y.begin(), ReduceType(0),
                                 std::plus<ReduceType>(), std::multiplies<ReduceType>());
}

template <int P>
ReduceType sp_async(std::span<NumberType> x, std::span<NumberType> y)
{
    assert(x.size() % P == 0);
    assert(x.size() == y.size());

    std::vector<std::future<ReduceType>> futures(P);
    std::vector<ReduceType> psums(P);

    size_t N = x.size();
    for (int rank = 0; rank < P; ++rank) {
        // build chunk
        size_t ibegin = N*rank/P;
        size_t iend = (N+1)*rank/P;

        auto fut = std::async(sp_seq, // default launch policy
                              std::span{x.begin()+ibegin, x.begin()+iend},
                              std::span{y.begin()+ibegin, y.begin()+iend});
        futures[rank] = std::move(fut);
    }
    for (int rank = 0; rank < P; ++rank) {
        psums[rank] = futures[rank].get();
    }
    return std::accumulate(psums.begin(), psums.end(), ReduceType(0));
}

template <int P>
ReduceType sp_async_unseq(std::span<NumberType> x, std::span<NumberType> y)
{
    assert(x.size() % P == 0);
    assert(x.size() == y.size());

    std::vector<std::future<ReduceType>> futures(P);
    std::vector<ReduceType> psums(P);

    size_t N = std::ssize(x);
    for (int rank = 0; rank < P; ++rank) {
        // build chunk
        size_t ibegin = N*rank/P;
        size_t iend = (N+1)*rank/P;

        auto fut = std::async(sp_unseq, // default launch policy
                              std::span{x.begin()+ibegin, x.begin()+iend},
                              std::span{y.begin()+ibegin, y.begin()+iend});
        futures[rank] = std::move(fut);
    }
    for (int rank = 0; rank < P; ++rank) {
        psums[rank] = futures[rank].get();
    }
    return std::reduce(std::execution::unseq, psums.begin(), psums.end(), ReduceType(0));
}

template <int P>
ReduceType sp_packaged_task(std::span<NumberType> x, std::span<NumberType> y)
{
    assert(x.size() % P == 0);
    assert(x.size() == y.size());

    std::vector<std::thread> threads(P);
    std::vector<std::future<ReduceType>> futures(P);
    std::vector<ReduceType> psums(P);

    size_t N = std::ssize(x);
    for (int rank = 0; rank < P; ++rank) {
        // build chunk
        size_t ibegin = N*rank/P;
        size_t iend = (N+1)*rank/P;

        using TaskType = decltype(sp_seq);
        std::packaged_task<TaskType> pt{sp_seq};
        std::future<ReduceType> fut{pt.get_future()};

        // explicitly launch thread
        threads[rank] = std::thread{std::move(pt),
                                    std::span{x.begin()+ibegin, x.begin()+iend},
                                    std::span{y.begin()+ibegin, y.begin()+iend}};
        futures[rank] = std::move(fut);
    }
    for (int rank = 0; rank < P; ++rank) {
        psums[rank] = futures[rank].get();
    }
    for (int rank = 0; rank < P; ++rank) {
        threads[rank].join();
    }
    return std::accumulate(psums.begin(), psums.end(), ReduceType(0));
}

ReduceType sp_openmp(std::span<NumberType> x, std::span<NumberType> y)
{
    ReduceType sum = 0;
    size_t N = x.size();

#pragma omp parallel for reduction(+:sum)
    for (size_t i = 0; i < N; ++i) {
        sum += x[i] * y[i];
    }
    return sum;
}

// Adapted from: Structured Parallel Programming, Listing 5.6
ReduceType sp_tbb(std::span<NumberType> x, std::span<NumberType> y)
{
    return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, x.size()), ReduceType(0),
        [=](tbb::blocked_range<size_t> &r, ReduceType in)
        {
            for (size_t i = r.begin(); i < r.end(); ++i) {
                in += x[i] * y[i];
            }
            return in;
        },
        std::plus<ReduceType>() // reducer for partial sums
        );
}

void sp_init(std::vector<NumberType> &x, std::vector<NumberType> &y, NumberType a, NumberType b)
{
    std::mt19937_64 gen{};
    std::uniform_real_distribution<NumberType> dist(a, b);

    for (size_t i = 0, n = x.size(); i < n; ++i) {
        x[i] = dist(gen);
    }
    for (size_t i = 0, n = y.size(); i < n; ++i) {
        y[i] = dist(gen);
    }
}

#endif // SCALAR_PRODUCT_ASYNC_HH
