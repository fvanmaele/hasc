#include <vectorclass/vectorclass.h>
#include <benchmark/benchmark.h>
#include <vector>

#include "functions.hh"
#include "midpoint_rule.hh"

using namespace hasc_ex01;

static void BM_midpoint_f1_seq(benchmark::State& state)
{
    double a = 0;
    double b = 1;
    Func1 f;
    
    for (auto _ : state)
        benchmark::DoNotOptimize(midpoint_rule_seq(a, b, state.range(0), f));
}

static void BM_midpoint_f2_seq(benchmark::State& state)
{
    double a = 0;
    double b = 1;
    Func2 f;
    
    for (auto _ : state)
        benchmark::DoNotOptimize(midpoint_rule_seq(a, b, state.range(0), f));
}

static void BM_midpoint_f1_vec(benchmark::State& state)
{
    double a = 0;
    double b = 1;
    Func1 f;
    
    for (auto _ : state)
        benchmark::DoNotOptimize(midpoint_rule_vec(a, b, state.range(0), f));
}

static void BM_midpoint_f2_vec(benchmark::State& state)
{
    double a = 0;
    double b = 1;
    Func2 f;
    
    for (auto _ : state)
        benchmark::DoNotOptimize(midpoint_rule_vec(a, b, state.range(0), f));
}

BENCHMARK(BM_midpoint_f1_seq)->Range(10, 10<<20);
BENCHMARK(BM_midpoint_f2_seq)->Range(10, 10<<20);
BENCHMARK(BM_midpoint_f1_vec)->Range(10, 10<<20);
BENCHMARK(BM_midpoint_f2_vec)->Range(10, 10<<20);

BENCHMARK_MAIN();