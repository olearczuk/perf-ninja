
#include <memory>

#include "benchmark/benchmark.h"
#include "solution.h"

static void bench_partial_sum(benchmark::State &state) {
  InputVector inA;
  init(inA);

  OutputVector outB;
  zero(outB, (int)inA.size());

  for (auto _ : state) {
    imageSmoothing(inA, radius, outB);
    benchmark::DoNotOptimize(outB);
  }
}

// Register the function as a benchmark
BENCHMARK(bench_partial_sum)->Unit(benchmark::kNanosecond);

// Run the benchmark
BENCHMARK_MAIN();
