## Memory Alignment

*Yes, it is matrix multiplication... again*

Contrary to what some people believe or may have heard somewhere, memory alignment is still required in some cases to achieve optimal performance. This lab assignment is one such case. First, a little introduction.

A typical case where data alignment is important is SIMD code, where loads and stores access large chunks of data with a single operation. In most processors, the L1 cache is designed to be able to read/write data at any alignment. Generally, even if a load/store is misaligned but does not cross the cache line boundary, it won't have any performance penalty. However, when a load or store crosses cache line boundary, such access requires two cache line reads (*split load/store*). It requires using a *split register*, which keeps the two parts and once both parts are fetched, they are combined into a single register. The number of split registers is limited. When executed sporadically, split accesses generally complete without any observable performance impact to overall execution. However, if that happens frequently, misaligned memory accesses will suffer delays.

Our simple matrix multiply generates SIMD instructions. For small-size matrices we use a regular version with loop interchange to achive cache-friendly accesses and vectorized code. For larger sizes, we rolled up a version that uses Loop Blocking. We provide all the boiler-plate code so that you only need to change a couple of lines.

When a matrix is misaligned, split loads happen very frequently which causes performance problems. In this lab you need to fix that. **Important**: It's not enough to only align the offset of a matrix, but also each row of the matrix has to be aligned. To do that you can insert dummy columns.

For AVX2 code, it is enough when each row is aligned at 32-byte boundary and for SSE and ARM Neon only 16-byte alignment is required. However, AVX-512 requires 64-byte alignment. To be on the safe side, you can align at the cacheline boundary. Keep in mind, in Apple processors (such as M1, M2 and later), L2 cache operates on 128-byte cache lines.

If you're actively using Intel's topdown methodology (TMA), then it will be reflected under `Memory_Bound -> L1_Bound -> Split Loads` category. When you improve the alignment, observe the change in the following events: `mem_inst_retired.split_loads`, and `mem_inst_retired.split_stores`.

Keep in mind, that performance penalty applies to both loads *and* stores, so all matrices should be aligned.

## Overview
Provided code allocates memory in a vector and performs matrix multiplication.<br/>
As mentioned in the problem statement, memory alignment is a very important performance detail.<br/>
The issue is that the code does not take it into account.
```c++
// Vector itself is not aligned to cache line size.
using Matrix = std::vector<float>;

// Rows are not aligned if N  * sizeof(float) is not divisible by cache line size.
int n_columns(int N) { return N; }
```

Amount of cases leading to performance cost due to misalignment can be measured using `mem_inst_retired.split_loads` and `mem_inst_retired.split_stores` counters.
```bash
$ perf stat -e mem_inst_retired.split_loads,mem_inst_retired.split_stores ./lab
 Performance counter stats for './lab':
     9,255,027,911      cpu_core/mem_inst_retired.split_loads/
     2,282,506,173      cpu_core/mem_inst_retired.split_stores/
```

## Solution
Let's address both issues independently.<br/>
To make vector itself aligned to cache line size, we can use a simple custom allocator.
```c++
template <typename T>
class CacheLineAlignedAllocator {
 public:
  using value_type = T;
  static std::align_val_t constexpr ALIGNMENT{CACHELINE_SIZE};
  [[nodiscard]] T* allocate(std::size_t N) {
    // Ensure that the allocation is aligned to the cache line size.
    return reinterpret_cast<T*>(::operator new[](N * sizeof(T), ALIGNMENT));
  }
  void deallocate(T* allocPtr, [[maybe_unused]] std::size_t N) {
    ::operator delete[](allocPtr, ALIGNMENT);
  }
};
template <typename T>
using AlignedVector = std::vector<T, CacheLineAlignedAllocator<T> >;

using Matrix = AlignedVector<float>;
```
To make each row aligned, all we need to do is add a bit of (unused) padding.
```c++
int n_columns(int N) {
  // Add dummy columns to ensure that each row is aligned to a cache line.
  constexpr int FLOATS_PER_CACHE_LINE = CACHELINE_SIZE / sizeof(float);
  int res = (N + FLOATS_PER_CACHE_LINE - 1) / FLOATS_PER_CACHE_LINE *
            FLOATS_PER_CACHE_LINE;
  assert(res >= N);
  assert((res * sizeof(float)) % CACHE_LINE_SIZE == 0);
  return res;
}
```

This results in few orders of magnitude less alignment-related performance issues.
```bash
$ getconf LEVEL1_DCACHE_LINESIZE
64
$ perf stat -e mem_inst_retired.split_loads,mem_inst_retired.split_stores ./lab
 Performance counter stats for './lab':
     32,151      cpu_core/mem_inst_retired.split_loads/
     47,909      cpu_core/mem_inst_retired.split_stores/
```

## Benchmark
Aligning entire vector as well as each row results in varying speedups, depeniding on row size.<br/>
Cache line on test Linux machine is 64 - which is probably why there is a minor speedup for rows smaller than 64.
```bash
Benchmark                      Time             CPU      Time Old      Time New       CPU Old       CPU New
-----------------------------------------------------------------------------------------------------------
bench1/_63                  -0.0312         -0.0313            33            32            33            32
bench1/_64                  -0.2869         -0.2874            15            10            15            10
bench1/_65                  -0.1897         -0.1904            15            12            15            12
bench1/_128                 -0.3015         -0.3021           118            83           118            83
bench1/_256                 -0.3229         -0.3235          1145           775          1144           774
bench1/_511                 -0.2487         -0.2489         10407          7819         10388          7802
bench1/_512                 -0.3391         -0.3397          9511          6286          9500          6273
bench1/_513                 -0.4056         -0.4059          9715          5774          9697          5761
bench1/_1024                -0.3739         -0.3744         82986         51954         82881         51850
Benchmark                      Time             CPU      Time Old      Time New       CPU Old       CPU New
-----------------------------------------------------------------------------------------------------------
bench1/_63                  -0.0435         -0.0435            37            35            37            35
bench1/_64                  -0.3285         -0.3284            16            10            16            10
bench1/_65                  -0.2329         -0.2330            16            12            16            12
bench1/_128                 -0.3198         -0.3199           122            83           122            83
bench1/_256                 -0.3537         -0.3537          1207           780          1205           779
bench1/_511                 -0.2842         -0.2842         10962          7847         10941          7832
bench1/_512                 -0.3684         -0.3684         10034          6338         10016          6326
bench1/_513                 -0.4243         -0.4244         10084          5805         10065          5793
bench1/_1024                -0.3927         -0.3927         88823         53945         88652         53842
Benchmark                      Time             CPU      Time Old      Time New       CPU Old       CPU New
-----------------------------------------------------------------------------------------------------------
bench1/_63                  -0.0356         -0.0356            36            34            36            34
bench1/_64                  -0.3084         -0.3084            15            11            15            11
bench1/_65                  -0.2384         -0.2383            16            12            16            12
bench1/_128                 -0.3326         -0.3325           125            83           125            83
bench1/_256                 -0.3575         -0.3576          1218           783          1216           781
bench1/_511                 -0.4171         -0.4171         13510          7875         13482          7858
bench1/_512                 -0.4201         -0.4202         10938          6343         10916          6329
bench1/_513                 -0.5372         -0.5371         12588          5826         12560          5814
bench1/_1024                -0.4380         -0.4379         96239         54086         96022         53972
```