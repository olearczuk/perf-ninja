This is a lab about [loop interchange](https://en.wikipedia.org/wiki/Loop_interchange), which is more advanced than the previous one.

In this lab assignment you will optimize [Gaussian blur](https://en.wikipedia.org/wiki/Gaussian_blur) algorithm applied to a grayscale image.
Modern cameras have good matrices and produce big files. How fast can modern CPU filter a camera shot?
Significant speedup has been already achieved by two passes of 1-dimensional digital filter instead of a plain 2D convolution.

## Solution
### Identifying the problem
Let's use [Top Down Analysis Method](https://perfwiki.github.io/main/top-down-analysis/) to find the problem

```bash
# Build
$ cmake -E make_directory build && cd build && cmake -DCMAKE_BUILD_TYPE=Release  -DCMAKE_CXX_FLAGS="-fno-inline -fno-omit-frame-pointer" .. && cmake --build . --config Release --parallel $(nproc)
$ perf stat ./validate ../pexels-pixabay-434334.pbm ../output-golden.pgm
Validation Successful

 Performance counter stats for './validate ../pexels-pixabay-434334.pbm ../output-golden.pgm':

             69.63 msec task-clock                       #    0.985 CPUs utilized             
                 2      context-switches                 #   28.724 /sec                      
                 2      cpu-migrations                   #   28.724 /sec                      
             2,385      page-faults                      #   34.253 K/sec                     
        77,038,339      cpu_atom/cycles/                 #    1.106 GHz                         (0.30%)
       288,790,264      cpu_core/cycles/                 #    4.148 GHz                         (99.63%)
        46,890,973      cpu_atom/instructions/           #    0.61  insn per cycle              (0.37%)
       542,252,396      cpu_core/instructions/           #    7.04  insn per cycle              (99.63%)
        10,014,152      cpu_atom/branches/               #  143.823 M/sec                       (0.37%)
        35,053,978      cpu_core/branches/               #  503.444 M/sec                       (99.63%)
           617,682      cpu_atom/branch-misses/          #    6.17% of all branches             (0.37%)
            39,389      cpu_core/branch-misses/          #    0.39% of all branches             (99.63%)
             TopdownL1 (cpu_core)                 #     67.5 %  tma_backend_bound  <------------------------------------
                                                  #      0.4 %  tma_bad_speculation    
                                                  #      2.0 %  tma_frontend_bound     
                                                  #     30.2 %  tma_retiring             (99.63%)
             TopdownL1 (cpu_atom)                 #     18.1 %  tma_bad_speculation    
                                                  #     17.0 %  tma_retiring             (0.37%)
                                                  #     19.9 %  tma_backend_bound      
                                                  #     45.0 %  tma_frontend_bound       (0.37%)
$ perf stat -M tma_backend_bound_group ./validate ../pexels-pixabay-434334.pbm ../output-golden.pgm
Validation Successful

 Performance counter stats for './validate ../pexels-pixabay-434334.pbm ../output-golden.pgm':

     1,866,177,110      cpu_core/TOPDOWN.SLOTS/          #     10.2 %  tma_core_bound         
                                                  #     60.0 %  tma_memory_bound         (99.55%) <------------------------
       519,602,253      cpu_core/topdown-retiring/                                              (99.55%)
     1,119,706,265      cpu_core/topdown-mem-bound/                                             (99.55%)
         7,318,341      cpu_core/topdown-bad-spec/                                              (99.55%)
        29,273,365      cpu_core/topdown-fe-bound/                                              (99.55%)
     1,309,983,147      cpu_core/topdown-be-bound/                                              (99.55%)
       231,009,732      cpu_atom/CPU_CLK_UNHALTED.CORE/  #      0.8 %  tma_core_bound         
                                                  #     28.8 %  tma_resource_bound       (0.45%)
         9,411,510      cpu_atom/TOPDOWN_BE_BOUND.ALLOC_RESTRICTIONS/                                        (0.45%)
       342,035,072      cpu_atom/TOPDOWN_BE_BOUND.ALL/                                          (0.45%)
```

As you can see, the program is backend-bound and the main culprit is the memory.

### perf mem-stores
```bash
$ perf mem record --call-graph=dwarf ./validate ../pexels-pixabay-434334.pbm ../output-golden.pgm 
$ perf report
```
For `mem-stores`, you can see both `filterVertically` and `filterHorizontally` are very similar
```bash
Samples: 283  of event 'cpu_core/mem-stores/P', Event count (approx.): 31053779
  Children      Self  Command   Shared Object        Symbol
+  100.00%     0.00%  validate  validate             [.] _start
+  100.00%     0.00%  validate  libc.so.6            [.] __libc_start_main@@GLIBC_2.34
+  100.00%     0.00%  validate  libc.so.6            [.] __libc_start_call_main
+  100.00%     0.00%  validate  validate             [.] main
+   98.33%     0.00%  validate  validate             [.] blur(unsigned char*, unsigned char const*, int, int, unsigned char*)
+   49.28%    48.83%  validate  validate             [.] filterVertically(unsigned char*, unsigned char const*, int, int, int const*, int, int) [clone .constprop.0]
+   49.05%    48.45%  validate  validate             [.] filterHorizontally(unsigned char*, unsigned char const*, int, int, int const*, int, int) [clone .constprop.0]
```
### perf mem-loads
However, for `mem-loads`, the difference is substantial
```bash
Samples: 215  of event 'cpu_core/mem-loads,ldlat=30/', Event count (approx.): 60850
  Children      Self  Command   Shared O  Symbol
+   81.16%    79.66%  validate  validate  [.] filterVertically(unsigned char*, unsigned char const*, int, int, int const*, int, int) [clone .constprop.0]
+    0.92%     0.92%  validate  validate  [.] filterHorizontally(unsigned char*, unsigned char const*, int, int, int const*, int, int) [clone .constprop.0]
```
### perf cache-misses
This is also reflected in cache misses
```bash
$ perf record -e cache-misses ./validate ../pexels-pixabay-434334.pbm ../output-golden.pgm 
$ perf report
Samples: 286  of event 'cpu_core/cache-misses/', Event count (approx.): 4772395
Overhead  Command   Shared Object      Symbol
  31.15%  validate  [kernel.kallsyms]  [k] clear_page_erms
  23.58%  validate  [kernel.kallsyms]  [k] _copy_to_iter
  19.28%  validate  validate           [.] filterVertically(unsigned char*, unsigned char const*, int, int, int const*, int, int) [clone .constprop.0]
  12.61%  validate  libc.so.6          [.] __memcmp_avx2_movbe
   3.93%  validate  validate           [.] filterHorizontally(unsigned char*, unsigned char const*, int, int, int const*, int, int) [clone .constprop.0]
```

This suggests we should take a deeper look at `filterVertically` function.

### `filterVertically` inefficiency
`filterVertically` code boils down to nested loops with strided memory accesses.<br/>
```c++
for (int c = 0; c < width; c++) {
  // Top part of line, partial kernel
  for (int r = 0; r < std::min(radius, height); r++) {
    // ...
    output[r * width + c] = static_cast<uint8_t>(value);
  }

  // Middle part of computations with full kernel
  for (int r = radius; r < height - radius; r++) {
    // ...
    output[r * width + c] = static_cast<uint8_t>(value);
  }

  // Bottom part of line, partial kernel
  for (int r = std::max(radius, height - radius); r < height; r++) {
    // ...
    output[r * width + c] = static_cast<uint8_t>(value);
  }
}
```
`output[r * width + c]` access pattern means that we are going to "jump" in memory by `width` because `r` is in the internal loop.<br/>
This results in inefficient memory accesses and does not utilize CPU cache as well as it could.<br/>
Switching order of `c` and `r` loops would result in more sequential memory access.
```c++
// Top part of line, partial kernel
for (int r = 0; r < std::min(radius, height); r++) {
  for (int c = 0; c < width; c++) {
    // ...
    output[r * width + c] = static_cast<uint8_t>(value);
  }
}

// Middle part of computations with full kernel
for (int r = radius; r < height - radius; r++) {
  for (int c = 0; c < width; c++) {
    // ...
    output[r * width + c] = static_cast<uint8_t>(value);
  }
}

// Bottom part of line, partial kernel
for (int r = std::max(radius, height - radius); r < height; r++) {
  for (int c = 0; c < width; c++) {
    // ...
    output[r * width + c] = static_cast<uint8_t>(value);
  }
}
```

## Benchmark
Interchanging loops results in 75% - 80% speedup.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.7465         -0.7467      42853770      10863025      42785426      10838043
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8111         -0.8112      52774192       9969937      52684018       9946041
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.7632         -0.7633      42472419      10057805      42406737      10039553
```