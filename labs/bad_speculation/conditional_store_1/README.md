This lab assignment focuses on improving performance by reducing the number of branch mispredictions. In this lab, we have a large collection of key-value pairs from which we select only useful items. Such algorithms are widely used in many real-world applications, for example, texture compression. Input collection contains random key values which makes it hard for a modern CPU to predict whether they should be selected or not. Your task here is to reduce the number of branch mispredictions.


## Overview
```c++
// Select items which have S.first in range [lower..upper]
std::size_t select(std::array<S, N> &output, const std::array<S, N> &input,
                   const std::uint32_t lower, const std::uint32_t upper) {
  std::size_t count = 0;
  for (const auto item : input) {
    if ((lower <= item.first) && (item.first <= upper)) {
      output[count++] = item;
    }
  }
  return count;
}
```

Let's use `perf` to see what's the problem with this function.
```bash
$ perf stat ./lab
 Performance counter stats for './lab':

          1,964.35 msec task-clock                       #    1.000 CPUs utilized             
                 4      context-switches                 #    2.036 /sec                      
                 1      cpu-migrations                   #    0.509 /sec                      
               395      page-faults                      #  201.085 /sec                      
     <not counted>      cpu_atom/cycles/                                                        (0.00%)
     9,125,570,609      cpu_core/cycles/                 #    4.646 GHz                       
     <not counted>      cpu_atom/instructions/                                                  (0.00%)
     6,847,475,961      cpu_core/instructions/                                                
     <not counted>      cpu_atom/branches/                                                      (0.00%)
       906,827,307      cpu_core/branches/               #  461.643 M/sec                     
     <not counted>      cpu_atom/branch-misses/                                                 (0.00%)
       277,881,962      cpu_core/branch-misses/                                               
             TopdownL1 (cpu_core)                 #      8.6 %  tma_backend_bound      
                                                  #     58.2 %  tma_bad_speculation      <---------- BAD          
                                                  #     20.0 %  tma_frontend_bound     
                                                  #     13.2 %  tma_retiring      
```
Clearly, there is a problem with speculative execution.
```bash
$ perf stat -M tma_bad_speculation_group ./lab
 Performance counter stats for './lab':

    54,636,429,073      cpu_core/TOPDOWN.SLOTS/          #     57.4 %  tma_branch_mispredicts     <---------- BAD       
                                                  #      1.0 %  tma_machine_clears       (99.88%)
     7,284,857,209      cpu_core/topdown-retiring/                                              (99.88%)
    31,496,294,406      cpu_core/topdown-bad-spec/                                              (99.88%)
    31,496,294,406      cpu_core/topdown-br-mispredict/                                         (99.88%)
    11,570,067,333      cpu_core/topdown-fe-bound/                                              (99.88%)
     4,499,470,629      cpu_core/topdown-be-bound/                                              (99.88%)
       538,240,419      cpu_core/INT_MISC.UOP_DROPPING/                                         (99.88%)
        86,016,787      cpu_atom/TOPDOWN_BAD_SPECULATION.MACHINE_CLEARS/ #      0.5 %  tma_machine_clears       (0.06%)
     3,746,485,138      cpu_atom/TOPDOWN_RETIRING.ALL/                                          (0.06%)
     3,548,413,831      cpu_atom/CPU_CLK_UNHALTED.CORE/                                         (0.06%)
     6,277,241,806      cpu_atom/TOPDOWN_FE_BOUND.ALL/                                          (0.06%)
     5,024,993,404      cpu_atom/TOPDOWN_BE_BOUND.ALL/                                          (0.06%)
     5,941,901,326      cpu_atom/TOPDOWN_RETIRING.ALL/   #     23.7 %  tma_branch_mispredicts   (0.06%)
     4,411,641,379      cpu_atom/TOPDOWN_BAD_SPECULATION.MISPREDICT/                                        (0.06%)
     3,723,367,765      cpu_atom/CPU_CLK_UNHALTED.CORE/                                         (0.06%)
     3,169,027,007      cpu_atom/TOPDOWN_FE_BOUND.ALL/                                          (0.06%)
     5,004,919,795      cpu_atom/TOPDOWN_BE_BOUND.ALL/                                          (0.06%)
```
The main issue seems to be branch misprediction, to confirm, we can gather branch-specific events.
```bash
$ perf stat -e branches,branch-misses ./lab
 Performance counter stats for './lab':
       906,738,050      cpu_core/branches/
       271,789,734      cpu_core/branch-misses/          #   29.97% of all branches
```
Quite a lot of branching in total, and **very** high branch miss rate.

## Solution
As presented above, the issue is clearly related to branch mispredictions, and there is only one conditional branch in the `select` function that we can eliminate.<br/>
Instead of branching, we can always overwrite `output` in current `count` and increase `count` when conditions are met.
```c++
// Select items which have S.first in range [lower..upper]
std::size_t select(std::array<S, N> &output, const std::array<S, N> &input,
                   const std::uint32_t lower, const std::uint32_t upper) {
  std::size_t count = 0;
  for (const auto item : input) {
    output[count] = item;
    // move to the new index when conditions are met
    count += (lower <= item.first) && (item.first <= upper);
  }
  return count;
}
```

To double check, let's see branch-related counters again.
```bash
$ perf stat -e branches,branch-misses ./lab 
 Performance counter stats for './lab':
       491,100,063      cpu_core/branches/                                                    
           443,144      cpu_core/branch-misses/          #    0.09% of all branches 
```
As you can see, there is quite a bit less branching and **significantly** less branch misses.

## Benchmark
Getting rid of branching results in ~85% speedup.
```bash
Benchmark                                 Time             CPU      Time Old      Time New       CPU Old       CPU New
----------------------------------------------------------------------------------------------------------------------
bench1/iterations:10000                -0.8486         -0.8486           192            29           191            29
Benchmark                                 Time             CPU      Time Old      Time New       CPU Old       CPU New
----------------------------------------------------------------------------------------------------------------------
bench1/iterations:10000                -0.8331         -0.8331           193            32           192            32
Benchmark                                 Time             CPU      Time Old      Time New       CPU Old       CPU New
----------------------------------------------------------------------------------------------------------------------
bench1/iterations:10000                -0.8560         -0.8562           191            27           190            27
```