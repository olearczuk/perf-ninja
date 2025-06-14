Welcome to the next lab assignment, where we will fight branch mispredictions by replacing them with lookup tables. The code in this lab assignment maps values from `[0;150]` into buckets, which involves a lot of comparisons, and so, branches. To solve this assignment you need to figure out a way how to replace *all* hard-to-predict branches.

## Overview
```c++
static std::size_t mapToBucket(std::size_t v) {
  //   size of a bucket
  if (v < 13)
    return 0;  //   13
  else if (v < 29)
    return 1;  //   16
  else if (v < 41)
    return 2;  //   12
  else if (v < 53)
    return 3;  //   12
  else if (v < 71)
    return 4;  //   18
  else if (v < 83)
    return 5;  //   12
  else if (v < 100)
    return 6;  //   17
  return DEFAULT_BUCKET;
}
```

Let's use `perf` to see what's the problem with this function.
```bash
$ perf stat ./lab
 Performance counter stats for './lab':

          1,224.73 msec task-clock                       #    0.999 CPUs utilized             
                10      context-switches                 #    8.165 /sec                      
                 3      cpu-migrations                   #    2.450 /sec                      
             1,167      page-faults                      #  952.863 /sec                      
     1,798,687,741      cpu_atom/cycles/                 #    1.469 GHz                         (0.06%)
     5,665,625,134      cpu_core/cycles/                 #    4.626 GHz                         (99.84%)
     4,793,508,321      cpu_atom/instructions/           #    2.67  insn per cycle              (0.15%)
     5,830,924,937      cpu_core/instructions/           #    3.24  insn per cycle              (99.84%)
       537,591,514      cpu_atom/branches/               #  438.947 M/sec                       (0.16%)
     1,737,124,752      cpu_core/branches/               #    1.418 G/sec                       (99.84%)
           426,582      cpu_atom/branch-misses/          #    0.08% of all branches             (0.16%)
       181,226,266      cpu_core/branch-misses/          #   33.71% of all branches             (99.84%)
             TopdownL1 (cpu_core)                 #      9.3 %  tma_backend_bound      
                                                  #     53.6 %  tma_bad_speculation      <---------- BAD
                                                  #     23.1 %  tma_frontend_bound     
                                                  #     14.0 %  tma_retiring             (99.84%)
             TopdownL1 (cpu_atom)                 #     -0.1 %  tma_bad_speculation    
                                                  #     59.1 %  tma_retiring             (0.16%)
                                                  #     26.1 %  tma_backend_bound      
                                                  #     14.9 %  tma_frontend_bound       (0.16%)
```
Clearly, there is a problem with speculative execution.
```bash
$ perf stat -M tma_bad_speculation_group ./lab
 Performance counter stats for './lab':

    34,107,761,311      cpu_core/TOPDOWN.SLOTS/          #     52.7 %  tma_branch_mispredicts     <---------- BAD       
                                                  #      1.0 %  tma_machine_clears       (99.43%)
     4,815,213,360      cpu_core/topdown-retiring/                                              (99.43%)
    18,057,050,105      cpu_core/topdown-bad-spec/                                              (99.43%)
    18,057,050,105      cpu_core/topdown-br-mispredict/                                         (99.43%)
     8,292,867,455      cpu_core/topdown-fe-bound/                                              (99.43%)
     3,076,386,313      cpu_core/topdown-be-bound/                                              (99.43%)
       349,442,614      cpu_core/INT_MISC.UOP_DROPPING/                                         (99.43%)
        32,313,976      cpu_atom/TOPDOWN_BAD_SPECULATION.MACHINE_CLEARS/ #      0.3 %  tma_machine_clears       (0.32%)
     5,845,032,401      cpu_atom/TOPDOWN_RETIRING.ALL/                                          (0.32%)
     2,440,428,993      cpu_atom/CPU_CLK_UNHALTED.CORE/                                         (0.32%)
     3,255,889,270      cpu_atom/TOPDOWN_FE_BOUND.ALL/                                          (0.32%)
     2,364,111,943      cpu_atom/TOPDOWN_BE_BOUND.ALL/                                          (0.32%)
     4,620,731,486      cpu_atom/TOPDOWN_RETIRING.ALL/   #      8.7 %  tma_branch_mispredicts   (0.25%)
       840,800,691      cpu_atom/TOPDOWN_BAD_SPECULATION.MISPREDICT/                                        (0.25%)
     1,926,601,303      cpu_atom/CPU_CLK_UNHALTED.CORE/                                         (0.25%)
     2,306,380,906      cpu_atom/TOPDOWN_FE_BOUND.ALL/                                          (0.25%)
     1,835,591,168      cpu_atom/TOPDOWN_BE_BOUND.ALL/                                          (0.25%)
```
The main issue seems to be branch misprediction, to confirm, we can gather branch-specific events.
```bash
$ perf stat -e branches,branch-misses ./lab
 Performance counter stats for './lab':
     1,741,203,776      cpu_core/branches/                                                    
       181,234,496      cpu_core/branch-misses/          #   10.41% of all branches   
```
Quite a lot of branching in total, and **very** high branch miss rate.

## Solution
The total universe of cases is relatively small, so we can eliminate branching by creating a lookup table with a fairly small memory cost.
```c++
// use short to decrease memory footprint
const std::array<short, 101> lookup_table {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                 // 0-12: bucket 0
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,        // 13-28: bucket 1
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,                    // 29-40: bucket 2
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,                    // 41-52: bucket 3
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 53-70: bucket 4
      5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,                    // 71-82: bucket 5
      6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,     // 83-99: bucket 6
      static_cast<short>(DEFAULT_BUCKET)                     // 100: default
};

static std::size_t mapToBucket(std::size_t v) {
  auto bucket = lookup_table[std::min(v, lookup_table.size() - 1)];
  return static_cast<std::size_t>(bucket);
}
```

As you can see, this change heavily reduces amount of branch mispredictions.
```bash
 Performance counter stats for './lab':
     1,491,129,065      cpu_core/branches/                                                    
            43,463      cpu_core/branch-misses/          #    0.00% of all branches 
```

## Benchmark
Removing branching results in ~85% speedup.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8774         -0.8774          4513           553          4507           553
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8783         -0.8783          4494           547          4487           546
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8787         -0.8787          4507           547          4501           546
```