This lab assignment focuses on improving performance by reducing the number of branch target mispredictions. In this lab,
we have a collection of objects of three distinct classes that inherit from a shared base class. An array of objects is
created in some arbitrary order and then repeatedly iterated. During the iteration, a virtual method is called on each
object. The CPU does not know which exact method will be called beforehand (it does not know the *target* of the call),
which slows its execution down.

Your task here is to reduce the number of branch target mispredictions by making the virtual method calls more predictable.

Authored-by: Jakub Beránek (@Kobzol)

## Overview
```c++
std::default_random_engine generator(0);
std::uniform_int_distribution<std::uint32_t> distribution(0, 2);

for (std::size_t i = 0; i < N; i++) {
    int value = distribution(generator);
    if (value == 0) {
        array.push_back(std::make_unique<ClassA>());
    } else if (value == 1) {
        array.push_back(std::make_unique<ClassB>());
    } else {
        array.push_back(std::make_unique<ClassC>());
    }
}
```
Let's use `perf` to see what's the problem with this function.
```bash
$ perf stat ./lab
 Performance counter stats for './lab':

            790.76 msec task-clock                       #    0.999 CPUs utilized             
                 3      context-switches                 #    3.794 /sec                      
                 2      cpu-migrations                   #    2.529 /sec                      
             2,749      page-faults                      #    3.476 K/sec                     
     1,509,070,037      cpu_atom/cycles/                 #    1.908 GHz                         (0.06%)
     3,664,985,927      cpu_core/cycles/                 #    4.635 GHz                         (99.93%)
     1,025,871,565      cpu_atom/instructions/           #    0.68  insn per cycle              (0.07%)
     6,609,054,609      cpu_core/instructions/           #    4.38  insn per cycle              (99.93%)
       215,865,146      cpu_atom/branches/               #  272.985 M/sec                       (0.07%)
     2,183,800,250      cpu_core/branches/               #    2.762 G/sec                       (99.93%)
        10,946,240      cpu_atom/branch-misses/          #    5.07% of all branches             (0.07%)
           362,175      cpu_core/branch-misses/          #    0.17% of all branches             (99.93%)
             TopdownL1 (cpu_core)                 #      0.4 %  tma_backend_bound      
                                                  #      0.4 %  tma_bad_speculation      <---------- BAD          
                                                  #     61.5 %  tma_frontend_bound     
                                                  #     37.6 %  tma_retiring             (99.93%)
             TopdownL1 (cpu_atom)                 #     18.0 %  tma_bad_speculation    
                                                  #     18.9 %  tma_retiring             (0.07%)
                                                  #     24.1 %  tma_backend_bound      
                                                  #     39.1 %  tma_frontend_bound       (0.07%)
```
Looks like there is a problem with speculative execution.
```bash
 Performance counter stats for './lab':

    31,137,583,044      cpu_core/TOPDOWN.SLOTS/          #     35.4 %  tma_branch_mispredicts      <---------- BAD       
                                                  #      0.6 %  tma_machine_clears     
     2,442,163,376      cpu_core/topdown-retiring/                                            
    11,111,843,360      cpu_core/topdown-bad-spec/                                            
    11,111,843,360      cpu_core/topdown-br-mispredict/                                       
    15,507,737,437      cpu_core/topdown-fe-bound/                                            
     2,320,055,207      cpu_core/topdown-be-bound/                                            
       173,753,746      cpu_core/INT_MISC.UOP_DROPPING/                                       
     <not counted>      cpu_atom/TOPDOWN_BAD_SPECULATION.MACHINE_CLEARS/                                        (0.00%)
     <not counted>      cpu_atom/TOPDOWN_RETIRING.ALL/                                          (0.00%)
     <not counted>      cpu_atom/CPU_CLK_UNHALTED.CORE/                                         (0.00%)
     <not counted>      cpu_atom/TOPDOWN_FE_BOUND.ALL/                                          (0.00%)
     <not counted>      cpu_atom/TOPDOWN_BE_BOUND.ALL/                                          (0.00%)
     <not counted>      cpu_atom/TOPDOWN_RETIRING.ALL/                                          (0.00%)
     <not counted>      cpu_atom/TOPDOWN_BAD_SPECULATION.MISPREDICT/                                        (0.00%)
     <not counted>      cpu_atom/CPU_CLK_UNHALTED.CORE/                                         (0.00%)
     <not counted>      cpu_atom/TOPDOWN_FE_BOUND.ALL/                                          (0.00%)
     <not counted>      cpu_atom/TOPDOWN_BE_BOUND.ALL/                                          (0.00%)
```
The main issue seems to be branch misprediction, to confirm, we can gather branch-specific events.
```bash
$ perf stat -e branches,branch-misses ./lab
Performance counter stats for './lab':
       641,943,029      cpu_core/branches/                                                    
       133,573,096      cpu_core/branch-misses/          #   20.81% of all branches     
```
That's a *very* high branch miss rate.

## Solution
There are no `if` or `switch` statements, however, virtual calls are special *indirect* case of branch prediction.<br/>
Branch predictor tries to spot a pattern in the source addresses, however, it struggles to because types are in a random order.<br/>
We can ensure to store the objects of the same type in a sequence to make the work of branch predictor much easier.
```c++
void generateObjects(InstanceArray& array) {
  std::default_random_engine generator(0);
  std::uniform_int_distribution<std::uint32_t> distribution(0, 2);

  std::size_t a_count = 0, b_count = 0, c_count = 0;
  for (std::size_t i = 0; i < N; i++) {
    int value = distribution(generator);
    if (value == 0) {
      ++a_count;
    } else if (value == 1) {
      ++b_count;
    } else {
      ++c_count;
    }
  }

  for (std::size_t i = 0; i < a_count; i++) {
    array.push_back(std::make_unique<ClassA>());
  }
  for (std::size_t i = 0; i < b_count; i++) {
    array.push_back(std::make_unique<ClassB>());
  }
  for (std::size_t i = 0; i < c_count; i++) {
    array.push_back(std::make_unique<ClassC>());
  }
}
```

As you can see, this change heavily reduces number of branch misses.
```bash
$ perf stat -e branches,branch-misses ./lab
Performance counter stats for './lab':
     2,112,619,010      cpu_core/branches/                                                    
           354,133      cpu_core/branch-misses/          #    0.02% of all branches   
```

## Benchmark
Resolving virtual calls misprediction rate results in ~80% speedup.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8030         -0.8030           360            71           360            71
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8027         -0.8027           362            71           361            71
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8034         -0.8034           359            71           359            71
```