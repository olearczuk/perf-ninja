# Software memory prefetching

When the CPU data prefetcher cannot figure out the memory access pattern, software prefetching comes in handy. The idea is to use special instructions that tell the CPU: "Hey, I plan to use this memory location a bit later, could you fetch it for me while I do other stuff so it waits for me when I am back".

In GCC and CLANG, you can use `__builtin_prefetch` to ask the CPU to prefetch data. Say, for example, that you are going to access an element of array `my_array[index]`, where `index` is some random number. To prefetch it, you will use `__builtin_prefetch(&my_array[index]);` or `__builtin_prefetch(&my_array + index);`.

Prefetching can benefit the performance, but it can also hurt the performance. It benefits it if the piece of data you are trying to access is not in the data cache. It hurts it if it is. So most of the time, it pays off when there are random memory accesses on a large data structure, such as a tree or a hash map.

An additional prerequisite for the speedup with prefetching is that between the time you request prefetching, and the time you actually access your data, some time needs to pass (known as "prefetching window"). Immediately accessing data that you want to prefetch will not give the expected results.

Authored-by: @ibogosavljevic


## Solution
The problem is a relatively problem that computes sum of values in a hash map for given collection of keys.
```c++
int solution(const hash_map_t *hash_map, const std::vector<int> &lookups) {
  int result = 0;
  for (int val : lookups) {
    if (hash_map->find(val)) {
      result += getSumOfDigits(val);
    }
  }
  return result;
}
```
We can quickly check what's the main bottleneck in the program using `perf`.
```bash
$ cmake -E make_directory build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build . --config Release --parallel $(nproc)
$ perf stat ./validate
Validation Successful

 Performance counter stats for './validate':

          3,370.17 msec task-clock                       #    0.997 CPUs utilized             
               300      context-switches                 #   89.016 /sec                      
                 6      cpu-migrations                   #    1.780 /sec                      
             1,212      page-faults                      #  359.626 /sec                      
     <not counted>      cpu_atom/cycles/                                                        (0.00%)
    13,785,814,156      cpu_core/cycles/                 #    4.091 GHz                       
     <not counted>      cpu_atom/instructions/                                                  (0.00%)
     7,912,608,024      cpu_core/instructions/                                                
     <not counted>      cpu_atom/branches/                                                      (0.00%)
       741,924,817      cpu_core/branches/               #  220.145 M/sec                     
     <not counted>      cpu_atom/branch-misses/                                                 (0.00%)
        47,554,806      cpu_core/branch-misses/                                               
             TopdownL1 (cpu_core)                 #     68.7 %  tma_backend_bound      <-----------------
                                                  #     18.1 %  tma_bad_speculation    
                                                  #      2.6 %  tma_frontend_bound     
                                                  #     10.5 %  tma_retiring   
```
As you can se, the program is heavily memory bound.<br/>
This makes sense because we are repeatedly fetching memory from a hash map which results in a lot of cache misses - there is no clear pattern in memory accesses so CPU is not able to prefetch it.<br/>
We can help it by explicitly fetching memory using software operations - we know what keys we are going to need in the future.<br/>
```c++
int solution(const hash_map_t *hash_map, const std::vector<int> &lookups) {
  constexpr int PREFETCH_DISTANCE = 16;

  int result = 0;
  const auto process_index = [&](int i) {
    int val = lookups[i];
    if (hash_map->find(val)) {
      result += getSumOfDigits(val);
    }
  };

  // Process elements while prefetching the next bucket.
  int i = 0;
  for (; i < lookups.size() - PREFETCH_DISTANCE; ++i) {
    process_index(i);
    __builtin_prefetch(
        hash_map->find_bucket_ptr(lookups[i + PREFETCH_DISTANCE]));
  }

  // Process the remaining elements without prefetching
  for (; i < lookups.size(); ++i) {
    process_index(i);
  }

  return result;
}
```

## Benchmark
Explicit software prefetching results in ~70% speedup.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.7134         -0.7136            66            19            66            19
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.7084         -0.7085            65            19            65            19
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.7141         -0.7141            66            19            66            19
```