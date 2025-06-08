This lab assignment focuses on improving performance by eliminating false sharing. In this lab, we
have several threads that modify data located close together in memory in parallel. This causes a lot
of overhead, because the individual cores must transfer cache lines containing the modified data amongst
themselves to satisfy cache coherence.

Your task here is to eliminate the false sharing by making sure that each thread will access a separate
cache line.

Expected speedup: at least 60%.

Authored-by: Jakub Beránek (@Kobzol)

## Solution
Program spawns multiple threads, each accessing a different index in provided vector.
```c++
struct Accumulator {
    std::atomic<uint32_t> value = 0;
};
std::vector<Accumulator> accumulators(thread_count);
```
On the surface everything looks good since each thread is accessing a different slot in the vector.<br/>
There is one problem - memory is fetched in cache lines - typicall 64 bytes in size.<br/>
However, `sizeof(std::atomic<uint32_t>) == 4` which means that 2 atomics are going to share a cache lines.<br/>
This means that any changes done to index `i` are going to invalidate this cache line, leading to `i + 1`, `i + 2`, ... being affected as well.<br/>
Because of that threads share memory even if they don't access exactly the same data - this results in expensive cache coherency work which heavily slows down the program.<br/>
The fix is very simple - enforce that each `Accumulator` does not share cache lines with a different value by aligning each to cache line size `64`.<br/>
```c++
struct alignas(64) Accumulator {
    std::atomic<uint32_t> value = 0;
};
```

## Benchmark
Getting rid of false sharing results in 85% - 90% speedup.
```bash
Benchmark                          Time             CPU      Time Old      Time New       CPU Old       CPU New
---------------------------------------------------------------------------------------------------------------
bench1/real_time                -0.8581         -0.8860           251            36           205            23
Benchmark                          Time             CPU      Time Old      Time New       CPU Old       CPU New
---------------------------------------------------------------------------------------------------------------
bench1/real_time                -0.8621         -0.8786           256            35           189            23
Benchmark                          Time             CPU      Time Old      Time New       CPU Old       CPU New
---------------------------------------------------------------------------------------------------------------
bench1/real_time                -0.8744         -0.8889           272            34           206            23
```