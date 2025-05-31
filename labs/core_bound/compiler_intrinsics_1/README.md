
This is a lab about using [compiler intrinsics](https://en.wikipedia.org/wiki/Intrinsic_function) to speed up parts of the code, where compilers fail to generate optimal code.

The kernel in this lab assignment is a part of the Average ImageSmoothing algorithm, which is reduced to 1 dimension and lacks division part. The algorithm uses sliding window approach to compute a sum in the subrange [-radius .. +radius]. It is a very fast approach compared to a classical Gaussian blur.

Author: @adamf88.

## Solution
Main hot loop basically performs of a difference between values starting in `pos + radius` and `pos - radius - 1`.<br/>
Prefix sums are good candidates for SIMD because of [Parallel Prefix Algorithm](https://en.wikipedia.org/wiki/Prefix_sum#Parallel_algorithms).<br/>
We can load chunks of corresponding arrays into SIMD registers, compute prefix sum and store values in `output`.
```cpp
  limit = size - radius;
  for (; pos < limit; ++pos) {
    currentSum -= input[pos - radius - 1];
    currentSum += input[pos + radius];
    output[pos] = currentSum;
  }
```

## Benchmark
### SSE2
Solution using 128-bit registers offers ~75% speedup.<br/>
```bash
Benchmark                           Time             CPU      Time Old      Time New       CPU Old       CPU New
----------------------------------------------------------------------------------------------------------------
bench_partial_sum                -0.7474         -0.7476         17148          4332         17133          4324
Benchmark                           Time             CPU      Time Old      Time New       CPU Old       CPU New
----------------------------------------------------------------------------------------------------------------
bench_partial_sum                -0.7570         -0.7570         17156          4169         17140          4165
Benchmark                           Time             CPU      Time Old      Time New       CPU Old       CPU New
----------------------------------------------------------------------------------------------------------------
bench_partial_sum                -0.7340         -0.7344         17137          4558         17123          4548
```

### AVX2
Solution using 256-bit registers is slightly faster than the previous one, offering ~78% speedup.<br/>
The reason for small gain is that in AVX2, shift operators on 256-bit registers are not performed on the entire register but independently on each half.<br.>
This means parallel prefix sum can not be achieved simply with bunch of shifts and adds.
```bash
Benchmark                           Time             CPU      Time Old      Time New       CPU Old       CPU New
----------------------------------------------------------------------------------------------------------------
bench_partial_sum                -0.7851         -0.7852         17181          3692         17163          3687
Benchmark                           Time             CPU      Time Old      Time New       CPU Old       CPU New
----------------------------------------------------------------------------------------------------------------
bench_partial_sum                -0.7740         -0.7742         17144          3875         17125          3867
Benchmark                           Time             CPU      Time Old      Time New       CPU Old       CPU New
----------------------------------------------------------------------------------------------------------------
bench_partial_sum                -0.7867         -0.7867         17125          3653         17111          3649
```