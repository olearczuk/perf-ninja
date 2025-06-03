This is a second lab about using [compiler intrinsics](https://en.wikipedia.org/wiki/Intrinsic_function) to speed up parts of the code, where compilers fail to generate optimal code.

The task of this lab assignment is to find the longest line in a file. There is a way to find end-of-line characters in a parallel way if you utilize compiler intrinsics.

Bonus exercise: whether solution that uses intrinsics is faster than the baseline is heavily affected by the input data. Run your solution on different input files to determine the speedup/slowdown.

The idea for this lab was proposed by Yuriy Lyfenko (@obender12).

Co-authored-by: Andrew Evstyukhin (@andrewevstyukhin)

Co-authored-by: Jakub Beránek (@Kobzol)


## Solution
The problem with naive approach is that it checks for newlines one byte at a time.<br/>
Instead, we can check multiple bytes at once for newlines, treat the result as a bitmap and calculate distances between consecutive newlines using the bitmap.<br/>
Using comipler intrinsics to trigger SIMD operations is the perfect candidate.

## Benchmark
Using AVX2-based approach results in ~8x speedup.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8878         -0.8878           306            34           306            34
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8845         -0.8845           306            35           306            35
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8836         -0.8837           306            36           306            36
```