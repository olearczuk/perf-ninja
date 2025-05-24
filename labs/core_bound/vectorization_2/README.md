This is a second lab about [auto vectorization](https://llvm.org/docs/Vectorizers.html). The subject of this lab assignment is a part of a checksum algorithm from the 80s, which has risen from the popularity of the Internet and [accompanying needs to validate transmitted packets](https://www.alpharithms.com/internet-checksum-calculation-steps-044921/). Even the problem is old, similar issues may exist nowadays in production code.

Modern compilers handle simple loops very well, including horizontal additions. In this lab computations inside the loop are slightly more difficult: we do an "add carry" operation. Some compilers recognize "add carry" and others don't. The [carry flag](https://en.wikipedia.org/wiki/Binary_number#Binary_arithmetic) is still a dark area in C++ while it exists more than 40 years. In this lab assignment, you will practice fixing auto-vectorization, which will improve performance significantly.

Hint: the [RFC 1071](http://www.faqs.org/rfcs/rfc1071.html) paper in the section "2. Calculating the Checksum" describes possible techniques to speed up this assignment. Also, clang can help to find [causes](https://llvm.org/docs/Vectorizers.html#diagnostics) of bad performance.


## Solution
The original `checksum` function computes a sum of `uint16_t` values, while also keeping track of amount of overflows.<br>
This extra conditional addition prevents CPU from using SIMD instructions.<br>
```c++
uint16_t acc = 0;
for (auto value : blob) {
    acc += value;
    acc += acc < value;  // add carry
}
```
This extra conditional addition `acc += acc < value` prevents the CPU from using SIMD instructions.<br/>
A better approach is to use a wider type like `uint32_t` and compute the sum without handling carry explicitly.<br/>
In this approach, the high 16 bits store the number of overflows (carry), and the low 16 bits store the raw sum.<br/>
The final result can be obtained by adding the upper and lower halves of the 32-bit value.

## Benchmark
This changes results in just over 10x speedup.
```bash
$ python ../../../tools/check_speedup.py -lab_path $(pwd) -num_runs 3 -v
-----------------------------------------------------
Benchmark           Time             CPU   Iterations
-----------------------------------------------------
bench1           1.76 us         1.75 us      1592778
solution: iteration 2 - Done
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.9299         -0.9299            28             2            28             2
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.9372         -0.9372            28             2            28             2
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.9377         -0.9377            28             2            28             2
```