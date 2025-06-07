This lab calculates an image of [Mandelbrot set](https://en.wikipedia.org/wiki/Mandelbrot_set). In this algorithm, you can process each pixel independently, however, the number of iterations for processing each pixel can be different. The issue you will face when you process two pixels in a SIMD fashion: what to do when the processing loop for one pixel finishes after 100 iterations while for the adjacent pixel it runs for 200 iterations?

The `-ffast-math` option is disabled for validation purposes. It doesn't help compiler to autovectorize the code.

Lab assignment developed by Oleg Makovski (@0legmak).

## Solution
As stated in the in the description above, the main problem is that different pixels require different amount of iterations until hardcoded limit is reached.<br/>
However, we can get around this problem by using a bit mask to keep track of pixels in SIMD chunk that have already exceeded the bound.

## Benchmark
AVX2-based solution results in 70% speedup.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.7017         -0.7018           769           229           769           229
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.7026         -0.7021           769           229           767           229
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.7003         -0.7003           764           229           763           229
```