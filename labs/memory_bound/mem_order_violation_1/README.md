This lab features an oldie-but-goodie Otsu's thresholding[^1] to convert a grayscale[^2] image to a binary[^3] image. One of the key parts in this algorithm is calculating a histogram of a grayscale image, i.e., calculate how many times a certain color appears in the image. Since the input image is 8-bit grayscale, there are only 256 different colors.

```cpp
std::array<uint32_t, 256> hist;
hist.fill(0);
for (int i = 0; i < image.width * image.height; ++i)
  hist[image.data[i]]++;
```

The implementation of the histogram algortihm is very simple but it has one nasty property. For each pixel on the image, you need to 1) read the current value of the corresponding color of the pixel, 2) increment it and 3) store it back.

When updates of the same color in the histogram occur at relatively high rates, the processor may not have completed updating pixel `i` prior to beginning pixel `i+1`. In such cases, a processor predicts whether the value loaded for the `i+1` update will come from memory or from the `i`'s store. If from memory, the two updates can be performed in parallel, otherwise the processor must serialize the updates.

Simple example: if you have the following pixels in the image:
```
0xFF 0xFF 0xFF 0xFF 0xFF 0xFF ...
```
Then all updates to `hist[0xFF]` will be serialized.

Think about how you can workaround this problem. Hint: you can use aditional memory.

Bonus exercise1: what would be the worst and the best cases for the original implementation and your solution?

Input images were taken from here: https://people.sc.fsu.edu/~jburkardt/data/pgmb/pgmb.html

[^1]: https://en.wikipedia.org/wiki/Otsu%27s_method
[^2]: https://en.wikipedia.org/wiki/Grayscale
[^3]: https://en.wikipedia.org/wiki/Binary_image

## Solution
As described in problem statement, the main problem is memory contention for subsequent pixels with the same value.<br/>
This can be solved by duplicating the array so that there are multiple slots in memory for each pixel value.<br/>
By rotating over the indices, we make sure memory location for subsequent pixels is never the same.<br/>
This approach introduces an overhead in scenarios of no memory contention but offers a meaningful speedup otherwise.
```c++
constexpr int UNROLL_FACTOR = 4;
std::array<std::array<uint32_t, UNROLL_FACTOR>, 256> hists;
for (int i = 0; i < hists.size(); ++i) {
  hists[i].fill(0);
}

int total_pixels = image.width * image.height;
// main "unrolling" loop
int i = 0;
for (; i < total_pixels - UNROLL_FACTOR; i += UNROLL_FACTOR) {
  for (int j = 0; j < UNROLL_FACTOR; ++j) {
    ++hists[image.data[i + j]][j];
  }
}

// finish the rest
for (; i < total_pixels; ++i) {
  ++hists[image.data[i]][0];
}

std::array<uint32_t, 256> hist;
for (int i = 0; i < hist.size(); ++i) {
  hist[i] = std::accumulate(hists[i].cbegin(), hists[i].cend(), 0);
}
return hist;
```

## Benchmark
This change results in varying speedups depending on the input file.<br/>
The optimization becomes less effective as the frequency of repeated pixel values decreases.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bird/0                   -0.2054         -0.2054            46            37            46            37
coins/1                  -0.2452         -0.2451            23            17            23            17
pepper/2                 -0.1176         -0.1178            16            14            16            14
pixabay/3                -0.0349         -0.0348             8             8             8             8
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bird/0                   -0.1553         -0.1559            46            39            46            39
coins/1                  -0.2334         -0.2326            25            19            25            19
pepper/2                 -0.1396         -0.1393            20            17            20            17
pixabay/3                -0.0964         -0.0965             9             8             9             8
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bird/0                   -0.2265         -0.2266            52            40            52            40
coins/1                  -0.2623         -0.2623            26            19            26            19
pepper/2                 -0.0314         -0.0298            18            17            18            17
pixabay/3                -0.0499         -0.0536             8             8             8             8
```