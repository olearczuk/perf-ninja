#include "solution.h"

#include <emmintrin.h>
#include <immintrin.h>

#include <memory>

std::pair<int, int> mainLoopSSE2(const InputVector& input, uint8_t radius,
                                 OutputVector& output, int current_sum,
                                 int pos) {
  int limit = input.size() - radius;

  __m128i current_sum_vec = _mm_set1_epi16(current_sum);
  for (; pos + 7 < limit; pos += 8) {
    __m128i add = _mm_cvtepu8_epi16(_mm_loadu_si64(&input[pos + radius]));
    __m128i sub = _mm_cvtepu8_epi16(_mm_loadu_si64(&input[pos - radius - 1]));
    __m128i diff = _mm_sub_epi16(add, sub);

    // Prefix sum for 8 elements
    diff = _mm_add_epi16(diff, _mm_slli_si128(diff, 1 * sizeof(uint16_t)));
    diff = _mm_add_epi16(diff, _mm_slli_si128(diff, 2 * sizeof(uint16_t)));
    diff = _mm_add_epi16(diff, _mm_slli_si128(diff, 4 * sizeof(uint16_t)));

    // Current prefix sum
    current_sum_vec = _mm_add_epi16(current_sum_vec, diff);

    // Store the result in the output vector
    _mm_storeu_si128((__m128i*)(&output[pos]), current_sum_vec);

    // Update the current sum and reset current_sum_vec to the last value in
    // prefix sum
    current_sum = (uint16_t)_mm_extract_epi16(current_sum_vec, 7);
    current_sum_vec = _mm_set1_epi16(current_sum);
  }

  for (; pos < limit; ++pos) {
    current_sum -= input[pos - radius - 1];
    current_sum += input[pos + radius];
    output[pos] = current_sum;
  }

  return {pos, current_sum};
}

std::pair<int, int> mainLoopAVX2(const InputVector& input, uint8_t radius,
                                 OutputVector& output, int current_sum,
                                 int pos) {
  int limit = input.size() - radius;

  auto current_sum_vec = _mm256_set1_epi16(current_sum);
  for (; pos + 15 < limit; pos += 16) {
    __m256i add = _mm256_cvtepu8_epi16(
        _mm_loadu_si128((const __m128i*)&input[pos + radius]));
    __m256i sub = _mm256_cvtepu8_epi16(
        _mm_loadu_si128((const __m128i*)&input[pos - radius - 1]));
    __m256i diff = _mm256_sub_epi16(add, sub);

    // Prefix sum for 16 elements
    diff =
        _mm256_add_epi16(diff, _mm256_slli_si256(diff, 1 * sizeof(uint16_t)));
    diff =
        _mm256_add_epi16(diff, _mm256_slli_si256(diff, 2 * sizeof(uint16_t)));
    diff =
        _mm256_add_epi16(diff, _mm256_slli_si256(diff, 4 * sizeof(uint16_t)));
    // All the above operations happened independently on the two halves.
    // Now we need to add the second half to the first half.
    __m256i shift_8 = _mm256_set_m128i(
        _mm_set1_epi16(_mm256_extract_epi16(diff, 7)), _mm_setzero_si128());
    diff = _mm256_add_epi16(diff, shift_8);

    // Current prefix sum
    current_sum_vec = _mm256_add_epi16(diff, current_sum_vec);

    // Store the result in the output vector
    _mm256_storeu_si256((__m256i*)&output[pos], current_sum_vec);

    // Update the current sum and reset current_sum_vec to the last value in
    // prefix sum
    current_sum = (uint16_t)_mm256_extract_epi16(current_sum_vec, 15);
    current_sum_vec = _mm256_set1_epi16(current_sum);
  }

  for (; pos < limit; ++pos) {
    current_sum -= input[pos - radius - 1];
    current_sum += input[pos + radius];
    output[pos] = current_sum;
  }

  return {pos, current_sum};
}

void imageSmoothing(const InputVector& input, uint8_t radius,
                    OutputVector& output) {
  int pos = 0;
  int currentSum = 0;
  int size = static_cast<int>(input.size());

  // 1. left border - time spend in this loop can be ignored, no need to
  // optimize it
  for (int i = 0; i < std::min<int>(size, radius); ++i) {
    currentSum += input[i];
  }

  int limit = std::min(radius + 1, size - radius);
  for (pos = 0; pos < limit; ++pos) {
    currentSum += input[pos + radius];
    output[pos] = currentSum;
  }

#ifdef SOLUTION
  std::tie(pos, currentSum) =
      mainLoopAVX2(input, radius, output, currentSum, pos);
#else
  // 2. main loop.
  limit = size - radius;
  for (; pos < limit; ++pos) {
    currentSum -= input[pos - radius - 1];
    currentSum += input[pos + radius];
    output[pos] = currentSum;
  }
#endif

  // 3. special case, executed only if size <= 2*radius + 1
  limit = std::min(radius + 1, size);
  for (; pos < limit; pos++) {
    output[pos] = currentSum;
  }

  // 4. right border - time spend in this loop can be ignored, no need to
  // optimize it
  for (; pos < size; ++pos) {
    currentSum -= input[pos - radius - 1];
    output[pos] = currentSum;
  }
}