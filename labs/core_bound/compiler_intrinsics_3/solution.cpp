#include "solution.hpp"

#include <immintrin.h>

#include <algorithm>
#include <cassert>

Position<std::uint32_t> solutionOld(
    std::vector<Position<std::uint32_t>> const &input) {
  std::uint64_t x = 0;
  std::uint64_t y = 0;
  std::uint64_t z = 0;

  for (auto pos : input) {
    x += pos.x;
    y += pos.y;
    z += pos.z;
  }

  return {
      static_cast<std::uint32_t>(x / std::max<std::uint64_t>(1, input.size())),
      static_cast<std::uint32_t>(y / std::max<std::uint64_t>(1, input.size())),
      static_cast<std::uint32_t>(z / std::max<std::uint64_t>(1, input.size())),
  };
}

Position<std::uint32_t> solution(
    std::vector<Position<std::uint32_t>> const &input) {
#ifdef SOLUTION
  __m256i sum_x_vec = _mm256_setzero_si256();
  __m256i sum_y_vec = _mm256_setzero_si256();
  __m256i sum_z_vec = _mm256_setzero_si256();

  size_t i = 0;
  __m256i xyzx_sum = _mm256_setzero_si256();
  __m256i yzxy_sum = _mm256_setzero_si256();
  __m256i zxyz_sum = _mm256_setzero_si256();

  for (; i + 3 < input.size(); i += 4) {
    uint32_t *p = (uint32_t *)&input[i];
    __m256i xyzx = _mm256_cvtepu32_epi64(_mm_loadu_si128((__m128i *)&p[0]));
    __m256i yzxy = _mm256_cvtepu32_epi64(_mm_loadu_si128((__m128i *)&p[4]));
    __m256i zxyz = _mm256_cvtepu32_epi64(_mm_loadu_si128((__m128i *)&p[8]));

    xyzx_sum = _mm256_add_epi64(xyzx_sum, xyzx);
    yzxy_sum = _mm256_add_epi64(yzxy_sum, yzxy);
    zxyz_sum = _mm256_add_epi64(zxyz_sum, zxyz);
  }

  std::uint64_t x = 0, y = 0, z = 0;
  x += _mm256_extract_epi64(xyzx_sum, 0);
  y += _mm256_extract_epi64(xyzx_sum, 1);
  z += _mm256_extract_epi64(xyzx_sum, 2);
  x += _mm256_extract_epi64(xyzx_sum, 3);

  y += _mm256_extract_epi64(yzxy_sum, 0);
  z += _mm256_extract_epi64(yzxy_sum, 1);
  x += _mm256_extract_epi64(yzxy_sum, 2);
  y += _mm256_extract_epi64(yzxy_sum, 3);

  z += _mm256_extract_epi64(zxyz_sum, 0);
  x += _mm256_extract_epi64(zxyz_sum, 1);
  y += _mm256_extract_epi64(zxyz_sum, 2);
  z += _mm256_extract_epi64(zxyz_sum, 3);

  // Process the rest
  for (; i < input.size(); ++i) {
    x += input[i].x;
    y += input[i].y;
    z += input[i].z;
  }

  return {
      static_cast<std::uint32_t>(x / std::max<std::uint64_t>(1, input.size())),
      static_cast<std::uint32_t>(y / std::max<std::uint64_t>(1, input.size())),
      static_cast<std::uint32_t>(z / std::max<std::uint64_t>(1, input.size())),
  };

#else
  return solutionOld(input);
#endif
}