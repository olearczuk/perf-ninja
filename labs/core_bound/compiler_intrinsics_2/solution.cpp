#include "solution.hpp"

#include <immintrin.h>

#include <cstdint>
#include <iostream>
// Find the longest line in a file.
// Implementation uses ternary operator with a hope that compiler will
// turn it into a CMOV instruction.
// The code inside the inner loop is equivalent to:
/*
if (s == '\n') {
  longestLine = std::max(curLineLength, longestLine);
  curLineLength = 0;
} else {
  curLineLength++;
}*/
unsigned solutionOld(const std::string& inputContents) {
  unsigned longestLine = 0;
  unsigned curLineLength = 0;

  for (auto s : inputContents) {
    curLineLength = (s == '\n') ? 0 : curLineLength + 1;
    longestLine = std::max(curLineLength, longestLine);
  }

  return longestLine;
}

unsigned solutionAVX2(const std::string& inputContents) {
  const __m256i newline = _mm256_set1_epi8('\n');

  int previous_newline = 0;
  int max_length = 0;
  size_t i = 0;
  int len = inputContents.size();
  for (; i + 32 <= len; i += 32) {
    __m256i chunk = _mm256_loadu_si256((const __m256i*)&inputContents[i]);
    // Check for newline characters in the chunk
    __m256i cmp = _mm256_cmpeq_epi8(chunk, newline);
    // Keep in mind that bit N is set to 1 if the N-th byte in the chunk is
    // equal to '\n'.
    // This means ordering is "flipped" compared to the original
    // string. For example, if the first byte is '\n', the least significant bit
    // of the mask will be set.
    // That's why __builtin_ctz is used (trailing 0s) to find positions of \n.
    unsigned int mask = _mm256_movemask_epi8(cmp);

    // Pop all 1s from the mask and calculate line lengths
    while (mask != 0) {
      int pos = __builtin_ctz(mask) + i;
      max_length = std::max(max_length, pos - previous_newline);
      previous_newline = pos + 1;

      mask &= mask - 1;
    }
  }

  // Process remaining bytes (if any) without SIMD
  for (; i < len; i++) {
    if (inputContents[i] == '\n') {
      int length = i - previous_newline;
      previous_newline = i + 1;
      max_length = std::max(max_length, length);
    }
  }

  max_length = std::max(max_length, len - previous_newline);
  return static_cast<unsigned>(max_length);
}

unsigned solution(const std::string& inputContents) {
#ifdef SOLUTION
  return solutionAVX2(inputContents);
#else
  return solutionOld(inputContents);
#endif
}
