#include "solution.hpp"

#include <array>

// use short to decrease memory footprint
const std::array<short, 101> lookup_table {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,                 // 0-12: bucket 0
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,        // 13-28: bucket 1
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,                    // 29-40: bucket 2
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,                    // 41-52: bucket 3
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 53-70: bucket 4
      5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,                    // 71-82: bucket 5
      6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,     // 83-99: bucket 6
      static_cast<short>(DEFAULT_BUCKET)                     // 100: default
};

static std::size_t mapToBucket(std::size_t v) {
#ifdef SOLUTION
  auto bucket = lookup_table[std::min(v, lookup_table.size() - 1)];
  return static_cast<std::size_t>(bucket);
#else
  //   size of a bucket
  if (v < 13)
    return 0;  //   13
  else if (v < 29)
    return 1;  //   16
  else if (v < 41)
    return 2;  //   12
  else if (v < 53)
    return 3;  //   12
  else if (v < 71)
    return 4;  //   18
  else if (v < 83)
    return 5;  //   12
  else if (v < 100)
    return 6;  //   17
  return DEFAULT_BUCKET;
#endif
}

std::array<std::size_t, NUM_BUCKETS> histogram(const std::vector<int> &values) {
  std::array<std::size_t, NUM_BUCKETS> retBuckets{0};
  for (auto v : values) {
    retBuckets[mapToBucket(v)]++;
  }
  return retBuckets;
}
