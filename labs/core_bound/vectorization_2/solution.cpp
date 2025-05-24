#include "solution.hpp"

uint16_t checksum_old(const Blob &blob) {
  uint16_t acc = 0;
  for (auto value : blob) {
    acc += value;
    acc += acc < value;  // add carry
  }
  return acc;
}

uint16_t checksum_new(const Blob &blob) {
  uint32_t acc = 0;
  for (auto value : blob) {
    acc += value;
  }

  // high 4 bytes store number of overflows
  // low 5 bytes store sum without taking extra carry into account
  const auto add_high_low_bytes = [](uint32_t v) {
    auto high = v >> 16;
    auto old = v & 0xffff;
    return high + old;
  };

  acc = add_high_low_bytes(acc);
  // compute again if previous call also resulted in an overflow
  acc = add_high_low_bytes(acc);
  return static_cast<uint16_t>(acc);
}

uint16_t checksum(const Blob &blob) {
#ifdef SOLUTION
  return checksum_new(blob);
#else
  return checksum_old(blob);
#endif
}