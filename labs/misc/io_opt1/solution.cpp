#include "solution.hpp"

#include <fstream>
#include <stdexcept>

#include "MappedFile.hpp"

uint32_t computeCrcNew(uint32_t crc, const char *file_name) {
  auto mapped_file = MappedFile{file_name};

  // Update the CRC32 value character by character
  for (char c : mapped_file.getContents()) {
    update_crc32(crc, static_cast<uint8_t>(c));
  }

  return crc;
}

uint32_t computeCrc(uint32_t crc, const char *file_name) {
  std::fstream file_stream{file_name};
  if (!file_stream.is_open())
    throw std::runtime_error{"The file could not be opened"};

  // Update the CRC32 value character by character
  char c;
  while (true) {
    file_stream.read(&c, 1);
    if (file_stream.eof()) break;
    update_crc32(crc, static_cast<uint8_t>(c));
  }

  return crc;
}

uint32_t solution(const char *file_name) {
  // Initial value has all bits set to 1
  uint32_t crc = 0xff'ff'ff'ff;
#ifdef SOLUTION
  crc = computeCrcNew(crc, file_name);
#else
  crc = computeCrc(crc, file_name);
#endif
  // Invert the bits
  return crc ^ 0xff'ff'ff'ff;
}
