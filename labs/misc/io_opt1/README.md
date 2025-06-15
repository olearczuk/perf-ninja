# Compute CRC32 checksum of a file

A [CRC](https://en.wikipedia.org/wiki/Cyclic_redundancy_check) (cyclic redundancy check) is often used to ensure that data transmitted over a network was not corrupted en route.
This operation involves iterating through the transmitted file and updating the result in 32-bit increments.
Since computing the CRC of a 32-bit chunk is very cheap[1], this problem is bound by the speed at which we can supply these 32-bit chunks to the CPU.

We don't recommend you trying to speed up the CRC computation using compiler intrinsics, focus on the IO instead. This lab is designed to teach you how to efficiently read data from a file.
Solution ideas include:

- Reduce the overall file reading overhead by processing large chunks of data. 
- Map the contents of a file into the address space (e.g. `mmap` on Linux). Take a look at `MappedFile.hpp`.

Authored-by: @kubagalecki

---

[1]: Hardware instructions dedicated to this purpose exist.
In the case of x86, the `crc32` instruction is part of the SSE 4.2 extension, and is accessible via a compiler intrinsic.
Its cost is comparable to a single signed integer multiplication.

## Overview
```c++
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
```
Intuitively, it looks like this program is very heavily IO-bound because. <br/>
`strace` evidence supports this suspicion.
```bash
$ strace -c ./lab
% time     seconds  usecs/call     calls    errors syscall
------ ----------- ----------- --------- --------- ----------------
 49.51    0.007336           0     70495           read
 26.25    0.003889           0     13628         2 openat
 16.93    0.002508           0     13626           close
  2.30    0.000341         341         1           execve
  2.03    0.000301          13        23           mmap
  1.68    0.000249           3        74           write
  0.40    0.000060          10         6           mprotect
  0.21    0.000031           3         8           fstat
  0.14    0.000021           1        18           clock_gettime
  0.10    0.000015          15         1           munmap
...
```

## Solution
`file_stream.read(&c, 1);` is the problematic part. Reading file byte by byte is very inefficient.<br/>
We can map file's content to process memory using `mmap` and then process it as a local memory (e.g. `stream_view`).

## Benchmark
mmap-based solution is 75% - 80% faster, depending on file size.
```bash
Benchmark                       Time             CPU      Time Old      Time New       CPU Old       CPU New
------------------------------------------------------------------------------------------------------------
Small file/0                 -0.7506         -0.7538            43            11            43            11
Medium file/1                -0.8126         -0.8125            21             4            21             4
Large file/2                 -0.8134         -0.8132             3             0             3             0
Benchmark                       Time             CPU      Time Old      Time New       CPU Old       CPU New
------------------------------------------------------------------------------------------------------------
Small file/0                 -0.7474         -0.7509            42            11            42            11
Medium file/1                -0.8128         -0.8126            21             4            21             4
Large file/2                 -0.8152         -0.8149             3             0             3             0
Benchmark                       Time             CPU      Time Old      Time New       CPU Old       CPU New
------------------------------------------------------------------------------------------------------------
Small file/0                 -0.7495         -0.7528            43            11            43            11
Medium file/1                -0.8131         -0.8130            21             4            21             4
Large file/2                 -0.8145         -0.8144             3             0             3             0
```