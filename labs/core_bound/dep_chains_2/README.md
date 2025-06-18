This small program simulates the random particle movement. We have 1000 particles moving on a 2D surface without constraints. It means that there are no bounds, they can move as far from their initial coordinates as they want. Each particle has initial x and y coordinates in the range [-1000,1000] and a constant speed in the range [0;1]. The program simulates 1000 steps of particle movement.

To validate the simulation (final positions of the particles), we use a deterministic random number generator (fake) that uses a global state and thus always generates the same sequence of numbers.

There is one very nasty performance problem that doesn't allow us to run the simulation fast. Can you find a dependency chain in the code and fix it?

Note: your solution is allowed to be not functionally equivalent to the original program if validation still passes. For example, if an RNG will generate a different sequence of random numbers than before - that's OK. Users would not be able to tell the difference since the motion of particles is random anyway.


## Overview
```c++
// Medium-quality random number generator.
// https://www.javamex.com/tutorials/random_numbers/xorshift.shtml
struct XorShift32 {
  uint32_t val;
  XorShift32(uint32_t seed) : val(seed) {}
  XorShift32() = default;

 public:
  uint32_t gen() {
    val ^= (val << 13);
    val ^= (val >> 17);
    val ^= (val << 5);
    return val;
  }
};

// Simulate the motion of the particles.
// For every particle, we generate a random angle and move the particle
// in the corresponding direction.
template <class RNG>
void randomParticleMotion(std::vector<Particle> &particles, uint32_t seed) {
  RNG rng(seed);
  for (int i = 0; i < STEPS; i++)
    for (auto &p : particles) {
      uint32_t angle = rng.gen();
      float angle_rad = angle * DEGREE_TO_RADIAN;
      p.x += cosine(angle_rad) * p.velocity;
      p.y += sine(angle_rad) * p.velocity;
    }
}
```
In the above implementation, each call to `XorShift32::gen()` depends on the previous one, as it mutates internal state (`val`).<br/>
This creates a dependency chain. No other dependency chains exist in the function.

## Solution
To eliminate the `RNG` dependency chain, multiple independent `RNG` instances can be used in parallel.<br/>
This introduces slight variability in results due to different `RNG` sequences but is acceptable given the stochastic nature of the simulation (the validate test still passes).
```c++
// Simulate the motion of the particles.
// For every particle, we generate a random angle and move the particle
// in the corresponding direction.
template <class RNG>
void randomParticleMotion(std::vector<Particle> &particles, uint32_t seed) {
  constexpr int RNG_COUNT = 8;
  assert(particles.size() % RNG_COUNT == 0 &&
         "The number of particles must be a multiple of RNG_COUNT");
  std::array<RNG, RNG_COUNT> rngs;
  for (int i = 0; i < RNG_COUNT; ++i) {
    rngs[i] = RNG(seed + i);
  }
  std::array<uint32_t, RNG_COUNT> angles;
  for (int step = 0; step < STEPS; step++) {
    int i = 0;
    for (; i < particles.size(); i += RNG_COUNT) {
      for (int j = 0; j < RNG_COUNT; ++j) {
        uint32_t angle = rngs[j].gen();
        float angle_rad = angle * DEGREE_TO_RADIAN;
        particles[i + j].x += cosine(angle_rad) * particles[i + j].velocity;
        particles[i + j].y += sine(angle_rad) * particles[i + j].velocity;
      }
    }
  }
}
```
Note: This requires `RNG` to be default-constructible or assignable into the `std::array`.


## Benchmark
This change produces varying results on the test machine
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.0408         -0.0398            13            13            13            13
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   +0.0019         +0.0019            13            13            13            13
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.0813         -0.0803            14            13            14            13
```
However, it results in [~40% speedup on Macbook M1](https://github.com/dendibakh/perf-ninja/actions/runs/15739331004/job/44360394284#step:3:119) and [~12% speedup on Linux with Alder Lake CPU](https://github.com/dendibakh/perf-ninja/actions/runs/15739331012/job/44360394338#step:3:132)
