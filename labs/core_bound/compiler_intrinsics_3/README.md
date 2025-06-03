This is a third lab about using [compiler intrinsics](https://en.wikipedia.org/wiki/Intrinsic_function) to speed up parts of the code, where compilers fail to generate optimal code.

The task of this lab is to find the average position of a series of positions in 3d space:
```
struct Position {
    u32 x;
    u32 y;
    u32 z;
}
...
    u64 x = 0;
    u64 y = 0;
    u64 z = 0;

    for (Position pos: input) {
        x += pos.x; // extend 32 bit to 64 bit and accumulate
        y += pos.y;
        z += pos.z;
    }
```

Currently, compilers fail to vectorize it for x86. However, for ARM they generate clever code using their magic LDM instruction that deinterleaves the XYZ components. [^1] However, we can still vectorize it more efficiently.

To vectorize this loop manually, think about the memory layout and operations applied to the elements.

Caveat: current benchmark processes 16k positions, which takes 192KB of space (`12 bytes * 16K positions`). As the input size grows, the benchmarks becomes more memory bound which reduces the gains. 

**The lab is contributed by Jonathan Hallström (@JonathanHallstrom).**

[^1]: Here is an example of the ARM's LDM instruction.
    ```asm
    // Assume input points to an array of structures {u32 x, u32 y, u32 z}
    LD3 {V0.4S, V1.4S, V2.4S}, [input] // Load 4 positions at once (each containing xyz elements) 
    // After that:
    //  V0.4S will have four X coordinates,
    //  V1.4S will have four Y coordinates,
    //  V2.4S will have four Z coordinates,
    // Reducing vectors with such a friendly layout is much easier.
    ```


## Solution
Actually, at the time of writing solution compilers manage to auto vectorize the code: `clang -O3 -c solution.cpp  -Rpass-analysis=loop-vectorize` doesn't show any loop not being vectorized.<br/>
However, it's still possible to write faster code using intrinsics.<br/>
Each `Position { uint32_t x, y, z; };` requires 96 bits - 4 such values can be loaded into 3 128-bit registers.<br/>
All we need to do is to compute total sum for each register and then extract all `x`, `y` and `z` values from it.

## Benchmark
AVX2 based solution results in ~50% speedup compared to auto-vectorization.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.4944         -0.4945          5272          2665          5268          2663
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.4978         -0.4978          5295          2659          5290          2657
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.4961         -0.4961          5268          2655          5264          2652
```