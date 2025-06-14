This is the good old game of life. The program takes 10 randomly generated 1024x1024 grids (boards) and simulates the next 10 game of life rounds. It first simulates 10 rounds for the first board, then goes off to the second board, and so on.

As written, the program experiences many branch mispredictions. Your job is to find where they happen and replace them with predicate instructions. On x86 you should see `cmov` instructions (hint: use `__builtin_unpredictable`, make sure you have Clang-17 or later version installed). On ARM you should see `csel` (conditional select) instructions. It is a good idea to experiment with the code in godbolt.org before you start modifying the actual code.

Bonus question: when you fix the main source of mispredictions, there are still many branches left in the hot loop nest. How can you get rid of them?

## Overview
As you can see in the provided code, the calculating logic is fairly simple, however, there are 2 main sources of branching:
- an `if` statement checking whether a potential neighbor is out of bounds
- `switch` statement to compute the future value
```c++
void simulateNextOld() {
    // printCurrentGrid();
    int M = current.size();
    int N = current[0].size();

    // Loop through every cell
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
        int aliveNeighbours = 0;
        // finding the number of neighbours that are alive
        for (int p = -1; p <= 1; p++) {    // row-offet (-1,0,1)
            for (int q = -1; q <= 1; q++) {  // col-offset (-1,0,1)
            if ((i + p < 0) ||      // if row offset less than UPPER boundary
                (i + p > M - 1) ||  // if row offset more than LOWER boundary
                (j + q < 0) ||      // if column offset less than LEFT boundary
                (j + q > N - 1))    // if column offset more than RIGHT boundary
                continue;
            aliveNeighbours += current[i + p][j + q];
            }
        }
        // The cell needs to be subtracted from
        // its neighbours as it was counted before
        aliveNeighbours -= current[i][j];

        // Implementing the Rules of Life:
        switch (aliveNeighbours) {
            // 1. Cell is lonely and dies
            case 0:
            case 1:
            future[i][j] = 0;
            break;
            // 2. Remains the same
            case 2:
            future[i][j] = current[i][j];
            break;
            // 3. A new cell is born
            case 3:
            future[i][j] = 1;
            break;
            // 4. Cell dies due to over population
            default:
            future[i][j] = 0;
        }
        }
    }
    std::swap(current, future);
}
```

As you can see there are quite a few branch mispredictions (and a lot of branching in general).
```bash
$ perf stat -e branches,branch-misses ./lab
 Performance counter stats for './lab':
     1,830,450,172      cpu_core/branches/
       116,375,207      cpu_core/branch-misses/          #    6.36% of all branches
```

## Solution
Let's remove switch statement first. `aliveNeighbours` value is pretty random and very hard to predict for branch predictor.<br/>
Intuitively, it means that using `cmov` instead of branching might be a good idea.
```c++
// __builtin_unpredictable is an extra suggestion for compiler to use cmov instructions.
future[i][j] = __builtin_unpredictable(aliveNeighbours == 2)
                    ? current[i][j]
                : __builtin_unpredictable(aliveNeighbours == 3) ? 1
                                                                : 0;
```

If statements in `aliveNeighbours` computation can be eliminated by adding extra padding for rows and columns.<br/>
By surrounding the grid with a border of zero-initialized cells, we eliminate the need for bounds checks in the neighbor computation.

As you can see, new code has significantly less branching overall and very small branch miss rate.
```bash
$ perf stat -e branches,branch-misses ./lab
 Performance counter stats for './lab':
       467,637,479      cpu_core/branches/
           455,814      cpu_core/branch-misses/          #    0.10% of all branches
```

## Benchmark
Removing both sources of branching results in ~95% speedup
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.9416         -0.9416          1088            64          1086            63
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.9503         -0.9503          1188            59          1185            59
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.9475         -0.9475          1194            63          1192            63
```