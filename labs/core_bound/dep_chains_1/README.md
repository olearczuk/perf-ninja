# Speed up data dependency chains #1.

Critical data dependency chains are increasingly becoming the [only thing that matters](https://easyperf.net/blog/2022/05/11/Visualizing-Performance-Critical-Dependency-Chains) for performance of a general-purpose application. That is why it is very important to identify those and know possible ways to make them run faster. On a SW level, you can sometimes occasionally introduce an artificial data dependency, which should not exist in the first place. Those cases are usually easy to find. In a contrast, some data dependency chains are inherent to a particular type of data structure.

Ahhhh, good old linked lists... Traversing a linked list is essentially a looooooooong data dependency chain. To get the node `N+1` you need to retrieve the node `N` first. Even if we set aside the problem with memory locality, a dependency chain will not go away. The data dependency effectively serializes the execution making your ILP (Instruction-Level Parallelism) be very low.

The task in this lab assignment is to look up all the values from linked list A in linked list B. This is an O(N^2) algorithm and involves a lot of pointer chasing. Both linked lists use an arena allocator to place individual nodes right next to each other, which improves memory locality. To improve performance of the benchmark in this lab assignment even further you need to overlap the execution of multiple dependency chains.

The idea for the lab was proposed by @ibogosavljevic.

## Solution
The problem with the naive "nested loops" approach is that for each value in `List1`, we iterate over `List2` and see if any value matches.<br/>
This means we iterate over `List2` `size(List1)` times.<br/>
We could make this more efficient if we processed `List1` in chunks. For each chunk of size `N`, we would iterate over `List2` and compare its values against all values in the chunk.<br/>
This way we iterate over `List2` `size(List1) / N` times instead.<br/>
This is very substantial since iterating over lists is very inefficient due to data dependency.

## Benchmark
This change results in ~80% speedup for `N = 8`.
```bash
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.8327         -0.8327           104            17           104            17
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.7910         -0.7909            86            18            86            18
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.7997         -0.7997            86            17            86            17
```