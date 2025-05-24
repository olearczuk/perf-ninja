This is a lab about [function inlining](https://en.wikipedia.org/wiki/Inline_expansion) to speed up sorting.

Function inlining is a transformation that replaces a call to a function `F` with the body for `F` specialized with the actual arguments of the call. Inlining is one of the most important compiler optimizations, not only because it eliminates the overhead of calling a function (prologue and epilogue), but also it enables other optimizations.

Whenever you find in a performance profile a function with hot prologue and epilogue, consider such function as one of the potential candidates for being inlined. In this lab assignment you will practice fixing such performance issues.

## Solution
First improvement would be to use `std::sort` instead of `qsort`.<br/>
Main advantage is that `std::sort` is implemented entirely in headers which allows the compiler to perform more optimizations.
This results in ~20% speedup.<br/>

However, there is one more problem.<br/>
When inspecting assembly code of `solution_free_function`, one can see following `call` instruction.<br/>
This means that `compare2` is passed via a function pointer, which makes it pretty hard for the compiler to inline this function.
```assembly
call  void std::__introsort_loop<S*, long, __gnu_cxx::__ops::_Iter_comp_iter<bool (*)(S const&, S const&)>>(S*, S*, long, __gnu_cxx::__ops::_Iter_comp_iter<bool (*)(S const&, S const&)>)
```

`solution_lambda` approach addresses this issue.<br>
Each lambda has its own distinct type which allows compiler to "look into" its code and perform more complex optimizations (like inlining).<br/>
When inspecting assembly code of `solution_lambda`, one can see following `call` instruction.<br/>
```assembly
call    void std::__introsort_loop<S*, long, __gnu_cxx::__ops::_Iter_comp_iter<solution_lambda(std::array<S, 10000ul>&)::$_0>>(S*, S*, long, __gnu_cxx::__ops::_Iter_comp_iter<solution_lambda(std::array<S, 10000ul>&)::$_0>)
```

## Benchmark
Lambda-based approach results in ~40% speedup.
```bash
$ python ../../../tools/check_speedup.py -lab_path $(pwd) -num_runs 3 -v
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.3867         -0.3869           622           381           621           381
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.3883         -0.3884           625           382           624           382
Benchmark                   Time             CPU      Time Old      Time New       CPU Old       CPU New
--------------------------------------------------------------------------------------------------------
bench1                   -0.3878         -0.3879           622           381           621           380
```