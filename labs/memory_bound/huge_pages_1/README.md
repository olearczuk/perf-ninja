# Finite element operator evaluation

## Motivation and problem description
The finite element method is a discretization technique used mainly for physics simulations.
The physical domain (usually 2D or 3D) where the physical problem is defined is broken up into a mesh of small, geometrically simple subdomains called elements (much like a surface is broken down into triangles in computer graphics).
The continuous fields representing physical quantities (e.g. density, pressure, temperature) are approximated using a finite number of values, which represent the values of the fields at the vertices (nodes) of the elements.
For this reason, these values are usually referred to as nodal values.
The partial differential equations describing the physics to be simulated are expressed as systems of linear algebraic equations, with the nodal values being the unknowns.
These algebraic systems are usually quite large, and therefore must be solved using (super-)computers.

Since the resulting system of algebraic equations is usually solved using iterative methods (e.g. the [conjugate gradient method](https://en.wikipedia.org/wiki/Conjugate_gradient_method)), it is not necessary to explicitly store the matrix describing it in memory.
It is sufficient to be able to evaluate the action of this matrix on a vector (compute the matrix-vector product).
This approach is usually referred to as matrix-free.

In this lab, you are meant to optimize a piece of code which evaluates the aforementioned matrix-vector product for a structural problem involving a truss.
The truss consists of many bars, each of which is represented using a single finite element.
To evaluate the matrix-free operator, we iterate over all elements, compute the local 4x4 matrix, gather the local 4-element left-hand side vector, perform the multiplication, and then scatter the result into the global right-hand side vector.
This example is rather simple (at the level of a 2nd year mechanical engineering course), but it illustrates the gather-scatter memory access pattern, which is ubiquitous in various simulation codes across the globe.

**Note**: The above description of the finite element method is extremely simplified, and is only meant to give a high level overview of the problem.

## Hint
Note that the compiler generates fairly optimal code for the floating-point computations, so there are little gains to be had by optimizing those.
Instead, focus on the memory access pattern, which is extremely random, meaning many distant memory addresses are accessed in rapid succession, which puts stress on the TLB.
This situation could be alleviated by allocating the memory which is accessed in a random fashion on huge pages (see [HugePagesSetupTips](HugePagesSetupTips.md)).
For the convenience of your solution, all such allocations are done using the `allocateDoublesArray` function.
In fact, this is the only place of the code which is modified in the suggested solution.

Authored-by: @kubagalecki

## Solution
### Code
Program allocates a lot of memory upfront using the function below and then does processing (which is not relevant to this lab).<br/>
```c++
inline auto allocateDoublesArray(size_t size) {
  return std::make_unique<double[]>(size);
}
```
Since the allocated memory is large, this puts quite a lot of pressure on Translation Lookaside Buffer.<br/>  
Using larger pages should significantly reduce the pressure.
```c++
inline auto allocateDoublesArray(size_t size) {
#define MAP_HUGE_2MB (21 << MAP_HUGE_SHIFT)
  // Allocate memory
  constexpr size_t huge_page_size = 2 * 1024 * 1024;  // 2MB huge pages
  size_t alloc_size = size * sizeof(double);
  // Round up to nearest huge page size
  alloc_size = (alloc_size + huge_page_size - 1) & ~(huge_page_size - 1);

  void* ptr =
      mmap(nullptr, alloc_size, PROT_READ | PROT_WRITE,
           MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_HUGE_2MB, -1, 0);
  if (ptr == MAP_FAILED) {
    throw std::runtime_error("Huge page allocation failed");
  }
  double* alloc = static_cast<double*>(ptr);

  // Use munmap to deallocate the memory
  auto deleter = [alloc_size](double* ptr) { munmap(ptr, alloc_size); };
  return std::unique_ptr<double[], decltype(deleter)>(alloc,
                                                      std::move(deleter));
}
```
### Reduce dTLB pressure
Let's compare normal vs Huge Pages solution in terms of dTLB pressure.<br/>
I had to explicitly turn Transparent Huge Pages off to showcase the advantage of using Huge Pages (otherwise, the default code would use it out of the box). Without it being turned off, the results are basically the same.<br/>>
```bash
echo never | sudo tee /sys/kernel/mm/transparent_hugepage/enabled
```
#### Normal 
```bash
$ cmake -E make_directory build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && cmake --build . --config Release --parallel $(nproc)
$ perf stat -e dTLB-loads,dTLB-load-misses,dTLB-store-misses ./lab
 Performance counter stats for './lab':

     3,061,260,865      cpu_atom/dTLB-loads/                                                    (0.02%)
     2,792,200,306      cpu_core/dTLB-loads/                                                    (99.98%)
         3,788,753      cpu_atom/dTLB-load-misses/       #    0.12% of all dTLB cache accesses  (0.02%)
       413,293,653      cpu_core/dTLB-load-misses/       #   14.80% of all dTLB cache accesses  (99.98%)
           328,776      cpu_atom/dTLB-store-misses/                                             (0.02%)
        66,286,685      cpu_core/dTLB-store-misses/                                             (99.98%)
```
### Huge Pages
As you can see, amount of misses is drastically smaller when using Huge Pages.
```bash
$ cmake -E make_directory build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-DSOLUTION .. && cmake --build . --config Release --parallel $(nproc)
$ perf stat -e dTLB-loads,dTLB-load-misses,dTLB-store-misses ./lab
 Performance counter stats for './lab':

     1,101,539,169      cpu_atom/dTLB-loads/                                                    (0.23%)
     2,549,148,429      cpu_core/dTLB-loads/                                                    (99.77%)
           361,590      cpu_atom/dTLB-load-misses/       #    0.03% of all dTLB cache accesses  (0.23%)
        49,193,675      cpu_core/dTLB-load-misses/       #    1.93% of all dTLB cache accesses  (99.77%)
             6,995      cpu_atom/dTLB-store-misses/                                             (0.23%)
         9,655,599      cpu_core/dTLB-store-misses/                                             (99.77%)
```

## Benchmark
Huge Pages result in ~20% speedup
```bash
Benchmark                                    Time             CPU      Time Old      Time New       CPU Old       CPU New
-------------------------------------------------------------------------------------------------------------------------
Apply matrix-free operator                -0.1616         -0.1617          4057          3402          4050          3395
Benchmark                                    Time             CPU      Time Old      Time New       CPU Old       CPU New
-------------------------------------------------------------------------------------------------------------------------
Apply matrix-free operator                -0.2354         -0.2351          4407          3369          4398          3364
Benchmark                                    Time             CPU      Time Old      Time New       CPU Old       CPU New
-------------------------------------------------------------------------------------------------------------------------
Apply matrix-free operator                -0.2002         -0.1984          4353          3481          4335          3475
```