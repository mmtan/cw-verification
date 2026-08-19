# Computational Verification for Circulant Weighing Matrices

This repository contains the C++17 verification programs accompanying the paper

**New Nonexistence Results for Circulant Weighing Matrices**.

The programs reproduce the finite computations used in the nonexistence proofs in the paper.

## Files

* `verify_lifting.cpp`
  Reproduces the computations involving the ternary ([35,12]) code, the Eisenstein-unit lifting search, the local order-4 character calculation, and the Gaussian orbit enumeration used in the paper.

* `verify_weight49.cpp`
  Reproduces the orbit data and exact filtering computations used for the weight-49 nonexistence results.

## Compilation and Execution

The computations reported in the paper were reproduced using GCC 14.2.0.

Compile the programs with

```bash
g++ -O3 -std=c++17 verify_lifting.cpp -o verify_lifting
g++ -O3 -std=c++17 verify_weight49.cpp -o verify_weight49
```

and run them with

```bash
./verify_lifting
./verify_weight49
```

The programs print the exact enumeration and data used in the proofs.

## Computational Details

The programs use integer arithmetic only. No floating-point comparisons, randomized searches, or heuristic acceptance criteria are used.

Pruning is limited to exact feasibility tests for the remaining augmentation and energy.

For the Eisenstein search, `verify_lifting.cpp` checks shifts $(1,\ldots,17)$; shifts $(18,\ldots,34)$ are their complex conjugates. The reciprocal ternary code is covered by the involution described in the paper.

For the Gaussian cases, the correlation conditions used in the program are necessary conditions. They eliminate every vector in the bounded search space considered in the paper, which contains all character images arising from four-entry fibers.
