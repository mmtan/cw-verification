# Verification code for circulant weighing matrices

This repository contains the exact finite computations accompanying the paper

**Ming Ming Tan, _Character and Multiplier Obstructions for Circulant Weighing Matrices_.**

The mathematical reductions are proved in the paper. The programs verify the remaining finite enumerations and correlation calculations. 
## Files

- `verify_cw.cpp` — exhaustive C++17 verification of the contracted order-35 cases, the ternary-code and Eisenstein-unit calculation for `CW(105,36)`, the Gaussian orbit calculations for weights 64, and the weight-49 orbit calculations.
- `verify_extensions.py` — short independent checks for the additional weight-64 arguments: the eight candidates for `CW(220,64)`, the arithmetic obstructions for orders 380 and 860, and the lifting obstruction for `CW(340,64)`.

## Run

Compile and run the C++ program:

```bash
g++ -O3 -std=c++17 -Wall -Wextra verify_cw.cpp -o verify_cw
./verify_cw
```

Run the Python checks:

```bash
python3 verify_extensions.py
```

A successful run ends with

```text
ALL C++ VERIFICATION CHECKS PASSED
ALL PYTHON VERIFICATION CHECKS PASSED
```

The programs use exact integer arithmetic only. They do not use floating-point tolerances, randomized searches, or external input files.

## Scope

The programs verify the finite statements explicitly used in the paper, including:

- the two contracted matrices of order 35 and weight 36;
- the ternary cyclic-code distribution and all normalized Eisenstein-unit lifts for `CW(105,36)`;
- the Gaussian multiplier-orbit calculations for moduli 35, 45, 49, 55, 85, 95, and 215;
- the order-340 lifting calculation;
- the multiplier-orbit calculations for `CW(116,49)`, `CW(120,49)`, and `CW(192,49)`.

For the hypotheses, notation, and proofs reducing the weighing-matrix questions to these finite computations, refer to the paper.
