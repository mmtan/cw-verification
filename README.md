# Computational Verification for Circulant Weighing Matrices

This repository contains the exact verification programs accompanying:

**Ming Ming Tan, _Character and Multiplier Obstructions for Circulant Weighing Matrices_**,  
[arXiv:2608.18468](https://arxiv.org/abs/2608.18468).

The revised manuscript proves nonexistence results for circulant weighing matrices using contraction, character evaluation, multiplier methods, and exact finite computations.

The parameters treated in the revised manuscript include

$CW(105,36), CW(140,36), CW(116,49), CW(120,49), CW(192,49),$ and 
$CW(v,64)$ for $v\in{ 140,180,196,220,340,380,860}.$

The paper also proves the nonexistence of $CW(20p,64)$ whenever $p$ is an odd prime for which the multiplicative order of $2$ modulo $p$ is congruent to $2 \pmod 4$. 

## Files

### `verify_cw.cpp`

This is the principal exact verification program for the finite computations in the paper. It reproduces the relevant contraction, multiplier-orbit, correlation, cyclic-code, and lifting calculations used in the nonexistence proofs.

The program uses exact arithmetic. It does not rely on floating-point approximations, randomized searches, or heuristic acceptance criteria.

### `verify_extensions.py`

This program contains the additional calculations used for the extended weight-\(64\) results. These include the finite orbit and coefficient calculations for the additional orders treated in the revised manuscript, including the extra lifting restrictions required in the order-\(340\) case.

The program is intended to be read together with the corresponding extension results in the paper.


## Compilation

Compile the C++ program with:

```bash
g++ -O3 -std=c++17 verify_cw.cpp -o verify_cw
```

## Execution

Run the principal verification program with:

```bash
./verify_cw
```

Run the extension calculations with:

```bash
python3 verify_extensions.py
```

The output should reproduce the orbit data, filtering counts, terminal candidates, and final nonexistence certificates stated in the paper.

## Nature of the Computations

The programs implement exhaustive finite searches arising after the mathematical reductions in the paper.

The reductions themselves are proved in the manuscript. The programs verify the remaining finite statements, including, as applicable:

- classification of bounded contracted group-ring elements;
- multiplier-orbit decompositions;
- exact periodic-correlation equations;
- enumeration of ternary cyclic-code words;
- testing of Eisenstein-unit lifts;
- Gaussian-integer orbit calculations;
- coefficient restrictions required to lift a character image to a circulant weighing matrix.

All acceptance and rejection conditions are exact.

