#!/usr/bin/env python3
"""Exact checks for the additional weight-64 cases.

Run run_all.sh first, or supply the Gaussian logs using --results-dir.
No third-party packages are required. All arithmetic uses Python integers.
The multiplier and lifting reductions are proved in Section 5 of manuscript.tex;
this program checks only the explicitly stated finite computations.
"""
from __future__ import annotations

import argparse
import re
import sys
from itertools import product
from pathlib import Path

Gaussian = tuple[int, int]


def require(condition: bool, message: str) -> None:
    """Keep verification checks active even when Python is run with -O."""
    if not condition:
        raise ValueError(message)


def doubling_orbits(modulus: int) -> list[list[int]]:
    require(modulus > 0 and modulus % 2 == 1, "The modulus must be positive and odd.")
    covered: set[int] = set()
    orbits: list[list[int]] = []
    for start in range(modulus):
        if start in covered:
            continue
        orbit: list[int] = []
        j = start
        while j not in covered:
            covered.add(j)
            orbit.append(j)
            j = 2 * j % modulus
        require(j == start, "An orbit did not return to its starting point.")
        orbits.append(orbit)
    return orbits


def read_gaussian_log(path: Path, modulus: int, expected: list[int]) -> list[list[Gaussian]]:
    """Check every printed cumulative count and read any full solutions."""
    lines = path.read_text(encoding="utf-8").splitlines()
    count_lines = [line for line in lines if line.startswith("counts ")]
    require(len(count_lines) == 1, f"Expected one count line in {path.name}.")
    counts = [int(value) for value in count_lines[0].split()[1:]]
    require(len(counts) == modulus // 2 + 1, f"Incomplete shift coverage in {path.name}.")
    require(counts == expected, f"Unexpected filtering counts in {path.name}: {counts}")
    rows = [
        [(int(a), int(b)) for a, b in re.findall(r"\((-?\d+),(-?\d+)\)", line)]
        for line in lines if line.startswith("(")
    ]
    require(len(rows) == counts[-1], f"Missing solution vectors in {path.name}.")
    return rows


def check_340(rows: list[list[Gaussian]]) -> None:
    orbits = doubling_orbits(85)
    require([min(o) for o in orbits] == [0, 1, 3, 5, 7, 9, 13, 15, 17, 21, 29, 37],
            "Unexpected orbit representatives modulo 85.")
    require([len(o) for o in orbits] == [1, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8],
            "Unexpected orbit sizes modulo 85.")
    require(len(rows) == 16 and len({tuple(row) for row in rows}) == 16,
            "Expected sixteen distinct Gaussian solutions.")
    T = [int(j % 17 == 0 and j != 0) for j in range(85)]
    counts = [0] * 44  # all assignments; augmentation zero; shifts 1,...,42
    for row in rows:
        require(len(row) == len(orbits), "Wrong number of Gaussian coefficients.")
        require(all(b == 0 for a, b in row), "A Gaussian solution is not real.")
        require(row[0] == (0, 0) and row[8] == (2, 0), "Wrong singleton or size-four coefficient.")
        require(sum(abs(a) == 1 for a, b in row) == 6, "Expected six unit-coefficient orbits.")
        require(all(a in (-1, 0, 1) for k, (a, b) in enumerate(row) if k != 8),
                "An unexpected coefficient occurs away from the size-four orbit.")
        # Check the actual group-ring equation independently of the C++ logs.
        E = [0] * 85
        for orbit, (a, _) in zip(orbits, row):
            for j in orbit:
                E[j] = a
        require(sum(E) == 8 and sum(a * a for a in E) == 64, "Incorrect scalar equations for E.")
        require(all(sum(E[(j + h) % 85] * E[j] for j in range(85)) == 0
                    for h in range(1, 85)), "A claimed E fails a correlation equation.")
        support = [k for k, (a, _) in enumerate(row) if abs(a) == 1]
        for signs in product((-1, 1), repeat=6):
            counts[0] += 1
            U = [0] * 85
            for k, sign in zip(support, signs):
                for j in orbits[k]:
                    U[j] = sign
            if sum(U) != 0:
                continue
            counts[1] += 1
            require(sum(a * a for a in U) == 48, "Wrong support size for U.")
            for h in range(1, 43):
                if sum(U[(j + h) % 85] * U[j] for j in range(85)) != -12 * T[h]:
                    break
                counts[h + 1] += 1
    print("CW(340,64): 16 real Gaussian solutions; all full correlations checked independently.")
    print("U assignments; augmentation-zero assignments; successive correlations h=1..42:")
    print(counts)
    require(counts == [1024, 320, 40, 40] + [0] * 40, "Unexpected lifting counts for U.")
    print("No lift: all candidates violate U U^(-1) = 48 - 12T.")


def check_arithmetic_obstructions() -> None:
    for d in (14, 18):
        candidates = []
        for ar, ai, br, bi in product(range(-2, 3), repeat=4):
            if (ar + 4 * br - 8) % d or (ai + 4 * bi) % d:
                continue
            squared_sum = ar * ar + ai * ai + 4 * (br * br + bi * bi)
            if squared_sum > 64 or (64 - squared_sum) % d:
                continue
            real_sum = (8 - ar - 4 * br) // d
            imag_sum = (-ai - 4 * bi) // d
            remaining_squares = (64 - squared_sum) // d
            # Integer squares have the same parity as the underlying integers.
            if (remaining_squares - real_sum - imag_sum) % 2:
                continue
            candidates.append((ar, ai, br, bi))
        require(not candidates, f"An arithmetic candidate remains for d={d}.")
        print(f"d={d}: no compatible singleton/order-five-orbit coefficients.")


def enumerate_55() -> None:
    """Independent direct recursion; no suffix dynamic programming is used."""
    orbits = doubling_orbits(55)
    weights = [len(orbit) for orbit in orbits]
    require(weights == [1, 20, 20, 10, 4], "Unexpected orbit sizes modulo 55.")
    candidates: list[list[Gaussian]] = []

    def visit(k: int, squares: int, real_sum: int, imag_sum: int, prefix: list[Gaussian]) -> None:
        if k == len(weights):
            if squares == 64 and real_sum == 8 and imag_sum == 0:
                candidates.append(prefix)
            return
        w = weights[k]
        for a, b in product(range(-2, 3), repeat=2):
            if squares + w * (a * a + b * b) <= 64:
                visit(k + 1, squares + w * (a * a + b * b), real_sum + w * a,
                      imag_sum + w * b, prefix + [(a, b)])

    visit(0, 0, 0, 0, [])
    require(len(candidates) == 8, "Expected eight scalar-equation solutions modulo 55.")
    for row in candidates:
        z = [(0, 0)] * 55
        for orbit, value in zip(orbits, row):
            for j in orbit:
                z[j] = value
        real_corr = imag_corr = 0
        for j in range(55):
            a, b = z[(j + 1) % 55]
            c, d = z[j]
            real_corr += a * c + b * d
            imag_corr += b * c - a * d
        require((real_corr, imag_corr) in {(-6, 3), (-6, -3), (-8, 1), (-8, -1)},
                "Unexpected first correlation modulo 55.")
        print("m=55 candidate:", row, "Gamma_1 =", (real_corr, imag_corr))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=Path(__file__).resolve().parent / "results")
    args = parser.parse_args()
    try:
        expected = {
            55: [8] + [0] * 27,
            85: [737310, 7141, 7141, 128, 128, 24, 24, 20, 20] + [16] * 34,
            95: [0] * 48,
            215: [0] * 108,
        }
        rows85: list[list[Gaussian]] = []
        for m, counts in expected.items():
            rows = read_gaussian_log(args.results_dir / f"extension_gauss_{m}.txt", m, counts)
            if m == 85:
                rows85 = rows
        print("All four Gaussian log files have the expected complete filtering counts.")
        check_340(rows85)
        check_arithmetic_obstructions()
        enumerate_55()
        print("ALL EXTENSION CHECKS PASSED")
    except (OSError, ValueError) as error:
        print(f"Verification failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
