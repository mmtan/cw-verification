
from __future__ import annotations

from itertools import product


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def doubling_orbits(modulus: int) -> list[list[int]]:
    covered: set[int] = set()
    orbits: list[list[int]] = []
    for start in range(modulus):
        if start in covered:
            continue
        orbit: list[int] = []
        x = start
        while x not in covered:
            covered.add(x)
            orbit.append(x)
            x = 2 * x % modulus
        require(x == start, "Invalid multiplication-by-2 orbit.")
        orbits.append(orbit)
    return orbits


# The sixteen normalized Gaussian solutions for m=85 listed in Appendix C.
M85_ROWS = [
    [0, -1, -1, 0, 1, 0, -1, 1, 2, 1, 0, 0],
    [0, -1, -1, 1, 1, 0, 1, 0, 2, 0, 0, -1],
    [0, -1, 0, 0, -1, -1, 1, 1, 2, 0, 0, 1],
    [0, -1, 0, 0, 0, 1, 0, 1, 2, -1, -1, 1],
    [0, -1, 0, 1, 1, 1, 0, 0, 2, 0, -1, -1],
    [0, -1, 1, 0, -1, 0, 1, 1, 2, -1, 0, 0],
    [0, 0, -1, 0, 0, -1, -1, 1, 2, 1, 1, 0],
    [0, 0, -1, 1, -1, 0, 1, 0, 2, -1, 1, 0],
    [0, 0, -1, 1, 0, -1, 0, 0, 2, 1, -1, 1],
    [0, 0, -1, 1, 0, 1, 0, 0, 2, -1, 1, -1],
    [0, 0, 1, 0, 0, 1, -1, 1, 2, -1, -1, 0],
    [0, 0, 1, 1, -1, 0, -1, 0, 2, 1, -1, 0],
    [0, 1, 0, 0, 0, -1, 0, 1, 2, -1, 1, -1],
    [0, 1, 0, 0, 1, -1, -1, 1, 2, 0, 0, -1],
    [0, 1, 0, 1, -1, -1, 0, 0, 2, 0, -1, 1],
    [0, 1, 1, 1, -1, 0, -1, 0, 2, 0, 0, -1],
]


def expand_orbit_values(modulus: int, values: list[int]) -> list[int]:
    orbits = doubling_orbits(modulus)
    require(len(values) == len(orbits), "Wrong number of orbit coefficients.")
    sequence = [0] * modulus
    for orbit, value in zip(orbits, values):
        for j in orbit:
            sequence[j] = value
    return sequence


def correlation(sequence: list[int], shift: int) -> int:
    n = len(sequence)
    return sum(sequence[(j + shift) % n] * sequence[j] for j in range(n))


def check_order_340() -> None:
    orbits = doubling_orbits(85)
    require([min(o) for o in orbits] == [0, 1, 3, 5, 7, 9, 13, 15, 17, 21, 29, 37],
            "Unexpected orbit representatives modulo 85.")
    require([len(o) for o in orbits] == [1, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8],
            "Unexpected orbit sizes modulo 85.")
    require(len(M85_ROWS) == 16 and len({tuple(row) for row in M85_ROWS}) == 16,
            "The Appendix C list must contain sixteen distinct rows.")

    counts = [0, 0, 0, 0, 0]  # all choices; augmentation; shifts 1, 2, 3
    for row in M85_ROWS:
        require(row[0] == 0 and row[8] == 2, "Unexpected special-orbit coefficients.")
        require(sum(abs(x) == 1 for x in row) == 6, "Expected six unit-coefficient orbits.")
        E = expand_orbit_values(85, row)
        require(sum(E) == 8 and sum(x * x for x in E) == 64,
                "A listed Gaussian image fails the scalar equations.")
        require(all(correlation(E, h) == 0 for h in range(1, 85)),
                "A listed Gaussian image fails a correlation equation.")

        support = [k for k, value in enumerate(row) if abs(value) == 1]
        for signs in product((-1, 1), repeat=6):
            counts[0] += 1
            U = [0] * 85
            for orbit_index, sign in zip(support, signs):
                for j in orbits[orbit_index]:
                    U[j] = sign
            if sum(U) != 0:
                continue
            counts[1] += 1
            require(sum(x * x for x in U) == 48, "Unexpected support size for U.")
            ok = True
            for h in (1, 2, 3):
                if correlation(U, h) != 0:
                    ok = False
                    break
                counts[h + 1] += 1
            if ok:
                raise ValueError("An order-340 lift survived the first three required correlations.")

    require(counts == [1024, 320, 40, 40, 0],
            f"Unexpected order-340 filtering counts: {counts}")
    print("CW(340,64): 16 Gaussian images checked; lifting counts 1024, 320, 40, 40, 0.")


def check_order_220() -> None:
    orbits = doubling_orbits(55)
    weights = [len(o) for o in orbits]
    require([min(o) for o in orbits] == [0, 1, 3, 5, 11],
            "Unexpected orbit representatives modulo 55.")
    require(weights == [1, 20, 20, 10, 4], "Unexpected orbit sizes modulo 55.")

    candidates: list[list[tuple[int, int]]] = []

    def visit(k: int, energy: int, real_sum: int, imag_sum: int,
              prefix: list[tuple[int, int]]) -> None:
        if k == len(weights):
            if energy == 64 and real_sum == 8 and imag_sum == 0:
                candidates.append(prefix)
            return
        w = weights[k]
        for a, b in product(range(-2, 3), repeat=2):
            new_energy = energy + w * (a * a + b * b)
            if new_energy <= 64:
                visit(k + 1, new_energy, real_sum + w * a,
                      imag_sum + w * b, prefix + [(a, b)])

    visit(0, 0, 0, 0, [])
    require(len(candidates) == 8, f"Expected 8 candidates modulo 55, found {len(candidates)}.")

    first_correlations: list[tuple[int, int]] = []
    for row in candidates:
        z = [(0, 0)] * 55
        for orbit, value in zip(orbits, row):
            for j in orbit:
                z[j] = value
        real = imag = 0
        for j in range(55):
            a, b = z[(j + 1) % 55]
            c, d = z[j]
            real += a * c + b * d
            imag += b * c - a * d
        first_correlations.append((real, imag))

    allowed = {(-6, 3), (-6, -3), (-8, 1), (-8, -1)}
    require(all(value in allowed for value in first_correlations),
            "An unexpected first correlation occurred modulo 55.")
    require(all(value != (0, 0) for value in first_correlations),
            "A candidate modulo 55 satisfies the first correlation.")
    print("CW(220,64): 8 scalar-equation candidates checked; none has zero first correlation.")


def arithmetic_candidates(d: int, use_parity: bool) -> list[tuple[int, int, int, int]]:
    candidates: list[tuple[int, int, int, int]] = []
    for ar, ai, br, bi in product(range(-2, 3), repeat=4):
        if (ar + 4 * br - 8) % d != 0 or (ai + 4 * bi) % d != 0:
            continue
        used = ar * ar + ai * ai + 4 * (br * br + bi * bi)
        if used > 64 or (64 - used) % d != 0:
            continue
        if use_parity:
            real_sum = (8 - ar - 4 * br) // d
            imag_sum = (-ai - 4 * bi) // d
            remaining_squares = (64 - used) // d
            if (remaining_squares - real_sum - imag_sum) % 2 != 0:
                continue
        candidates.append((ar, ai, br, bi))
    return candidates


def check_orders_380_and_860() -> None:
    require(not arithmetic_candidates(18, use_parity=False),
            "An arithmetic candidate remains for CW(380,64).")
    require(not arithmetic_candidates(14, use_parity=True),
            "An arithmetic candidate remains for CW(860,64).")
    print("CW(380,64) and CW(860,64): arithmetic coefficient checks passed.")


def main() -> None:
    check_order_220()
    check_orders_380_and_860()
    check_order_340()
    print("ALL PYTHON VERIFICATION CHECKS PASSED")


if __name__ == "__main__":
    main()
