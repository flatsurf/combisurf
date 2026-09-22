r"""
Measure the crossover between ``PartialSumsNaive`` and ``PartialSumsFenwick``.

This is the measurement behind ``_PARTIAL_SUMS_NAIVE_THRESHOLD`` in
``combisurf/partial_sums.pyx``. Each step does what one step of the
startpoint term of a geometric intersection does: two ``update`` at random
positions and one ``partial_sum(0, i)``. The two classes are timed in the
interleaved order A B B A, and each size keeps the minimum over the rounds,
because a non-interleaved run drifts by several percent (CPU frequency),
which is the size of the difference being measured.

Usage::

    python extra/partial_sum_threshold.py [n ...]

It prints the time in ns per step for each ``n``; the threshold is the
largest ``n`` at which the naive class is still at least as fast.
"""
import random
import sys
import time

from combisurf.partial_sums import PartialSumsNaive, PartialSumsFenwick

STEPS = 20000
ROUNDS = 7


def ops(n, seed):
    rng = random.Random(seed)
    return [(rng.randrange(n), rng.randrange(n), rng.randrange(n + 1))
            for _ in range(STEPS)]


def run(cls, n, steps):
    P = cls(n)
    update = P.update
    partial_sum = P.partial_sum
    t0 = time.perf_counter()
    for i, j, k in steps:
        update(i, 1)
        update(j, -1)
        partial_sum(0, k)
    return (time.perf_counter() - t0) / len(steps) * 1e9


def measure(n):
    best = {PartialSumsNaive: float('inf'), PartialSumsFenwick: float('inf')}
    for r in range(ROUNDS):
        steps = ops(n, r)
        for cls in (PartialSumsNaive, PartialSumsFenwick,
                    PartialSumsFenwick, PartialSumsNaive):
            best[cls] = min(best[cls], run(cls, n, steps))
    return best[PartialSumsNaive], best[PartialSumsFenwick]


if __name__ == '__main__':
    sizes = [int(a) for a in sys.argv[1:]] or [64, 128, 256, 512, 1024]
    print(f"{'n':>6} {'naive':>8} {'Fenwick':>8}   (ns per step)")
    for n in sizes:
        a, b = measure(n)
        print(f"{n:>6} {a:8.0f} {b:8.0f}")
