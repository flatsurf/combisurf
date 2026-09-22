import importlib

import pytest


@pytest.mark.parametrize("module", ["combisurf.partial_sums", "combisurf.partial_sums_naive"])
@pytest.mark.parametrize("n, repeat", [(2, 20), (3, 20), (4, 20), (5, 20), (6, 20), (7, 20), (8, 20), (9, 20),
                                       (100, 10), (255, 100), (256, 100), (257, 100), (3333, 100)])
def test_partial_sums(module, n, repeat):
    # PartialSumsNaive and PartialSumsFenwick must agree with each other
    # under random operations, both as the pure Python oracle
    # (combisurf.partial_sums_naive) and as the Cython structures actually
    # used by the package (combisurf.partial_sums).
    from random import randrange
    mod = importlib.import_module(module)
    P0 = mod.PartialSumsNaive(n)
    P1 = mod.PartialSumsFenwick(n)
    for _ in range(repeat):
        i = randrange(n)
        x = randrange(-10, 10)
        P0.update(i, x)
        P1.update(i, x)
        for _ in range(10):
            start = randrange(n - 1)
            end = randrange(start + 1, n)
            s0 = P0.partial_sum(start, end)
            s1 = P1.partial_sum(start, end)
            assert s0 == s1, (module, n, start, end, s0, s1)


@pytest.mark.parametrize("cls_name", ["PartialSumsNaive", "PartialSumsFenwick"])
@pytest.mark.parametrize("n, repeat", [(1, 20), (2, 20), (5, 20), (9, 20), (100, 20), (255, 20), (256, 20), (257, 20), (1024, 20)])
def test_partial_sums_cython_matches_oracle(cls_name, n, repeat):
    # The Cython port must agree with the pure Python reference
    # implementation it was tested against while it was being written, not
    # merely with its sibling algorithm.
    from random import randrange
    from combisurf import partial_sums as fast
    from combisurf import partial_sums_naive as oracle
    Pf = getattr(fast, cls_name)(n)
    Po = getattr(oracle, cls_name)(n)
    for _ in range(repeat):
        i = randrange(n)
        x = randrange(-10, 10)
        Pf.update(i, x)
        Po.update(i, x)
        for _ in range(10):
            start = randrange(n - 1) if n > 1 else 0
            end = randrange(start + 1, n) if n > 1 else 1
            assert Pf.partial_sum(start, end) == Po.partial_sum(start, end), (cls_name, n, start, end)
    Pf.reset()
    Po.reset()
    assert Pf.partial_sum(0, n) == 0 == Po.partial_sum(0, n)


def test_partial_sums_factory():
    # PartialSums(n) picks PartialSumsNaive up to the threshold and
    # PartialSumsFenwick above it, and either way answers queries correctly.
    from combisurf.partial_sums import PartialSums, PartialSumsNaive, PartialSumsFenwick

    assert type(PartialSums(1)) is PartialSumsNaive
    assert type(PartialSums(256)) is PartialSumsNaive
    assert type(PartialSums(257)) is PartialSumsFenwick
    assert type(PartialSums(100000)) is PartialSumsFenwick

    for n in (1, 256, 257, 5000):
        P = PartialSums(n)
        P.update(0, 3)
        P.update(n - 1, 4)
        assert P.partial_sum(0, n) == 7
