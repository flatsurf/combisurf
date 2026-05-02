import pytest


@pytest.mark.parametrize("n, repeat", [(2, 20), (3, 20), (4, 20), (5, 20), (6, 20), (7, 20), (8, 20), (9, 20),
                                       (100, 10), (127, 100), (128, 100), (129, 100), (3333, 100)])
def test_partial_sums(n, repeat):
    from random import randrange
    from combisurf.partial_sums import PartialSumsNaive, PartialSumsBinarySplitting
    P0 = PartialSumsNaive(n)
    P1 = PartialSumsBinarySplitting(n)
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
            assert s0 == s1, (s0, s1, P0._values, P1._values)



