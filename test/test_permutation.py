import pytest

def test_str_to_cycles():
    from combisurf.permutation import str_to_cycles

    assert str_to_cycles("(0,1)") == [[0, 1]]
    assert str_to_cycles("(0,1)(3,2)") == [[0, 1], [3, 2]]
    assert str_to_cycles("()(0,1)()(2,3)") == [[0, 1], [2, 3]]
    assert str_to_cycles("(0,1,2)(~0,~1,~2)") == [[0, 1, 2], [-1, -2, -3]]

    with pytest.raises(TypeError):
        str_to_cycles(2)


def test_order():
    from combisurf.permutation import perm_init, perm_order

    p = perm_init("()")
    assert perm_order(p) == 1

    p = perm_init("(1)")
    assert perm_order(p) == 1

    p = perm_init("(1,2)")
    assert perm_order(p) == 2

    p = perm_init("(1,2)(3,4)(6,7,8)")
    assert perm_order(p) == 6


def test_cycles():
    from combisurf.permutation import perm_random, perm_cycles, perm_are_in_same_orbit

    for cycle_length in [1, 2, 3, 5, 10, 50]:
        p = perm_random(10)
        cycles = perm_cycles(p)
        for ci in range(len(cycles)):
            for cj in range(len(cycles)):
                for i in cycles[ci]:
                    for j in cycles[cj]:
                        assert perm_are_in_same_orbit(p, i, j) == (ci == cj)
