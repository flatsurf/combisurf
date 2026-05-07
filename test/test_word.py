import pytest

def test_word_is_reduced():
    from array import array
    from combisurf.word import word_is_reduced

    assert word_is_reduced(array('i', [0, 2]))
    assert word_is_reduced(array('i', [2, 0]))
    assert word_is_reduced(array('i', [0, 2, 1, 3]))
    assert word_is_reduced(array('i', [2, 1, 3, 0]))
    assert word_is_reduced(array('i', [0, 2, 1]))
    assert not word_is_reduced(array('i', [0, 1]))
    assert not word_is_reduced(array('i', [1, 0]))
    assert not word_is_reduced(array('i', [2, 3]))
    assert not word_is_reduced(array('i', [3, 2]))
    assert not word_is_reduced(array('i', [0, 2, 3, 0]))
    assert not word_is_reduced(array('i', [0, 3, 2, 0]))


def test_word_reduce():
    from array import array
    from combisurf.word import word_reduce

    for w in [[0, 2], [2, 0], [0, 2, 1, 3], [2, 1, 3, 0], [0, 2, 1]]:
        w = array('i', w)
        assert word_reduce(w) == w

    assert word_reduce(array('i', [1, 0])) == array('i', [])
    assert word_reduce(array('i', [1, 2, 3, 0, 0])) == array('i', [0])


def test_word_is_cyclically_reduced():
    from array import array
    from combisurf.word import word_is_cyclically_reduced

    assert word_is_cyclically_reduced(array('i', [0, 2]))
    assert word_is_cyclically_reduced(array('i', [2, 0]))
    assert word_is_cyclically_reduced(array('i', [0, 2, 1, 3]))
    assert word_is_cyclically_reduced(array('i', [2, 1, 3, 0]))
    assert not word_is_cyclically_reduced(array('i', [0, 2, 1]))
    assert not word_is_cyclically_reduced(array('i', [0, 1]))
    assert not word_is_cyclically_reduced(array('i', [1, 0]))
    assert not word_is_cyclically_reduced(array('i', [2, 3]))
    assert not word_is_cyclically_reduced(array('i', [3, 2]))
    assert not word_is_cyclically_reduced(array('i', [0, 2, 3, 0]))
    assert not word_is_cyclically_reduced(array('i', [0, 3, 2, 0]))


def test_word_free_group_inverse():
    from array import array
    from combisurf.word import word_free_group_inverse, word_free_group_mul

    for w in [[], [0], [0, 2], [0, 3], [2, 0, 3]]:
        w = array('i', w)
        inv = word_free_group_inverse(w)
        assert not word_free_group_mul(w, inv)
        assert not word_free_group_mul(inv, w)
