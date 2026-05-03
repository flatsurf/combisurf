import pytest

def test_fg_word_is_reduced():
    from array import array
    from combisurf.word import fg_word_is_reduced

    assert fg_word_is_reduced(array('i', [0, 2]))
    assert fg_word_is_reduced(array('i', [2, 0]))
    assert fg_word_is_reduced(array('i', [0, 2, 1, 3]))
    assert fg_word_is_reduced(array('i', [2, 1, 3, 0]))
    assert fg_word_is_reduced(array('i', [0, 2, 1]))
    assert not fg_word_is_reduced(array('i', [0, 1]))
    assert not fg_word_is_reduced(array('i', [1, 0]))
    assert not fg_word_is_reduced(array('i', [2, 3]))
    assert not fg_word_is_reduced(array('i', [3, 2]))
    assert not fg_word_is_reduced(array('i', [0, 2, 3, 0]))
    assert not fg_word_is_reduced(array('i', [0, 3, 2, 0]))


def test_fg_word_reduce():
    from array import array
    from combisurf.word import fg_word_reduce

    for w in [[0, 2], [2, 0], [0, 2, 1, 3], [2, 1, 3, 0], [0, 2, 1]]:
        w = array('i', w)
        assert fg_word_reduce(w) == w

    assert fg_word_reduce(array('i', [1, 0])) == array('i', [])
    assert fg_word_reduce(array('i', [1, 2, 3, 0, 0])) == array('i', [0])


def test_fg_word_is_cyclically_reduced():
    from array import array
    from combisurf.word import fg_word_is_cyclically_reduced

    assert fg_word_is_cyclically_reduced(array('i', [0, 2]))
    assert fg_word_is_cyclically_reduced(array('i', [2, 0]))
    assert fg_word_is_cyclically_reduced(array('i', [0, 2, 1, 3]))
    assert fg_word_is_cyclically_reduced(array('i', [2, 1, 3, 0]))
    assert not fg_word_is_cyclically_reduced(array('i', [0, 2, 1]))
    assert not fg_word_is_cyclically_reduced(array('i', [0, 1]))
    assert not fg_word_is_cyclically_reduced(array('i', [1, 0]))
    assert not fg_word_is_cyclically_reduced(array('i', [2, 3]))
    assert not fg_word_is_cyclically_reduced(array('i', [3, 2]))
    assert not fg_word_is_cyclically_reduced(array('i', [0, 2, 3, 0]))
    assert not fg_word_is_cyclically_reduced(array('i', [0, 3, 2, 0]))


def test_fg_word_inverse():
    from array import array
    from combisurf.word import fg_word_inverse, fg_word_mul

    for w in [[], [0], [0, 2], [0, 3], [2, 0, 3]]:
        w = array('i', w)
        inv = fg_word_inverse(w)
        assert not fg_word_mul(w, inv)
        assert not fg_word_mul(inv, w)
