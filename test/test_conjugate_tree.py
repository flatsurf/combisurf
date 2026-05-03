import itertools
import pytest

def small_binary_lyndon_words():
    return ((0,), (1,),
            (0,1), (0,0,1), (0,1,1),
            (0,0,0,1), (0,0,1,1), (0,1,1,1),
            (0,0,0,0,1), (0,0,0,1,1), (0,0,1,0,1), (0,0,1,1,1), (0,1,0,1,1), (0,1,1,1,1),
            (0,0,0,0,0,1), (0,0,0,0,1,1), (0,0,0,1,0,1), (0,0,0,1,1,1), (0,0,1,0,1,1),
            (0,0,1,1,0,1), (0,0,1,1,1,1), (0,1,0,1,1,1), (0,1,1,1,1,1))


def small_ternary_lyndon_words():
    return ((0,), (1,), (2,),
            (0,1), (0,2), (1,2),
            (0,0,1), (0,0,2), (0,1,1), (0,1,2), (0,2,1), (0,2,2), (1,1,2), (1,2,2),
            (0,0,0,1), (0,0,0,2), (0,0,1,1), (0,0,1,2), (0,0,2,1), (0,0,2,2), (0,1,0,2),
            (0,1,1,1), (0,1,1,2), (0,1,2,1), (0,1,2,2), (0,2,1,1), (0,2,1,2), (0,2,2,1),
            (0,2,2,2), (1,1,1,2), (1,1,2,2), (1,2,2,2))


def test_process():
    from combisurf.conjugate_tree import ConjugateTree

    T = ConjugateTree()

    # a third power of a primitive word
    assert T.process([0,0,1,0,0,1,0,0,1]) == 3
    # identical root
    assert T.process([0,0,1,0,0,1]) == 0

    # a primitive word
    assert T.process([0,0,1,0]) == 1
    # identical root
    assert T.process([1,0,0,0,1,0,0,0]) == -1


def test_leaf_as_conjugate():
    from combisurf.conjugate_tree import ConjugateTree

    for W in [small_binary_lyndon_words(), small_ternary_lyndon_words()]:
        for k in range(1, 4):
            for words in itertools.combinations(W, k):
                for swords in itertools.permutations(words):
                    T = ConjugateTree()
                    for word in swords:
                        assert T.process(word) == 1
                    leaves = T.leaves()
                    assert len(leaves) == sum(map(len, swords))
                    leaves_by_words = []
                    c = 0
                    for w in swords:
                        leaves_by_words.append(leaves[c:c+len(w)])
                        c += len(w)
                    for i, w_leaves in enumerate(leaves_by_words):
                        for k, leaf in enumerate(w_leaves):
                            assert T.leaf_as_conjugate(leaf) == (i, k)
                            assert T._leaf_shift(leaf) == w_leaves[(k + 1) % len(w_leaves)]
