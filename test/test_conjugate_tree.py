import itertools
import random
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


def test_hard_check():
    # exercises the three possible outcomes of process() with hard_check=True:
    # a new primitive word, a new non-primitive word, and a word that is a
    # conjugate (possibly of a power) of an already registered word.
    from combisurf.conjugate_tree import ConjugateTree

    T = ConjugateTree()
    assert T.process([0], hard_check=True) == 1
    assert T.process([0, 1, 0, 0, 1], hard_check=True) == 1
    assert T.process([1, 0, 1, 0], hard_check=True) == 2
    assert T.process([0, 1], hard_check=True) == -2
    T._check()


def test_hard_check_random():
    from combisurf.conjugate_tree import ConjugateTree

    rng = random.Random(0)
    outcomes = {"primitive": 0, "non_primitive": 0, "conjugate": 0}
    for _ in range(50):
        T = ConjugateTree()
        alphabet = rng.randint(2, 5)
        for _ in range(20):
            length = rng.randint(1, 8)
            base = [rng.randrange(alphabet) for _ in range(length)]
            w = base * rng.choice([1, 1, 1, 2, 3])
            ans = T.process(w, hard_check=True)
            if ans > 1:
                outcomes["non_primitive"] += 1
            elif ans == 1:
                outcomes["primitive"] += 1
            else:
                outcomes["conjugate"] += 1
        T._check()

    # all three outcomes of process() must be exercised, otherwise the
    # bijection check at the end of process() and the structural check
    # during process() would not both be covered
    assert all(outcomes.values())


def cyclic_order_key(conjugate, angles, depth):
    r"""
    Return the sequence that the cyclic order at infinity compares: the angle
    of the first letter, then at each further step the angle from the reverse
    of the previous letter to the current one.

    This is an independent reimplementation of the ordering that
    :meth:`~combisurf.conjugate_tree.ConjugateTree.cyclically_sorted_leaves`
    realizes through the tree.
    """
    n = len(angles)
    l = len(conjugate)
    ans = [angles[conjugate[0]]]
    for d in range(1, depth):
        ans.append((angles[conjugate[d % l]] - angles[conjugate[(d - 1) % l] ^ 1]) % n)
    return ans


def check_cyclically_sorted_leaves(T, angles):
    words = T.words()
    depth = 2 * sum(len(w) for w in words) + 4
    expected = sorted(((i, k) for i, w in enumerate(words) for k in range(len(w))),
                      key=lambda ik: cyclic_order_key(list(words[ik[0]][ik[1]:]) + list(words[ik[0]][:ik[1]]),
                                                      angles, depth))
    leaves = T.cyclically_sorted_leaves(angles)
    assert sorted(leaves) == T.leaves()
    assert [T.leaf_as_conjugate(s) for s in leaves] == expected


def test_cyclically_sorted_leaves():
    from combisurf.conjugate_tree import ConjugateTree

    T = ConjugateTree()
    assert T.process([0, 1, 1]) == 1
    assert T.process([0, 1]) == 1
    assert T.cyclically_sorted_leaves([0, 1]) == [6, 1, 8, 4, 2]
    check_cyclically_sorted_leaves(T, [0, 1])

    # the torus and the octagon, with the angle tables that
    # GeometricIntersection builds on them
    T = ConjugateTree()
    for w in [[0, 0, 2, 2], [0, 2, 0, 0, 3], [1, 3, 1, 2]]:
        T.process(w)
    check_cyclically_sorted_leaves(T, [0, 2, 1, 3])

    T = ConjugateTree()
    for w in [[0, 2, 2, 5, 2, 2, 5], [0, 3, 6], [1, 4, 7, 0]]:
        T.process(w)
    check_cyclically_sorted_leaves(T, [0, 2, 4, 6, 1, 3, 5, 7])


def test_cyclically_sorted_leaves_random():
    import random
    from combisurf.conjugate_tree import ConjugateTree

    rng = random.Random(20260922)
    for _ in range(60):
        n = 2 * rng.randint(1, 4)
        angles = list(range(n))
        rng.shuffle(angles)
        T = ConjugateTree()
        for _ in range(rng.randint(1, 5)):
            T.process([rng.randrange(n) for _ in range(rng.randint(1, 9))])
        if not T.words():
            continue
        check_cyclically_sorted_leaves(T, angles)
