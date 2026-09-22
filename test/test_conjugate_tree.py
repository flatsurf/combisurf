import itertools
import random
import pytest

# the conjugate tree comes in two implementations and the fast one in two
# transition representations; every test below runs against all of them
#
# - "naive"   the pure Python reference, combisurf.conjugate_tree_naive
# - "dense"   the Cython one with an alphabet-indexed transition table
# - "sparse"  the Cython one walking a linked list of siblings
# - "auto"    the Cython one picking between the two on the alphabet
# - "unknown" the Cython one that was not told the alphabet
KINDS = ["naive", "dense", "sparse", "auto", "unknown"]


def make_tree(kind, alphabet, reserve=0):
    from combisurf.conjugate_tree import ConjugateTree
    from combisurf.conjugate_tree_naive import ConjugateTreeNaive

    if kind == "naive":
        return ConjugateTreeNaive()
    if kind == "unknown":
        return ConjugateTree(reserve=reserve)
    if kind == "auto":
        return ConjugateTree(alphabet, reserve=reserve)
    return ConjugateTree(alphabet, reserve=reserve, algorithm=kind)


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


@pytest.mark.parametrize("kind", KINDS)
def test_constructor(kind):
    T = make_tree(kind, 4)
    assert T.num_states() == 1
    assert T.num_words() == 0
    assert T.words() == []
    assert T.leaves() == []
    assert T.internal_states() == []
    assert T.transitions(0) == {}


def test_alphabet_is_checked():
    from combisurf.conjugate_tree import ConjugateTree

    T = ConjugateTree(3)
    assert T.alphabet() == 3
    with pytest.raises(ValueError):
        T.process([0, 1, 3])
    # the rejected word left nothing behind
    assert T.num_words() == 0
    assert T.process([0, 1, 2]) == 1

    with pytest.raises(ValueError):
        ConjugateTree(-1)
    with pytest.raises(ValueError):
        ConjugateTree(4, reserve=-1)
    with pytest.raises(ValueError):
        ConjugateTree(algorithm='dense')
    with pytest.raises(ValueError):
        ConjugateTree(4, algorithm='triangular')


def test_algorithm_dispatch():
    from combisurf.conjugate_tree import ConjugateTree

    # the threshold is on the alphabet alone, so that the choice survives not
    # knowing how long the words will be
    assert ConjugateTree(4).algorithm() == 'dense'
    assert ConjugateTree(32).algorithm() == 'dense'
    assert ConjugateTree(33).algorithm() == 'sparse'
    assert ConjugateTree(256).algorithm() == 'sparse'
    assert ConjugateTree().algorithm() == 'sparse'
    assert ConjugateTree(256, algorithm='dense').algorithm() == 'dense'
    assert ConjugateTree(4, algorithm='sparse').algorithm() == 'sparse'


@pytest.mark.parametrize("kind", KINDS)
def test_process(kind):
    T = make_tree(kind, 2)

    # a third power of a primitive word
    assert T.process([0,0,1,0,0,1,0,0,1]) == 3
    # identical root
    assert T.process([0,0,1,0,0,1]) == 0

    # a primitive word
    assert T.process([0,0,1,0]) == 1
    # identical root
    assert T.process([1,0,0,0,1,0,0,0]) == -1


@pytest.mark.parametrize("kind", KINDS)
def test_process_errors(kind):
    T = make_tree(kind, 2)
    with pytest.raises(ValueError):
        T.process([])
    with pytest.raises(ValueError):
        T.process([0, -1])
    assert T.num_words() == 0


@pytest.mark.parametrize("kind", KINDS)
@pytest.mark.parametrize("reserve", [0, 1, 1000])
def test_leaf_as_conjugate(kind, reserve):
    for W, alphabet in [(small_binary_lyndon_words(), 2), (small_ternary_lyndon_words(), 3)]:
        for k in range(1, 4):
            for words in itertools.combinations(W, k):
                for swords in itertools.permutations(words):
                    T = make_tree(kind, alphabet, reserve=reserve)
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


@pytest.mark.parametrize("kind", KINDS)
def test_hard_check(kind):
    # exercises the three possible outcomes of process() with hard_check=True:
    # a new primitive word, a new non-primitive word, and a word that is a
    # conjugate (possibly of a power) of an already registered word.
    T = make_tree(kind, 2)
    assert T.process([0], hard_check=True) == 1
    assert T.process([0, 1, 0, 0, 1], hard_check=True) == 1
    assert T.process([1, 0, 1, 0], hard_check=True) == 2
    assert T.process([0, 1], hard_check=True) == -2
    T._check()


@pytest.mark.parametrize("kind", KINDS)
def test_hard_check_random(kind):
    rng = random.Random(0)
    outcomes = {"primitive": 0, "non_primitive": 0, "conjugate": 0}
    for _ in range(50):
        alphabet = rng.randint(2, 5)
        T = make_tree(kind, alphabet)
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


@pytest.mark.parametrize("kind", KINDS)
def test_cyclically_sorted_leaves(kind):
    T = make_tree(kind, 2)
    assert T.process([0, 1, 1]) == 1
    assert T.process([0, 1]) == 1
    assert T.cyclically_sorted_leaves([0, 1]) == [6, 1, 8, 4, 2]
    check_cyclically_sorted_leaves(T, [0, 1])

    # the torus and the octagon, with the angle tables that
    # GeometricIntersection builds on them
    T = make_tree(kind, 4)
    for w in [[0, 0, 2, 2], [0, 2, 0, 0, 3], [1, 3, 1, 2]]:
        T.process(w)
    check_cyclically_sorted_leaves(T, [0, 2, 1, 3])

    T = make_tree(kind, 8)
    for w in [[0, 2, 2, 5, 2, 2, 5], [0, 3, 6], [1, 4, 7, 0]]:
        T.process(w)
    check_cyclically_sorted_leaves(T, [0, 2, 4, 6, 1, 3, 5, 7])


@pytest.mark.parametrize("kind", KINDS)
def test_cyclically_sorted_leaves_random(kind):
    rng = random.Random(20260922)
    for _ in range(60):
        n = 2 * rng.randint(1, 4)
        angles = list(range(n))
        rng.shuffle(angles)
        T = make_tree(kind, n)
        for _ in range(rng.randint(1, 5)):
            T.process([rng.randrange(n) for _ in range(rng.randint(1, 9))])
        if not T.words():
            continue
        check_cyclically_sorted_leaves(T, angles)



@pytest.mark.parametrize("kind", ["dense", "sparse", "auto", "unknown"])
def test_cyclically_sorted_leaf_arcs_random(kind):
    # the leaves of cyclically_sorted_leaf_arcs are the ones of
    # cyclically_sorted_leaves, described through leaf_as_conjugate
    rng = random.Random(20260923)
    for _ in range(60):
        n = 2 * rng.randint(1, 4)
        angles = list(range(n))
        rng.shuffle(angles)
        T = make_tree(kind, n)
        for _ in range(rng.randint(1, 5)):
            T.process([rng.randrange(n) for _ in range(rng.randint(1, 9))])
        words = T.words()
        expected = []
        for s in T.cyclically_sorted_leaves(angles):
            i, k = T.leaf_as_conjugate(s)
            w = words[i]
            expected.append((i, w[k], (angles[w[k - 1] ^ 1] - angles[w[k]]) % n - 1))
        word_index, firsts, turns = T.cyclically_sorted_leaf_arcs(angles)
        assert all(a.typecode == 'q' for a in (word_index, firsts, turns))
        assert list(zip(word_index, firsts, turns)) == expected, (words, angles)

def assert_same_tree(T0, T1):
    r"""
    Check that two conjugate trees are the same down to the numbering of
    their nodes.
    """
    assert T0.num_states() == T1.num_states()
    assert T0.num_words() == T1.num_words()
    assert T0.words() == T1.words()
    assert T0.size() == T1.size()
    assert T0.leaves() == T1.leaves()
    assert T0.internal_states() == T1.internal_states()
    for s in range(T0.num_states()):
        assert T0.transitions(s) == T1.transitions(s), s
    for s in T0.leaves():
        assert T0.leaf_as_conjugate(s) == T1.leaf_as_conjugate(s), s
        assert T0._leaf_shift(s) == T1._leaf_shift(s), s
    for s in T0.internal_states():
        assert T0.internal_state_word(s) == T1.internal_state_word(s), s


# 4 and 8 are below the dense/sparse threshold, 34 and 64 above it
@pytest.mark.parametrize("alphabet", [2, 4, 8, 34, 64])
def test_against_naive(alphabet):
    r"""
    Run the Cython conjugate tree and the pure Python one side by side.
    """
    from combisurf.conjugate_tree import ConjugateTree
    from combisurf.conjugate_tree_naive import ConjugateTreeNaive

    rng = random.Random(1000 + alphabet)
    for trial in range(40):
        reserve = rng.choice([0, 1, 5, 500])
        algorithm = rng.choice([None, 'dense', 'sparse'])
        T0 = ConjugateTree(alphabet, reserve=reserve, algorithm=algorithm)
        T1 = ConjugateTreeNaive()
        history = []
        for _ in range(rng.randint(1, 8)):
            base = [rng.randrange(alphabet) for _ in range(rng.randint(1, 12))]
            w = base * rng.choice([1, 1, 1, 2, 3])
            history.append(w)
            hard = (trial % 8 == 0)
            assert T0.process(list(w), hard_check=hard) == T1.process(list(w), hard_check=hard), history
            assert_same_tree(T0, T1)
        T0._check()
        T1._check()
        if T0.num_words():
            angles = list(range(alphabet))
            rng.shuffle(angles)
            assert T0.cyclically_sorted_leaves(angles) == T1.cyclically_sorted_leaves(angles), (history, angles)


def test_against_naive_long_words():
    r"""
    The same, on words long enough that the node arrays are reallocated many
    times over.
    """
    from combisurf.conjugate_tree import ConjugateTree
    from combisurf.conjugate_tree_naive import ConjugateTreeNaive

    rng = random.Random(4242)
    for alphabet, length in [(2, 4000), (4, 2000), (64, 2000)]:
        T0 = ConjugateTree(alphabet)
        T1 = ConjugateTreeNaive()
        for _ in range(3):
            w = [rng.randrange(alphabet) for _ in range(length)]
            assert T0.process(list(w)) == T1.process(list(w))
        assert_same_tree(T0, T1)
        angles = list(range(alphabet))
        rng.shuffle(angles)
        assert T0.cyclically_sorted_leaves(angles) == T1.cyclically_sorted_leaves(angles)


def test_reserve_is_only_a_hint():
    r"""
    A tree given the exact number of nodes it needs never reallocates, and a
    tree given a wrong hint still answers the same.
    """
    from combisurf.conjugate_tree import ConjugateTree

    rng = random.Random(7)
    words = [[rng.randrange(6) for _ in range(30)] for _ in range(5)]
    total = sum(len(w) for w in words)

    exact = ConjugateTree(6, reserve=2 * total + 1)
    lean = ConjugateTree(6)
    for w in words:
        assert exact.process(list(w)) == lean.process(list(w))
    assert_same_tree(exact, lean)
    # 2 T + 1 is an upper bound on the number of nodes, never reached here
    assert exact.num_states() <= 2 * total + 1
