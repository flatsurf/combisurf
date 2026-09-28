import itertools
import random
import pytest

# the conjugate tree comes in two implementations and the fast one in three
# transition representations; every test below runs against all of them
#
# - "naive"   the pure Python reference, combisurf.conjugate_tree_naive
# - "dense"   the Cython one with an alphabet-indexed transition table
# - "sparse"  the Cython one walking a linked list of siblings
# - "rows"    the Cython one with sibling lists and dense rows for the nodes
#             with many children
# - "auto"    the Cython one picking between them on the alphabet
# - "unknown" the Cython one that was not told the alphabet
KINDS = ["naive", "dense", "sparse", "rows", "auto", "unknown"]


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


def process_and_check(T, w):
    r"""
    Process ``w`` in ``T`` and check all invariants of ``T``; the naive tree
    also checks its structure at each step of the insertion.
    """
    from combisurf.conjugate_tree_naive import ConjugateTreeNaive

    if isinstance(T, ConjugateTreeNaive):
        ans = T.process(w, hard_check=True)
    else:
        ans = T.process(w)
    T._check()
    return ans


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
    assert ConjugateTree(33).algorithm() == 'rows'
    assert ConjugateTree(256).algorithm() == 'rows'
    assert ConjugateTree().algorithm() == 'rows'
    assert ConjugateTree(256, algorithm='dense').algorithm() == 'dense'
    assert ConjugateTree(4, algorithm='sparse').algorithm() == 'sparse'
    assert ConjugateTree(4, algorithm='rows').algorithm() == 'rows'


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
@pytest.mark.parametrize("check", [True, False])
def test_process_copies_its_input(kind, check):
    # the tree keeps the primitive root of a power, which must not truncate
    # the list of the caller
    T = make_tree(kind, 2)
    w = [0, 1, 0, 1]
    assert T.process(w, check=check) == 2
    assert w == [0, 1, 0, 1]
    assert list(T.words()[0]) == [0, 1]
    assert T.process((0, 0, 1), check=check) == 1
    assert [list(u) for u in T.words()] == [[0, 1], [0, 0, 1]]


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
def test_checked_process(kind):
    # exercises the three possible outcomes of process() with the invariants
    # checked after each of them: a new primitive word, a new non-primitive
    # word, and a word that is a conjugate (possibly of a power) of an already
    # registered word.
    T = make_tree(kind, 2)
    assert process_and_check(T, [0]) == 1
    assert process_and_check(T, [0, 1, 0, 0, 1]) == 1
    assert process_and_check(T, [1, 0, 1, 0]) == 2
    assert process_and_check(T, [0, 1]) == -2
    T._check()


@pytest.mark.parametrize("kind", KINDS)
def test_checked_process_random(kind):
    rng = random.Random(0)
    outcomes = {"primitive": 0, "non_primitive": 0, "conjugate": 0}
    for _ in range(50):
        alphabet = rng.randint(2, 5)
        T = make_tree(kind, alphabet)
        for _ in range(20):
            length = rng.randint(1, 8)
            base = [rng.randrange(alphabet) for _ in range(length)]
            w = base * rng.choice([1, 1, 1, 2, 3])
            ans = process_and_check(T, w)
            if ans > 1:
                outcomes["non_primitive"] += 1
            elif ans == 1:
                outcomes["primitive"] += 1
            else:
                outcomes["conjugate"] += 1
        T._check()

    # all three outcomes of process() must be exercised, otherwise the
    # checks after each of them would not all be covered
    assert all(outcomes.values())


def cyclic_order_key(conjugate, order, pivot, depth):
    r"""
    Return the sequence that the order of the leaves compares: ``order`` of
    the first letter, then at each further letter ``c`` preceded by ``b`` the
    value ``(order[c] - pivot[b]) % n``.

    This is an independent reimplementation of the ordering that
    :meth:`~combisurf.conjugate_tree.ConjugateTree.cyclically_sorted_leaves`
    realizes through the tree.
    """
    n = len(order)
    l = len(conjugate)
    ans = [order[conjugate[0]]]
    for d in range(1, depth):
        ans.append((order[conjugate[d % l]] - pivot[conjugate[(d - 1) % l]]) % n)
    return ans


def check_cyclically_sorted_leaves(T, order, pivot):
    words = T.words()
    depth = 2 * sum(len(w) for w in words) + 4
    expected = sorted(((i, k) for i, w in enumerate(words) for k in range(len(w))),
                      key=lambda ik: cyclic_order_key(list(words[ik[0]][ik[1]:]) + list(words[ik[0]][:ik[1]]),
                                                      order, pivot, depth))
    leaves = T.cyclically_sorted_leaves(order, pivot)
    assert sorted(leaves) == T.leaves()
    assert [T.leaf_as_conjugate(s) for s in leaves] == expected


def random_order_and_pivot(rng, n):
    order = list(range(n))
    rng.shuffle(order)
    pivot = [rng.randrange(n) for _ in range(n)]
    return order, pivot


@pytest.mark.parametrize("kind", KINDS)
def test_cyclically_sorted_leaves(kind):
    T = make_tree(kind, 2)
    assert T.process([0, 1, 1]) == 1
    assert T.process([0, 1]) == 1
    assert T.cyclically_sorted_leaves([0, 1], [1, 0]) == [6, 1, 8, 4, 2]
    for pivot in ([1, 0], [0, 1], [0, 0], [1, 1]):
        check_cyclically_sorted_leaves(T, [0, 1], pivot)
        check_cyclically_sorted_leaves(T, [1, 0], pivot)

    rng = random.Random(20260923)
    T = make_tree(kind, 4)
    for w in [[0, 0, 2, 2], [0, 2, 0, 0, 3], [1, 3, 1, 2]]:
        T.process(w)
    check_cyclically_sorted_leaves(T, [0, 2, 1, 3], [2, 0, 3, 1])
    for _ in range(10):
        check_cyclically_sorted_leaves(T, *random_order_and_pivot(rng, 4))

    T = make_tree(kind, 8)
    for w in [[0, 2, 2, 5, 2, 2, 5], [0, 3, 6], [1, 4, 7, 0]]:
        T.process(w)
    check_cyclically_sorted_leaves(T, [0, 2, 4, 6, 1, 3, 5, 7], [2, 0, 6, 4, 3, 1, 7, 5])
    for _ in range(10):
        check_cyclically_sorted_leaves(T, *random_order_and_pivot(rng, 8))


@pytest.mark.parametrize("kind", KINDS)
def test_cyclically_sorted_leaves_pivot_from_order(kind):
    # pivot[b] = order[b ^ 1] gives the cyclic order at infinity of curves
    # on a surface; the leaves are hard-coded so that any change of that
    # order is caught
    for n, words, order, expected in [
            (2, [[0, 1, 1], [0, 1]], [1, 0], [8, 4, 2, 6, 1]),
            (4, [[0, 0, 2, 2], [0, 2, 0, 0, 3], [1, 3, 1, 2]], [0, 2, 1, 3],
             [13, 12, 1, 8, 3, 10, 6, 4, 20, 19, 15, 17, 14]),
            (8, [[0, 2, 2, 5, 2, 2, 5], [0, 3, 6], [1, 4, 7, 0]], [0, 2, 4, 6, 1, 3, 5, 7],
             [19, 1, 13, 17, 16, 5, 11, 4, 9, 2, 7, 15, 14, 18]),
            (6, [[0, 2, 4, 1, 5, 3, 3], [5, 5, 2, 0], [1, 2, 3, 4, 0]], [3, 0, 5, 1, 4, 2],
             [17, 4, 6, 8, 19, 11, 5, 10, 22, 15, 1, 21, 3, 18, 13, 2])]:
        T = make_tree(kind, n)
        for w in words:
            T.process(w)
        pivot = [order[b ^ 1] for b in range(n)]
        assert T.cyclically_sorted_leaves(order, pivot) == expected


@pytest.mark.parametrize("kind", KINDS)
def test_cyclically_sorted_leaves_random(kind):
    rng = random.Random(20260922)
    for _ in range(60):
        n = rng.randint(1, 8)
        T = make_tree(kind, n)
        for _ in range(rng.randint(1, 5)):
            T.process([rng.randrange(n) for _ in range(rng.randint(1, 9))])
        if not T.words():
            continue
        check_cyclically_sorted_leaves(T, *random_order_and_pivot(rng, n))


@pytest.mark.parametrize("kind", ["dense", "sparse", "rows", "auto", "unknown"])
def test_sorted_leaves_as_conjugates(kind):
    rng = random.Random(20260925)
    for _ in range(60):
        n = rng.randint(1, 8)
        T = make_tree(kind, n)
        for _ in range(rng.randint(1, 5)):
            T.process([rng.randrange(n) for _ in range(rng.randint(1, 9))])
        order, pivot = random_order_and_pivot(rng, n)
        words, shifts = T.sorted_leaves_as_conjugates(order, pivot)
        assert words.typecode == shifts.typecode == 'i'
        expected = [T.leaf_as_conjugate(s) for s in T.cyclically_sorted_leaves(order, pivot)]
        assert list(zip(words, shifts)) == expected


def test_sorted_leaves_as_conjugates_against_naive():
    from combisurf.conjugate_tree import ConjugateTree
    from combisurf.conjugate_tree_naive import ConjugateTreeNaive

    def outcome(T, order, pivot):
        try:
            return T.sorted_leaves_as_conjugates(order, pivot)
        except ValueError as e:
            return str(e)

    rng = random.Random(20260924)
    for _ in range(100):
        n = rng.randint(1, 8)
        T0 = ConjugateTree(n)
        T1 = ConjugateTreeNaive()
        for _ in range(rng.randint(0, 5)):
            w = [rng.randrange(n) for _ in range(rng.randint(1, 9))]
            assert T0.process(list(w)) == T1.process(list(w))
        order, pivot = random_order_and_pivot(rng, n)
        assert outcome(T0, order, pivot) == outcome(T1, order, pivot)
        # invalid arguments give the same error
        m = rng.randint(1, n)
        order, pivot = random_order_and_pivot(rng, m)
        if rng.random() < 0.5:
            pivot[rng.randrange(m)] = rng.choice([-1, m])
        if rng.random() < 0.3:
            order[rng.randrange(m)] = order[rng.randrange(m)]
        if rng.random() < 0.3:
            pivot = pivot[:-1]
        assert outcome(T0, order, pivot) == outcome(T1, order, pivot), (T0.words(), order, pivot)


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
    assert sorted(T0.graph().edges()) == sorted(T1.graph().edges())
    for s in range(T0.num_states()):
        assert T0.transitions(s) == T1.transitions(s), s
    for s in T0.leaves():
        assert T0.leaf_as_conjugate(s) == T1.leaf_as_conjugate(s), s
        assert T0._leaf_shift(s) == T1._leaf_shift(s), s
    for s in T0.internal_states():
        assert T0.internal_state_word(s) == T1.internal_state_word(s), s


# 4 and 8 are below the threshold of the dense layout, 34 and 64 above it
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
        algorithm = rng.choice([None, 'dense', 'sparse', 'rows'])
        T0 = ConjugateTree(alphabet, reserve=reserve, algorithm=algorithm)
        T1 = ConjugateTreeNaive()
        history = []
        for _ in range(rng.randint(1, 8)):
            base = [rng.randrange(alphabet) for _ in range(rng.randint(1, 12))]
            w = base * rng.choice([1, 1, 1, 2, 3])
            history.append(w)
            hard = (trial % 8 == 0)
            assert T0.process(list(w)) == T1.process(list(w), hard_check=hard), history
            if hard:
                T0._check()
            assert_same_tree(T0, T1)
        T0._check()
        T1._check()
        if T0.num_words():
            order, pivot = random_order_and_pivot(rng, alphabet)
            assert T0.cyclically_sorted_leaves(order, pivot) == T1.cyclically_sorted_leaves(order, pivot), (history, order, pivot)


def test_against_naive_long_words():
    r"""
    The same, on words long enough that the node arrays are reallocated many
    times over.
    """
    from combisurf.conjugate_tree import ConjugateTree
    from combisurf.conjugate_tree_naive import ConjugateTreeNaive

    rng = random.Random(4242)
    for alphabet, length in [(2, 2200), (5, 1300), (70, 1300)]:
        T0 = ConjugateTree(alphabet)
        T1 = ConjugateTreeNaive()
        for _ in range(3):
            w = [rng.randrange(alphabet) for _ in range(length)]
            assert T0.process(list(w)) == T1.process(list(w))
        assert_same_tree(T0, T1)
        order, pivot = random_order_and_pivot(rng, alphabet)
        assert T0.cyclically_sorted_leaves(order, pivot) == T1.cyclically_sorted_leaves(order, pivot)


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


def test_pprint_against_naive(capsys):
    r"""
    The two implementations print the same transitions, including the labels
    that wrap around the end of their word.
    """
    from combisurf.conjugate_tree import ConjugateTree
    from combisurf.conjugate_tree_naive import ConjugateTreeNaive

    rng = random.Random(1072)
    for trial in range(40):
        words = [[1, 1, 2, 1, 0, 0, 1, 0], [1, 0]]
        words += [[rng.randrange(3) for _ in range(rng.randint(1, 9))] for _ in range(rng.randint(0, 4))]
        rng.shuffle(words)
        T0 = ConjugateTree(3)
        T1 = ConjugateTreeNaive()
        for w in words:
            assert T0.process(list(w)) == T1.process(list(w))
        T0._pprint()
        out0 = capsys.readouterr().out
        T1._pprint()
        out1 = capsys.readouterr().out
        assert out0 == out1, words
        assert "array" not in out0
