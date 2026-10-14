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


def test_word_reduce_random():
    # exhaustive agreement with a reduction by repeated cancellation of one
    # pair of adjacent inverse letters, freely and then cyclically
    from array import array
    from itertools import product
    from combisurf.word import word_reduce, word_cyclically_reduce, word_cancel_ends

    def naive_reduce(w, cyclic):
        w = list(w)
        while True:
            for i in range(len(w) - 1):
                if w[i] ^ 1 == w[i + 1]:
                    del w[i:i + 2]
                    break
            else:
                if cyclic and len(w) > 1 and w[0] ^ 1 == w[-1]:
                    del w[-1]
                    del w[0]
                else:
                    return w

    for n in range(8):
        for ww in product(range(4), repeat=n):
            w = array('i', ww)
            assert list(word_reduce(w)) == naive_reduce(w, False), ww
            assert list(word_cyclically_reduce(w)) == naive_reduce(w, True), ww
            assert word_cancel_ends(word_reduce(w)) == word_cyclically_reduce(w), ww
            assert w == array('i', ww)  # the input is untouched


def test_word_reduce_typecode():
    # other typecodes are read too, and the output always has typecode 'i'
    # (arrays of distinct typecodes compare equal letter by letter, hence the
    # explicit check of the typecode)
    from array import array
    from combisurf.word import word_reduce, word_cyclically_reduce, word_cancel_ends, word_is_cyclically_reduced

    for typecode in 'lqhI':
        for f, w, ans in [(word_reduce, [], []),
                          (word_reduce, [4], [4]),
                          (word_reduce, [4, 2, 3, 0], [4, 0]),
                          (word_cyclically_reduce, [], []),
                          (word_cyclically_reduce, [4], [4]),
                          (word_cyclically_reduce, [0, 2, 3, 0, 1], [0]),
                          (word_cyclically_reduce, [0, 2, 3, 2, 1], [2]),
                          (word_cancel_ends, [], []),
                          (word_cancel_ends, [4], [4]),
                          (word_cancel_ends, [0, 2], [0, 2]),
                          (word_cancel_ends, [0, 4, 1], [4])]:
            v = f(array(typecode, w))
            assert v.typecode == 'i' and list(v) == ans, (f.__name__, typecode, w)

        assert word_is_cyclically_reduced(array(typecode, [0, 2, 4]))
        assert not word_is_cyclically_reduced(array(typecode, [0, 2, 3]))
        assert not word_is_cyclically_reduced(array(typecode, [0, 2, 1]))


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


def test_word_find():
    from itertools import product
    from combisurf.word import word_init, word_failure_table, word_find

    # exhaustive agreement with str.find, for every start
    for m in range(5):
        for uu in product([0, 1], repeat=m):
            u = word_init(uu)
            su = ''.join('ab'[x] for x in uu)
            t = word_failure_table(u)
            for n in range(7):
                for vv in product([0, 1], repeat=n):
                    v = word_init(vv)
                    sv = ''.join('ab'[x] for x in vv)
                    for start in range(-n - 2, n + 3):
                        expected = sv.find(su, start)
                        assert word_find(u, v, start) == expected, (uu, vv, start)
                        assert word_find(u, v, start, t) == expected, (uu, vv, start)


def test_word_is_factor():
    import random
    from itertools import product
    from combisurf.word import word_init, word_is_factor

    def naive(u, v):
        return any(v[i:i + len(u)] == u for i in range(len(v) - len(u) + 1))

    for m in range(5):
        for uu in product([0, 1], repeat=m):
            u = word_init(uu)
            for n in range(7):
                for vv in product([0, 1], repeat=n):
                    v = word_init(vv)
                    assert word_is_factor(u, v) == naive(u, v), (uu, vv)

    # larger alphabets and longer words
    rng = random.Random(0)
    for _ in range(20000):
        a = rng.choice([1, 2, 3, 9])
        u = word_init([rng.randrange(a) for _ in range(rng.randint(0, 10))])
        v = word_init([rng.randrange(a) for _ in range(rng.randint(0, 30))])
        assert word_is_factor(u, v) == naive(u, v), (list(u), list(v))


def test_word_border_table():
    from itertools import product
    from combisurf.word import word_init, word_border_table, word_failure_table

    for n in range(9):
        for uu in product([0, 1], repeat=n):
            u = word_init(uu)
            b = word_border_table(u)
            assert len(b) == n + 1
            # b[i] is the length of the longest proper border of u[:i]
            for i in range(n + 1):
                pre = uu[:i]
                expected = max([k for k in range(i) if pre[:k] == pre[i - k:]] or [0])
                assert b[i] == expected, (uu, i, list(b))
            # the failure table is the strong variant of the border table
            t = word_failure_table(u)
            ref = [-1] if n else []
            for i in range(1, n):
                ref.append(b[i] if u[i] != u[b[i]] else ref[b[i]])
            assert list(t) == ref, (uu, list(t), ref)


def test_word_apply_morphism_check():
    from array import array
    from combisurf.word import word_init, word_apply_morphism

    mor = [word_init([0, 2, 1, 2]), word_init([1]), word_init([3, 0]), word_init([0, 0])]
    w = word_init([0, 2, 0, 3, 1, 2, 1, 3])
    for reduce in [False, True]:
        assert word_apply_morphism(w, mor, reduce=reduce, check=False) == word_apply_morphism(w, mor, reduce=reduce)

    # a bad entry is rejected by check=True only
    for bad in [[0], array('l', [0]), None]:
        bad_mor = [word_init([0]), bad]
        with pytest.raises(ValueError, match="invalid morphism"):
            word_apply_morphism(word_init([0]), bad_mor)
        with pytest.raises(ValueError, match="invalid morphism"):
            word_apply_morphism(word_init([0]), bad_mor, check=True)
    assert word_apply_morphism(word_init([0]), [word_init([0]), [0]], check=False) == word_init([0])

    # the letters are always checked, a negative one does not wrap around
    for letter in [-1, 4]:
        for check in [False, True]:
            with pytest.raises(ValueError, match="invalid morphism"):
                word_apply_morphism(word_init([0, letter]), mor, check=check)


def check_walk_reference(w, vp, h2v, closed):
    n = len(vp)
    hprev = None
    for i, h in enumerate(w):
        if h < 0 or h >= n:
            return i
        if vp[h] == -1:
            return i
        if hprev is not None and h2v[hprev ^ 1] != h2v[h]:
            return i
        hprev = h
    if w and closed and h2v[w[0]] != h2v[w[-1] ^ 1]:
        return len(w)
    return -1


def test_word_check_walk():
    import random
    from array import array
    from combisurf import OrientedMap
    from combisurf.word import word_check_walk
    from test_barycentric_subdivision import random_maps
    from test_geometric_intersection_oracles import random_closed_walk

    rng = random.Random(5)
    maps = random_maps() + [OrientedMap(vp="(0,3,6,~3)(~0,1,~6,~1)"),
                            OrientedMap(vp="(0,1,4,5,~0,~1,~4,~5)")]
    for m in maps:
        vp = m.vertex_permutation()
        h2v = m.half_edge_to_vertex()
        n = len(vp)
        active = [h for h in range(n) if vp[h] != -1]
        inactive = [h for h in range(n) if vp[h] == -1]
        if not active:
            continue

        def compare(w):
            w = array('i', w)
            for closed in [False, True]:
                assert word_check_walk(w, vp, h2v, closed) == check_walk_reference(w, vp, h2v, closed), (list(w), closed)

        compare([])
        compare([rng.choice(active)])
        for _ in range(5):
            try:
                w = list(random_closed_walk(m, rng, rng.randint(1, 30)))
            except IndexError:
                # the helper needs a half-edge to leave every vertex
                break
            compare(w)
            compare(w[:rng.randint(1, len(w))])
            for bad in [-1, n, n + 7, rng.choice(active), rng.choice(active)] + inactive[:2]:
                v = list(w)
                v[rng.randrange(len(v))] = bad
                compare(v)
            v = list(w)
            v[0], v[-1] = v[-1] ^ 1, v[0] ^ 1
            compare(v)
