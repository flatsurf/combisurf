import pytest


def test_geometric_intersection():
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection

    torus = OrientedMap(fp="(0,1,~0,~1)")
    octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")

    torus_gi = GeometricIntersection(torus)
    octagon_gi = GeometricIntersection(octagon)

    for (cmap, u, v, expected_intersection) in [
        (torus_gi, [0], [2], 1),
        (torus_gi, [0], [0,2], 1),
        (torus_gi, [0,2], [0,2,0], 1),
        (torus_gi, [0,2,0], [0,2,0,0,2], 1),
        (torus_gi, [0,3,1,2], [0,2,0,0,2], 0),
        (torus_gi, [0,0,2,2], [1], 2),
        (torus_gi, [0,2,0,0,3], [0], 0),
        (torus_gi, [0,2,0,0,3], [2], 3),
        (torus_gi, [0], None, 0),
        (torus_gi, [0,2,0,0,2], None, 0),
        (torus_gi, [0,0,2,2], None, 1),
        (torus_gi, [0,2,2,0], None, 1),
        (torus_gi, [0,2,0,2,0,0], None, 1),
        (torus_gi, [0,2,2,0,2,2,2,2], None, 1),
        (torus_gi, [0,2,0,0,2,0,0,2,0,2], None, 1),
        (torus_gi, [0,2,0,0,2,0,2,0], None, 0),
        (torus_gi, [0,2,0,3], None, 1),
        (torus_gi, [1,3,1,2], None, 1),
        (torus_gi, [0,0,2,2,2], None, 2),
        (torus_gi, [0,0,2,2,2,2], None, 3),
        (torus_gi, [0,0,3,3,0,3], None, 2),
        (octagon_gi, [0,6], [6,0,2], 0)]:
        if v is None:
            computed_intersection = cmap.geometric_intersection([u], None)
        else:
            computed_intersection = cmap.geometric_intersection([u], [v])
        assert computed_intersection == expected_intersection, (cmap, u, v, expected_intersection, computed_intersection)


def test_torus_mcg():
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection

    torus = OrientedMap(fp="(0,1,~0,~1)")
    gi = GeometricIntersection(torus)

    f0 = [[2], [3], [0], [1]]
    f1 = [[1], [0], [2], [3]]
    f2 = [[0, 2], [3, 1], [2], [3]]
    f3 = [[0], [1], [2, 0], [1, 3]]
    f = [f0, f1, f2, f3]
    def apply_mcg(f, w):
        ww = []
        for i in w:
            ww.extend(f[i])
        return ww

    for w in [[0], [0, 0, 2, 2], [0, 0, 0, 2, 2], [0, 0, 2, 2, 0, 0, 2, 2], [0, 3, 1, 2], [0, 2, 0, 0, 3]]:
        intersection = gi.geometric_intersection([w])
        for s in [[0], [1], [2], [3], [0, 1], [0, 2], [1, 2], [2, 2, 3, 0], [3, 2, 1]]:
            ww = w[:]
            for i in s:
                ww = apply_mcg(f[i], ww)
            assert gi.geometric_intersection([ww]) == intersection


def test_geometric_intersection_multilinearity():
    import itertools
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection

    torus = OrientedMap(fp="(0,1,~0,~1)")
    gi = GeometricIntersection(torus)

    # non-primitivity self-intersections
    for u in [[0, 0, 2], [0, 2, 1, 3], [0, 0, 2, 2], [0, 0, 2, 2, 1, 1, 3]]:
        intersection = gi.geometric_intersection([u])
        assert gi.geometric_intersection([u * 2]) == 4 * intersection + 1
        assert gi.geometric_intersection([u, u]) == 4 * intersection
        assert gi.geometric_intersection([u * 2, u]) == 9 * intersection + 1
        assert gi.geometric_intersection([u * 3]) == 9 * intersection + 2

    wlist = [[0], [1], [0, 1], [0, 0, 1], [0, 1, 1], [0, 0, 1, 1]]
    Q = [[gi.geometric_intersection([w0], [w1]) for w1 in wlist] for w0 in wlist]
    for ucoeffs in [[3,5,1,0,2,4], [1,3,0,1,2,1], [5,1,2,3,0,3]]:
        ulist = [w * mult for w, mult in zip(wlist, ucoeffs)]
        for vcoeffs in [[0,1,0,2,0,3], [1,1,1,1,0,2], [2,6,4,1,5,3]]:
            vlist = [w * mult for w, mult in zip(wlist, vcoeffs)]
            ans0 = sum(ucoeffs[i] * vcoeffs[j] * Q[i][j] for i in range(6) for j in range(6))
            ans1 = gi.geometric_intersection(ulist, vlist)
            assert ans0 == ans1, gi


def polygon_4g(g):
    r"""
    Return the one vertex one face map obtained by identifying the sides of a
    ``4g``-gon, that is a surface of genus ``g`` with ``4 * g`` half-edges.
    """
    from combisurf import OrientedMap
    sides = [str(i) for i in range(2 * g)] + ["~%d" % i for i in range(2 * g)]
    return OrientedMap(fp="(" + ",".join(sides) + ")")


def random_primitive_curves(n, length, num, rng):
    r"""
    Return ``num`` distinct primitive cyclically reduced words of length
    ``length`` on the ``n`` half-edges ``0``, ..., ``n - 1``.
    """
    from combisurf.conjugate_tree import ConjugateTree
    from combisurf.word import word_init, word_cyclically_reduce

    curves = []
    seen = set()
    while len(curves) < num:
        w = [rng.randrange(n)]
        while len(w) < length:
            letter = rng.randrange(n)
            if letter != w[-1] ^ 1:
                w.append(letter)
        w = word_cyclically_reduce(word_init(w))
        if len(w) != length or tuple(w) in seen:
            continue
        # a conjugate of something already kept would share a slot, which is
        # fine, but a power would be rejected by GeometricIntersectionMatrix
        if ConjugateTree().process(w[:]) != 1:
            continue
        seen.add(tuple(w))
        curves.append(list(w))
    return curves


def test_intersection_matrix_torus_benchmark():
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection
    from combisurf.lyndon_word_family import cyclically_reduced_lyndon_words

    torus = OrientedMap(fp="(0,1,~0,~1)")
    gi = GeometricIntersection(torus)
    curves = [list(w) for w in cyclically_reduced_lyndon_words(torus.num_edges(), 1, 7, up_to_inverse=True)]
    assert len(curves) == 99

    I = gi.intersection_matrix(curves)
    for x, u in enumerate(curves):
        for y, v in enumerate(curves):
            assert I.entry(x, y) == gi.geometric_intersection([u], [v]), (x, y)

    mat = I.matrix()
    assert sum(sum(row) for row in mat.rows()) == 62902


def test_intersection_matrix_genus_and_length():
    import random
    from combisurf.geometric_intersection import GeometricIntersection

    rng = random.Random(20260922)
    for g in [1, 2, 4, 8, 16]:
        m = polygon_4g(g)
        gi = GeometricIntersection(m)
        for length in [3, 8, 40, 200]:
            curves = random_primitive_curves(4 * g, length, 5, rng)
            I = gi.intersection_matrix(curves)
            # both sides of the conjugate tree dense/sparse threshold are
            # covered by this range of genera
            assert I._tree.algorithm() == ('dense' if 4 * g <= 32 else 'sparse')
            for x, u in enumerate(curves):
                for y, v in enumerate(curves):
                    if y < x:
                        continue
                    assert I.entry(x, y) == gi.geometric_intersection([u], [v]), (g, length, x, y)


def test_intersection_matrix_large_alphabet():
    # A large alphabet, well beyond the tree's dense/sparse threshold.
    # GeometricIntersectionMatrix and the unit computation must pick the same
    # PartialSums structure for this map size, since both now go through the
    # shared combisurf.partial_sums.PartialSums factory.
    import random
    from combisurf.geometric_intersection import GeometricIntersection
    from combisurf.partial_sums import PartialSums

    rng = random.Random(1234)
    g = 70
    m = polygon_4g(g)
    gi = GeometricIntersection(m)
    for length in [5, 30]:
        curves = random_primitive_curves(4 * g, length, 4, rng)
        I = gi.intersection_matrix(curves)
        assert type(I._Nu) is type(PartialSums(4 * g - 1))
        assert type(I._Nv) is type(PartialSums(4 * g - 1))
        assert I._tree.algorithm() == 'sparse'
        for x, u in enumerate(curves):
            for y, v in enumerate(curves):
                assert I.entry(x, y) == gi.geometric_intersection([u], [v]), (length, x, y)


def test_intersection_matrix_conjugates_and_inverses():
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection
    from combisurf.word import word_init, word_free_group_inverse

    torus = OrientedMap(fp="(0,1,~0,~1)")
    gi = GeometricIntersection(torus)

    w = word_init([0, 0, 2, 0, 3])
    curves = [w, w[2:] + w[:2], word_free_group_inverse(w), word_init([0, 2]), word_init([0, 2, 2, 0, 3])]
    I = gi.intersection_matrix(curves)

    # a word, one of its conjugates and its inverse share a slot
    assert I._slot == [0, 0, 0, 1, 2]

    for x, u in enumerate(curves):
        for y, v in enumerate(curves):
            assert I.entry(x, y) == gi.geometric_intersection([list(u)], [list(v)]), (x, y)


def test_intersection_matrix_non_primitive():
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection
    from combisurf.word import word_init, word_free_group_inverse

    torus = OrientedMap(fp="(0,1,~0,~1)")
    gi = GeometricIntersection(torus)

    # a power of a curve that is not in the list
    with pytest.raises(NotImplementedError):
        gi.intersection_matrix([[0, 2, 0, 2]])
    # a power of a curve that is already in the list
    with pytest.raises(NotImplementedError):
        gi.intersection_matrix([[0, 2], [0, 2, 0, 2]])
    # a power of the inverse of a curve that is already in the list
    inverse_square = list(word_free_group_inverse(word_init([0, 2, 0, 2])))
    with pytest.raises(NotImplementedError):
        gi.intersection_matrix([[0, 2], inverse_square])
    # a curve that is trivial in the free group
    with pytest.raises(ValueError):
        gi.intersection_matrix([[0, 1]])


def test_intersection_matrix_row_and_matrix():
    from sage.rings.integer_ring import ZZ
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection

    octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")
    gi = GeometricIntersection(octagon)
    curves = [[0], [3], [0, 3, 6], [0, 2, 2, 5, 2, 2, 5], [0, 4, 1, 5]]
    I = gi.intersection_matrix(curves)

    mat = I.matrix()
    assert mat.is_symmetric()
    assert mat.base_ring() is ZZ
    assert mat.nrows() == mat.ncols() == len(curves) == len(I)
    for x in range(len(curves)):
        assert I.row(x) == list(mat.row(x)), x


def naive_double_sum_matrix_class():
    r"""
    Return a subclass of ``GeometricIntersectionMatrix`` whose crossing arcs
    term goes through the pure Python oracle of the sorted sweep, so that it
    can be compared against the Cython one.
    """
    from combisurf.crossing_arcs_naive import crossing_arcs_sweep_sorted
    from combisurf.geometric_intersection import GeometricIntersectionMatrix

    class NaiveDoubleSum(GeometricIntersectionMatrix):
        def _double_sum(self, sx, sy):
            keys = self._arc_keys
            weights = self._arc_weights
            if sx == sy:
                return crossing_arcs_sweep_sorted(self._n, keys[sx], weights[sx])
            return crossing_arcs_sweep_sorted(self._n, keys[sx], weights[sx], keys[sy], weights[sy])

    return NaiveDoubleSum


def test_intersection_matrix_double_sum_paths():
    # the Cython sweep against its pure Python oracle, and matrix(), row() and
    # entry() against each other and against geometric_intersection
    import random
    from combisurf.geometric_intersection import GeometricIntersection, GeometricIntersectionMatrix

    naive = naive_double_sum_matrix_class()
    rng = random.Random(20260923)
    for g, length in [(1, 8), (2, 8), (4, 8), (8, 8), (16, 8), (32, 8), (8, 100)]:
        m = polygon_4g(g)
        gi = GeometricIntersection(m)
        curves = random_primitive_curves(4 * g, length, 12, rng)
        fast = GeometricIntersectionMatrix(gi, curves)
        slow = naive(gi, curves)

        mat = fast.matrix()
        assert mat == slow.matrix(), (g, length)
        for x in range(len(curves)):
            row = [fast.entry(x, y) for y in range(len(curves))]
            assert row == [slow.entry(x, y) for y in range(len(curves))], (g, length, x)
            assert fast.row(x) == row, (g, length, x)
            assert slow.row(x) == row, (g, length, x)
            assert list(mat.row(x)) == row, (g, length, x)
        for x, y in [(0, 0), (0, 1), (3, 7), (11, 5)]:
            assert mat[x, y] == gi.geometric_intersection([curves[x]], [curves[y]]), (g, length, x, y)


def test_intersection_matrix_double_sum_identity():
    # the crossing arcs term of an entry is the double sum over the arcs of
    # the two slots, with u-weights from the first and v-weights from the
    # second, as geometric_intersection builds them
    import random
    from combisurf.crossing_arcs_naive import crossing_arcs_double_sum
    from combisurf.geometric_intersection import GeometricIntersection

    rng = random.Random(20260926)
    for g, length in [(1, 8), (2, 8), (4, 8), (8, 8), (8, 100)]:
        gi = GeometricIntersection(polygon_4g(g))
        n = 4 * g
        I = gi.intersection_matrix(random_primitive_curves(n, length, 8, rng))
        words = I._tree.words()
        num_slots = len(words) // 2
        for sx in range(num_slots):
            for sy in range(num_slots):
                arcs = {}
                for s, side in [(sx, 0), (sy, 1)]:
                    w = words[2 * s]
                    for p in range(len(w)):
                        first, last = sorted([gi._angles[w[p]], gi._angles[w[p - 1] ^ 1]])
                        weights = arcs.setdefault(last * n + first, [0, 0])
                        weights[side] += 1
                if sx == sy:
                    for weights in arcs.values():
                        weights[1] = weights[0]
                assert I._double_sum(sx, sy) == crossing_arcs_double_sum(n, arcs), (g, length, sx, sy)


def test_crossing_arcs_sweep_sorted_random():
    # the sorted sweep against its pure Python oracle and against the dict
    # sweep, exactly, on random weighted arcs, with and without a scratch
    import random
    from array import array
    from combisurf.crossing_arcs import crossing_arcs_sweep, crossing_arcs_sweep_sorted
    from combisurf import crossing_arcs_naive

    rng = random.Random(20260927)
    for n in list(range(1, 20)) + [64, 257]:
        scratch = array('q', [0]) * (2 * (n + 1))
        for num in [0, 1, 2, 5, 30, 200]:
            if n < 2 and num:
                continue
            sides = []
            for _ in range(2):
                weights = {}
                for _ in range(num):
                    first, last = sorted(rng.sample(range(n), 2))
                    weights[last * n + first] = rng.randrange(1, 4)
                keys = sorted(weights)
                sides.append((keys, [weights[k] for k in keys]))
            (ukeys, uweights), (vkeys, vweights) = sides
            arcs = {}
            for side, (keys, weights) in enumerate(sides):
                for k, x in zip(keys, weights):
                    arcs.setdefault(k, [0, 0])[side] = x
            expected = crossing_arcs_sweep(n, arcs)
            q = [array('q', x) for x in (ukeys, uweights, vkeys, vweights)]
            assert crossing_arcs_naive.crossing_arcs_sweep_sorted(n, ukeys, uweights, vkeys, vweights) == expected
            assert crossing_arcs_sweep_sorted(n, *q) == expected
            assert crossing_arcs_sweep_sorted(n, *q, scratch) == expected
            assert not any(scratch)

            expected = crossing_arcs_sweep(n, {k: [x, x] for k, x in zip(ukeys, uweights)}, True)
            assert crossing_arcs_naive.crossing_arcs_sweep_sorted(n, ukeys, uweights) == expected
            assert crossing_arcs_sweep_sorted(n, q[0], q[1]) == expected
            assert crossing_arcs_sweep_sorted(n, q[0], q[1], None, None, scratch) == expected
            assert crossing_arcs_sweep_sorted(n, q[0], q[1], q[0], q[1], scratch) == expected
            assert not any(scratch)


def crossing_arcs_implementations():
    r"""
    Return the three functions computing the crossing arcs term of
    ``GeometricIntersection.geometric_intersection``: the Cython sweep, its
    pure Python version and the `O(n^2)` double sum.
    """
    from combisurf.crossing_arcs import crossing_arcs_sweep
    from combisurf import crossing_arcs_naive

    return [crossing_arcs_sweep, crossing_arcs_naive.crossing_arcs_sweep,
            crossing_arcs_naive.crossing_arcs_double_sum]


def test_crossing_arcs_random():
    # the sweeps against the double sum, exactly, on random weighted arcs,
    # including arcs sharing one endpoint with many others
    import random

    rng = random.Random(20260924)
    implementations = crossing_arcs_implementations()
    for n in list(range(1, 20)) + [64, 257]:
        for num in [0, 1, 2, 5, 30, 200]:
            if n < 2 and num:
                continue
            arcs = {}
            for _ in range(num):
                first, last = sorted(rng.sample(range(n), 2))
                arcs[last * n + first] = [rng.randrange(4), rng.randrange(4)]
            answers = [f(n, arcs) for f in implementations]
            assert answers.count(answers[0]) == 3, (n, arcs, answers)
            for weights in arcs.values():
                weights[1] = weights[0]
            answers = [f(n, arcs, True) for f in implementations]
            answers.append(implementations[0](n, arcs, False))
            assert answers.count(answers[0]) == 4, (n, arcs, answers)


def test_geometric_intersection_crossing_arcs_paths(monkeypatch):
    # geometric_intersection with the crossing arcs term computed by each of
    # the three implementations, which must agree exactly on every call
    import random
    import combisurf.geometric_intersection as geometric_intersection
    from combisurf.geometric_intersection import GeometricIntersection
    from combisurf.word import word_init, word_free_group_inverse

    implementations = crossing_arcs_implementations()
    calls = []

    def all_paths(n, arcs, symmetric):
        answers = [f(n, arcs, symmetric) for f in implementations]
        assert answers.count(answers[0]) == 3, (n, arcs, symmetric, answers)
        calls.append(answers[0])
        return answers[0]

    def forced(f):
        def call(ulist, vlist=None):
            monkeypatch.setattr(geometric_intersection, "crossing_arcs_sweep", f)
            return gi.geometric_intersection(ulist, vlist)
        return call

    paths = [forced(f) for f in implementations + [all_paths]]

    rng = random.Random(20260925)
    for g in [1, 2, 4, 8, 16, 32]:
        gi = GeometricIntersection(polygon_4g(g))
        for length in [1, 2, 8, 100]:
            c0, c1, c2, c3 = random_primitive_curves(4 * g, length, 4, rng)
            conj = c1[1:] + c1[:1]
            inv = list(word_free_group_inverse(word_init(c2)))
            inputs = [
                ([c0], [c1]),
                ([c0], None),
                ([c0, c0, c1 * 2, conj, inv], [c1, c2 * 3, c0, c3]),
                ([c0, c0, c1 * 2, conj, inv], None),
                ([c2 * 2, c3, c3], [c3 * 2, inv]),
                ([c2 * 3, c3], None)]
            for ulist, vlist in inputs:
                ncalls = len(calls)
                answers = [path(ulist, vlist) for path in paths]
                assert len(calls) == ncalls + 1
                assert answers.count(answers[0]) == 4, (g, length, ulist, vlist, answers)
