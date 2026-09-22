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


def forced_matrix_classes():
    r"""
    Return the two subclasses of ``GeometricIntersectionMatrix`` that pin the
    evaluation of the `O(n^2)` term to one of its two paths, so that they can
    be compared against each other.
    """
    from combisurf.geometric_intersection import GeometricIntersectionMatrix

    class DoubleSumPerPair(GeometricIntersectionMatrix):
        _DOT_SETUP_WORK = float('inf')

    class DoubleSumByProduct(GeometricIntersectionMatrix):
        _DOT_SETUP_WORK = 0

    return DoubleSumPerPair, DoubleSumByProduct


def test_intersection_matrix_double_sum_paths():
    # the matrix product path against the entry-by-entry double sum, exactly
    import random
    from combisurf.geometric_intersection import GeometricIntersection

    per_pair, by_product = forced_matrix_classes()
    rng = random.Random(20260923)
    for g, length in [(1, 8), (2, 8), (4, 8), (8, 8), (16, 8), (32, 8), (8, 100)]:
        m = polygon_4g(g)
        gi = GeometricIntersection(m)
        curves = random_primitive_curves(4 * g, length, 12, rng)
        slow = per_pair(gi, curves)
        fast = by_product(gi, curves)
        assert not slow._dot_arrays(1)
        assert fast._dot_arrays(0)

        assert slow.matrix() == fast.matrix(), (g, length)
        for x in range(len(curves)):
            row = [slow.entry(x, y) for y in range(len(curves))]
            assert slow.row(x) == row, (g, length, x)
            assert fast.row(x) == row, (g, length, x)


def test_intersection_matrix_double_sum_gate():
    # the matrix product is set up once the calls add up to more than it costs,
    # and both sides of that switch give the same answers
    import random
    from combisurf.geometric_intersection import GeometricIntersection

    rng = random.Random(4242)
    g = 16
    gi = GeometricIntersection(polygon_4g(g))
    curves = random_primitive_curves(4 * g, 8, 200, rng)
    I = gi.intersection_matrix(curves)
    work = len(curves) * I._K
    assert work < I._DOT_SETUP_WORK < 2 * work

    expected = [I.entry(0, y) for y in range(len(curves))]
    assert I.row(0) == expected
    assert not I._dot                      # first row, still the double sum
    assert I.row(0) == expected
    assert I._dot                          # second row, now the product
    assert I.row(0) == expected


def test_intersection_matrix_float64_bound():
    # when the products could leave the exactly representable integers the
    # class must fall back to the Python double sum rather than round
    import random
    from combisurf.geometric_intersection import GeometricIntersection

    _, by_product = forced_matrix_classes()
    rng = random.Random(55)
    g = 8
    gi = GeometricIntersection(polygon_4g(g))
    curves = random_primitive_curves(4 * g, 8, 10, rng)

    I = by_product(gi, curves)
    assert I._dot_arrays(0)                # K * L^2 is nowhere near 2^53 here
    expected = I.matrix()

    J = by_product(gi, curves)
    J._K = 2 ** 53                         # as if the curves were enormous
    assert not J._dot_arrays(0)
    assert J._double_sum_row(0) is None
    assert J._double_sum_table() is None
    assert J.matrix() == expected
