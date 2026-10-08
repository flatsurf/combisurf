import pytest

def test_constructor():
    from combisurf import OrientedMap

    OrientedMap(None, None)
    OrientedMap([], None)
    OrientedMap(None, [])

    OrientedMap([2, 1, 3, 0], None)
    OrientedMap(None, [2, 1, 3, 0])
    OrientedMap([0, -1, 3, 4, 5, 2])
    OrientedMap(None, [0, -1, 3, 4, 5, 2])

    OrientedMap([0, -1, 3, 4, 5, 2], None)
    OrientedMap(None, [0, -1, 3, 4, 5, 2])

    OrientedMap("(0,1,2)(~0,~1,~2)", None)
    OrientedMap(None, "(0,1,2)(~0,~1,~2)")

    with pytest.raises(ValueError, match="vp is not a permutation"):
        OrientedMap([0, 0], None)

    with pytest.raises(ValueError, match="fp is not a permutation"):
        OrientedMap(None, [0, 0])

    with pytest.raises(ValueError, match="different lengths"):
        OrientedMap([1, 0], [1, 0, 3, 2])
    with pytest.raises(ValueError, match="different lengths"):
        OrientedMap([1, 0, 3, 2], [1, 0])

    with pytest.raises(ValueError, match="different domains"):
        OrientedMap([1, 0, 3, 2, 5, 4], [1, 0, -1, -1, 5, 4])

    OrientedMap([0, -1], [0, -1])
    with pytest.raises(ValueError, match="is active but its twin"):
        OrientedMap([-1, 1], [-1, 1])
    with pytest.raises(ValueError, match="is active but its twin"):
        OrientedMap([-1,1,2,-1],[-1,1,2,-1])

    with pytest.raises(ValueError, match="trailing inactive edges"):
        OrientedMap([0, -1, -1, -1], [0, -1, -1, -1])

    with pytest.raises(ValueError, match="fev relation not satisfied"):
        OrientedMap([0, 1], [0, 1])


def test_check_half_edge():
    from combisurf import OrientedMap

    assert OrientedMap(vp="(0,~0)")._check_half_edge(0) == 0

    with pytest.raises(TypeError):
        OrientedMap(vp="(0,~0)")._check_half_edge("1")

    with pytest.raises(ValueError):
        OrientedMap(vp="(0,~0)")._check_half_edge(2)

    assert OrientedMap(vp="(0)")._check_half_edge(0) == 0

    with pytest.raises(ValueError):
        OrientedMap(vp="(0")._check_half_edge(1)

    with pytest.raises(ValueError):
        OrientedMap(vp="(0,2)")._check_half_edge(2)


def test_check_edge():
    from combisurf import OrientedMap

    assert OrientedMap(fp="(0,1,~1)")._check_edge(0) == 0
    assert OrientedMap(fp="(0,1,~1)")._check_edge(1) == 1

    with pytest.raises(ValueError):
        OrientedMap(fp="(0,1,~1)")._check_edge(-1)

    with pytest.raises(ValueError):
        OrientedMap("(0,1,~1)")._check_edge(2)

    with pytest.raises(ValueError):
        OrientedMap(fp="(0,2)")._check_edge(1)


def test_cmp():
    from combisurf import OrientedMap

    assert OrientedMap("(0,1,2)(~0,~1,~2)", None) == OrientedMap("(0,1,2)(~0,~1,~2)", None)
    assert not (OrientedMap("(0,1,2)(~0,~1,~2)", None) == OrientedMap("(0,~0,1)(~1,2,~2)", None))
    assert OrientedMap("(0,1,2,~0,~1,~2)", None) == OrientedMap("(0,1,2,~0,~1,~2)", None)
    assert not (OrientedMap("(0,1,2,~0,~1,~2)", None) == OrientedMap("(0,1,~2,~0,~1,2)", None))

    assert not (OrientedMap("(0,1,2)(~0,~1,~2)", None) != OrientedMap("(0,1,2)(~0,~1,~2)", None))
    assert OrientedMap("(0,1,2)(~0,~1,~2)", None) != OrientedMap("(0,~0,1)(~1,2,~2)", None)
    assert not (OrientedMap("(0,1,2,~0,~1,~2)", None) != OrientedMap("(0,1,2,~0,~1,~2)", None))
    assert OrientedMap("(0,1,2,~0,~1,~2)", None) != OrientedMap("(0,1,~2,~0,~1,2)", None)


def test_pickling():
    from combisurf import OrientedMap
    from pickle import loads, dumps

    t = OrientedMap(fp="(0,1,2)")
    assert loads(dumps(t)) == t
    t = OrientedMap("(0,1,2)(~0,3,4)")
    assert loads(dumps(t)) == t

    t0 = OrientedMap("(0,1,2)", mutable=False)
    t1 = OrientedMap("(0,1,2)", mutable=True)
    s0 = loads(dumps(t0))  # indirect doctest
    assert s0 == t0 and not s0.is_mutable()
    s0._check()
    s1 = loads(dumps(t1))  # indirect doctest
    assert s1 == t1 and s1.is_mutable()
    s1._check()


def test_hash():
    from itertools import permutations, combinations
    from combisurf import OrientedMap


    m0 = OrientedMap(vp="(0,~0)")
    m1 = OrientedMap(vp="(0,~0)", mutable=True)
    hash(m0)
    with pytest.raises(ValueError):
        hash(m1)


    maps = []
    maps.append(OrientedMap(fp="(0,1,2)"))
    for p in permutations(["1", "~1", "2", "~2"]):
        maps.append(OrientedMap(fp="(0,{},{})(~0,{},{})".format(*p)))

    for i, j in combinations([0, 1, 2, 3], 2):
        for k, l in permutations(set(range(4)).difference([i, j])):
            vars = {'i': i, 'j': j, 'k': k, 'l': l}
            m = OrientedMap(fp="({i},{j},{k})(~{i},~{j},{l})".format(**vars))
            maps.append(m)
            m = OrientedMap(fp="({i},~{j},{k})(~{i},{j},{l})".format(**vars))
            maps.append(m)
            m = OrientedMap(fp="({i},{k},{j})(~{i},~{j},{l})".format(**vars))
            maps.append(m)
            m = OrientedMap(fp="({i},{k},~{j})(~{i},{j},{l})".format(**vars))
            maps.append(m)
            m = OrientedMap(fp="({i},{j},{k})(~{i},{l},~{j})".format(**vars))
            maps.append(m)
            m = OrientedMap(fp="({i},~{j},{k})(~{i},{l},{j})".format(**vars))
            maps.append(m)
            m = OrientedMap(fp="({i},{k},{j})(~{i},{l},~{j})".format(**vars))
            maps.append(m)
            m = OrientedMap(fp="({i},{k},~{j})(~{i},{l},{j})".format(**vars))
            maps.append(m)

    for i in range(len(maps)):
        for j in range(len(maps)):
            assert (maps[i] == maps[j]) == (i == j), (i, j)
            assert (maps[i] != maps[j]) == (i != j), (i, j)

    hashes = {}
    for t in maps:
        h = hash(t)
        assert h not in hashes
        hashes[h] = t
    assert len(hashes) == len(maps)

    hashes1 = {}
    hashes2 = {}
    for m in maps:
        h1 = hash(m) % (2 ** 16)
        h2 = (hash(m) >> 16) % (2 ** 16)
        assert h1 not in hashes1
        hashes1[h1] = m
        assert h2 not in hashes2
        hashes2[h2] = m
    assert len(hashes1) == len(hashes2) == len(maps)

def test_copy():
    from combisurf import OrientedMap

    for mutable in [True, False]:
        m = OrientedMap(fp="(0,1,2)(~0,~1,~2)", mutable=mutable)
        m1 = m.copy()
        m2 = m.copy(mutable=True)
        m3 = m.copy(mutable=False)
        assert m == m1
        assert m == m2
        assert m == m3
        assert m1.is_mutable() == m.is_mutable()
        assert m2.is_mutable()
        assert not m3.is_mutable()


def small_maps(folded=True):
    from combisurf import OrientedMap

    if folded:
        yield OrientedMap("(0)")

    yield OrientedMap("(0,~0)")
    yield OrientedMap("(0)(~0)")

    if folded:
        yield OrientedMap("(0,1)")
        yield OrientedMap("(0)(~0,1)")
        yield OrientedMap("(~0)(0,1)")
        yield OrientedMap("(1)(0,~1)")
        yield OrientedMap("(~1)(0,1)")
        yield OrientedMap("(0,1,~1)")
        yield OrientedMap("(0,~1,1)")
        yield OrientedMap("(0,~0,1)")
        yield OrientedMap("(0,1,~0)")

    yield OrientedMap("(0,1)(~0)(~1)")
    yield OrientedMap("(~0,1)(0)(~1)")
    yield OrientedMap("(0,~1)(~0)(1)")
    yield OrientedMap("(~0,~1)(0)(1)")

    yield OrientedMap("(0,1,~1)(~0)")

    yield OrientedMap("(0,~0,1,~1)")
    yield OrientedMap("(0,~0,~1,1)")
    yield OrientedMap("(0,1,~0,~1)")
    yield OrientedMap("(0,1,~1,~0)")
    yield OrientedMap("(0,~1,~0,1)")
    yield OrientedMap("(0,~1,1,~0)")

    # the same shapes with edge 0 inactive: the labels of a map need not start
    # at 0, and neither the mutations nor the vertex and face indexing may
    # assume that they do
    if folded:
        yield OrientedMap("(1)")

    yield OrientedMap("(1,~1)")
    yield OrientedMap("(1)(~1)")
    yield OrientedMap("(1,2)(~1)(~2)")
    yield OrientedMap("(1,~1,2,~2)")


def test_inactive_edge_zero():
    # vertex and face indices come from perm_dense_cycles, which numbers only
    # the active cycles, so index 0 is a real vertex even here
    from combisurf import OrientedMap
    from pickle import loads, dumps

    m = OrientedMap(vp="(2,1,5)(~1,~2,~5)")
    m._check()
    assert m.vertex_permutation(copy=False)[0] == -1
    assert list(m.half_edges()) == [2, 3, 4, 5, 10, 11]
    assert m.num_vertices() == 2 and m.num_faces() == 3

    assert loads(dumps(m)) == m
    assert m.copy() == m
    assert hash(m) == hash(m.copy())

    # the default roots are vertex 0 and face 0, which exist
    assert m.forest_coforest_decomposition() == m.forest_coforest_decomposition((0,), (0,))

    r = m.copy(mutable=True)
    r.relabel()
    r._check()
    assert list(r.half_edges()) == [0, 1, 2, 3, 4, 5]


def test_reverse_orientation():
    for m0 in small_maps(folded=False):
        m1 = m0.copy(mutable=True)
        for e in m1.edge_indices():
            m1.reverse_orientation(e)
            m1._check()
            assert m0.vertex_profile(sort=True) == m1.vertex_profile(sort=True), m0
            assert m1.face_profile(sort=True) == m1.face_profile(sort=True), m0
            assert m0.euler_characteristic() == m1.euler_characteristic(), m0
            m1.reverse_orientation(e)
            m1._check()
            assert m1 == m0


def test_insert_edge_contract_edge():
    from combisurf import OrientedMap

    for m in small_maps():
        for e in m.edge_indices():
            mm = m.copy(mutable=True)
            mm.contract_edge(e)
            mm._check()
            assert mm.num_edges() == m.num_edges() - 1

    m = OrientedMap("", "", mutable=True)

    m.insert_edge(-1, -1)
    m._check()
    assert m == OrientedMap("(0)(~0)", "(0,~0)")

    m.insert_edge(-1, -2)
    m._check()
    assert m == OrientedMap("(0)(~0)(1,~1)", "(0,~0)(1)(~1)")

    m.insert_edge(0, -1)
    m.insert_edge(-1, 2)
    m._check()
    assert m == OrientedMap("(0)(~0,~2,2)(1,3,~3,~1)", "(0,2,~0)(1,~3)(~1)(~2)(3)")

    m.insert_edge(0, 6)
    m._check()
    assert m == OrientedMap("(0)(~0,~2,2,~4,~3,~1,1,3,4)", "(0,4,2,~0)(1,~3)(~1)(~2)(3,~4)")

    m.contract_edge(2)
    m._check()
    assert m == OrientedMap("(0)(~0,~4,~3,~1,1,3,4)", "(0,4,~0)(1,~3)(~1)(3,~4)")

    m.contract_edge(0)
    m._check()
    assert m == OrientedMap("(1,3,4,~4,~3,~1)", "(1,~3)(~1)(3,~4)(4)")


def test_add_edge_delete_edge():
    from combisurf import OrientedMap

    m = OrientedMap("", "", mutable=True)

    m.add_edge(-1, -1)
    m._check()
    assert m == OrientedMap("(0,~0)", "(0)(~0)")

    m.add_edge(-1, -2)
    m._check()
    assert m == OrientedMap("(0,~0)(1)(~1)", "(0)(~0)(1,~1)")

    m.add_edge(0, -1)
    m._check()
    assert m == OrientedMap("(0,2,~0)(1)(~1)(~2)", "(0,2,~2)(~0)(1,~1)")

    m.add_edge(-1, 0)
    m._check()
    assert m == OrientedMap("(0,~3,2,~0)(1)(~1)(~2)(3)", "(0,2,~2,~3,3)(~0)(1,~1)")

    m.add_edge(1, 1)
    m._check()
    assert m == OrientedMap("(0,~3,2,~0,4,~4)(1)(~1)(~2)(3)", "(0,2,~2,~3,3)(~0,~4)(1,~1)(4)")

    m.add_edge(0, 2)
    m._check()
    assert m == OrientedMap("(0,5,~3,2,~0,4,~4)(1,~5)(~1)(~2)(3)", "(0,2,~2,~3,3,5,1,~1,~5)(~0,~4)(4)")

    m.delete_edge(5)
    m._check()
    assert m == OrientedMap("(0,~3,2,~0,4,~4)(1)(~1)(~2)(3)", "(0,2,~2,~3,3)(~0,~4)(1,~1)(4)")

    m.delete_edge(4)
    m._check()
    assert m == OrientedMap("(0,~3,2,~0)(1)(~1)(~2)(3)", "(0,2,~2,~3,3)(~0)(1,~1)")

    m.delete_edge(3)
    m._check()
    assert m == OrientedMap("(0,2,~0)(1)(~1)(~2)", "(0,2,~2)(~0)(1,~1)")

    m.delete_edge(2)
    m._check()
    assert m == OrientedMap("(0,~0)(1)(~1)", "(0)(~0)(1,~1)")

    m.delete_edge(1)
    m._check()
    assert m == OrientedMap("(0,~0)", "(0)(~0)")

    m.delete_edge(0)
    m._check()
    assert m == OrientedMap("", "")


def test_relabel():
    from combisurf import OrientedMap
    from combisurf.permutation import perm_compose, perm_random_centralizer

    m = OrientedMap(fp="(0,1,2)(~0,~1,~2)", mutable=True)
    for _ in range(10):
        r = perm_random_centralizer(m.edge_permutation())
        m.relabel(r)
        m._check()

    fp = "(0,16,~15)(1,19,~18)(2,22,~21)(3,21,~20)(4,20,~19)(5,23,~22)(6,18,~17)(7,17,~16)(8,~1,~23)(9,~2,~8)(10,~3,~9)(11,~4,~10)(12,~5,~11)(13,~6,~12)(14,~7,~13)(15,~0,~14)"
    m = OrientedMap(fp=fp)
    ep = m.edge_permutation()
    for _ in range(10):
        p1 = perm_random_centralizer(ep)
        p2 = perm_random_centralizer(ep)
        m1 = m.copy(mutable=True)
        m1.relabel(p1)
        m1.relabel(p2)
        m2 = m.copy(mutable=True)
        m2.relabel(perm_compose(p1, p2))
        assert m1  == m2


def half_edge_to_cell_corpus():
    # the small maps, plus maps with folded edges, with inactive half-edges,
    # disconnected ones and the empty map
    from combisurf import OrientedMap
    from test_fold import sample_maps, maps_with_a_folded_edge

    out = list(small_maps()) + sample_maps() + maps_with_a_folded_edge()
    out += [
        OrientedMap(),                                   # empty
        OrientedMap(vp="(2,1,5)(~1,~2,~5)"),             # edge 0 inactive
        OrientedMap(fp="(0,1,3)(~0,~1,~3)(2,4,5)(~2,~4,~5)"),   # disconnected
        OrientedMap("(0,2,~2)"),                         # a folded edge
    ]
    return out


def test_half_edge_to_vertex_and_face():
    # the arrays are the inverse of vertices() and faces(): entry h is the
    # index of the cell containing h. This alignment is what
    # forest_coforest_decomposition relies on to go from a half-edge to the
    # cell it belongs to.
    for m in half_edge_to_cell_corpus():
        for cells, h2c in ((m.vertices(), m.half_edge_to_vertex()),
                           (m.faces(), m.half_edge_to_face())):
            assert len(h2c) == len(m.vertex_permutation(copy=False)), m

            seen = set()
            for i, c in enumerate(cells):
                for h in c:
                    assert h2c[h] == i, (m, h, i, list(h2c))
                    seen.add(h)

            # exactly the active half-edges are labelled, and with -1 elsewhere
            assert seen == set(m.half_edges()), m
            for h in range(len(h2c)):
                if h in seen:
                    assert h2c[h] != -1, (m, h)
                else:
                    assert h2c[h] == -1, (m, h)


def test_half_edge_to_cell_labels_are_consecutive():
    # the labels run 0, 1, ..., and there are as many as cells
    for m in half_edge_to_cell_corpus():
        for n, h2c in ((m.num_vertices(), m.half_edge_to_vertex()),
                       (m.num_faces(), m.half_edge_to_face())):
            labels = sorted({x for x in h2c if x != -1})
            if not list(m.half_edges()):
                # no half-edge to label, though the map still has one cell
                assert labels == [], m
                assert n == 1, m
            else:
                assert labels == list(range(n)), (m, labels, n)


def test_half_edge_to_cell_follows_the_permutations():
    # a half-edge and its image under vp are on the same vertex, and likewise
    # for fp and the faces
    for m in half_edge_to_cell_corpus():
        vp = m.vertex_permutation(copy=False)
        fp = m.face_permutation(copy=False)
        h2v = m.half_edge_to_vertex()
        h2f = m.half_edge_to_face()
        for h in m.half_edges():
            assert h2v[vp[h]] == h2v[h], (m, h)
            assert h2f[fp[h]] == h2f[h], (m, h)


def test_half_edge_to_cell_on_the_empty_map():
    # the empty map has one vertex and one face but no half-edge, so the two
    # arrays are empty rather than of length one
    from combisurf import OrientedMap

    z = OrientedMap()
    assert z.num_vertices() == 1 and z.num_faces() == 1
    assert z.vertices() == [[]] and z.faces() == [[]]
    assert list(z.half_edge_to_vertex()) == []
    assert list(z.half_edge_to_face()) == []


def test_half_edge_to_cell_with_inactive_edge_zero():
    # the labels come from the active cycles only, so index 0 is a real cell
    # even when half-edge 0 is inactive
    from combisurf import OrientedMap

    m = OrientedMap(vp="(2,1,5)(~1,~2,~5)")
    assert list(m.half_edge_to_vertex()) == [-1, -1, 0, 1, 0, 1, -1, -1, -1, -1, 0, 1]
    assert list(m.half_edge_to_face()) == [-1, -1, 0, 1, 1, 2, -1, -1, -1, -1, 2, 0]


def radial_maps():
    from test_fold import sample_maps
    from test_quad_system import degenerate_maps
    return list(small_maps(folded=False)) + sample_maps() + degenerate_maps()


def test_radial_map_punctured_shape():
    # one bigon per punctured face, on the edge 2 n + 2 k, 2 n + 2 k + 1; the
    # quadrilaterals and the walks of length two stay
    import itertools

    for m in radial_maps():
        n = len(m.vertex_permutation(copy=False))
        faces = range(m.num_faces())
        subsets = [list(P) for k in range(len(faces) + 1) for P in itertools.combinations(faces, k)][:16]
        for P in subsets:
            b = len(P)
            r, mor, bigons = m.radial_map(mapping=True, punctured_faces=P)
            assert len(bigons) == b, (m, P)
            assert r.num_vertices() == m.num_vertices() + m.num_faces(), (m, P)
            assert r.num_edges() == 2 * m.num_edges() + b, (m, P)
            assert r.num_faces() == m.num_edges() + b, (m, P)
            assert r.genus(connected=False) == m.genus(connected=False), (m, P)
            assert list(r.half_edges()) == [y for y in range(2 * n + 2 * b) if y >= 2 * n or y // 2 in m.half_edges()]
            for i, face in enumerate(r.faces()):
                assert len(face) == (2 if i in bigons else 4), (m, P, i)
            # each bigon has a side among the new edges
            h2f = r.half_edge_to_face()
            assert sorted(h2f[2 * n + 2 * k] for k in range(b)) == bigons, (m, P)
            # odd half-edges have their tail at a face of m, even ones at a
            # vertex of m
            h2v = r.half_edge_to_vertex()
            white = {h2v[y] for y in r.half_edges() if y % 2}
            black = {h2v[y] for y in r.half_edges() if y % 2 == 0}
            assert not white & black and len(white) == m.num_faces(), (m, P)
            for h in m.half_edges():
                assert len(mor[h]) == 2, (m, P, h)
                r._check_walk(mor[h])
                assert list(mor[h ^ 1]) == [x ^ 1 for x in reversed(mor[h])], (m, P, h)


def test_radial_map_punctured_unchanged():
    from combisurf import OrientedMap

    for m in radial_maps():
        r0, mor0 = m.radial_map(mapping=True)
        assert m.radial_map(mapping=True, punctured_faces=False) == (r0, mor0, [])
        assert m.radial_map(punctured_faces=[]) == (r0, [])
        r1, mor1, _ = m.radial_map(mapping=True, punctured_faces=True)
        assert not r1.is_mutable()
        assert m.radial_map(punctured_faces=True, mutable=True)[0].is_mutable()
        # the walks are unchanged, but for those that take a new edge
        n = len(m.vertex_permutation(copy=False))
        assert all(mor1[h] == mor0[h] or max(mor1[h]) >= 2 * n for h in m.half_edges())

    m = OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~4,~5)(~2,5,~3)")
    with pytest.raises(ValueError, match="out of range"):
        m.radial_map(punctured_faces=[-1])
    with pytest.raises(TypeError):
        m.radial_map(punctured_faces=0)
    with pytest.raises(NotImplementedError):
        OrientedMap("(0)").radial_map(punctured_faces=True)


def test_half_edge_to_cell_with_a_folded_edge():
    # only the even half-edge of a folded edge is active, so the odd one is
    # labelled -1 while the even one belongs to a genuine cell
    from combisurf import OrientedMap

    m = OrientedMap("(0,2,~2)")
    assert m.has_folded_edge()
    vp = m.vertex_permutation(copy=False)
    h2v = m.half_edge_to_vertex()
    h2f = m.half_edge_to_face()
    folded = [e for e in m.edge_indices() if vp[(2 * e) ^ 1] == -1]
    assert folded
    for e in folded:
        assert h2v[2 * e] != -1 and h2f[2 * e] != -1, (m, e)
        assert h2v[2 * e + 1] == -1 and h2f[2 * e + 1] == -1, (m, e)
