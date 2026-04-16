import pytest

def test_constructor():
    from combisurf import OrientedMap

    OrientedMap([2, 1, 3, 0], None)
    OrientedMap(None, [2, 1, 3, 0])
    OrientedMap([0, -1, 3, 4, 5, 2])
    OrientedMap(None, [0, -1, 3, 4, 5, 2])

    OrientedMap([0, -1, 3, 4, 5, 2], None)
    OrientedMap(None, [0, -1, 3, 4, 5, 2])

    OrientedMap("(0,1,2)(~0,~1,~2)", None)
    OrientedMap(None, "(0,1,2)(~0,~1,~2)")

    with pytest.raises(ValueError):
        OrientedMap([0, 0], None)

    with pytest.raises(ValueError):
        OrientedMap(None, [0, 0])


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
