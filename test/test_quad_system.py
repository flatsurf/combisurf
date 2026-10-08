import pytest

from test_fold import find_isomorphism, sample_maps


def all_isomorphisms(A, B):
    avp, afp = A.vertex_permutation(copy=False), A.face_permutation(copy=False)
    bvp, bfp = B.vertex_permutation(copy=False), B.face_permutation(copy=False)
    AH, BH = list(A.half_edges()), list(B.half_edges())
    out = []
    if len(AH) != len(BH):
        return out
    for k0 in BH:
        h0 = AH[0]
        f = {h0: k0}
        todo = [h0]
        ok = True
        while todo and ok:
            x = todo.pop()
            for xx, kk in ((avp[x], bvp[f[x]]), (afp[x], bfp[f[x]]), (x ^ 1, f[x] ^ 1)):
                if xx in f:
                    if f[xx] != kk:
                        ok = False
                        break
                else:
                    f[xx] = kk
                    todo.append(xx)
        if ok and len(f) == len(AH) and len(set(f.values())) == len(AH):
            out.append(f)
    return out


def degenerate_maps():
    # maps having a quadrilateral two of whose sides are carried by the same
    # edge of the radial map. This happens exactly at a leaf edge, where
    # fp[h] == h ^ 1, and at a monogon face, where fp[h] == h; none of the maps
    # of sample_maps() has either.
    from combisurf import OrientedMap
    return [
        OrientedMap(vp="(0,1,~0,~1,2)(~2)"),       # a leaf edge
        OrientedMap(vp="(0,1,~0,2,~2,~1)"),        # a monogon face
        OrientedMap(vp="(0,1,~0,~1,2,3,~3)(~2)"),  # both
    ]


def quad_system_maps():
    from combisurf import OrientedMap
    # the last one has a deleted edge whose image is empty, so that a
    # contracted edge and a monogon face are not the only ways to get one
    return sample_maps() + degenerate_maps() + [
        OrientedMap("(0,~2,4,~1,~0,3)(1,~4)(2,~3)",
                    "(0,~1,~4,~2,~3,~0,3,2)(1,4)"),
    ]


def test_degenerate_maps_are_degenerate():
    # guards degenerate_maps() against silently losing what it is there for
    for m in degenerate_maps():
        fp = m.face_permutation(copy=False)
        assert m.is_connected() and m.genus() >= 1, m
        assert any(fp[h] == h or fp[h] == h ^ 1 for h in m.half_edges()), m


def test_quad_system_shape():
    # a spanning forest and coforest give a quadrangulation with two vertices
    # of degree 4g, 4g edges and 2g quadrilateral faces
    for m in quad_system_maps():
        if not m.is_connected():
            continue
        g = m.genus()
        if g < 1:
            continue
        q, _ = m.quad_system()
        assert q.num_vertices() == 2, m
        assert q.num_edges() == 4 * g, m
        assert q.num_faces() == 2 * g, m
        assert q.vertex_profile() == [4 * g, 4 * g], m
        assert set(q.face_profile()) == {4}, m
        assert q.genus() == g, m


def test_quad_system_projection():
    for m in quad_system_maps():
        if not m.is_connected() or m.genus() < 1:
            continue
        fp = m.face_permutation(copy=False)
        forest, coforest, _ = m.forest_coforest_decomposition()
        contracted = {h // 2 for h in forest if h != -1}
        q, proj, _ = m.quad_system(forest, coforest, mapping=True)

        for h in m.half_edges():
            p = proj[h]
            # the image is a walk of length two unless the folds identify its
            # two half-edges. That always happens for a contracted edge and
            # for a monogon face, and it happens for some deleted edges too,
            # so the converse does not hold.
            assert len(p) in (0, 2), (m, h)
            if h // 2 in contracted or fp[h] == h or fp[h ^ 1] == h ^ 1:
                assert len(p) == 0, (m, h)
            # the walk lives in the quad system
            assert all(x in list(q.half_edges()) for x in p), (m, h)
        # the two half-edges of an edge give reverse walks
        for h in m.half_edges():
            p, pp = proj[h], proj[h ^ 1]
            if len(p) == 2:
                assert list(pp) == [p[1] ^ 1, p[0] ^ 1], (m, h)


def test_quad_system_matches_QuadSystem():
    # fed the decomposition that QuadSystem uses, the two agree, projection
    # included
    from combisurf import QuadSystem
    from combisurf.quad_systems import tree_co_tree

    tested = 0
    for m in sample_maps():
        if not m.is_connected() or m.genus() < 2:
            continue
        # tree_co_tree does not support maps with inactive half-edges
        if list(m.half_edges()) != list(range(2 * m.num_edges())):
            continue
        code = tree_co_tree(m)
        forest = [2 * e for e in m.edge_indices() if code[e] == 0]
        coforest = [2 * e for e in m.edge_indices() if code[e] == 1]
        q, proj, _ = m.quad_system(forest, coforest, mapping=True)
        Q = QuadSystem(m)

        # QuadSystem._proj is not orientation consistent: _proj[ep(h)] is not
        # the reverse of _proj[h]. The two agree on the even half-edges, which
        # already pins the quadrangulation and the projection.
        carried = [f for f in all_isomorphisms(q, Q._quad)
                   if all((not proj[h] and not Q._proj[h])
                          or [f[x] for x in proj[h]] == list(Q._proj[h])
                          for h in m.half_edges() if h % 2 == 0)]
        assert carried, m
        tested += 1
    assert tested


def test_quad_system_genus_zero():
    # in genus zero the forest and the coforest are spanning, since Euler
    # gives ne = (nv - 1) + (nf - 1) and leaves no complementary edge, so
    # everything collapses: the quad system is the empty map and every
    # half-edge projects to the empty walk
    from combisurf import OrientedMap

    maps = [m for m in quad_system_maps() if m.is_connected() and m.genus() == 0]
    maps += [OrientedMap(vp="(0,~0)"),
             OrientedMap(vp="(0,1,~1,~0)"),
             OrientedMap(vp="(0,1,2)(~0)(~1)(~2)")]
    assert maps
    for m in maps:
        forest, coforest, comp = m.forest_coforest_decomposition()
        assert not list(comp), m
        q, proj, _ = m.quad_system(forest, coforest, mapping=True)
        assert q == OrientedMap("", ""), m
        assert q.num_edges() == 0, m
        assert all(p is None or len(p) == 0 for p in proj), m

    # fed a coforest that is not spanning, the same map gives a genuine
    # quadrangulation instead
    m = OrientedMap(vp="(0,~0)")
    forest, coforest, _ = m.forest_coforest_decomposition()
    q, _ = m.quad_system(forest, [])
    assert q == m.radial_map()[0], m
    assert set(q.face_profile()) == {4}, m
    assert q.genus() == 0, m


def test_quad_system_arguments():
    from combisurf import OrientedMap

    m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
    forest, coforest, _ = m.forest_coforest_decomposition()
    with pytest.raises(ValueError, match="together"):
        m.quad_system(forest)
    with pytest.raises(ValueError, match="together"):
        m.quad_system(coforest=coforest)

    assert not m.quad_system()[0].is_mutable()
    assert m.quad_system(mutable=True)[0].is_mutable()

    h = next(x for x in forest if x != -1)
    with pytest.raises(ValueError, match="listed twice"):
        m.quad_system(list(forest) + [h], coforest)
    # an edge counts as listed whichever of its half-edges names it
    with pytest.raises(ValueError, match="listed twice"):
        m.quad_system(forest, list(coforest) + [h ^ 1])

    f = OrientedMap("(0,2,~2)")
    assert f.has_folded_edge()
    with pytest.raises(NotImplementedError):
        f.quad_system()


def test_quad_system_relabel():
    for m in quad_system_maps():
        if not m.is_connected() or m.genus() < 1:
            continue
        forest, coforest, _ = m.forest_coforest_decomposition()
        plain, proj, _ = m.quad_system(forest, coforest, mapping=True)
        compact, proj2, _ = m.quad_system(forest, coforest, relabel=True, mapping=True)

        # no inactive half-edge left
        assert list(compact.half_edges()) == list(range(2 * compact.num_edges())), m

        # the two describe the same quadrangulation
        assert find_isomorphism(plain, compact) is not None, m

        # the two projections differ exactly by the compaction, which goes
        # edge by edge and keeps the parity inside each edge
        survivors = sorted({y // 2 for y in plain.half_edges()})
        relabelling = {}
        for i, e in enumerate(survivors):
            relabelling[2 * e] = 2 * i
            relabelling[2 * e + 1] = 2 * i + 1
        for h in m.half_edges():
            assert list(proj2[h]) == [relabelling[x] for x in proj[h]], (m, h)

        # the parity of the labels still gives the bipartition of the vertices
        h2v = compact.half_edge_to_vertex()
        assert len({h2v[y] for y in compact.half_edges() if y % 2 == 0}) == 1, m
        assert len({h2v[y] for y in compact.half_edges() if y % 2 == 1}) == 1, m

        # and the walks still reverse properly
        for h in m.half_edges():
            if proj2[h]:
                assert list(proj2[h ^ 1]) == [proj2[h][1] ^ 1, proj2[h][0] ^ 1], (m, h)


def punctured_maps():
    r"""
    Return the list of triples ``(name, map, sets of punctured faces)``.

    Each set ``P`` satisfies ``χ - |P| < 0``, so that the pairing applies on
    the punctured surface; the first one is the set of all faces.
    """
    import json
    import os
    from combisurf import OrientedMap

    data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
    out = [
        ("sphere4", OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~4,~5)(~2,5,~3)"), [[0, 1, 2, 3], [1, 2, 3], [0, 1, 3]]),
        ("torus", OrientedMap(vp=[[0, 2, 1, 3]]), [[0]]),
        ("genus1_two_faces", OrientedMap(vp=[[0, 2, 4, 1, 3, 5]]), [[0, 1], [0], [1]]),
        ("genus2", OrientedMap(vp=[[0, 2, 1, 3, 4, 6, 5, 7]]), [[0]]),
        ("pants", OrientedMap(fp="(0)(1)(~0,~1)"), [[0, 1, 2]]),
    ]
    for name, subsets in (("lshape", [[1], [0, 4], [2, 3, 5]]),
                          ("emmw_1117", [[0], [3, 7, 11, 19]])):
        with open(os.path.join(data, f"cylinders_{name}.json")) as f:
            m = OrientedMap(vp=json.load(f)["vp"])
        out.append((name, m, [list(range(m.num_faces()))] + subsets))
    return out


def test_punctured_maps():
    # guards punctured_maps() against an instance the pairing would reject
    for name, m, subsets in punctured_maps():
        assert subsets[0] == list(range(m.num_faces())), name
        for P in subsets:
            assert m.euler_characteristic() - len(P) < 0, (name, P)


def test_quad_system_punctured_shape():
    # spanning forest and coforest: 2 vertices, 2g quadrilaterals and one
    # bigon per punctured face, hence 4g + b edges
    for name, m, subsets in punctured_maps():
        g = m.genus()
        for P in subsets:
            b = len(P)
            for relabel in (False, True):
                q, proj, bigons = m.quad_system(relabel=relabel, mapping=True, punctured_faces=P)
                assert (q.num_vertices(), q.num_edges(), q.num_faces()) == (2, 4 * g + b, 2 * g + b), (name, P)
                assert q.genus() == g, (name, P)
                assert len(bigons) == b, (name, P)
                faces = q.faces()
                for i, face in enumerate(faces):
                    assert len(face) == (2 if i in bigons else 4), (name, P, i)
                for h in m.half_edges():
                    assert len(proj[h]) in (0, 2), (name, P, h)
                    q._check_walk(proj[h])
                    assert list(proj[h ^ 1]) == [x ^ 1 for x in reversed(proj[h])], (name, P, h)
                if relabel:
                    # the parity of the labels gives the bipartition
                    h2v = q.half_edge_to_vertex()
                    assert len({h2v[y] for y in q.half_edges() if y % 2 == 0}) == 1, (name, P)
                    assert len({h2v[y] for y in q.half_edges() if y % 2 == 1}) == 1, (name, P)


def test_quad_system_no_puncture():
    # no punctured face is the default, and gives no bigon
    for name, m, subsets in punctured_maps():
        forest, coforest, _ = m.forest_coforest_decomposition()
        q0, proj0, bigons = m.quad_system(forest, coforest, mapping=True)
        assert bigons == [], name
        assert m.quad_system(forest, coforest, mapping=True, punctured_faces=False) == (q0, proj0, []), name
        assert m.quad_system(punctured_faces=[]) == m.quad_system(), name


def test_quad_system_punctured_arguments():
    from combisurf import OrientedMap

    m = OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~4,~5)(~2,5,~3)")
    assert m.quad_system(punctured_faces=[2, 0, 2])[1] == m.quad_system(punctured_faces=[0, 2])[1]
    with pytest.raises(ValueError, match="out of range"):
        m.quad_system(punctured_faces=[4])
    with pytest.raises(TypeError):
        m.quad_system(punctured_faces=1)
    with pytest.raises(NotImplementedError):
        OrientedMap("(0,2,~2)").quad_system(punctured_faces=True)


def image(mor, w):
    return [x for h in w for x in mor[h]]


def free_word(r, bigons, w):
    r"""
    Return the cyclically reduced word of the closed walk ``w`` of ``r`` in the
    free fundamental group of ``r`` punctured at the faces ``bigons``.
    """
    from combisurf.word import word_init, word_cyclically_reduce
    forest, coforest, _ = r.forest_coforest_decomposition(root_faces=bigons)
    _, proj = r.reduced_map(forest, coforest, mapping=True)
    return list(word_cyclically_reduce(word_init(image(proj, w))))


def assert_peripheral(m, P, r, mor, bigons):
    # the boundary of each punctured face, which has the face on its left,
    # goes to a curve freely homotopic to the boundary of a bigon, which has
    # the bigon on its left: the puncture stays on the same side
    from test_geometric_intersection_oracles import is_cyclic_conjugate

    peripheral = [free_word(r, bigons, r.faces()[i]) for i in bigons]
    faces = m.faces()
    for f in P:
        u = free_word(r, bigons, image(mor, faces[f]))
        assert u, (m, P, f)
        assert any(is_cyclic_conjugate(u, v) for v in peripheral), (m, P, f)


def test_punctured_peripheral():
    for name, m, subsets in punctured_maps():
        for P in subsets:
            r, mor, bigons = m.radial_map(mapping=True, punctured_faces=P)
            assert_peripheral(m, P, r, mor, bigons)
            q, proj, bigons = m.quad_system(mapping=True, punctured_faces=P)
            assert_peripheral(m, P, q, proj, bigons)


def pairing_instances():
    r"""
    Return the list of triples ``(map, sets of punctured faces, curves)`` on
    which the pairing is compared before and after the reduction.
    """
    import random
    from combisurf import OrientedMap
    from test_geometric_intersection_oracles import (COHEN_LUSTIG_ROWS, cylinder_data, genus1_two_faces, genus2,
                                                     pants, random_closed_walk, random_word, sphere4, torus)

    def rows(surface):
        # the words V and W of the rows of Cohen–Lustig on that surface
        return [w for row in COHEN_LUSTIG_ROWS if row[1] is surface for w in row[2:4] if w is not None]

    rng = random.Random(0)
    out = [(torus(), [[0]], rows(torus)),
           (genus2(), [[0]], rows(genus2)),
           (genus1_two_faces(), [[0, 1], [0], [1]], rows(genus1_two_faces))]
    out.append((pants(), [[0, 1, 2]], [random_word(rng, 2, 8) for _ in range(6)]))
    m = sphere4()
    eight01 = [0, 2, 4, 0, 9, 7]
    out.append((m, [[0, 1, 2, 3], [0, 1, 3]],
                [eight01] + [random_closed_walk(m, rng, length) for length in (6, 9, 12, 15)]))
    m = OrientedMap(vp=cylinder_data("lshape")["vp"])
    out.append((m, [[0, 1, 2, 3, 4, 5], [1], [0, 4], [2, 3, 5]],
                [random_closed_walk(m, rng, length) for length in (3, 5, 7, 9, 12)]))
    return out


def test_pairing_instances():
    # guards pairing_instances() against an instance that tests nothing
    instances = pairing_instances()
    assert [m.genus() for m, _, _ in instances] == [1, 2, 1, 0, 0, 2]
    for m, subsets, curves in instances:
        assert len(curves) >= 3
        for w in curves:
            m._check_walk(w, closed=True)


@pytest.mark.parametrize("reduction", ["radial_map", "quad_system"])
def test_punctured_geometric_intersection(reduction):
    # the reduction to a map with bigons preserves the geometric intersection
    # numbers of the punctured surface
    from combisurf.geometric_intersection import GeometricIntersectionPairing

    # GeometricIntersectionPairing rejects a map with inactive half-edges,
    # which the quad system has unless relabelled
    kwds = {"relabel": True} if reduction == "quad_system" else {}
    for m, subsets, curves in pairing_instances():
        for P in subsets:
            r, mor, bigons = getattr(m, reduction)(mapping=True, punctured_faces=P, **kwds)
            gi = GeometricIntersectionPairing(m, P)
            gr = GeometricIntersectionPairing(r, bigons)
            images = [image(mor, w) for w in curves]
            for u, v in zip(curves, images):
                assert gi([u]) == gr([v]), (m, P, u)
            for i in range(len(curves)):
                for j in range(i, len(curves)):
                    assert gi([curves[i]], [curves[j]]) == gr([images[i]], [images[j]]), (m, P, i, j)
