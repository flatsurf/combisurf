import pytest

from test_fold import sample_maps


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
        q = m.quad_system()
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
        q, proj = m.quad_system(forest, coforest, mapping=True)

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
        q, proj = m.quad_system(forest, coforest, mapping=True)
        Q = QuadSystem(m)

        # QuadSystem._proj is not orientation consistent: _proj[ep(h)] is not
        # the reverse of _proj[h]. The two agree on the even half-edges, which
        # already pins the quadrangulation and the projection.
        # every isomorphism from q to Q._quad is one of them composed with an
        # automorphism of q
        ok, r = q.is_isomorphic(Q._quad, certificate=True)
        isomorphisms = [[r[a[x]] for x in range(len(r))] for a in q.automorphisms()] if ok else []
        carried = [f for f in isomorphisms
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
        q, proj = m.quad_system(forest, coforest, mapping=True)
        assert q == OrientedMap("", ""), m
        assert q.num_edges() == 0, m
        assert all(p is None or len(p) == 0 for p in proj), m

    # fed a coforest that is not spanning, the same map gives a genuine
    # quadrangulation instead
    m = OrientedMap(vp="(0,~0)")
    forest, coforest, _ = m.forest_coforest_decomposition()
    q = m.quad_system(forest, [])
    assert q == m.radial_map(), m
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

    assert not m.quad_system().is_mutable()
    assert m.quad_system(mutable=True).is_mutable()

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
        plain, proj = m.quad_system(forest, coforest, mapping=True)
        compact, proj2 = m.quad_system(forest, coforest, relabel=True, mapping=True)

        # no inactive half-edge left
        assert list(compact.half_edges()) == list(range(2 * compact.num_edges())), m

        # the two describe the same quadrangulation
        assert plain.is_isomorphic(compact), m

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
