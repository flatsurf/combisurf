import pytest


def sample_maps():
    from combisurf import OrientedMap
    return [
        OrientedMap(vp="(0,1,~0,2)(~1,~2)"),
        OrientedMap(vp="(0,1,2)(~0,~1,~2)"),
        OrientedMap(vp="(0,~0,1,~1)"),
        OrientedMap(vp="(0,1,~1,~0)"),
        OrientedMap(vp="(2,3,1,5,6,~5)(~1,4,~2,~6)(~3,7,~4,8)(~7,9,10,11,~9)(~8,~10,~11)"),
        OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]]),
        OrientedMap(vp=[[0, 2, 4, 6], [7, 8, 5], [9, 10, 12, 11], [3, 15, 1, 13, 14]]),
    ]



def ep(m, h):
    r"""
    Return the image of ``h`` under the edge permutation of ``m``.

    This is ``h ^ 1``, except on a folded edge where it is ``h`` itself.
    """
    vp = m.vertex_permutation(copy=False)
    return h if vp[h ^ 1] == -1 else h ^ 1


def check_fold_corner_counts(m, h, r):
    r"""
    Check the four counts of ``r``, the result of ``m.fold_corner(h)``.

    What happens is decided by the two ends the fold merges, the head of
    ``fp[h]`` and the tail of ``h``. Note that the head of a half-edge is the
    tail of its image under ``ep``, which is not ``^ 1`` on a folded edge.
    """
    fp = m.face_permutation(copy=False)
    h2v = m.half_edge_to_vertex()
    de = r.num_edges() - m.num_edges()
    dv = r.num_vertices() - m.num_vertices()
    df = r.num_faces() - m.num_faces()
    dfolded = r.num_folded_edges() - m.num_folded_edges()

    if fp[h] == h:
        # a monogon glues the edge to itself, so it survives as a folded edge
        # and its face closes up
        assert (de, dfolded) == (0, 1), (m, h)
        assert (dv, df) == (0, -1), (m, h)
        return

    assert (de, dfolded) == (-1, 0), (m, h)

    if fp[h] == ep(m, h):
        # a leaf half-edge: the edge is pruned and its head, a vertex of
        # degree one, goes with it
        assert (dv, df) == (-1, 0), (m, h)
    elif h2v[h] != h2v[ep(m, fp[h])]:
        # two distinct vertices, which merge
        assert (dv, df) == (-1, 0), (m, h)
    elif fp[fp[h]] == h:
        # the same vertex, and the face of h is the bigon (h, fp[h]), which
        # closes up
        assert (dv, df) == (0, -1), (m, h)
    else:
        # the same vertex twice over: it is pinched and splits in two. As for
        # contract_edge and delete_edge the degenerate case is performed
        # rather than refused, so the euler characteristic jumps by two
        assert (dv, df) == (1, 0), (m, h)
        assert r.euler_characteristic() == m.euler_characteristic() + 2, (m, h)


def check_fold_half_edge_counts(m, h, r):
    r"""
    Check the four counts of ``r``, the result of ``m.fold_half_edge(h)``.
    """
    fp = m.face_permutation(copy=False)
    h2v = m.half_edge_to_vertex()
    de = r.num_edges() - m.num_edges()
    dv = r.num_vertices() - m.num_vertices()
    df = r.num_faces() - m.num_faces()
    dfolded = r.num_folded_edges() - m.num_folded_edges()

    # the edge survives, folded onto itself, keeping its even half-edge
    assert (de, dfolded) == (0, 1), (m, h)
    assert 2 * (h // 2) in list(r.folded_half_edges()), (m, h)

    if h2v[h] != h2v[ep(m, h)]:
        # the two ends of the edge are distinct and merge
        assert (dv, df) == (-1, 0), (m, h)
    elif fp[h] == h:
        # a loop already bounding a monogon: that face closes up
        assert (dv, df) == (0, -1), (m, h)
    else:
        # a loop otherwise: its vertex is pinched and splits in two, and the
        # degenerate case is performed rather than refused
        assert (dv, df) == (1, 0), (m, h)
        assert r.euler_characteristic() == m.euler_characteristic() + 2, (m, h)

    # chi = F - E + (V + folded)
    assert (r.euler_characteristic() - m.euler_characteristic()
            == df - de + dv + dfolded), (m, h)


def test_fold_corner_twice_is_contract_edge():
    # on a radial map, folding the two corners of the quadrilateral of an edge
    # along the odd half-edge contracts that edge
    done = [0]
    for m in sample_maps():
        R = m.radial_map()
        rfp = R.face_permutation(copy=False)
        for e in m.edge_indices():
            h = 2 * e
            if h not in list(m.half_edges()):
                continue
            mm = m.copy(mutable=True)
            mm.contract_edge(e)
            x0 = 2 * h + 1
            x1 = rfp[x0]
            x3 = rfp[rfp[x1]]
            r = R.copy(mutable=True)
            try:
                r.fold_corner(x3)
                r.fold_corner(x1)
            except ValueError:
                continue          # degenerate quadrilateral
            assert mm.radial_map().is_isomorphic(r), (m, e)
            done[0] += 1
    assert done[0] > 20


def test_fold_corner_twice_is_delete_edge():
    # ... and along the even half-edge it deletes it
    done = [0]
    for m in sample_maps():
        R = m.radial_map()
        fp = m.face_permutation(copy=False)
        rfp = R.face_permutation(copy=False)
        for e in m.edge_indices():
            h = 2 * e
            if h not in list(m.half_edges()):
                continue
            mm = m.copy(mutable=True)
            mm.delete_edge(e)
            x0 = 2 * fp[h]
            x1 = rfp[x0]
            x3 = rfp[rfp[x1]]
            r = R.copy(mutable=True)
            try:
                r.fold_corner(x3)
                r.fold_corner(x1)
            except ValueError:
                continue          # degenerate quadrilateral
            assert mm.radial_map().is_isomorphic(r), (m, e)
            done[0] += 1
    assert done[0] > 20


def test_fold_corner_effect():
    for m in sample_maps():
        for h in m.half_edges():
            r = m.copy(mutable=True)
            try:
                r.fold_corner(h)
            except ValueError:
                continue
            r._check()
            check_fold_corner_counts(m, h, r)


def test_fold_half_edge_effect():
    for m in sample_maps():
        assert m.num_folded_edges() == 0, m
        for h in m.half_edges():
            r = m.copy(mutable=True)
            r.fold_half_edge(h)
            r._check()
            assert r.num_folded_edges() == 1, (m, h)
            check_fold_half_edge_counts(m, h, r)


def test_fold_corner_prunes_a_leaf():
    from combisurf import OrientedMap

    # the opposite half-edge follows h in its face: the corner is bounded by
    # the edge of h on both sides and folding prunes it
    for m in sample_maps():
        R = m.radial_map()
        rvp = R.vertex_permutation(copy=False)
        rfp = R.face_permutation(copy=False)
        leaves = [h for h in R.half_edges() if rfp[h] == h ^ 1]
        for h in leaves:
            # the fev relation makes the head of h a vertex of degree one
            assert rvp[h ^ 1] == h ^ 1, (m, h)
            r = R.copy(mutable=True)
            r.fold_corner(h)
            r._check()
            # the same effect as any other fold: one edge less, the face of h
            # two shorter, and the two ends of the edge merged
            assert r.num_edges() == R.num_edges() - 1, (m, h)
            assert r.num_vertices() == R.num_vertices() - 1, (m, h)
            assert r.num_faces() == R.num_faces(), (m, h)
            assert r.euler_characteristic() == R.euler_characteristic(), (m, h)
            assert not r.has_folded_edge(), (m, h)

    # folding the two corners of the single quadrilateral of a one edge map
    # empties it
    r = OrientedMap(fp="(0,~0,1,~1)", mutable=True)
    r.fold_corner(1)
    r.fold_corner(3)
    r._check()
    assert r == OrientedMap("", "")


def test_fold_corner_folds_a_monogon():
    from combisurf import OrientedMap

    # the face of h is a monogon: the rule identifies h with ep(h), so the
    # edge is glued to itself and stays on as a folded edge
    tested = 0
    for m in sample_maps():
        vp = m.vertex_permutation(copy=False)
        fp = m.face_permutation(copy=False)
        for h in m.half_edges():
            if fp[h] != h:
                continue
            tested += 1
            # the fev relation puts h and ep(h) next to each other at the tail
            assert vp[h] == h ^ 1, (m, h)
            r = m.copy(mutable=True)
            r.fold_corner(h)
            r._check()
            assert r.num_edges() == m.num_edges(), (m, h)
            assert r.num_folded_edges() == m.num_folded_edges() + 1, (m, h)
            assert r.num_faces() == m.num_faces() - 1, (m, h)
            assert r.num_vertices() == m.num_vertices(), (m, h)
            assert r.euler_characteristic() == m.euler_characteristic(), (m, h)
            # which is exactly folding the edge outright
            s = m.copy(mutable=True)
            s.fold_half_edge(h)
            assert r == s, (m, h)
    assert tested


def test_fold_corner_errors():
    from combisurf import OrientedMap

    # folded edges are not supported by fold_corner
    f = OrientedMap("(0,2,~2)", mutable=True)
    assert f.has_folded_edge()
    with pytest.raises(NotImplementedError):
        f.fold_corner(0)


def test_fold_half_edge_twice():
    from combisurf import OrientedMap

    # an edge can not be folded twice
    for m in sample_maps():
        for h in m.half_edges():
            r = m.copy(mutable=True)
            r.fold_half_edge(h)
            for x in (h, h ^ 1):
                if x in list(r.half_edges()):
                    with pytest.raises(ValueError, match="already folded"):
                        r.fold_half_edge(x)


def test_fold_requires_mutable():
    for m in sample_maps():
        h = list(m.half_edges())[0]
        with pytest.raises(ValueError):
            m.fold_corner(h)
        with pytest.raises(ValueError):
            m.fold_half_edge(h)


def maps_with_a_folded_edge():
    # folding one edge of each sample map gives maps in which the folding
    # primitives meet a neighbour lying on a folded edge, where ep(x) is x and
    # not x ^ 1
    out = []
    for m in sample_maps():
        for h in m.half_edges():
            r = m.copy(mutable=True)
            r.fold_half_edge(h)
            out.append(r)
    return out


def test_fold_corner_next_to_a_folded_edge():
    from combisurf import OrientedMap

    m = OrientedMap("(0,~0,3,~1,2)(1,~3)", "(0)(~0,2,~1,~3)(1,3)")
    r = m.copy(mutable=True)
    r.fold_corner(2)
    r._check()
    check_fold_corner_counts(m, 2, r)

    tested = 0
    for m in maps_with_a_folded_edge():
        assert m.has_folded_edge(), m
        for h in list(m.half_edges()):
            r = m.copy(mutable=True)
            try:
                r.fold_corner(h)
            except (ValueError, NotImplementedError):
                continue
            r._check()
            check_fold_corner_counts(m, h, r)
            tested += 1
    assert tested > 500, tested


def test_fold_half_edge_next_to_a_folded_edge():
    from combisurf import OrientedMap

    m = OrientedMap("(0,~0,2)(1,~1)", "(0)(~0,2)(1)(~1)")
    assert m.num_folded_edges() == 1
    r = m.copy(mutable=True)
    r.fold_half_edge(0)
    r._check()
    check_fold_half_edge_counts(m, 0, r)

    tested = 0
    for m in maps_with_a_folded_edge():
        assert m.has_folded_edge(), m
        for h in list(m.half_edges()):
            r = m.copy(mutable=True)
            try:
                r.fold_half_edge(h)
            except (ValueError, NotImplementedError):
                continue
            r._check()
            check_fold_half_edge_counts(m, h, r)
            tested += 1
    assert tested > 500, tested
