r"""
Triangle systems: strands, sizes, colours and images of the half-edges, on the
maps of genus 2 to 4 of the barycentric subdivision tests, each with several
random choices of tree, cotree and cocotree. The homotopy classes of closed
walks are compared through the hyperbolic quadrangulations of the quad systems
of the map and of its triangle system. With punctured faces, on these maps and
on the maps of genus 0 and 1 of the same tests, the images of the face
boundaries and of random closed walks are compared through the strongly
reduced closed walks of the triangle system.
"""
import random

import pytest

from combisurf import OrientedMap
from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation

from test_barycentric_subdivision import is_rotation, one_vertex_map, random_maps

CHOICES = 3


def hyperbolic_maps():
    r"""
    Return the list of the pairs ``(name, m)`` of the maps of genus 2 to 4.
    """
    maps = [(f"one_vertex_{g}", one_vertex_map(g).copy(mutable=False)) for g in range(2, 5)]
    maps.extend((f"random_{i}", m) for i, m in enumerate(random_maps()) if m.genus() >= 2)
    return maps


MAPS = hyperbolic_maps()


def random_strand_tree(strands, rng, first=None):
    r"""
    Return a random spanning tree of the graph of strands, as the list of the
    first half-edges of its strands, holding the strand starting at ``first``
    if given and not a loop of the graph.
    """
    strands = list(strands)
    rng.shuffle(strands)
    if first is not None:
        strands.sort(key=lambda s: s[0] != first)
    root = {}

    def find(e):
        root.setdefault(e, e)
        while root[e] != e:
            e = root[e]
        return e

    tree = []
    for h, hh in strands:
        u = find(h // 2)
        v = find(hh // 2)
        if u != v:
            root[u] = v
            tree.append(h)
    return tree


def choices(m, seed, punctured_faces=False):
    r"""
    Return ``CHOICES`` triples ``(tree, cotree, cocotree)`` for ``m``, from
    random roots and random spanning trees of the graph of strands, the one
    of the map with a loop in each face of ``punctured_faces``.
    """
    rng = random.Random(seed)
    result = []
    for _ in range(CHOICES):
        tree, cotree, _ = m.forest_coforest_decomposition((rng.randrange(m.num_vertices()),),
                                                          (rng.randrange(m.num_faces()),))
        cocotree = random_strand_tree(m.strands(tree, cotree, punctured_faces=punctured_faces), rng)
        result.append((tree, cotree, cocotree))
    return result


CASES = [(name, m, c) for i, (name, m) in enumerate(MAPS) for c in choices(m, i)]
IDS = [f"{name}-{j}" for name, _ in MAPS for j in range(CHOICES)]


def punctured_sets(m):
    r"""
    Return the sets of punctured faces tested on ``m``: the face ``0``, every
    other face and every face, without repetition, keeping those that leave a
    surface of negative Euler characteristic once removed.
    """
    nf = m.num_faces()
    g = m.genus()
    result = []
    for s in ([0], list(range(0, nf, 2)), list(range(nf))):
        if s not in result and 2 - 2 * g - len(s) < 0:
            result.append(s)
    return result


TETRAHEDRON = OrientedMap(fp="(0,1,2)(~0,4,~3)(~1,3,~5)(~2,5,~4)")
TORUS = OrientedMap(vp="(0,1,~0,~1)")


def punctured_maps():
    r"""
    Return the list of the triples ``(name, m, sets)`` of the maps tested with
    punctured faces and their sets of punctured faces: the maps of genus 2 to
    4, the random maps of genus 0 and 1, the tetrahedron and the torus with
    one face.
    """
    maps = [(name, m, punctured_sets(m)) for name, m in MAPS]
    maps.extend((f"random_{i}", m, punctured_sets(m)) for i, m in enumerate(random_maps()) if m.genus() <= 1)
    maps.append(("tetrahedron", TETRAHEDRON, [[0, 1, 2], [0, 1, 2, 3]]))
    maps.append(("torus", TORUS, [[0]]))
    return [(name, m, sets) for name, m, sets in maps if sets]


PUNCTURED_MAPS = punctured_maps()

# the sets of punctured faces of a map have distinct sizes, which name them
PUNCTURED_CASES = [(f"{name}-b{len(p)}", m, p, c)
                   for i, (name, m, sets) in enumerate(PUNCTURED_MAPS) for p in sets
                   for c in choices(m, i, p)]
PUNCTURED_IDS = [f"{name}-{j % CHOICES}" for j, (name, _, _, _) in enumerate(PUNCTURED_CASES)]


def test_number_of_maps():
    assert len(MAPS) == 3 + 18
    sizes = {g: [] for g in range(5)}
    for name, m, sets in PUNCTURED_MAPS:
        if name.startswith("random_") and m.genus() <= 1:
            sizes[m.genus()].append([len(p) for p in sets])
    assert sizes[0] == [[12, 23], [3, 6], [3], [9, 17], [4]]
    assert sizes[1] == [[1, 5, 10], [1, 5, 10], [1, 8, 15], [1], [1, 6, 12], [1, 2, 4]]
    assert sum(len(p) for p in sizes[0] + sizes[1]) == 8 + 16
    assert sum(m.genus() <= 1 for _, m, _, _ in PUNCTURED_CASES) == CHOICES * (8 + 16 + 3)


def random_closed_walk(m, rng, length):
    r"""
    Return a closed walk on ``m`` made of a random walk of ``length``
    half-edges, which backtracks only at a vertex of degree one, followed by
    a shortest path back to its start.
    """
    vp = m.vertex_permutation()
    vertex = m.half_edge_to_vertex()
    half_edges = [h for h in range(len(vp)) if vp[h] != -1]
    out = {}
    for h in half_edges:
        out.setdefault(vertex[h], []).append(h)

    h = rng.choice(half_edges)
    w = [h]
    for _ in range(length - 1):
        h = rng.choice([x for x in out[vertex[h ^ 1]] if x != h ^ 1] or [h ^ 1])
        w.append(h)

    start = vertex[w[0]]
    parent = {vertex[w[-1] ^ 1]: None}
    queue = [vertex[w[-1] ^ 1]]
    for v in queue:
        if v == start:
            break
        for x in out[v]:
            u = vertex[x ^ 1]
            if u not in parent:
                parent[u] = x
                queue.append(u)
    path = []
    v = start
    while parent[v] is not None:
        path.append(parent[v])
        v = vertex[parent[v]]
    w.extend(reversed(path))
    return w


@pytest.mark.parametrize("name, m, choice", CASES, ids=IDS)
def test_strands(name, m, choice):
    tree, cotree, _ = choice
    g = m.genus()
    strands = m.strands(tree, cotree)
    assert len(strands) == 4 * g
    removed = {h // 2 for h in list(tree) + list(cotree) if h != -1}
    remaining = [h for h in m.half_edges() if h // 2 not in removed]
    assert len(remaining) == 4 * g
    assert sorted(h for h, _ in strands) == sorted(remaining)
    assert sorted(hh for _, hh in strands) == sorted(remaining)
    # the vertex permutation of the map with the tree contracted and the
    # cotree deleted, whose labels are those of the remaining edges
    r = m.reduced_map(tree, cotree)
    assert r.num_vertices() == r.num_faces() == 1
    rvp = r.vertex_permutation()
    assert all(rvp[h] == hh for h, hh in strands)


@pytest.mark.parametrize("name, m, choice", CASES, ids=IDS)
def test_sizes(name, m, choice):
    g = m.genus()
    for relabel in (False, True):
        t, punctured = m.triangle_system(*choice, relabel=relabel)
        assert punctured == []
        assert t.num_vertices() == 3
        assert t.num_edges() == 6 * g + 3
        assert t.num_faces() == 4 * g + 2
        assert t.vertex_profile() == [4 * g + 2] * 3
        assert t.face_profile() == [3] * (4 * g + 2)
        assert t.genus() == g
        assert not t.is_mutable()
        BipartiteHyperbolicTriangulation(t)
        if relabel:
            assert list(t.half_edges()) == list(range(2 * (6 * g + 3)))


def image(mor, w):
    return [x for h in w for x in mor[h]]


@pytest.mark.parametrize("name, m, choice", CASES, ids=IDS)
def test_images(name, m, choice):
    tree, cotree, _ = choice
    t, mor, _ = m.triangle_system(*choice, mapping=True)
    tree = {h // 2 for h in tree if h != -1}
    cotree = {h // 2 for h in cotree if h != -1}
    vp = m.vertex_permutation()
    assert len(mor) == len(vp)
    for h in range(len(vp)):
        if vp[h] == -1:
            assert mor[h] is None
            continue
        w = list(mor[h])
        if h // 2 in tree:
            assert w == []
        elif h // 2 in cotree:
            # empty when the edge closed by the tree is contractible
            assert len(w) <= 4
        else:
            assert len(w) == 2
        if w:
            t._check_walk(w)
        assert all(w[i + 1] != w[i] ^ 1 for i in range(len(w) - 1))
        assert list(mor[h ^ 1]) == [x ^ 1 for x in reversed(w)]
    for f in m.faces():
        w = image(mor, f)
        if w:
            t._check_walk(w, closed=True)


@pytest.mark.parametrize("name, m, choice", CASES, ids=IDS)
def test_homotopy(name, m, choice):
    t, mor, _ = m.triangle_system(*choice, mapping=True, relabel=True)
    q, proj, _ = t.quad_system(relabel=True, mapping=True)
    Q = HyperbolicQuadrangulation(q)
    mq, mproj, _ = m.quad_system(relabel=True, mapping=True)
    mQ = HyperbolicQuadrangulation(mq)

    def geodesic(w, canonical=None):
        w = image(proj, image(mor, w))
        return Q.geodesic(w, canonical=canonical) if w else []

    for f in m.faces():
        assert geodesic(f) == []

    # each half-edge followed by the path of the tree back to its tail, a
    # closed walk which is contractible exactly when its image, that of the
    # half-edge, is
    to_root = {v: x for v, x in enumerate(choice[0]) if x != -1}
    h2v = m.half_edge_to_vertex()
    for h in m.half_edges():
        path_up = []
        v = h2v[h ^ 1]
        while v in to_root:
            path_up.append(to_root[v])
            v = h2v[to_root[v] ^ 1]
        path_down = []
        v = h2v[h]
        while v in to_root:
            path_down.append(to_root[v] ^ 1)
            v = h2v[to_root[v] ^ 1]
        loop = [h] + path_up + path_down[::-1]
        m._check_walk(loop, closed=True)
        loop = image(mproj, loop)
        assert (not loop or not mQ.geodesic(loop)) == (not geodesic([h]))

    rng = random.Random(len(m.vertex_permutation()))
    h2f = m.half_edge_to_face()
    faces = m.faces()
    for length in (1, 3, 8, 15):
        w = random_closed_walk(m, rng, length)
        wq = image(mproj, w)
        if not wq or not mQ.geodesic(wq):
            continue
        right = geodesic(w, "right")
        assert right
        # a detour around the face on the left of w[i], inserted before w[i]
        i = rng.randrange(len(w))
        face = faces[h2f[w[i]]]
        j = face.index(w[i])
        detour = w[:i] + face[j:] + face[:j] + w[i:]
        assert is_rotation(geodesic(detour, "right"), right)


def reverse(w):
    return [h ^ 1 for h in reversed(w)]


@pytest.mark.parametrize("name, m, punctured, choice", PUNCTURED_CASES, ids=PUNCTURED_IDS)
def test_punctured_strands(name, m, punctured, choice):
    tree, cotree, _ = choice
    g = m.genus()
    b = len(punctured)
    n = len(m.vertex_permutation())
    strands = m.strands(tree, cotree, punctured_faces=punctured)
    full = m.strands(tree, cotree, punctured_faces=punctured, punctured_flag=False)
    assert len(strands) == 4 * g + b - 1
    assert len(full) == 4 * g + b
    # the default leaves out exactly the strand starting at the first loop
    assert [s for s in full if s[0] != n] == strands
    assert sum(h == n for h, _ in full) == 1
    removed = {h // 2 for h in list(tree) + list(cotree) if h != -1}
    remaining = [h for h in m.half_edges() if h // 2 not in removed]
    loops = [n + 2 * k for k in range(b)]
    assert sorted(h for h, _ in full) == sorted(remaining + loops)
    # each corner (n + 2k + 1, n + 2k) around a monogon is left out, so the
    # half-edges n + 2k end no strand
    assert sorted(hh for _, hh in full) == sorted(remaining + [x + 1 for x in loops])
    # a strand that neither starts nor ends at a loop is a strand without
    # punctured face, and a strand of the closed map is cut by the loops of
    # its corners into strands chained through the monogon corners
    closed = dict(m.strands(tree, cotree))
    assert all(closed[h] == hh for h, hh in full if h < n and hh < n)
    chain = dict(full)
    chain.update((x + 1, x) for x in loops)
    for h, hh in closed.items():
        y = chain[h]
        while y >= n:
            y = chain[y]
        assert y == hh
    # the loop of the k-th punctured face sits at the tail of the first
    # half-edge r of the face: n + 2k + 1 follows r around the vertex
    faces = m.faces()
    for k, f in enumerate(punctured):
        r = faces[f][0]
        if r in remaining:
            assert (r, n + 2 * k + 1) in full


@pytest.mark.parametrize("name, m, punctured, choice", PUNCTURED_CASES, ids=PUNCTURED_IDS)
def test_punctured_sizes(name, m, punctured, choice):
    g = m.genus()
    b = len(punctured)
    for relabel in (False, True):
        t, t_punctured = m.triangle_system(*choice, relabel=relabel, punctured_faces=punctured)
        assert t.num_vertices() == 3
        assert t.num_edges() == 6 * g + 2 + b
        assert t.num_faces() == 4 * g + 1 + b
        assert len(t_punctured) == b and t_punctured == sorted(t_punctured)
        faces = t.faces()
        # b - 1 punctured bigons and one punctured triangle
        assert sorted(len(faces[i]) for i in t_punctured) == [2] * (b - 1) + [3]
        triangles = [f for i, f in enumerate(faces) if i not in t_punctured]
        assert all(len(f) == 3 and f[0] % 2 == f[1] % 2 == f[2] % 2 for f in triangles)
        assert sorted(t.vertex_profile()) == sorted([4 * g + 1 + b, 4 * g + 1 + b, 4 * g + 2])
        assert t.genus() == g
        assert not t.is_mutable()
        BipartiteHyperbolicTriangulation(t, t_punctured)
        if relabel:
            assert list(t.half_edges()) == list(range(2 * (6 * g + 2 + b)))
            assert m.triangle_system(*choice, relabel=True, check=False, punctured_faces=punctured) == \
                (t, t_punctured)


@pytest.mark.parametrize("name, m, punctured, choice", PUNCTURED_CASES, ids=PUNCTURED_IDS)
def test_punctured_images(name, m, punctured, choice):
    tree = {h // 2 for h in choice[0] if h != -1}
    t, mor, _ = m.triangle_system(*choice, mapping=True, punctured_faces=punctured)
    # The docstring bounds the images by four. The bound three in genus 0 is
    # observed, not proved: no image of length four occurs on any genus 0
    # case of PUNCTURED_CASES.
    longest = 3 if m.genus() == 0 else 4
    vp = m.vertex_permutation()
    assert len(mor) == len(vp)
    for h in range(len(vp)):
        if vp[h] == -1:
            assert mor[h] is None
            continue
        w = list(mor[h])
        if h // 2 in tree:
            assert w == []
        assert len(w) <= longest
        if w:
            t._check_walk(w)
        assert all(w[i + 1] != w[i] ^ 1 for i in range(len(w) - 1))
        assert list(mor[h ^ 1]) == reverse(w)
    for f in m.faces():
        w = image(mor, f)
        if w:
            t._check_walk(w, closed=True)


def check_punctured_homotopy(m, punctured, t, mor, t_punctured):
    r"""
    Check that the images under ``mor`` of the boundaries of the filled faces
    of ``m`` are contractible in ``t``, and that those of the punctured faces
    are freely homotopic to the boundaries of distinct punctured faces of
    ``t``.
    """
    T = BipartiteHyperbolicTriangulation(t, t_punctured)
    faces = t.faces()
    around = {i: T.strongly_reduced(faces[i]) for i in t_punctured}
    assert all(around.values())
    matched = []
    for i, f in enumerate(m.faces()):
        r = T.strongly_reduced(image(mor, f))
        if i not in punctured:
            assert r == [], (i, r)
            continue
        # the puncture lies on the left of the boundary of the face, as the
        # punctured face of t does of its own boundary
        found = [j for j in t_punctured if is_rotation(r, around[j])]
        assert len(found) == 1, (i, r)
        matched.append(found[0])
    assert sorted(matched) == t_punctured
    return T


@pytest.mark.parametrize("name, m, punctured, choice", PUNCTURED_CASES, ids=PUNCTURED_IDS)
def test_punctured_homotopy(name, m, punctured, choice):
    t, mor, t_punctured = m.triangle_system(*choice, mapping=True, relabel=True, punctured_faces=punctured)
    T = check_punctured_homotopy(m, punctured, t, mor, t_punctured)

    def reduced(w):
        w = image(mor, w)
        return T.strongly_reduced(w) if w else []

    # beyond face boundaries: the strongly reduced image of a closed walk is
    # unchanged, up to rotation, by a rotation of the walk and by a detour
    # around a filled face, and that of the reverse walk is the reverse
    rng = random.Random(len(m.vertex_permutation()) + len(punctured))
    h2f = m.half_edge_to_face()
    faces = m.faces()
    filled = [f for f in range(m.num_faces()) if f not in punctured]
    for length in (3, 8, 15):
        w = random_closed_walk(m, rng, length)
        r = reduced(w)
        i = rng.randrange(len(w))
        assert is_rotation(reduced(w[i:] + w[:i]), r)
        assert is_rotation(reduced(reverse(w)), reverse(r))
        candidates = [i for i in range(len(w)) if h2f[w[i]] in filled]
        if candidates:
            i = rng.choice(candidates)
            face = faces[h2f[w[i]]]
            j = face.index(w[i])
            assert is_rotation(reduced(w[:i] + face[j:] + face[:j] + w[i:]), r)


@pytest.mark.parametrize("name, m, choice", CASES, ids=IDS)
def test_no_punctured_face(name, m, choice):
    for relabel in (False, True):
        for mapping in (False, True):
            reference = m.triangle_system(*choice, mapping=mapping, relabel=relabel)
            for punctured in (False, []):
                assert m.triangle_system(*choice, mapping=mapping, relabel=relabel,
                                         punctured_faces=punctured) == reference
    assert m.strands(*choice[:2], punctured_faces=[]) == m.strands(*choice[:2])
    assert m.strands(*choice[:2], punctured_flag=False) == m.strands(*choice[:2])


def test_punctured_mutable():
    name, m, punctured, choice = PUNCTURED_CASES[0]
    t, t_punctured = m.triangle_system(*choice, punctured_faces=punctured, mutable=True)
    assert t.is_mutable()
    reference = m.triangle_system(*choice, punctured_faces=punctured)
    assert (t.copy(mutable=False), t_punctured) == reference


def test_punctured_loop_strand_in_cocotree():
    # a cocotree holding the strand of the loop of the second punctured face,
    # 14 + 2, with relabel
    m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
    n = len(m.vertex_permutation())
    assert n == 14
    tree, cotree, _ = m.tree_cotree_decomposition()
    strands = m.strands(tree, cotree, punctured_faces=True)
    cocotree = random_strand_tree(strands, random.Random(0), first=n + 2)
    assert len(cocotree) == 5
    assert n + 2 in cocotree
    for c in [cocotree, random_strand_tree(strands, random.Random(1))]:
        t, mor, t_punctured = m.triangle_system(tree, cotree, c, relabel=True, mapping=True, punctured_faces=True)
        assert list(t.half_edges()) == list(range(2 * 16))
        assert sorted(len(t.faces()[i]) for i in t_punctured) == [2, 3]
        check_punctured_homotopy(m, [0, 1], t, mor, t_punctured)


def test_punctured_errors():
    TORUS.triangle_system()
    TORUS.triangle_system(punctured_faces=True)
    for m, punctured in [(OrientedMap(fp="(0,1,2)(~0,~2,~1)"), True),
                         (TETRAHEDRON, [0]), (TETRAHEDRON, [1, 3]), (TETRAHEDRON, False)]:
        with pytest.raises(NotImplementedError, match="genus zero|negative Euler characteristic"):
            m.triangle_system(punctured_faces=punctured)
    with pytest.raises(NotImplementedError, match=r"negative Euler characteristic .* \(got 0\)"):
        TETRAHEDRON.triangle_system(punctured_faces=[1, 3])
    with pytest.raises(NotImplementedError, match=r"negative Euler characteristic .* \(got 1\)"):
        TETRAHEDRON.triangle_system(punctured_faces=[2])
    TETRAHEDRON.triangle_system(punctured_faces=[0, 2, 3])
    m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
    tree, cotree, _ = m.tree_cotree_decomposition()
    cocotree = random_strand_tree(m.strands(tree, cotree, punctured_faces=True), random.Random(0))
    assert len(cocotree) == 5
    m.triangle_system(tree, cotree, cocotree, punctured_faces=True)
    with pytest.raises(ValueError, match="cocotree must have 5 strands"):
        m.triangle_system(tree, cotree, cocotree[:4], punctured_faces=True)
    # the strand around the monogon of the loop 14, 15 is left out, and so is
    # the strand starting at the first loop 14
    with pytest.raises(ValueError, match="the half-edge 15 starts no strand"):
        m.triangle_system(tree, cotree, cocotree[:4] + [15], punctured_faces=True)
    with pytest.raises(ValueError, match="the half-edge 14 starts no strand"):
        m.triangle_system(tree, cotree, cocotree[:4] + [14], punctured_faces=True)
    with pytest.raises(ValueError, match="face number out of range"):
        m.triangle_system(punctured_faces=[2])


@pytest.mark.parametrize("name, m, punctured, choice", PUNCTURED_CASES[::CHOICES], ids=PUNCTURED_IDS[::CHOICES])
def test_punctured_first_loop_in_cocotree(name, m, punctured, choice):
    # a cocotree holding the strand of the first loop n is refused; it is a
    # spanning tree of the graph of strands with punctured_flag=False unless
    # that strand is a loop of this graph
    tree, cotree, cocotree = choice
    n = len(m.vertex_permutation())
    full = m.strands(tree, cotree, punctured_faces=punctured, punctured_flag=False)
    (start, end), = [s for s in full if s[0] == n]
    if start // 2 != end // 2:
        cocotree = random_strand_tree(full, random.Random(0), first=n)
        assert n in cocotree
    else:
        cocotree = list(cocotree[:-1]) + [n]
    with pytest.raises(ValueError, match=f"the half-edge {n} starts no strand"):
        m.triangle_system(tree, cotree, cocotree, punctured_faces=punctured)


def test_errors():
    with pytest.raises(NotImplementedError):
        OrientedMap(vp="(0,1,~1)").triangle_system()
    with pytest.raises(NotImplementedError):
        OrientedMap(fp="(0,1,2)(~0,~2,~1)").triangle_system()
    with pytest.raises(ValueError):
        OrientedMap(fp="(0,1,~0,~1)(2,3,~2,~3)").triangle_system()
    m = one_vertex_map(2)
    strands = m.strands([], [])
    # two strands between the same two edges make a cycle
    pairs = {}
    for h, hh in strands:
        pairs.setdefault(frozenset((h // 2, hh // 2)), []).append(h)
    edges, cycle = next((e, p[:2]) for e, p in pairs.items() if len(p) >= 2)
    other = next(h for h, hh in strands if not {h // 2, hh // 2} <= edges)
    with pytest.raises(ValueError, match="cycle"):
        m.triangle_system([], [], cycle + [other])
    with pytest.raises(ValueError, match="cocotree must have 3 strands"):
        m.triangle_system([], [], cycle[:1])
    with pytest.raises(ValueError, match="cocotree must have 3 strands"):
        m.triangle_system([], [], [h for h, _ in strands][:4])
    m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
    tree, cotree, _ = m.tree_cotree_decomposition()
    edge = next(h for h in tree if h != -1)
    with pytest.raises(ValueError, match="listed twice"):
        m.triangle_system(tree, list(cotree) + [edge])
    with pytest.raises(ValueError, match="listed twice"):
        m.triangle_system(tree, list(cotree) + [edge ^ 1])
