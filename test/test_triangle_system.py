r"""
Triangle systems: strands, sizes, colours and images of the half-edges, on the
maps of genus 2 to 4 of the barycentric subdivision tests, each with several
random choices of tree, cotree and cocotree. The homotopy classes of closed
walks are compared through the hyperbolic quadrangulations of the quad systems
of the map and of its triangle system.
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


def random_strand_tree(strands, rng):
    r"""
    Return a random spanning tree of the graph of strands, as the list of the
    first half-edges of its strands.
    """
    strands = list(strands)
    rng.shuffle(strands)
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


def choices(m, seed):
    r"""
    Return ``CHOICES`` triples ``(tree, cotree, cocotree)`` for ``m``, from
    random roots and random spanning trees of the graph of strands.
    """
    rng = random.Random(seed)
    result = []
    for _ in range(CHOICES):
        tree, cotree, _ = m.forest_coforest_decomposition((rng.randrange(m.num_vertices()),),
                                                          (rng.randrange(m.num_faces()),))
        cocotree = random_strand_tree(m.strands(tree, cotree), rng)
        result.append((tree, cotree, cocotree))
    return result


CASES = [(name, m, c) for i, (name, m) in enumerate(MAPS) for c in choices(m, i)]
IDS = [f"{name}-{j}" for name, _ in MAPS for j in range(CHOICES)]


def test_number_of_maps():
    assert len(MAPS) == 3 + 18


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
        t = m.triangle_system(*choice, relabel=relabel)
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
    t, mor = m.triangle_system(*choice, mapping=True)
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
    t, mor = m.triangle_system(*choice, mapping=True, relabel=True)
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
