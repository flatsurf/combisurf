r"""
Barycentric subdivisions: sizes, images of the half-edges and colours of the
triangles, on the maps of the hyperbolic quadrangulation tests and on random
connected maps of genus 0 to 4, with no, one, some or every face punctured.
"""
import random

import pytest

from combisurf import OrientedMap

from test_geometric_intersection_oracles import cylinder_data, genus2, random_closed_walk
from test_hyperbolic_quadrangulation import INSTANCES, build, genus2_pentagons, genus2_two_hexagons

VERTEX, EDGE, FACE = "vertex", "edge", "face"

MAX_EDGES = 30


def one_vertex_map(genus):
    r"""
    Return the map with one vertex and one face of genus ``genus``, or a single
    edge between two vertices for genus ``0``.
    """
    if genus == 0:
        return OrientedMap(vp=[[0], [1]], mutable=True)
    return OrientedMap(vp=[[x + 4 * i for i in range(genus) for x in (0, 2, 1, 3)]], mutable=True)


def random_map(rng, genus, num_edges):
    r"""
    Return a random connected map of genus ``genus`` with ``num_edges`` edges.

    Starting from :func:`one_vertex_map`, each step adds either a pendant edge
    at a new vertex or an edge between two corners of the same face, which
    splits that face; neither changes the genus.
    """
    if num_edges > MAX_EDGES:
        raise ValueError(f"at most {MAX_EDGES} edges")
    m = one_vertex_map(genus)
    while m.num_edges() < num_edges:
        h1 = rng.choice(list(m.half_edges()))
        if rng.random() < 0.3:
            m.add_edge(-1, h1)
        else:
            face = [h1]
            h = m.next_in_face(h1)
            while h != h1:
                face.append(h)
                h = m.next_in_face(h)
            m.add_edge(rng.choice(face), h1)
    assert m.is_connected() and m.genus() == genus
    return m.copy(mutable=False)


def random_maps():
    rng = random.Random(0)
    maps = []
    for genus in range(5):
        for _ in range(6):
            maps.append(random_map(rng, genus, rng.randint(max(1, 2 * genus), MAX_EDGES)))
    return maps


def example_maps():
    r"""
    Return the list of the pairs ``(name, m)`` of the tested maps.
    """
    maps = [("genus2", genus2()), ("hexagons", genus2_two_hexagons()), ("pentagons", genus2_pentagons()),
            ("lshape", OrientedMap(vp=cylinder_data("lshape")["vp"])),
            ("emmw_1117", OrientedMap(vp=cylinder_data("emmw_1117")["vp"])),
            ("inactive", OrientedMap(vp="(0,3,6,~3)(~0,1,~6,~1)"))]
    maps.extend((name, build(name)[1]) for name in INSTANCES)
    maps.extend((f"random_{i}", m) for i, m in enumerate(random_maps()))
    return maps


MAPS = example_maps()

PUNCTURES = ["none", "one", "some", "all"]


def punctured_faces(m, choice, rng):
    r"""
    Return the sorted list of the punctured faces of ``m`` for ``choice``.
    """
    F = m.num_faces()
    if choice == "none":
        return []
    if choice == "one":
        return [rng.randrange(F)]
    if choice == "some":
        return sorted(f for f in range(F) if rng.random() < 0.5)
    return list(range(F))


def cases():
    rng = random.Random(1)
    return [(name, m, punctured_faces(m, choice, rng)) for name, m in MAPS for choice in PUNCTURES]


CASES = cases()
IDS = [f"{name}-{choice}" for name, _ in MAPS for choice in PUNCTURES]


def test_number_of_maps():
    assert len(MAPS) == 6 + len(INSTANCES) + 30


def vertex_kinds(m, t):
    r"""
    Return the list of the kinds of the vertices of ``t``: ``VERTEX`` for a
    vertex of ``m``, ``EDGE`` for an edge vertex and ``FACE`` for a face vertex.
    """
    n = len(m.vertex_permutation())
    h2v = t.half_edge_to_vertex()
    kinds = [None] * t.num_vertices()

    def mark(x, kind):
        if x < len(h2v) and h2v[x] != -1:
            v = h2v[x]
            assert kinds[v] in (None, kind), (x, kinds[v], kind)
            kinds[v] = kind

    for h in m.half_edges():
        mark(2 * h, VERTEX)
        mark(2 * h + 1, EDGE)
        mark(2 * n + 2 * h, EDGE)
        mark(2 * n + 2 * h + 1, FACE)
        mark(4 * n + 2 * h, FACE)
        mark(4 * n + 2 * h + 1, VERTEX)
    assert None not in kinds
    return kinds


@pytest.mark.parametrize("name, m, punctured", CASES, ids=IDS)
def test_sizes(name, m, punctured):
    t, t_punctured = m.barycentric_subdivision(punctured_faces=punctured)
    faces = m.faces()
    filled = [f for i, f in enumerate(faces) if i not in punctured]
    assert t.num_vertices() == m.num_vertices() + m.num_edges() + len(filled)
    assert t.genus() == m.genus()
    assert t.num_faces() == 2 * sum(len(f) for f in filled) + len(punctured)
    assert t_punctured == sorted(t_punctured)
    assert len(t_punctured) == len(punctured)
    t_faces = t.faces()
    t_h2f = t.half_edge_to_face()
    for f in punctured:
        g = t_h2f[2 * faces[f][0]]
        assert g in t_punctured
        assert len(t_faces[g]) == 2 * len(faces[f])
    for i, f in enumerate(t_faces):
        if i not in t_punctured:
            assert len(f) == 3


@pytest.mark.parametrize("name, m, punctured", CASES, ids=IDS)
def test_images(name, m, punctured):
    t, mor, _ = m.barycentric_subdivision(mapping=True, punctured_faces=punctured)
    assert len(mor) == len(m.vertex_permutation())
    for h in range(len(mor)):
        if m.vertex_permutation()[h] == -1:
            assert mor[h] is None
            continue
        assert len(mor[h]) == 2
        t._check_walk(mor[h])
        assert list(mor[h ^ 1]) == [mor[h][1] ^ 1, mor[h][0] ^ 1]
    walks = list(m.faces())
    if min(m.vertex_profile()) > 1:
        # random_closed_walk avoids backtracking, which a vertex of degree 1 forces
        rng = random.Random(len(m.vertex_permutation()))
        walks.extend(random_closed_walk(m, rng, length) for length in (1, 3, 8))
    for w in walks:
        t._check_walk([x for h in w for x in mor[h]], closed=True)


def is_rotation(u, v):
    return len(u) == len(v) and any(u[i:] + u[:i] == v for i in range(len(u)))


@pytest.mark.parametrize("name, m, punctured", CASES, ids=IDS)
def test_colours(name, m, punctured):
    t, t_punctured = m.barycentric_subdivision(punctured_faces=punctured)
    kinds = vertex_kinds(m, t)
    filled = m.num_faces() - len(punctured)
    assert kinds.count(VERTEX) == m.num_vertices()
    assert kinds.count(EDGE) == m.num_edges()
    assert kinds.count(FACE) == filled
    h2v = t.half_edge_to_vertex()
    for i, f in enumerate(t.faces()):
        corners = [kinds[h2v[x]] for x in f]
        if i in t_punctured:
            assert all(f[j] % 2 != f[(j + 1) % len(f)] % 2 for j in range(len(f)))
            assert FACE not in corners
            continue
        parities = {x % 2 for x in f}
        assert len(parities) == 1, f
        # the face cycle runs counterclockwise around its face
        ccw = is_rotation(corners, [VERTEX, EDGE, FACE])
        assert ccw or is_rotation(corners, [VERTEX, FACE, EDGE]), corners
        assert ccw == (parities == {0}), (f, corners)


def test_mutable():
    m = genus2()
    assert not m.barycentric_subdivision()[0].is_mutable()
    assert m.barycentric_subdivision(mutable=True)[0].is_mutable()


def test_folded():
    with pytest.raises(NotImplementedError):
        OrientedMap(vp="(0,1,~1)").barycentric_subdivision()
