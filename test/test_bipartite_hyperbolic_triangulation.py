r"""
Checks of the constructor of bipartite hyperbolic triangulations, on the
barycentric subdivisions of the maps of ``test_barycentric_subdivision``.
"""
import random

import pytest

from combisurf import OrientedMap
from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation

from test_barycentric_subdivision import CASES, IDS


def expected_accepted(m, punctured):
    r"""
    Return whether the barycentric subdivision of ``m`` with the faces
    ``punctured`` punctured is a bipartite hyperbolic triangulation.

    Its triangles are 2-coloured by construction. An edge vertex has degree `4`
    unless its edge touches a punctured face, a face vertex has twice the
    degree of its face, and a vertex of ``m`` with no punctured corner has
    twice its degree.
    """
    is_punctured = [f in punctured for f in range(m.num_faces())]
    h2f = m.half_edge_to_face()
    if m.euler_characteristic() - len(punctured) >= 0:
        return False
    if any(not is_punctured[h2f[h]] and not is_punctured[h2f[h ^ 1]] for h in m.half_edges()):
        return False
    if any(not is_punctured[i] and len(f) < 4 for i, f in enumerate(m.faces())):
        return False
    return all(len(v) >= 4 or any(is_punctured[h2f[h]] for h in v) for v in m.vertices())


@pytest.mark.parametrize("name, m, punctured", CASES, ids=IDS)
def test_barycentric_subdivision(name, m, punctured):
    t, t_punctured = m.barycentric_subdivision(punctured_faces=punctured)
    if expected_accepted(m, punctured):
        T = BipartiteHyperbolicTriangulation(t, t_punctured)
        assert T.map() == t
        assert T.punctured_faces() == t_punctured
    else:
        with pytest.raises(ValueError):
            BipartiteHyperbolicTriangulation(t, t_punctured)
    BipartiteHyperbolicTriangulation(t, t_punctured, check=False)


def test_accepted_cases():
    accepted = sum(expected_accepted(m, punctured) for _, m, punctured in CASES)
    assert 0 < accepted < len(CASES)


def reversed_edge(t, t_punctured, e):
    r"""
    Return a copy of ``t`` with the orientation of the edge ``e`` reversed,
    and the indices of its faces that were ``t_punctured`` in ``t``: the
    faces are numbered from their labels, which the reversal changes.
    """
    u = t.copy(mutable=True)
    u.reverse_orientation(e)
    faces = t.faces()
    h2f = u.half_edge_to_face()
    u_punctured = sorted(h2f[x ^ 1 if x // 2 == e else x] for x in (faces[g][0] for g in t_punctured))
    return u, u_punctured


@pytest.mark.parametrize("name, m, punctured", CASES, ids=IDS)
def test_parity(name, m, punctured):
    if not expected_accepted(m, punctured):
        return
    t, t_punctured = m.barycentric_subdivision(punctured_faces=punctured)
    h2f = t.half_edge_to_face()
    rng = random.Random(t.num_half_edges())
    edges = list(t.edge_indices())
    rng.shuffle(edges)
    for e in edges[:10]:
        u, u_punctured = reversed_edge(t, t_punctured, e)
        if h2f[2 * e] in t_punctured and h2f[2 * e + 1] in t_punctured:
            # no parity constraint between two punctured faces
            BipartiteHyperbolicTriangulation(u, u_punctured)
        else:
            with pytest.raises(ValueError, match="mixed parities"):
                BipartiteHyperbolicTriangulation(u, u_punctured)


def test_between_punctured_faces():
    # every edge of the torus lies between its two punctured faces
    m = OrientedMap(fp="(0,1,2)(~0,~1,~2)")
    t, t_punctured = m.barycentric_subdivision(punctured_faces=True)
    BipartiteHyperbolicTriangulation(t, t_punctured)
    for e in t.edge_indices():
        BipartiteHyperbolicTriangulation(*reversed_edge(t, t_punctured, e))


def test_errors():
    with pytest.raises(TypeError, match="oriented map"):
        BipartiteHyperbolicTriangulation("(0,1,2)(~0,~1,~2)")
    with pytest.raises(NotImplementedError, match="folded"):
        BipartiteHyperbolicTriangulation(OrientedMap(vp="(0,1,~1)"))
    with pytest.raises(ValueError, match="connected"):
        BipartiteHyperbolicTriangulation(OrientedMap(fp="(0,1,2)(~0,~1,~2)(3,4,5)(~3,~4,~5)"), True)
    m = OrientedMap(fp="(0,1,2,3)(~0,~1,~2,~3)")
    with pytest.raises(ValueError, match="negative Euler characteristic"):
        BipartiteHyperbolicTriangulation(m.barycentric_subdivision()[0])
    with pytest.raises(ValueError, match="degree 4 instead of 3"):
        BipartiteHyperbolicTriangulation(m, [1])
    t, t_punctured = m.barycentric_subdivision(punctured_faces=[1])
    BipartiteHyperbolicTriangulation(t, t_punctured)
    with pytest.raises(ValueError, match="mixed parities"):
        BipartiteHyperbolicTriangulation(*reversed_edge(t, t_punctured, 8))
    t, _ = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)").barycentric_subdivision()
    with pytest.raises(ValueError, match="instead of at least 8"):
        BipartiteHyperbolicTriangulation(t)
    with pytest.raises(ValueError, match="out of range"):
        BipartiteHyperbolicTriangulation(t, [16])


def test_check_false():
    m = OrientedMap(fp="(0,1,2,3)(~0,~1,~2,~3)")
    T = BipartiteHyperbolicTriangulation(m, check=False)
    assert T.map() == m
    assert T.punctured_faces() == []
