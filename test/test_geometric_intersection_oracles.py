r"""
Reference values and invariants of geometric intersection numbers.

Letters are the half-edges of a one-vertex map: ``a = 0``, ``a⁻¹ = 1``,
``b = 2``, ``b⁻¹ = 3``, ``c = 4``, ``c⁻¹ = 5``, ``d = 6``, ``d⁻¹ = 7``; the
inverse of ``h`` is ``h ^ 1``.

The tests on closed surfaces go through :func:`pairing`, which skips them as
long as :class:`GeometricIntersectionPairing` raises ``NotImplementedError``
there.
"""
import itertools
import json
import os
import random

import pytest

from combisurf import OrientedMap
from combisurf.geometric_intersection import GeometricIntersectionPairing
from combisurf.word import word_init, word_reduce, word_cyclically_reduce

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

LETTERS = {"a": 0, "b": 2, "c": 4, "d": 6}


def _skip_not_implemented(f, *args, **kwds):
    try:
        return f(*args, **kwds)
    except NotImplementedError as e:
        pytest.skip(str(e))


class pairing:
    r"""
    A :class:`GeometricIntersectionPairing` on ``m`` with the given punctured
    faces, whose ``NotImplementedError``, at construction or at a call, skips
    the test.
    """
    def __init__(self, m, punctured_faces):
        self._gi = _skip_not_implemented(GeometricIntersectionPairing, m, punctured_faces=punctured_faces)

    def __call__(self, u, v=None):
        return _skip_not_implemented(self._gi, u, v)


def word(*pairs):
    r"""
    Return the word of the ``(letter, exponent)`` pairs, a negative exponent
    standing for a power of the inverse letter.
    """
    w = []
    for x, e in pairs:
        if isinstance(x, str):
            h = LETTERS[x]
            w.extend([h if e > 0 else h ^ 1] * abs(e))
        else:
            # x is itself a word
            w.extend((x if e > 0 else inverse(x)) * abs(e))
    return w


def inverse(w):
    return [h ^ 1 for h in reversed(w)]


def rotations(w):
    return [w[i:] + w[:i] for i in range(len(w))]


def is_cyclic_conjugate(u, v):
    return len(u) == len(v) and u in rotations(v)


def determinant(M):
    r"""
    Return the determinant of the small integer matrix ``M`` (Laplace
    expansion along the first row).
    """
    if len(M) == 1:
        return M[0][0]
    return sum((-1) ** j * M[0][j] * determinant([row[:j] + row[j + 1:] for row in M[1:]])
               for j in range(len(M)) if M[0][j])


# NOTE: the punctured torus and the genus 2 surface are the one-vertex maps
# whose vertex sees the alphabets a, b, a⁻¹, b⁻¹ and a, b, a⁻¹, b⁻¹, c, d,
# c⁻¹, d⁻¹ of Cohen–Lustig 1987, §5, in this cyclic order; the third map sees
# a, b, c, a⁻¹, b⁻¹, c⁻¹.
def torus():
    return OrientedMap(vp=[[0, 2, 1, 3]])


def genus2():
    return OrientedMap(vp=[[0, 2, 1, 3, 4, 6, 5, 7]])


def genus1_two_faces():
    return OrientedMap(vp=[[0, 2, 4, 1, 3, 5]])


def pants():
    return OrientedMap(fp="(0)(1)(~0,~1)")


def sphere4():
    return OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~4,~5)(~2,5,~3)")


X = word(("b", -1), ("a", 1), ("b", 1), ("a", -1))
Y = word(("c", 1), ("d", -1), ("c", -1), ("d", 1))

B1_V = word(("c", 1), (X, 1), ("c", 2), ("d", 2), (Y, 1), ("c", 1))
B1_W = word(("c", 1), (X, 1), ("c", 2), ("d", 2), (X, 1), ("c", 1))
B2_V = word(("a", 70), ("b", -3), ("c", 17))
B3_V = word(("a", 3), ("b", 1), ("a", -1), ("c", 2))
B3_W = word(("b", 2), ("c", 1))
C1_V = word(("a", 5), ("b", 1), ("a", -1), ("c", 28), ("b", -9))
C2_V = B3_V
C2_W = B3_W

# Cohen–Lustig 1987, §5, p. 498: the self-intersections |V| and |W| and the
# intersection |V ∩ W|, all faces punctured; None where the article gives no
# value.
COHEN_LUSTIG_ROWS = [
    ("A1", torus, word(("a", 1), ("b", 2), ("a", -1), ("b", -3)),
     word(("b", 1), ("a", 2), ("b", -1), ("a", -3)), 3, 3, 13),
    ("A2", torus, word(("a", 3), ("b", 7)), word(("a", 5), ("b", 13)), 12, 48, 54),
    ("A3", torus, word(("a", 6), ("b", 15)), word(("a", 10), ("b", 8)), 70, 63, 172),
    ("B1", genus2, B1_V, B1_W, 11, 9, 20),
    ("B2", genus2, B2_V, None, 292, None, None),
    ("B3", genus2, B3_V, B3_W, 2, 1, 6),
    ("C1", genus1_two_faces, C1_V, None, 408, None, None),
    ("C2", genus1_two_faces, C2_V, C2_W, 8, 0, 11),
]


@pytest.mark.parametrize("row", COHEN_LUSTIG_ROWS, ids=[row[0] for row in COHEN_LUSTIG_ROWS])
def test_cohen_lustig_tables(row):
    _, surface, V, W, iV, iW, iVW = row
    gi = GeometricIntersectionPairing(surface(), punctured_faces=True)
    assert gi([V]) == iV
    if W is not None:
        assert gi([W]) == iW
        assert gi([V], [W]) == iVW
        assert gi([W], [V]) == iVW


def cohen_lustig_remark1(p, q, r, s):
    r"""
    Return ``|aᵖ bᵠ ∩ aʳ bˢ|`` on the punctured torus for ``p > r``, from
    Cohen–Lustig 1987, §5 Remark 1, pp. 498–499.
    """
    if s >= q:
        delta = p - r + s - q
    else:
        delta = p - r + q - s - 2
    return (p - 1) * (s - 1) + (q - 1) * (r - 1) + delta


def test_cohen_lustig_remark1_example():
    # the row A2 of the table: p, q, r, s = 5, 13, 3, 7
    assert cohen_lustig_remark1(5, 13, 3, 7) == 54


@pytest.mark.parametrize("p", range(2, 9))
def test_cohen_lustig_remark1(p):
    gi = GeometricIntersectionPairing(torus(), punctured_faces=True)
    for r in range(1, p):
        for q in range(1, 9):
            for s in range(1, 9):
                V = word(("a", p), ("b", q))
                W = word(("a", r), ("b", s))
                assert gi([V], [W]) == cohen_lustig_remark1(p, q, r, s), (p, q, r, s)


# Cohen–Lustig 1987, §5 Remark 2, p. 499: with its face filled, the genus 2
# map has the relator b⁻¹ a b a⁻¹ d⁻¹ c d c⁻¹, that is X = Y, so that the two
# words V, W of the row B1 are the same closed curve.
GENUS2_RELATOR = word(("b", -1), ("a", 1), ("b", 1), ("a", -1), ("d", -1), ("c", 1), ("d", 1), ("c", -1))


def test_genus2_relator():
    faces = genus2().faces()
    assert len(faces) == 1
    assert is_cyclic_conjugate(faces[0], GENUS2_RELATOR)
    assert GENUS2_RELATOR == X + inverse(Y)


def test_cohen_lustig_remark2():
    gi = pairing(genus2(), punctured_faces=False)
    assert gi([B1_V]) == gi([B1_W])
    assert gi([B1_V], [B1_W]) == gi([B1_V], [B1_V]) == 2 * gi([B1_V])
    for U in [B2_V, B3_V, B3_W, C1_V]:
        assert gi([B1_V], [U]) == gi([B1_W], [U])
    # filling the face does not increase the intersection, and |W| = 9 on the
    # punctured surface
    assert gi([B1_V]) <= 9


# On a translation surface the core curves of two cylinders cross with the
# same sign at each crossing, so that their geometric intersection is the
# absolute value of their algebraic intersection, and core curves are simple.
# The walks follow the left or the right side of the boundary saddle
# connections of each cylinder, on the triangulation map whose vertices are
# the singularities and whose faces are the triangles.
CYLINDER_EXAMPLES = [
    # name, number of curves, vertices, faces, genus
    ("lshape", 26, 1, 6, 2),
    ("emmw_1117", 50, 4, 20, 4),
]


def cylinder_data(name):
    with open(os.path.join(DATA, f"cylinders_{name}.json")) as f:
        return json.load(f)


@pytest.mark.parametrize("name, num_curves, nv, nf, g", CYLINDER_EXAMPLES, ids=[e[0] for e in CYLINDER_EXAMPLES])
def test_cylinder_core_curves_data(name, num_curves, nv, nf, g):
    data = cylinder_data(name)
    assert data["name"] == name
    m = OrientedMap(vp=data["vp"])
    assert (m.num_vertices(), m.num_faces(), m.genus()) == (nv, nf, g)
    assert len(data["left"]) == len(data["right"]) == num_curves
    for w in data["left"] + data["right"]:
        m._check_walk(w, closed=True)
    M = data["algebraic_intersection"]
    assert len(M) == num_curves
    assert all(len(row) == num_curves for row in M)
    assert all(M[i][j] == M[j][i] for i in range(num_curves) for j in range(num_curves))
    assert all(M[i][i] == 0 for i in range(num_curves))


@pytest.mark.parametrize("name", [e[0] for e in CYLINDER_EXAMPLES])
def test_cylinder_core_curves(name):
    data = cylinder_data(name)
    m = OrientedMap(vp=data["vp"])
    left = data["left"]
    right = data["right"]
    M = data["algebraic_intersection"]
    gi = pairing(m, punctured_faces=False)
    N = len(left)
    for i in range(N):
        assert gi([left[i]]) == 0, i
        assert gi([right[i]]) == 0, i
        # the two sides of a cylinder are homotopic simple curves
        assert gi([left[i]], [right[i]]) == 0, i
    for i in range(N):
        for j in range(i + 1, N):
            assert gi([left[i]], [left[j]]) == M[i][j], (i, j)
            assert gi([left[i]], [right[j]]) == M[i][j], (i, j)


def random_closed_walk(m, rng, length):
    r"""
    Return a closed walk on ``m`` that starts with a random non-backtracking
    walk of ``length`` half-edges and goes back to its start along a shortest
    path.
    """
    vp = m.vertex_permutation()
    vertex = m.half_edge_to_vertex()
    n = len(vp)
    half_edges = [h for h in range(n) if vertex[h] != -1]
    out = {}
    for h in half_edges:
        out.setdefault(vertex[h], []).append(h)

    h = rng.choice(half_edges)
    w = [h]
    for _ in range(length - 1):
        choices = [x for x in out[vertex[h ^ 1]] if x != h ^ 1]
        h = rng.choice(choices)
        w.append(h)

    # breadth-first search from the end of w back to its start
    start = vertex[w[0]]
    end = vertex[w[-1] ^ 1]
    parent = {end: None}
    queue = [end]
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
        x = parent[v]
        path.append(x)
        v = vertex[x]
    w.extend(reversed(path))
    m._check_walk(w, closed=True)
    return w


def filling_instances():
    rng = random.Random(0)
    m = sphere4()
    sphere_curves = [list(word_init("0,1,2,0,~4,~3"))]
    sphere_curves.extend(random_closed_walk(m, rng, length) for length in (2, 3, 4, 5, 6, 8))
    lshape = cylinder_data("lshape")
    return [
        ("sphere4", m, sphere_curves),
        ("C", genus1_two_faces(), [C1_V, C2_V, C2_W]),
        ("lshape", OrientedMap(vp=lshape["vp"]), lshape["left"][:6]),
    ]


def intersection_table(gi, curves):
    N = len(curves)
    return ([gi([u]) for u in curves],
            [[gi([curves[i]], [curves[j]]) for j in range(i, N)] for i in range(N)])


def assert_not_larger(table0, table1):
    self0, pairs0 = table0
    self1, pairs1 = table1
    assert all(x <= y for x, y in zip(self0, self1)), (self0, self1)
    for row0, row1 in zip(pairs0, pairs1):
        assert all(x <= y for x, y in zip(row0, row1)), (row0, row1)


def min_punctures(m):
    r"""
    Return the least positive number of punctures for which ``m`` has
    negative Euler characteristic.

    :class:`GeometricIntersectionPairing` raises ``NotImplementedError`` on
    the disk, the annulus, the sphere and the torus.
    """
    return max(1, m.euler_characteristic() + 1)


@pytest.mark.parametrize("instance", filling_instances(), ids=lambda instance: instance[0])
def test_filling_punctures(instance):
    # For P' ⊂ P, i_{P'} ≤ i_P; by transitivity it is enough to compare P'
    # with the sets P' ∪ {f}.
    _, m, curves = instance
    faces = range(m.num_faces())
    tables = {}
    for k in range(min_punctures(m), len(faces) + 1):
        for P in itertools.combinations(faces, k):
            gi = GeometricIntersectionPairing(m, punctured_faces=list(P))
            tables[P] = intersection_table(gi, curves)
    for P in tables:
        for f in faces:
            if f not in P:
                assert_not_larger(tables[P], tables[tuple(sorted(P + (f,)))])


@pytest.mark.parametrize("instance", [instance for instance in filling_instances() if instance[1].euler_characteristic() < 0],
                         ids=lambda instance: instance[0])
def test_filling_all_punctures(instance):
    _, m, curves = instance
    table = intersection_table(pairing(m, punctured_faces=False), curves)
    for P in itertools.combinations(range(m.num_faces()), 1):
        gi = GeometricIntersectionPairing(m, punctured_faces=list(P))
        assert_not_larger(table, intersection_table(gi, curves))


def morphism(images, rank):
    r"""
    Return the morphism of the free group on the ``rank`` first letters that
    sends the letter ``h`` to ``images[h]`` if ``h`` is in ``images`` and
    fixes it otherwise, as the list of the images of all half-edges.
    """
    mor = [None] * (2 * rank)
    for h in range(0, 2 * rank, 2):
        mor[h] = list(word_reduce(word_init(images.get(h, [h]))))
        mor[h ^ 1] = inverse(mor[h])
    return mor


def apply(mor, w):
    return list(word_cyclically_reduce(word_init(sum((mor[h] for h in w), []))))


def homology_action(mor, rank):
    r"""
    Return the matrix of the action of ``mor`` on the abelianisation of the
    free group on the ``rank`` first letters; its column ``j`` is the image
    of the letter ``2 j``.
    """
    M = [[0] * rank for _ in range(rank)]
    for j in range(rank):
        for h in mor[2 * j]:
            M[h >> 1][j] += -1 if h & 1 else 1
    return M


# Each surface: its name, its map, the punctured faces (True for all of them,
# False for none), the number of letters, and mapping classes as (name, images
# of the letters that are not fixed, whether it is a Dehn twist). A Dehn twist
# acts on homology by a transvection, of determinant 1. The twist along a
# simple curve γ fixes γ and the letters disjoint from it.
_sep = X  # b⁻¹ a b a⁻¹, separates the two handles of the genus 2 surface
MAPPING_CLASSES = [
    ("pants", pants, True, 2, [
        ("exchange", {0: [2], 2: [0]}, False),
        ("rotation", {0: [2], 2: [1, 3]}, False),
    ]),
    ("torus", torus, True, 2, [
        ("twist_a", {2: word(("b", 1), ("a", 1))}, True),
        ("twist_a_inverse", {2: word(("b", 1), ("a", -1))}, True),
        ("twist_b", {0: word(("a", 1), ("b", 1))}, True),
        ("twist_b_inverse", {0: word(("a", 1), ("b", -1))}, True),
    ]),
    ("genus2", genus2, False, 4, [
        ("twist_a", {2: word(("a", -1), ("b", 1))}, True),
        ("twist_a_inverse", {2: word(("a", 1), ("b", 1))}, True),
        ("twist_b", {0: word(("a", 1), ("b", -1))}, True),
        ("twist_b_inverse", {0: word(("a", 1), ("b", 1))}, True),
        ("twist_c", {6: word(("c", -1), ("d", 1))}, True),
        ("twist_c_inverse", {6: word(("c", 1), ("d", 1))}, True),
        ("twist_d", {4: word(("c", 1), ("d", -1))}, True),
        ("twist_d_inverse", {4: word(("c", 1), ("d", 1))}, True),
        # the non-separating curve b c meets both handles
        ("twist_bc", {0: word(("b", 1), ("c", 1), ("a", 1)), 6: word(("d", 1), ("c", -1), ("b", -1))}, True),
        ("twist_bc_inverse", {0: word(("c", -1), ("b", -1), ("a", 1)), 6: word(("d", 1), ("b", 1), ("c", 1))}, True),
        # the twist along the separating curve b⁻¹ a b a⁻¹ conjugates the
        # first handle
        ("twist_separating", {0: word((_sep, 1), ("a", 1), (_sep, -1)),
                              2: word((_sep, 1), ("b", 1), (_sep, -1))}, True),
    ]),
]
MAPPING_CLASS_CASES = [(surface, punctured, rank, name, images, is_twist)
                       for surface, _, punctured, rank, classes in MAPPING_CLASSES
                       for name, images, is_twist in classes]
_SURFACES = {surface: m for surface, m, _, _, _ in MAPPING_CLASSES}


def boundary_words(surface):
    if surface == "genus2":
        return [GENUS2_RELATOR]
    return _SURFACES[surface]().faces()


@pytest.mark.parametrize("surface, punctured, rank, name, images, is_twist", MAPPING_CLASS_CASES,
                         ids=[f"{case[0]}-{case[3]}" for case in MAPPING_CLASS_CASES])
def test_mapping_class(surface, punctured, rank, name, images, is_twist):
    mor = morphism(images, rank)
    boundaries = boundary_words(surface)
    for w in boundaries:
        image = apply(mor, w)
        assert any(is_cyclic_conjugate(image, u) or is_cyclic_conjugate(image, inverse(u)) for u in boundaries), (w, image)
    M = homology_action(mor, rank)
    det = determinant(M)
    assert det in (1, -1)
    if is_twist:
        assert det == 1
        # a transvection v ↦ v + (v · γ) γ: M - 1 has rank at most one
        N = [[M[i][j] - (i == j) for j in range(rank)] for i in range(rank)]
        assert all(N[i][j] * N[k][jj] == N[i][jj] * N[k][j]
                   for i in range(rank) for k in range(rank) for j in range(rank) for jj in range(rank))


def test_mapping_class_inverses():
    # the pairs (twist, twist_inverse) compose to the identity, up to
    # conjugation of the images of the curves
    rng = random.Random(1)
    for surface, _, _, rank, classes in MAPPING_CLASSES:
        named = {name: morphism(images, rank) for name, images, _ in classes}
        for name in named:
            if name + "_inverse" not in named:
                continue
            mor = named[name]
            mor_inv = named[name + "_inverse"]
            for _ in range(20):
                u = random_word(rng, rank, 12)
                image = apply(mor_inv, apply(mor, u))
                assert is_cyclic_conjugate(image, u), (surface, name, u, image)


def random_word(rng, rank, max_length):
    r"""
    Return a random non-empty cyclically reduced word of length at most
    ``max_length`` in the ``rank`` first letters.
    """
    while True:
        length = rng.randint(1, max_length)
        w = [rng.randrange(2 * rank)]
        for _ in range(length - 1):
            w.append(rng.choice([h for h in range(2 * rank) if h != w[-1] ^ 1]))
        w = list(word_cyclically_reduce(word_init(w)))
        if w:
            return w


@pytest.mark.parametrize("surface, punctured, rank, name, images, is_twist", MAPPING_CLASS_CASES,
                         ids=[f"{case[0]}-{case[3]}" for case in MAPPING_CLASS_CASES])
def test_mapping_class_invariance(surface, punctured, rank, name, images, is_twist):
    rng = random.Random(2)
    mor = morphism(images, rank)
    curves = [random_word(rng, rank, 12) for _ in range(30)]
    images = [apply(mor, u) for u in curves]
    gi = pairing(_SURFACES[surface](), punctured_faces=punctured)
    for u, phi_u in zip(curves, images):
        assert gi([phi_u]) == gi([u]), u
    for (u, phi_u), (v, phi_v) in itertools.combinations(zip(curves, images), 2):
        assert gi([phi_u], [phi_v]) == gi([u], [v]), (u, v)
