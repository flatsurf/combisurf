r"""
Weighted turns and geodesics on hyperbolic punctured quadrangulations.

The instances are quad systems and radial maps, closed or punctured at their
bigons. On a punctured surface the fundamental group is free, and two closed
walks are freely homotopic exactly when their cyclically reduced words in the
reduced map are cyclic conjugates. On a closed surface only the homology
class is compared.
"""
import random

import pytest

from combisurf import OrientedMap
from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation, WeightedTurnRuns
from combisurf.word import word_init, word_cyclically_reduce

from test_geometric_intersection_oracles import (B1_V, B1_W, cylinder_data, genus2, inverse, is_cyclic_conjugate,
                                                 random_closed_walk)


def image(mor, w):
    return [x for h in w for x in mor[h]]


def genus2_two_hexagons():
    # two vertices of degree 6 and two faces of degree 6
    return OrientedMap(vp=[[0, 9, 10, 7, 8, 11], [2, 1, 4, 3, 6, 5]])


def genus2_pentagons():
    # four vertices of degree 5 and four faces of degree 5; its radial map
    # has eight vertices of degree 5
    return OrientedMap(vp=[[1, 2, 13, 16, 6], [9, 10, 15, 3, 11], [17, 8, 12, 7, 4], [19, 0, 14, 18, 5]])


def build(name):
    r"""
    Return ``(m, q, mor, bigons)`` for the instance ``name``: the walks of
    ``m`` map to walks of ``q`` through ``mor``.
    """
    surface, construction, punctured = INSTANCES[name][:3]
    if surface == "genus2":
        m = genus2()
    elif surface == "hexagons":
        m = genus2_two_hexagons()
    elif surface == "pentagons":
        m = genus2_pentagons()
    else:
        m = OrientedMap(vp=cylinder_data(surface)["vp"])
    if construction == "quad_system":
        q, mor, bigons = m.quad_system(mapping=True, punctured_faces=punctured)
    else:
        q, mor, bigons = m.radial_map(mapping=True, punctured_faces=punctured)
    return m, q, mor, bigons


# name: surface, construction, punctured faces, (V, E, F), sorted vertex
# degrees, number of bigons
INSTANCES = {
    "genus2": ("genus2", "quad_system", False, (2, 8, 4), [8, 8], 0),
    "genus2_punctured": ("genus2", "quad_system", True, (2, 9, 5), [9, 9], 1),
    "lshape": ("lshape", "quad_system", False, (2, 8, 4), [8, 8], 0),
    "lshape_01": ("lshape", "quad_system", [0, 1], (2, 10, 6), [10, 10], 2),
    "lshape_punctured": ("lshape", "quad_system", True, (2, 14, 10), [14, 14], 6),
    "emmw_1117": ("emmw_1117", "quad_system", False, (2, 16, 8), [16, 16], 0),
    "emmw_1117_punctured": ("emmw_1117", "quad_system", True, (2, 36, 28), [36, 36], 20),
    "hexagons_radial": ("hexagons", "radial_map", False, (4, 12, 6), [6, 6, 6, 6], 0),
    "hexagons_radial_0": ("hexagons", "radial_map", [0], (4, 13, 7), [6, 6, 7, 7], 1),
    "lshape_radial_punctured": ("lshape", "radial_map", True, (7, 24, 15), [4] * 6 + [24], 6),
    "pentagons_radial": ("pentagons", "radial_map", False, (8, 20, 10), [5] * 8, 0),
    "pentagons_radial_0": ("pentagons", "radial_map", [0], (8, 21, 11), [5] * 6 + [6, 6], 1),
}
NAMES = list(INSTANCES)
CLOSED = [name for name in NAMES if INSTANCES[name][5] == 0]
PUNCTURED = [name for name in NAMES if INSTANCES[name][5] > 0]


class Instance:
    r"""
    An instance with the data used by the checks: the punctured corners,
    the weighted degrees, and the projection to the reduced map that gives
    free words (punctured) or homology classes (closed).
    """
    def __init__(self, name):
        self.name = name
        self.m, self.q, self.mor, self.bigons = build(name)
        self.Q = HyperbolicQuadrangulation(self.q, self.bigons)
        q = self.q
        self.vp = q.vertex_permutation()
        self.fp = q.face_permutation()
        self.h2v = q.half_edge_to_vertex()
        self.h2f = q.half_edge_to_face()
        self.half_edges = list(q.half_edges())
        self.out = {}
        for h in self.half_edges:
            self.out.setdefault(self.h2v[h], []).append(h)
        self.punctured_corner = {h: self.h2f[h] in self.bigons for h in self.half_edges}
        self.weighted_degree = {v: sum(4 if self.punctured_corner[h] else 1 for h in hs) for v, hs in self.out.items()}
        if self.bigons:
            forest, coforest, _ = q.forest_coforest_decomposition(root_faces=self.bigons)
        else:
            forest, coforest, _ = q.forest_coforest_decomposition()
        _, self.proj = q.reduced_map(forest, coforest, mapping=True)

    def free_word(self, w):
        return list(word_cyclically_reduce(word_init(image(self.proj, w))))

    def homology(self, w):
        c = {}
        for h in image(self.proj, w):
            c[h // 2] = c.get(h // 2, 0) + (-1 if h % 2 else 1)
        return {e: x for e, x in c.items() if x}

    def homotopic(self, u, v):
        r"""
        Whether the closed walks ``u`` and ``v`` are freely homotopic
        (punctured), or homologous (closed).
        """
        if self.bigons:
            return is_rotation(self.free_word(u), self.free_word(v))
        return self.homology(u) == self.homology(v)

    def link(self, x):
        r"""
        Return the closed walk around the vertex of ``x``, which leaves the
        vertex on its left and starts at the end of ``x``.

        For each corner, from the corner between ``x`` and ``vp[x]`` on, it
        goes along the sides of the face of the corner that do not touch the
        vertex.
        """
        fp = self.fp
        vp = self.vp
        w = []
        y = x
        while True:
            z = fp[y]
            while z != vp[y] ^ 1:
                w.append(z)
                z = fp[z]
            y = vp[y]
            if y == x:
                return w

    def random_walks(self, rng, count, max_length):
        q = self.q
        walks = [random_closed_walk(q, rng, rng.randint(1, max_length)) for _ in range(count)]
        walks += [image(self.mor, random_closed_walk(self.m, rng, rng.randint(1, max_length // 2)))
                  for _ in range(count)]
        return [w for w in walks if w]


_INSTANCE_CACHE = {}


def instance(name):
    if name not in _INSTANCE_CACHE:
        _INSTANCE_CACHE[name] = Instance(name)
    return _INSTANCE_CACHE[name]


@pytest.mark.parametrize("name", NAMES)
def test_instances(name):
    # guards the instances against a change of the constructions
    I = instance(name)
    q = I.q
    _, _, _, counts, degrees, num_bigons = INSTANCES[name]
    assert (q.num_vertices(), q.num_edges(), q.num_faces()) == counts
    assert sorted(q.vertex_profile()) == degrees
    assert len(I.bigons) == num_bigons
    assert all(len(q.faces()[i]) == 2 for i in I.bigons)
    touching = sorted(v for v, hs in I.out.items() if any(I.punctured_corner[h] for h in hs))
    if name == "hexagons_radial_0":
        assert sorted(I.weighted_degree.values()) == [6, 6, 10, 10]
        assert len(touching) == 2
    if name == "lshape_radial_punctured":
        assert len(touching) == 7


def test_rejected():
    from test_geometric_intersection_oracles import torus

    q, _ = torus().quad_system(punctured_faces=False)
    assert (q.num_vertices(), q.num_edges(), q.num_faces()) == (2, 4, 2)
    with pytest.raises(ValueError, match="Euler characteristic"):
        HyperbolicQuadrangulation(q)

    tetrahedron = OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~4,~5)(~2,5,~3)")
    assert tetrahedron.euler_characteristic() == 2
    with pytest.raises(ValueError):
        HyperbolicQuadrangulation(tetrahedron)
    with pytest.raises(ValueError, match="degree 3 instead of 4"):
        HyperbolicQuadrangulation(tetrahedron, [0, 1, 2], check=True)

    square = OrientedMap(fp="(0,1,2,3)(~0,~3,~2,~1)")
    with pytest.raises(ValueError, match="Euler characteristic"):
        HyperbolicQuadrangulation(square, [1])  # a disk
    with pytest.raises(ValueError, match="Euler characteristic"):
        HyperbolicQuadrangulation(square, True)  # an annulus

    # a vertex of degree four that touches no punctured face
    m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
    r, _ = m.radial_map()
    with pytest.raises(ValueError, match="degree 4 instead of at least 5"):
        HyperbolicQuadrangulation(r)

    with pytest.raises(NotImplementedError):
        HyperbolicQuadrangulation(OrientedMap(vp="(0,1,~1)"), True)
    with pytest.raises(TypeError):
        HyperbolicQuadrangulation(q.vertex_permutation())


@pytest.mark.parametrize("name", NAMES)
def test_weighted_turn_definition(name):
    I = instance(name)
    Q = I.Q
    vp = I.vp
    for v, hs in I.out.items():
        d = len(hs)
        wd = I.weighted_degree[v]
        for g in hs:
            h1 = g ^ 1
            ring = [g]
            while vp[ring[-1]] != g:
                ring.append(vp[ring[-1]])
            assert len(ring) == d
            for k, h2 in enumerate(ring):
                # the corners swept counterclockwise from g to h2, and the
                # others; True for a punctured corner
                ccw = [I.punctured_corner[x] for x in ring[:k]]
                cw = [I.punctured_corner[x] for x in ring[k:]]
                assert Q.turn(h1, h2) == k
                weight_ccw = sum(4 if p else 1 for p in ccw)
                weight_cw = sum(4 if p else 1 for p in cw)
                assert weight_ccw + weight_cw == wd
                t = weight_ccw % wd
                if 2 * t > wd:
                    t -= wd
                wt = Q.weighted_turn(h1, h2)
                assert wt == t
                assert -wd < 2 * wt <= wd
                # the integer turn of Despré-Lazarus counts the corners on the
                # side with no punctured corner; at a vertex without punctured
                # corner of even degree, the two sides tie at half the degree
                if weight_ccw == weight_cw:
                    assert wt == weight_ccw
                elif not any(ccw) and len(ccw) <= 3 and weight_ccw < weight_cw:
                    assert wt == len(ccw), (h1, h2)
                elif not any(cw) and len(cw) <= 3 and weight_cw < weight_ccw:
                    assert wt == -len(cw), (h1, h2)
                if any(ccw) and any(cw):
                    assert abs(wt) >= 4, (h1, h2)
                assert (wt == 0) == (h2 == g)

    with pytest.raises(ValueError, match="does not end where"):
        h = I.half_edges[0]
        Q.turn(h, h)


def rebuild(I, h, turns):
    r"""
    Return the walk that starts with ``h`` and has the weighted turns
    ``turns``, by search in the half-edges of each vertex.
    """
    w = [h]
    for t in turns[:-1]:
        g = [x for x in I.out[I.h2v[h ^ 1]] if I.Q.weighted_turn(h, x) == t]
        assert len(g) == 1
        h = g[0]
        w.append(h)
    return w


@pytest.mark.parametrize("name", NAMES)
def test_weighted_turns_rebuild(name):
    I = instance(name)
    rng = random.Random(10)
    for w in I.random_walks(rng, 10, 30):
        t = I.Q.weighted_turns(w)
        assert len(t) == len(w)
        assert rebuild(I, w[0], t) == w


def is_rotation(u, v):
    return u == v or is_cyclic_conjugate(u, v)


def check_geodesic(I, w, g):
    r"""
    Check that ``g`` is a geodesic homotopic to the closed walk ``w``.
    """
    Q = I.Q
    I.q._check_walk(g, closed=True)
    assert Q.is_geodesic(g), (w, g)
    assert len(g) <= len(w)
    assert is_rotation(Q.geodesic(g), g), g
    assert I.homotopic(w, g), (w, g)
    if I.bigons:
        assert (not g) == (not I.free_word(w)), (w, g)
    elif not g:
        assert not I.homology(w)


@pytest.mark.parametrize("name", NAMES)
def test_geodesic_random(name):
    I = instance(name)
    rng = random.Random(0)
    for w in I.random_walks(rng, 20, 40):
        check_geodesic(I, w, I.Q.geodesic(w))


def detours(I, rng, w):
    r"""
    Return closed walks freely homotopic to ``w`` by construction: a
    rotation of ``w``, and ``w`` with a spur, the boundary of a filled face
    and the link of a vertex without punctured corner inserted at a random
    position.
    """
    i = rng.randrange(len(w))
    out = I.out[I.h2v[w[i]]]
    variants = [w[i:] + w[:i]]
    x = rng.choice(out)
    variants.append(w[:i] + [x, x ^ 1] + w[i:])
    filled = [x for x in out if I.h2f[x] not in I.bigons]
    if filled:
        x = rng.choice(filled)
        face = [x]
        while I.fp[face[-1]] != x:
            face.append(I.fp[face[-1]])
        variants.append(w[:i] + face + w[i:])
    # the link of a neighbour without punctured corner, through the
    # starting vertex of w[i]
    neighbours = [x for x in out if not any(I.punctured_corner[y] for y in I.out[I.h2v[x ^ 1]])]
    if neighbours:
        x = rng.choice(neighbours)
        variants.append(w[:i] + I.link(x ^ 1) + w[i:])
    return variants


def random_homotopic(I, rng, w, moves):
    r"""
    Return a random closed walk freely homotopic to ``w``, obtained by
    ``moves`` successive elementary moves of :func:`detours`.
    """
    u = w
    for _ in range(moves):
        u = rng.choice(detours(I, rng, u))
    return u


@pytest.mark.parametrize("name", NAMES)
def test_geodesic_detours(name):
    # homotopic geodesics have the same length: Despré-Lazarus, Corollary 11
    I = instance(name)
    Q = I.Q
    rng = random.Random(1)
    for w in I.random_walks(rng, 8, 30):
        g = Q.geodesic(w)
        for u in detours(I, rng, w):
            I.q._check_walk(u, closed=True)
            h = Q.geodesic(u)
            assert len(h) == len(g), (w, u, g, h)
            check_geodesic(I, u, h)
            assert I.homotopic(g, h)


@pytest.mark.parametrize("name", NAMES)
def test_geodesic_powers(name):
    I = instance(name)
    Q = I.Q
    rng = random.Random(2)
    for w in I.random_walks(rng, 5, 20):
        g = Q.geodesic(w)
        for k in (2, 3):
            assert len(Q.geodesic(w * k)) == k * len(g)
            assert Q.is_geodesic(g * k)


@pytest.mark.parametrize("name", NAMES)
def test_geodesic_corner_cases(name):
    I = instance(name)
    Q = I.Q
    faces = I.q.faces()
    for i, face in enumerate(faces):
        if i in I.bigons:
            # peripheral, and a geodesic: its weighted turns are at least 4
            # in absolute value
            assert Q.is_geodesic(face)
            assert Q.geodesic(face) == face
        else:
            assert Q.geodesic(face) == []
            assert Q.geodesic(inverse(face)) == []
    for h in I.half_edges:
        assert Q.geodesic([h, h ^ 1]) == []
        if I.h2v[h] == I.h2v[h ^ 1]:
            assert Q.geodesic([h]) == [h]
    for v, hs in I.out.items():
        corners = [x for x in hs if I.punctured_corner[x]]
        link = I.link(hs[0])
        g = Q.geodesic(link)
        if not corners:
            # contractible
            assert g == []
        elif len(corners) == 1:
            # around the vertex and the puncture of its punctured corner:
            # homotopic to the boundary of that face
            face = faces[I.h2f[corners[0]]]
            assert len(g) == len(face)
            assert I.homotopic(g, face)


@pytest.mark.parametrize("name", [name for name in PUNCTURED if INSTANCES[name][5] >= 2])
def test_both_sides_punctured(name):
    # a weighted turn that sweeps a punctured corner on each side is never
    # reduced
    I = instance(name)
    Q = I.Q
    vp = I.vp

    def both_sides(h1, h2):
        g = h1 ^ 1
        ccw = []
        x = g
        while x != h2:
            ccw.append(I.punctured_corner[x])
            x = vp[x]
        cw = []
        while x != g:
            cw.append(I.punctured_corner[x])
            x = vp[x]
        return any(ccw) and any(cw)

    for i in I.bigons:
        face = I.q.faces()[i]
        assert both_sides(face[0], face[1])
        assert Q.geodesic(face) == face

    found = 0
    rng = random.Random(3)
    for w in I.random_walks(rng, 20, 40):
        g = Q.geodesic(w)
        n = len(g)
        if any(both_sides(g[j], g[(j + 1) % n]) for j in range(n)):
            found += 1
            assert is_rotation(Q.geodesic(g), g)
    assert found


@pytest.mark.parametrize("surface, num_curves", [("lshape", 26), ("emmw_1117", 50)])
def test_cylinder_sides(surface, num_curves):
    # O4: the two sides of a cylinder are homotopic on the closed surface
    data = cylinder_data(surface)
    m = OrientedMap(vp=data["vp"])
    q, mor, bigons = m.quad_system(mapping=True)
    Q = HyperbolicQuadrangulation(q, bigons)
    assert len(data["left"]) == len(data["right"]) == num_curves
    for u, v in zip(data["left"], data["right"]):
        gu = Q.geodesic(image(mor, u))
        gv = Q.geodesic(image(mor, v))
        assert gu and len(gu) == len(gv), (u, v, gu, gv)


def test_cohen_lustig_remark2():
    # O3: the words V and W of the row B1 of Cohen-Lustig are the same closed
    # curve on the closed genus two surface
    q, mor, bigons = genus2().quad_system(mapping=True)
    Q = HyperbolicQuadrangulation(q, bigons)
    V = image(mor, B1_V)
    W = image(mor, B1_W)
    assert len(V) == len(W) == 28
    gV = Q.geodesic(V)
    gW = Q.geodesic(W)
    assert gV and len(gV) == len(gW)
    assert Q.is_geodesic(gV) and Q.is_geodesic(gW)


def step(I, h, t):
    r"""
    Return the half-edge after ``h`` with weighted turn ``t``, or ``None``
    if there is none.
    """
    g = [x for x in I.out[I.h2v[h ^ 1]] if I.Q.weighted_turn(h, x) == t]
    assert len(g) <= 1, (h, t)
    return g[0] if g else None


def cyclic_brackets(I):
    r"""
    Return the closed walks of ``I`` with weighted turns ``s, 2s, ..., 2s``
    for ``s = ±1``, one per starting half-edge.

    Weighted turns `2s` from a given half-edge determine the walk, and two
    half-edges are never sent to the same one, so that following them from
    ``step(h, s)`` either comes back to ``h`` or never meets it.
    """
    walks = []
    for s in (1, -1):
        for h in I.half_edges:
            w = [h]
            x = step(I, h, s)
            while x is not None and x != h and len(w) <= len(I.half_edges):
                w.append(x)
                x = step(I, x, 2 * s)
            if x == h:
                walks.append((s, w))
    return walks


@pytest.mark.parametrize("name", NAMES)
def test_cyclic_brackets(name):
    I = instance(name)
    Q = I.Q
    for s, w in cyclic_brackets(I):
        t = Q.weighted_turns(w)
        assert t == [s] + [2 * s] * (len(w) - 1)
        assert len(w) >= 3
        assert not Q.is_geodesic(w)
        g = Q.geodesic(w)
        check_geodesic(I, w, g)


def test_cyclic_brackets_exist():
    assert any(cyclic_brackets(instance(name)) for name in NAMES)


def near_cyclic_brackets(I, max_r):
    r"""
    Return the closed walks of ``I`` with weighted turns ``s, 2s, ..., 2s,
    s, x`` for ``s = ±1``, at most ``max_r`` turns `2s` and ``x`` neither
    `0` nor ``s``, one per starting half-edge and number of turns `2s`.

    A bracket then covers all the turns but ``x``: this is the near-cyclic
    case of the reduction. Excluding ``x = 0`` and ``x = s`` keeps a spur
    or a shorter bracket from being found first.
    """
    walks = []
    for s in (1, -1):
        for h in I.half_edges:
            w = [h]
            x = step(I, h, s)
            for r in range(max_r + 1):
                if x is None:
                    break
                y = step(I, x, s)
                if y is not None and I.h2v[y ^ 1] == I.h2v[h] and I.Q.weighted_turn(y, h) not in (0, s):
                    walks.append((s, r, w + [x, y]))
                w.append(x)
                x = step(I, x, 2 * s)
    return walks


@pytest.mark.parametrize("name", NAMES)
def test_near_cyclic_brackets(name):
    I = instance(name)
    Q = I.Q
    for s, r, w in near_cyclic_brackets(I, 6):
        t = Q.weighted_turns(w)
        assert t[:-1] == [s] + [2 * s] * r + [s]
        assert not Q.is_geodesic(w)
        g = Q.geodesic(w)
        assert len(g) <= len(w) - 2
        check_geodesic(I, w, g)


def test_near_cyclic_brackets_exist():
    assert any(near_cyclic_brackets(instance(name), 6) for name in NAMES)


def test_reductions_exercised():
    # every kind of elementary reduction happens on the random walks
    counts = {"_spur": 0, "_bracket": 0, "_cyclic_bracket": 0, "_replace_all": 0, "_splice": 0}
    saved = {f: getattr(WeightedTurnRuns, f) for f in counts}

    def counting(f):
        def wrapper(self, *args):
            counts[f] += 1
            return saved[f](self, *args)
        return wrapper

    try:
        for f in counts:
            setattr(WeightedTurnRuns, f, counting(f))
        for name in NAMES:
            I = instance(name)
            rng = random.Random(4)
            for w in I.random_walks(rng, 20, 40):
                I.Q.geodesic(w)
                for v in I.out:
                    I.Q.geodesic(I.link(I.out[v][0]))
            for _, w in cyclic_brackets(I):
                I.Q.geodesic(w)
    finally:
        for f in counts:
            setattr(WeightedTurnRuns, f, saved[f])
    assert all(counts.values()), counts


def linear_is_geodesic(t):
    r"""
    Whether the weighted turns ``t`` of a path have no ``0`` and no window
    ``s, 2s, ..., 2s, s`` for ``s = ±1``, by a scan of the windows.
    """
    if 0 in t:
        return False
    n = len(t)
    for s in (1, -1):
        for i in range(n):
            if t[i] != s:
                continue
            j = i + 1
            while j < n and t[j] == 2 * s:
                j += 1
            if j < n and t[j] == s:
                return False
    return True


@pytest.mark.parametrize("name", NAMES)
def test_paths(name):
    I = instance(name)
    Q = I.Q
    rng = random.Random(5)
    for w in I.random_walks(rng, 10, 30):
        t = Q.weighted_turns(w)
        for k in range(len(w) + 1):
            p = w[:k]
            tp = Q.weighted_turns(p, closed=False)
            assert tp == t[:max(k - 1, 0)], (w, k)
            assert Q.is_geodesic(p, closed=False) == linear_is_geodesic(tp), (w, k)
        # a path is checked as a walk, not as a closed walk
        if w[:-1] and I.h2v[w[0]] != I.h2v[w[-2] ^ 1]:
            assert Q.weighted_turns(w[:-1], closed=False) == t[:len(w) - 2]
            with pytest.raises(ValueError):
                Q.weighted_turns(w[:-1])


def test_paths_corner_cases():
    I = instance("genus2")
    Q = I.Q
    assert Q.weighted_turns([], closed=False) == []
    assert Q.is_geodesic([], closed=False)
    for h in I.half_edges:
        assert Q.weighted_turns([h], closed=False) == []
        assert Q.is_geodesic([h], closed=False)
        assert not Q.is_geodesic([h, h ^ 1], closed=False)
    with pytest.raises(NotImplementedError):
        Q.geodesic(I.half_edges[:1], closed=False)


# the weighted turn that a right-most (left-most) geodesic never has
SIDES = {"right": 1, "left": -1}


@pytest.mark.parametrize("name", NAMES)
def test_canonical_uniqueness(name):
    # one canonical form per free homotopy class: Erickson-Whittlesey, Lemma
    # 5.2; Despré-Lazarus, Theorem 8
    I = instance(name)
    Q = I.Q
    rng = random.Random(6)
    for w in I.random_walks(rng, 8, 30):
        variants = detours(I, rng, w) + [random_homotopic(I, rng, w, 10)]
        for side in SIDES:
            c = Q.geodesic(w, canonical=side)
            for u in variants:
                assert is_rotation(Q.geodesic(u, canonical=side), c), (side, w, u)


@pytest.mark.parametrize("name", NAMES)
def test_canonical_properties(name):
    I = instance(name)
    Q = I.Q
    rng = random.Random(7)
    for w in I.random_walks(rng, 20, 40):
        n = len(Q.geodesic(w))
        for side, s in SIDES.items():
            c = Q.geodesic(w, canonical=side)
            check_geodesic(I, w, c)
            assert len(c) == n
            t = Q.weighted_turns(c)
            assert s not in t, (side, w, c, t)
            assert not c or any(x != 2 * s for x in t), (side, w, c)
            assert Q.geodesic(c, canonical=side) == c


@pytest.mark.parametrize("name", NAMES)
def test_canonical_powers(name):
    # powers of a canonical closed geodesic are canonical: Despré-Lazarus,
    # Remark 2
    I = instance(name)
    Q = I.Q
    rng = random.Random(8)
    for w in I.random_walks(rng, 5, 20):
        for side in SIDES:
            c = Q.geodesic(w, canonical=side)
            for k in (2, 3):
                assert is_rotation(Q.geodesic(w * k, canonical=side), c * k)


def test_canonical_corner_cases():
    I = instance("genus2_punctured")
    Q = I.Q
    for side in SIDES:
        assert Q.geodesic([], canonical=side) == []
        for i in I.bigons:
            face = I.q.faces()[i]
            assert is_rotation(Q.geodesic(face, canonical=side), face)
            assert is_rotation(Q.geodesic(inverse(face), canonical=side), inverse(face))
        with pytest.raises(NotImplementedError):
            Q.geodesic(I.half_edges[:1], closed=False, canonical=side)
    for h in I.half_edges:
        assert Q.geodesic([h, h ^ 1], canonical="right") == []
    with pytest.raises(ValueError, match="canonical must be"):
        Q.geodesic([], canonical="up")


def shift_row(R, a):
    r"""
    Return the row of the table of right shifts (see
    ``WeightedTurnRuns._right_shift``) that applies to the run ``a`` of
    value `1`: 1 to 3 for a bend that does not cover the walk, according to
    ``(s > 0, t > 0)``, ``s = 0`` and ``t = 0``, 4 to 6 likewise for a bend
    that covers it, and 7 when it covers it with ``x = 3``.
    """
    turns = []
    u = a
    for _ in range(R._size):
        turns += [R._val[u]] * R._len[u]
        u = R._next[u]
    # the turns from the turn 1 of a
    n = len(turns)
    t = 0
    while turns[1 + t] == 2:
        t += 1
    s = 0
    while turns[-1 - s] == 2:
        s += 1
    cyclic = s + t + 2 == n
    if cyclic:
        # the vertex where the bend starts, at the weighted turn x
        v = R._h2v[R._first[a if s == 0 else R._prev[a]]]
        if (turns[1 + t] - 3) % R._wdeg[v] == 0:
            return 7
    row = 1 if s and t else 2 if not s else 3
    return row + 3 if cyclic else row


def test_right_shifts_exercised():
    # every row of the table of right shifts happens
    counts = {row: 0 for row in range(1, 9)}
    saved = WeightedTurnRuns._right_shift, WeightedTurnRuns._right_shift_cylinder

    def right_shift(self, a):
        counts[shift_row(self, a)] += 1
        return saved[0](self, a)

    def right_shift_cylinder(self, u):
        counts[8] += 1
        return saved[1](self, u)

    try:
        WeightedTurnRuns._right_shift = right_shift
        WeightedTurnRuns._right_shift_cylinder = right_shift_cylinder
        for name in NAMES:
            I = instance(name)
            rng = random.Random(9)
            for w in I.random_walks(rng, 20, 40):
                for side, s in SIDES.items():
                    assert s not in I.Q.weighted_turns(I.Q.geodesic(w, canonical=side))
    finally:
        WeightedTurnRuns._right_shift, WeightedTurnRuns._right_shift_cylinder = saved
    assert all(counts.values()), counts


def test_right_shift_sign():
    # a standard right shift x, 2, 1, 2, y on the radial map of the genus 2
    # surface with four pentagons, all of whose vertices have degree 5
    I = instance("pentagons_radial")
    Q = I.Q
    g = [0, 25, 34, 23, 6, 39]
    assert Q.weighted_turns(g) == [2, -2, -2, -1, 2, 1]
    assert Q.is_geodesic(g)
    c = Q.geodesic(g, canonical="right")
    assert is_rotation(c, [23, 30, 5, 32, 15, 34])
    # x = -1 from 23 to the bend 6, 39, 0, 25 and y = -2 after it; x - 1 =
    # -2 merges with the turn -2 from 34 to 23, and -3 and y - 1 = -3 are 2
    # at vertices of degree 5
    assert Q.weighted_turns([23, 30, 5, 32, 15, 34]) == [-2, -1, 2, -1, 2, -2]
    # the faces on the right of the bend: R_1 along 6, the corner face C
    # along 39 and 0, R'_1 along 25
    R1, C, R1p = I.h2f[6 ^ 1], I.h2f[39 ^ 1], I.h2f[25 ^ 1]
    assert I.h2f[0 ^ 1] == C and len({R1, C, R1p}) == 3
    # the new walk 30, 5, 32, 15 goes around R_1 and R'_1, on its left,
    # through the corner of C opposite to the turn 1
    assert [I.h2f[h] for h in [30, 5, 32, 15]] == [R1, R1, R1p, R1p]
    assert I.h2v[5 ^ 1] == I.h2v[32] and C in [I.h2f[h] for h in I.out[I.h2v[32]]]
    # the left-most geodesic shifts the other way
    assert Q.weighted_turns(Q.geodesic(g, canonical="left")).count(-1) == 0
