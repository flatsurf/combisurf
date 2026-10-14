r"""
Strongly reduced closed walks of bipartite hyperbolic triangulations, on the
triangle systems of the maps of genus 2 to 4 of the triangle system tests,
closed or with one to three punctured triangles, and on the barycentric
subdivisions of the barycentric subdivision tests that are hyperbolic. Free
homotopy is decided by the right-most geodesics of the hyperbolic
quadrangulation of the quad system of each triangulation. The Cython
reduction and its pure Python reference return the same walks, rotation
included.
"""
import random
from array import array

import pytest

from combisurf import OrientedMap
from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
from combisurf.bipartite_hyperbolic_triangulation_reduction import strongly_reduced as cython_strongly_reduced
from combisurf.bipartite_hyperbolic_triangulation_reduction_naive import TurnRuns
from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation

import test_barycentric_subdivision
import test_triangle_system
from test_triangle_system import random_closed_walk

LENGTHS = (1, 3, 8, 20, 60)


def instances():
    r"""
    Return the list of the triples ``(name, t, punctured)``.
    """
    result = []
    for i, (name, m, choice) in enumerate(test_triangle_system.CASES):
        t, _ = m.triangle_system(*choice, relabel=True)
        result.append((f"closed-{name}-{i % 3}", t, []))
        rng = random.Random(i)
        punctured = sorted(rng.sample(range(t.num_faces()), rng.randint(1, 3)))
        result.append((f"punctured-{name}-{i % 3}", t, punctured))
    for i, (name, m, punctured) in enumerate(test_barycentric_subdivision.CASES):
        t, t_punctured = m.barycentric_subdivision(punctured_faces=punctured)
        try:
            BipartiteHyperbolicTriangulation(t, t_punctured)
        except ValueError:
            continue
        result.append((f"barycentric-{test_barycentric_subdivision.IDS[i]}", t, t_punctured))
    return result


INSTANCES = instances()
IDS = [name for name, _, _ in INSTANCES]


def test_number_of_instances():
    assert sum(name.startswith("closed") for name in IDS) == 63
    assert sum(name.startswith("punctured") for name in IDS) == 63
    assert sum(name.startswith("barycentric") for name in IDS) == 56


def is_rotation(u, v):
    r"""
    Return whether the lists ``u`` and ``v`` are cyclic rotations of each
    other, two empty lists included.
    """
    return len(u) == len(v) and (not u or any(u[i:] + u[:i] == v for i in range(len(u))))


def reverse(w):
    return [h ^ 1 for h in reversed(w)]


def naive_reduced(T, w):
    r"""
    Return the strongly reduced closed walk of ``w`` computed by the pure
    Python reference.
    """
    return TurnRuns(T, w).strongly_reduced() if w else []


def check_same(T, w):
    r"""
    Check that the Cython reduction of ``w``, through
    :meth:`BipartiteHyperbolicTriangulation.strongly_reduced` and directly on
    the tables, equals the one of the reference, and return it.
    """
    r = T.strongly_reduced(w)
    assert r == naive_reduced(T, w)
    assert cython_strongly_reduced(T._h2v, T._wpos, T._wdeg, T._offset, T._at, array('i', w)) == r
    return r


class Oracle:
    r"""
    Free homotopy of closed walks of ``t`` through the right-most geodesics
    of the hyperbolic quadrangulation of its quad system.
    """
    def __init__(self, t, punctured):
        q, self.proj, bigons = t.quad_system(relabel=True, mapping=True, punctured_faces=punctured)
        self.Q = HyperbolicQuadrangulation(q, bigons)

    def __call__(self, w):
        w = [x for h in w for x in self.proj[h]]
        return self.Q.geodesic(w, canonical="right") if w else []

    def homotopic(self, u, v):
        return is_rotation(self(u), self(v))


def filled_half_edges_at(t, punctured):
    r"""
    Return a dictionary from the vertices of ``t`` to the half-edges that
    start there and whose face is filled.
    """
    h2v = t.half_edge_to_vertex()
    h2f = t.half_edge_to_face()
    out = {}
    for h in t.half_edges():
        if h2f[h] not in punctured:
            out.setdefault(h2v[h], []).append(h)
    return out


def face_from(t, h):
    r"""
    Return the boundary of the face of ``h``, starting with ``h``.
    """
    f = [h]
    x = t.next_in_face(h)
    while x != h:
        f.append(x)
        x = t.next_in_face(x)
    return f


def variants(t, punctured, w, rng):
    r"""
    Return closed walks freely homotopic to ``w``: a rotation, a detour around
    a filled triangle, a backtrack and both.
    """
    h2v = t.half_edge_to_vertex()
    filled = filled_half_edges_at(t, punctured)
    n = len(w)
    out = []
    i = rng.randrange(n)
    out.append(w[i:] + w[:i])
    v = h2v[w[i]]
    if v in filled:
        detour = face_from(t, rng.choice(filled[v]))
        out.append(w[:i] + detour + w[i:])
    h = rng.choice([x for x in t.half_edges() if h2v[x] == v])
    backtrack = w[:i] + [h, h ^ 1] + w[i:]
    out.append(backtrack)
    j = rng.randrange(len(backtrack))
    out.append(backtrack[j:] + backtrack[:j])
    return out


@pytest.mark.parametrize("name, t, punctured", INSTANCES, ids=IDS)
def test_strongly_reduced(name, t, punctured):
    T = BipartiteHyperbolicTriangulation(t, punctured)
    oracle = Oracle(t, punctured)
    rng = random.Random(len(t.vertex_permutation()) + len(punctured))
    walks = [random_closed_walk(t, rng, length) for length in LENGTHS]
    walks.append(walks[-2] + walks[-2])
    for w in walks:
        r = check_same(T, w)
        assert T.is_strongly_reduced(r)
        assert len(r) <= len(w)
        assert oracle.homotopic(r, w)
        assert (r == []) == (oracle(w) == [])
        assert is_rotation(check_same(T, r), r)
        assert is_rotation(check_same(T, reverse(w)), reverse(r))
        for u in variants(t, punctured, w, rng):
            assert is_rotation(check_same(T, u), r)


@pytest.mark.parametrize("name, t, punctured", INSTANCES, ids=IDS)
def test_contractible(name, t, punctured):
    T = BipartiteHyperbolicTriangulation(t, punctured)
    h2v = t.half_edge_to_vertex()
    filled = filled_half_edges_at(t, punctured)
    for i, f in enumerate(t.faces()):
        if i not in punctured:
            assert check_same(T, f) == []
            assert check_same(T, reverse(f)) == []
        else:
            check_same(T, f)
            check_same(T, reverse(f))
    rng = random.Random(t.num_faces())
    for length in LENGTHS:
        c = random_closed_walk(t, rng, length)
        v = h2v[c[0]]
        if v not in filled:
            continue
        f = face_from(t, rng.choice(filled[v]))
        g = face_from(t, rng.choice(filled[v]))
        w = c + f + reverse(c)
        assert check_same(T, w) == []
        # a product of two conjugates of boundaries of triangles, all based
        # at the start of c
        u = w + reverse(c) + reverse(g) + c
        assert check_same(T, u) == []


def weighted_turns_with_weight(t, punctured, w, weight):
    r"""
    Return the weighted turns of the closed walk ``w`` when a punctured corner
    weighs ``weight``.
    """
    h2f = t.half_edge_to_face()
    pos = {}
    wdeg = {}
    for i, vertex in enumerate(t.vertices()):
        p = 0
        for h in vertex:
            pos[h] = p
            p += weight if h2f[h] in punctured else 1
        for h in vertex:
            wdeg[h] = p
    out = []
    for i, h in enumerate(w):
        g = h ^ 1
        x = (pos[w[(i + 1) % len(w)]] - pos[g]) % wdeg[g]
        out.append(x - wdeg[g] if 2 * x > wdeg[g] else x)
    return out


def test_weight_three_fails():
    r"""
    The boundary of a white punctured triangle, with the triangle on its left.

    At each corner the walk turns around the punctured corner alone, on its
    left, the vertex having all its other corners, filled, on its right. With
    the weight `4` of a punctured corner the weighted turns are `-4`: the
    walk is strongly reduced, and it is not contractible since it goes
    around the puncture. With the weight `3`, they would read as white
    weighted turns `-3`, the turns `3_r` of the article, so that the walk
    would not be strongly reduced: the last move of the reduction would
    replace it by the other side of a cylinder of triangles on its left,
    which is the punctured triangle itself, pushing the walk across the
    puncture (the positions inside the punctured corner hold no half-edge).
    """
    name, m, choice = test_triangle_system.CASES[0]
    t, _ = m.triangle_system(*choice, relabel=True)
    punctured = [0]
    w = t.faces()[0]
    assert all(h % 2 == 0 for h in w)
    T = BipartiteHyperbolicTriangulation(t, punctured)
    assert T.weighted_turns(w) == [-4, -4, -4]
    assert weighted_turns_with_weight(t, punctured, w, 4) == [-4, -4, -4]
    assert weighted_turns_with_weight(t, punctured, w, 3) == [-3, -3, -3]
    assert T.is_strongly_reduced(w)
    assert is_rotation(T.strongly_reduced(w), w)
    oracle = Oracle(t, punctured)
    assert oracle(w)
    assert oracle.homotopic(T.strongly_reduced(w), w)
    # the reverse walk, with the triangle on its right, has the black
    # weighted turns 4, which the weight 3 would read as the black weighted
    # turns 3, the turns -3_b of the article
    assert T.weighted_turns(reverse(w)) == [4, 4, 4]
    assert weighted_turns_with_weight(t, punctured, reverse(w), 3) == [3, 3, 3]
    assert is_rotation(T.strongly_reduced(reverse(w)), reverse(w))


MAX_STEPS = 10000
MAX_TRIES = 200


class MoveCounter:
    r"""
    Record the cylinder moves and the flips ``(s, k, l)`` of the reductions,
    by wrapping the methods of :class:`TurnRuns` through ``monkeypatch``.
    """
    def __init__(self, monkeypatch):
        self.cylinders = 0
        self.flips = []
        flip = TurnRuns._flip
        cylinder = TurnRuns._cylinder
        counter = self

        def counted_flip(runs, u, s):
            val, length = runs._val, runs._len
            p, q = runs._prev[u], runs._next[u]
            k = length[p] if val[p] == 3 * s else 0
            l = length[q] if val[q] == 3 * s else 0
            r = flip(runs, u, s)
            if r != -1:
                counter.flips.append((s, k, l))
            return r

        def counted_cylinder(runs, u):
            counter.cylinders += 1
            return cylinder(runs, u)

        monkeypatch.setattr(TurnRuns, "_flip", counted_flip)
        monkeypatch.setattr(TurnRuns, "_cylinder", counted_cylinder)


def half_edges_at(t):
    r"""
    Return a dictionary from the vertices of ``t`` to the half-edges that
    start there.
    """
    h2v = t.half_edge_to_vertex()
    out = {}
    for h in t.half_edges():
        out.setdefault(h2v[h], []).append(h)
    return out


def next_at_turn(T, out, h, turn):
    r"""
    Return the half-edge ``h2`` such that the weighted turn from ``h`` to
    ``h2`` is ``turn``, or ``None`` if there is none.
    """
    h2v = T.map().half_edge_to_vertex()
    for x in out[h2v[h ^ 1]]:
        if T.weighted_turn(h, x) == turn:
            return x
    return None


def cylinder_walk(T, out, h, s):
    r"""
    Return the closed walk from ``h`` whose weighted turns are all ``3 s``,
    or ``None`` if the walk meets a punctured corner.

    The map from a half-edge to the next one is injective, so that the walk
    comes back to ``h``; a weighted turn `\pm 3` keeps the parity of the
    half-edges, hence the colour of the turns.
    """
    w = [h]
    for _ in range(MAX_STEPS):
        x = next_at_turn(T, out, w[-1], 3 * s)
        if x is None:
            return None
        if x == h:
            return w
        w.append(x)
    raise RuntimeError(f"no closed walk from {h} in {MAX_STEPS} steps")


@pytest.mark.parametrize("name, t, punctured", INSTANCES[::6], ids=IDS[::6])
def test_cylinder_move(name, t, punctured, monkeypatch):
    # Fig. 3.7 of Colin de Verdière, Despré and Dubois, "Untangling graphs
    # on surfaces" (2024): a closed walk whose weighted turns are all
    # white -3 or all black 3 goes around a cylinder of triangles and is
    # replaced by the other side of the cylinder. It is reached directly and
    # from homotopic walks with a detour or a backtrack.
    T = BipartiteHyperbolicTriangulation(t, punctured)
    oracle = Oracle(t, punctured)
    counter = MoveCounter(monkeypatch)
    out = half_edges_at(t)
    rng = random.Random(len(t.vertex_permutation()))
    num_walks = 0
    for s, parity in ((-1, 0), (1, 1)):
        starts = [h for h in t.half_edges() if h % 2 == parity]
        for h in rng.sample(starts, min(3, len(starts))):
            w = cylinder_walk(T, out, h, s)
            if w is None:
                continue
            num_walks += 1
            assert all(x % 2 == parity for x in w)
            assert T.weighted_turns(w) == [3 * s] * len(w)
            assert not T.is_strongly_reduced(w)
            before = counter.cylinders
            # the moves are counted on the reference, which the Cython
            # reduction follows move for move
            r = TurnRuns(T, w).strongly_reduced()
            assert counter.cylinders == before + 1
            assert check_same(T, w) == r
            assert r and T.is_strongly_reduced(r)
            assert T.weighted_turns(r) == [-3 * s] * len(r)
            assert oracle.homotopic(r, w)
            assert is_rotation(T.strongly_reduced(reverse(w)), reverse(r))
            for u in variants(t, punctured, w, rng):
                assert is_rotation(T.strongly_reduced(u), r)
    if not punctured:
        # a closed surface has no punctured corner to stop a cylinder walk
        assert num_walks == 6


def flip_walk(T, out, h, s, k, l, rng):
    r"""
    Return a closed walk whose weighted turns read cyclically ``x, (3s)^k, 2s,
    (3s)^l, y`` with a white ``2s`` and ``s x``, ``s y`` not in `\{-1, 0, 1,
    2, 3\}`, from the even half-edge ``h``, or ``None`` if none is found.

    The pattern is followed from ``h``, then a random walk goes back to the
    start of ``h``.
    """
    path = [h]
    for turn in [3 * s] * k + [2 * s] + [3 * s] * l:
        x = next_at_turn(T, out, path[-1], turn)
        if x is None:
            return None
        path.append(x)
    t = T.map()
    h2v = t.half_edge_to_vertex()
    start = h2v[h]
    for _ in range(MAX_TRIES):
        w = list(path)
        for _ in range(rng.randint(1, 6)):
            w.append(rng.choice([x for x in out[h2v[w[-1] ^ 1]] if x != w[-1] ^ 1]))
        closing = [x for x in out[h2v[w[-1] ^ 1]] if h2v[x ^ 1] == start and x != w[-1] ^ 1]
        if not closing:
            continue
        w.append(rng.choice(closing))
        turns = T.weighted_turns(w)
        x = turns[-1]
        y = turns[k + l + 1]
        if s * x not in (-1, 0, 1, 2, 3) and s * y not in (-1, 0, 1, 2, 3):
            assert turns[:k + l + 1] == [3 * s] * k + [2 * s] + [3 * s] * l
            return w
    return None


def test_long_flips(monkeypatch):
    # flips x, (3s)^k, 2s, (3s)^l, y with k, l >= 2, on both sides s, from
    # walks whose weighted turns contain the pattern
    counter = MoveCounter(monkeypatch)
    rng = random.Random(3)
    num_walks = 0
    for name, t, punctured in INSTANCES[::4]:
        T = BipartiteHyperbolicTriangulation(t, punctured)
        oracle = Oracle(t, punctured)
        out = half_edges_at(t)
        evens = [h for h in t.half_edges() if h % 2 == 0]
        for s in (-1, 1):
            for k, l in ((2, 2), (2, 3), (3, 2)):
                w = flip_walk(T, out, rng.choice(evens), s, k, l, rng)
                if w is None:
                    continue
                num_walks += 1
                before = len(counter.flips)
                r = TurnRuns(T, w).strongly_reduced()
                # the reduction of w itself takes a flip with k, l >= 2
                assert any(kk >= 2 and ll >= 2 for _, kk, ll in counter.flips[before:])
                assert check_same(T, w) == r
                assert T.is_strongly_reduced(r)
                assert oracle.homotopic(r, w)
                assert is_rotation(T.strongly_reduced(reverse(w)), reverse(r))
                for u in variants(t, punctured, w, rng):
                    assert is_rotation(T.strongly_reduced(u), r)
    assert num_walks >= 100
    long_flips = [(s, k, l) for s, k, l in counter.flips if k >= 2 and l >= 2]
    assert any(s == -1 for s, _, _ in long_flips)
    assert any(s == 1 for s, _, _ in long_flips)


@pytest.mark.parametrize("name, t, punctured", INSTANCES[::40], ids=IDS[::40])
def test_long_walks(name, t, punctured):
    # walks of 10^4 turns, their reverses and their squares, on which the
    # run arrays of the Cython reduction are grown many times
    T = BipartiteHyperbolicTriangulation(t, punctured)
    rng = random.Random(len(t.vertex_permutation()))
    for _ in range(2):
        w = random_closed_walk(t, rng, 10000)
        r = check_same(T, w)
        assert T.is_strongly_reduced(r)
        check_same(T, reverse(w))
        check_same(T, w + w)


# a closed walk on the triangle system of a torus, whose vertices have degree
# 6 with no punctured corner, on which the moves do not terminate
TORUS_WALK = [12, 15, 7, 5, 17, 12, 15, 11, 6, 12, 17, 9, 10, 9, 5, 2, 1, 17, 11, 8, 16, 2, 6, 9, 1, 13, 14, 17]


def test_non_termination():
    t, _ = OrientedMap(fp="(0,1,~0,~1)").triangle_system(relabel=True)
    assert t.vertex_profile() == [6, 6, 6]
    T = BipartiteHyperbolicTriangulation(t, check=False)
    message = "the reduction does not terminate; the triangulation is not a bipartite hyperbolic triangulation"
    with pytest.raises(RuntimeError, match=message):
        TurnRuns(T, TORUS_WALK).strongly_reduced()
    with pytest.raises(RuntimeError, match=message):
        T.strongly_reduced(TORUS_WALK)
    # the error leaves no state behind: a second call raises it again
    with pytest.raises(RuntimeError, match=message):
        T.strongly_reduced(TORUS_WALK)


def test_errors():
    name, m, choice = test_triangle_system.CASES[0]
    T = BipartiteHyperbolicTriangulation(*m.triangle_system(*choice, relabel=True))
    with pytest.raises(NotImplementedError):
        T.strongly_reduced([0, 1], closed=False)
    with pytest.raises(ValueError):
        T.strongly_reduced([0])
    assert T.strongly_reduced([]) == []
    assert T.is_strongly_reduced([])
