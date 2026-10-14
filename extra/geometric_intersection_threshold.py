r"""
Threshold between the two reductions of the default algorithm of
:class:`~combisurf.geometric_intersection.GeometricIntersectionPairing`.

With a punctured face the default algorithm builds the reduced map (rm) and
keeps it when its longest image of a half-edge, the rm L_inf, is at most
``_MAX_IMAGE_LENGTH``; otherwise it builds the triangle system (ts). This
script prints the data behind that constant.

It measures the strips ``torus_grid_rect(a, 8)`` (torus grids of 8 columns and
``a`` rows) for ``a`` in 12, 16, 20, 24, 28, 32, 40, each with three sets of
punctured faces (``half``, ``all_but_col0``, ``odd_cols``). For each group
(``a``, punctured set) it prints the rm L_inf and the ratio ts/rm of the pair
time, the geometric mean over random pairs of closed walks of length 1000 and
10000 of the best time of one call ``g([u], [v])``. Both algorithms must give
the same intersection number on every pair, otherwise it raises.

The spanning forest and coforest of the rm are depth first here. Those of the
package are breadth first, which keeps the rm L_inf of these strips at 1 to 4
and so never reaches the range where the two algorithms cross; the depth first
ones give an L_inf of 10 to 49.

The constant ``C`` is the largest measured rm L_inf ``x`` such that every
group with rm L_inf at most ``x`` has ts/rm at least 1. The groups on the wrong
side of ``C`` (ts/rm below 1 at L_inf at most ``C``, or above 1 beyond) are
printed as exceptions.

Usage::

    sage -python extra/geometric_intersection_threshold.py [--repeat K] [--pairs N]

``--repeat`` is the number of timed runs of each pair, the best being kept
(default 5), and ``--pairs`` the number of pairs for each of the two lengths
(default 3). ``_MAX_IMAGE_LENGTH`` of ``combisurf/geometric_intersection.py``
quotes this output.
"""
import math
import random
import sys
import time
from array import array
from collections import deque

from combisurf.oriented_map import OrientedMap
from combisurf.geometric_intersection import GeometricIntersectionPairing, _morphism_arrays
from combisurf.permutation import perm_dense_cycle_positions

ROWS = (12, 16, 20, 24, 28, 32, 40)
COLUMNS = 8
LENGTHS = (1000, 10000)
SETS = ("half", "all_but_col0", "odd_cols")


def torus_grid_rect(a, length):
    r"""
    Return the torus grid with ``length`` columns and ``a`` rows, and the
    function ``E(i, j, d)`` giving the even half-edge leaving the vertex
    ``(i, j)`` to the right (``d = 0``) or up (``d = 1``).
    """
    def E(i, j, d):
        return 2 * (2 * ((i % length) * a + (j % a)) + d)

    vp = []
    for i in range(length):
        for j in range(a):
            vp.append([E(i, j, 0), E(i, j, 1), E(i - 1, j, 0) + 1, E(i, j - 1, 1) + 1])
    return OrientedMap(vp=vp), E


def punctured_sets(m, E, a):
    r"""
    Return the dictionary of the punctured sets of the strip: ``half`` (the
    first half of the faces), ``all_but_col0`` (the faces of the columns other
    than the first) and ``odd_cols`` (the faces of the odd columns), each as
    a sorted list of faces.
    """
    h2f = m.half_edge_to_face()
    nf = m.num_faces()
    col = {}
    for i in range(COLUMNS):
        for j in range(a):
            col[h2f[E(i, j, 0)]] = i
    return {
        "half": list(range(nf // 2)),
        "all_but_col0": sorted(f for f in range(nf) if col[f] != 0),
        "odd_cols": sorted(f for f in range(nf) if col[f] % 2 == 1),
    }


def depth_first_forest(m, cycles, h2c, roots, used):
    r"""
    Return a depth first spanning forest of the cells ``cycles`` rooted at
    ``roots``, in the format of ``OrientedMap._spanning_forest``: ``-1`` at a
    root, the half-edge to the parent at the other cells, ``-2`` at a cell
    not reached. The edges flagged in ``used`` are skipped, and the edges taken
    are flagged.
    """
    vp = m._vp
    forest = array('i', [-2] * len(cycles))
    todo = []
    for c in roots:
        forest[c] = -1
        todo.append(c)
    while todo:
        c = todo.pop()
        for h in cycles[c]:
            if used[h // 2]:
                continue
            if vp[h ^ 1] == -1:
                continue
            h = h ^ 1
            cc = h2c[h]
            if forest[cc] == -2:
                forest[cc] = h
                used[h // 2] = 1
                todo.append(cc)
    return forest


def depth_first_pairing(m, punctured, algorithm):
    r"""
    Return the pairing of ``m`` with the given punctured faces and algorithm.
    For ``"reduced_map"`` its reduced map is the one of depth first forest and
    coforest.
    """
    g = GeometricIntersectionPairing(m, punctured, algorithm=algorithm)
    if algorithm == "reduced_map":
        # the class builds its reduced map from breadth first forests and
        # takes no forest as argument: replace the reduction by the one of
        # depth first forests, setting the attributes as _build_reduced_map
        # and __init__ do
        cm = g._cm
        used = array('b', [0] * cm.num_edges())
        f = depth_first_forest(cm, cm.vertices(), cm.half_edge_to_vertex(), [0], used)
        c = depth_first_forest(cm, cm.faces(), cm.half_edge_to_face(), g._punctured_faces, used)
        r, mor = cm.reduced_map(forest=f, coforest=c, relabel=True, mapping=True, mutable=False)
        g._r = r
        g._mor = _morphism_arrays(mor)
        g._T = None
        g._angles = perm_dense_cycle_positions(g._r._vp)
    return g


def random_closed_walk(m, length, rng):
    r"""
    Return a closed walk of ``m``: a random walk of ``length`` half-edges
    that does not backtrack unless forced, closed by a shortest path.
    """
    h2v = m.half_edge_to_vertex()
    out = m.vertices()
    start = rng.randrange(m.num_vertices())
    v = start
    walk = []
    prev = None
    for _ in range(length):
        choices = [h for h in out[v] if prev is None or h != prev ^ 1]
        if not choices:
            choices = list(out[v])
        h = rng.choice(choices)
        walk.append(h)
        prev = h
        v = h2v[h ^ 1]
    # shortest path back to start, by breadth first search from start
    # backwards: parent[w] is the half-edge leaving w on a shortest path
    parent = {start: None}
    queue = deque([start])
    while queue:
        w = queue.popleft()
        for h in out[w]:
            x = h2v[h ^ 1]
            if x not in parent:
                parent[x] = h ^ 1
                queue.append(x)
    while v != start:
        h = parent[v]
        walk.append(h)
        v = h2v[h ^ 1]
    return walk


def best_time(g, u, v, repeat):
    r"""
    Return the intersection of ``u`` and ``v`` in ``g`` and the best time of
    ``repeat`` runs.
    """
    best = None
    for _ in range(repeat):
        t0 = time.perf_counter()
        x = g([u], [v])
        dt = time.perf_counter() - t0
        if best is None or dt < best:
            best = dt
    return x, best


def measure_group(m, a, name, punctured, repeat, npairs):
    r"""
    Return the rm L_inf and the ts/rm ratio of the group ``(a, name)``.
    """
    rm = depth_first_pairing(m, punctured, "reduced_map")
    ts = depth_first_pairing(m, punctured, "triangle_system")
    logs = []
    for length in LENGTHS:
        for k in range(npairs):
            rng = random.Random(length * 1000 + k)
            u = random_closed_walk(m, length, rng)
            v = random_closed_walk(m, length, rng)
            x_rm, t_rm = best_time(rm, u, v, repeat)
            x_ts, t_ts = best_time(ts, u, v, repeat)
            if x_rm != x_ts:
                raise RuntimeError(f"a={a} {name}: reduced map gives {x_rm} and triangle system {x_ts} "
                                   f"(length {length}, pair {k})")
            logs.append(math.log(t_ts / t_rm))
    return rm._max_image_length(), math.exp(sum(logs) / len(logs))


def threshold(groups):
    r"""
    Return ``C`` and the exceptions for the list of tuples
    ``(name, a, linf, ratio)``.
    """
    C = None
    for x in sorted({g[2] for g in groups}):
        if all(g[3] >= 1 for g in groups if g[2] <= x):
            C = x
        else:
            break
    if C is None:
        raise RuntimeError("the triangle system is faster on the group of smallest L_inf")
    exceptions = [g for g in groups if (g[2] <= C and g[3] < 1) or (g[2] > C and g[3] > 1)]
    return C, exceptions


def main(argv):
    repeat = 5
    npairs = 3
    args = list(argv)
    while args:
        opt = args.pop(0)
        if opt == "--repeat":
            repeat = int(args.pop(0))
        elif opt == "--pairs":
            npairs = int(args.pop(0))
        else:
            raise SystemExit(f"unknown argument {opt}")

    groups = []
    for a in ROWS:
        m, E = torus_grid_rect(a, COLUMNS)
        sets = punctured_sets(m, E, a)
        for name in SETS:
            linf, ratio = measure_group(m, a, name, sets[name], repeat, npairs)
            groups.append((name, a, linf, ratio))
    groups.sort(key=lambda g: (SETS.index(g[0]), g[1]))

    print(f"{'punctured':<14}{'a':>4}{'rm L_inf':>10}{'ts/rm':>8}")
    for name, a, linf, ratio in groups:
        print(f"{name:<14}{a:>4}{linf:>10}{ratio:>8.2f}")
    C, exceptions = threshold(groups)
    print(f"C = {C}")
    if exceptions:
        print("exceptions: " + ", ".join(f"{name} a={a} (L_inf {linf}, ts/rm {ratio:.2f})"
                                         for name, a, linf, ratio in exceptions))
    else:
        print("exceptions: none")


if __name__ == "__main__":
    main(sys.argv[1:])
