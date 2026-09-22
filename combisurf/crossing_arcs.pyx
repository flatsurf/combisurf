r"""
Counting crossing chords of a circle

The chords here are arcs between the positions ``0``, ..., ``n - 1`` of a
circle. Two of them, `A = (i_0, j_0)` and `B = (i_1, j_1)` with `i_0 < j_0`
and `i_1 < j_1`, cross when their endpoints strictly interleave, that is
`i_0 < i_1 < j_0 < j_1` or `i_1 < i_0 < j_1 < j_0`; arcs that share an
endpoint do not cross. Each arc carries a `u`-weight and a `v`-weight, and
:func:`crossing_arcs_sweep` returns

.. MATH::

    S = \sum_{A, B} u(A) v(B) + v(A) u(B)

over the pairs `(A, B)` of crossing arcs with `i_0 < i_1`.

This is the `O(n^2)` term of
:meth:`~combisurf.geometric_intersection.GeometricIntersection.geometric_intersection`,
where the positions are the angles at the vertex and the arcs the
consecutive pairs of letters of the curves.

EXAMPLES::

    sage: from combisurf.crossing_arcs import crossing_arcs_sweep
    sage: n = 4
    sage: crossing_arcs_sweep(n, {2 * n + 0: [1, 0], 3 * n + 1: [0, 1]})
    1

.. SEEALSO::

    :mod:`combisurf.crossing_arcs_naive` holds the same sweep in pure Python,
    and :func:`~combisurf.crossing_arcs_naive.crossing_arcs_double_sum` the
    `O(n^2)` double sum. They are the reference this is tested against.
"""
# ****************************************************************************
#  This file is part of combisurf
#
#       Copyright (C) 2026 Vincent Delecroix
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# ****************************************************************************

from libc.stdlib cimport calloc, free


cdef inline void _fenwick_add(long long *tree, Py_ssize_t size, Py_ssize_t i, long long x) noexcept:
    # add x at position i (0-based) of the Fenwick tree tree[1..size]
    i += 1
    while i <= size:
        tree[i] += x
        i += i & (-i)


cdef inline long long _fenwick_prefix(long long *tree, Py_ssize_t i) noexcept:
    # sum of the positions 0, ..., i - 1 of the Fenwick tree
    cdef long long s = 0
    while i > 0:
        s += tree[i]
        i -= i & (-i)
    return s


def crossing_arcs_sweep(int n, dict arcs, bint symmetric=False):
    r"""
    Return the weighted number of pairs of crossing arcs in ``arcs``.

    INPUT:

    - ``n`` -- positive integer, the number of positions on the circle

    - ``arcs`` -- dictionary; the arc from ``first`` to ``last`` (with
      ``0 <= first < last < n``) is the key ``last * n + first`` and its value
      is the pair ``[u, v]`` of its weights, both non-negative

    - ``symmetric`` -- boolean (default: ``False``); whether the `u`-weight
      and the `v`-weight are equal on every arc, in which case the
      `v`-weights are not read

    OUTPUT: the integer `S` of the module documentation

    ALGORITHM:

    The arcs `B = (i_1, j_1)` are visited by increasing right endpoint `j_1`.
    When `B` is visited, the arcs already inserted are exactly the arcs
    `A = (i_0, j_0)` with `j_0 < j_1`, and among those `A` crosses `B` from
    the left if and only if `i_0 < i_1` and not `j_0 \le i_1` (the second
    condition implies the first since `i_0 < j_0`). So `A` is inserted as
    `+u(A)` at position `i_0 + 1` and `-u(A)` at position `j_0` of a Fenwick
    tree (:ref:`fenwick1994`), and its prefix sum up to `i_1` is the total
    `u`-weight of the arcs crossing `B` from the left. The same goes for `v`.
    The arcs sharing the right endpoint `j_1` do not cross each other, so
    they are all visited before any of them is inserted.

    This is the easy case of counting segment intersections by a sweep
    (:ref:`chazelle1986`): on a circle the cyclic order of the endpoints
    already is the sweep order. It costs `O(|\text{arcs}| \log(n))` after
    sorting the keys.

    The weights and the result are held in C ``long long``. The `u`-weights
    add up to the total length of the curves they come from, and the same
    for `v`, so the result is at most twice the product of these two lengths
    and does not overflow for any input that fits in memory.

    EXAMPLES::

        sage: from combisurf.crossing_arcs import crossing_arcs_sweep
        sage: n = 4
        sage: crossing_arcs_sweep(n, {2 * n + 0: [1, 0], 3 * n + 1: [0, 1]})
        1
        sage: crossing_arcs_sweep(n, {2 * n + 0: [1, 1], 3 * n + 1: [1, 1]}, True)
        2

    Arcs sharing an endpoint do not cross::

        sage: crossing_arcs_sweep(n, {2 * n + 0: [1, 1], 3 * n + 2: [1, 1]})
        0
        sage: crossing_arcs_sweep(n, {2 * n + 0: [1, 1], 2 * n + 1: [1, 1]})
        0

    It agrees with the double sum::

        sage: from combisurf.crossing_arcs_naive import crossing_arcs_double_sum
        sage: n = 9
        sage: arcs = {}
        sage: for _ in range(20):
        ....:     first, last = sorted(sample(range(n), 2))
        ....:     arcs[last * n + first] = [randint(0, 3), randint(0, 3)]
        sage: crossing_arcs_sweep(n, arcs) == crossing_arcs_double_sum(n, arcs)
        True
        sage: for weights in arcs.values():
        ....:     weights[1] = weights[0]
        sage: crossing_arcs_sweep(n, arcs, True) == crossing_arcs_double_sum(n, arcs)
        True

    TESTS::

        sage: crossing_arcs_sweep(1, {})
        0
        sage: crossing_arcs_sweep(0, {})
        Traceback (most recent call last):
        ...
        ValueError: n must be positive
    """
    if n <= 0:
        raise ValueError("n must be positive")

    cdef list keys = sorted(arcs)
    cdef Py_ssize_t num = len(keys)
    cdef Py_ssize_t k, k0, t
    cdef long long key, last, first, u, v
    cdef long long S = 0

    # the positions go from 0 to n - 1, stored at 1..n
    cdef long long *fu = <long long *> calloc(n + 1, sizeof(long long))
    cdef long long *fv = <long long *> calloc(n + 1, sizeof(long long))
    if fu == NULL or fv == NULL:
        free(fu)
        free(fv)
        raise MemoryError

    try:
        k = 0
        while k < num:
            last = (<long long> keys[k]) // n
            k0 = k
            while k < num:
                key = keys[k]
                if key // n != last:
                    break
                first = key - last * n
                weights = arcs[keys[k]]
                u = weights[0]
                if symmetric:
                    S += 2 * u * _fenwick_prefix(fu, first + 1)
                else:
                    v = weights[1]
                    S += v * _fenwick_prefix(fu, first + 1) + u * _fenwick_prefix(fv, first + 1)
                k += 1
            for t in range(k0, k):
                key = keys[t]
                first = key - last * n
                weights = arcs[keys[t]]
                u = weights[0]
                _fenwick_add(fu, n, first + 1, u)
                _fenwick_add(fu, n, last, -u)
                if not symmetric:
                    v = weights[1]
                    _fenwick_add(fv, n, first + 1, v)
                    _fenwick_add(fv, n, last, -v)
    finally:
        free(fu)
        free(fv)

    return S
