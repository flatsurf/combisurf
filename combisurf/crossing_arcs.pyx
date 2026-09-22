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
:func:`crossing_arcs_sweep_sorted` computes the same number from two sorted
arrays of arcs, one per curve, which is how
:class:`~combisurf.geometric_intersection.GeometricIntersectionMatrix` keeps
them.

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

from cpython cimport array
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


cdef inline void _fenwick_clear(long long *tree, Py_ssize_t size, Py_ssize_t i) noexcept:
    # zero the cells of the Fenwick tree that _fenwick_add(tree, size, i, x)
    # touches, so that a tree can be restored to zero in the time it took to
    # fill it rather than in O(size)
    i += 1
    while i <= size:
        tree[i] = 0
        i += i & (-i)


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


cdef long long *_as_longlongs(array.array a, str name) except? NULL:
    # the data pointer of an empty array may be NULL
    if a.ob_descr.typecode != b'q':
        raise TypeError(f"{name} must be an array of typecode 'q'")
    return a.data.as_longlongs


cdef int _check_arcs(int n, long long *keys, Py_ssize_t num, str name) except -1:
    cdef Py_ssize_t k
    cdef long long first, last
    for k in range(num):
        if keys[k] < 0:
            raise ValueError(f"invalid arc in {name}")
        last = keys[k] // n
        first = keys[k] - last * n
        if not first < last < n:
            raise ValueError(f"invalid arc in {name}")
        if k and keys[k] <= keys[k - 1]:
            raise ValueError(f"{name} must be strictly increasing")
    return 0


def crossing_arcs_sweep_sorted(int n, array.array ukeys not None, array.array uweights not None,
                               array.array vkeys=None, array.array vweights=None,
                               array.array scratch=None, bint check=True):
    r"""
    Return the weighted number of pairs of crossing arcs, the arcs being given
    as two sorted arrays, one for the `u`-weights and one for the `v`-weights.

    This is :func:`crossing_arcs_sweep` on the arcs ``arcs`` with
    ``arcs[ukeys[k]][0] = uweights[k]``, ``arcs[vkeys[k]][1] = vweights[k]``
    and the weights missing from these two arrays equal to zero.

    INPUT:

    - ``n`` -- positive integer, the number of positions on the circle

    - ``ukeys``, ``uweights`` -- two arrays of typecode ``'q'`` and of the
      same length; ``ukeys`` holds strictly increasing keys ``last * n +
      first`` (with ``0 <= first < last < n``) and ``uweights`` the
      `u`-weights of these arcs, non-negative

    - ``vkeys``, ``vweights`` -- (default: ``None``) the same for the
      `v`-weights; when they are ``None``, the `v`-weights are taken equal to
      the `u`-weights, which is the ``symmetric`` case of
      :func:`crossing_arcs_sweep`

    - ``scratch`` -- (default: ``None``) an array of typecode ``'q'``, of
      length at least ``2 * (n + 1)`` and filled with zeros, used as the
      memory of the sweep and filled with zeros again on return; it saves an
      allocation for each call when many of them are made with the same ``n``

    - ``check`` -- boolean (default: ``True``); whether to check that the
      keys are sorted and valid

    OUTPUT: the integer `S` of the module documentation

    ALGORITHM:

    The sweep of :func:`crossing_arcs_sweep`, with the arcs visited by
    merging the two arrays group by group, a group being the arcs with a
    given right endpoint. The `u`-arcs of a group are queried against the
    `v`-weights inserted so far, the `v`-arcs against the `u`-weights, and
    both are inserted afterwards. So there is no dictionary and no sort, and
    the cost is `O((|ukeys| + |vkeys|) \log(n))`. The cells of the Fenwick
    trees touched by the insertions are set back to zero at the end, in the
    same time.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.crossing_arcs import crossing_arcs_sweep, crossing_arcs_sweep_sorted
        sage: n = 4
        sage: crossing_arcs_sweep_sorted(n, array('q', [2 * n + 0]), array('q', [1]),
        ....:                               array('q', [3 * n + 1]), array('q', [1]))
        1
        sage: crossing_arcs_sweep_sorted(n, array('q', [2 * n + 0, 3 * n + 1]), array('q', [1, 1]))
        2

    It agrees with :func:`crossing_arcs_sweep`, and leaves ``scratch`` as it
    found it::

        sage: n = 9
        sage: ukeys = sorted(set(b * n + a for a, b in (sorted(sample(range(n), 2)) for _ in range(15))))
        sage: vkeys = sorted(set(b * n + a for a, b in (sorted(sample(range(n), 2)) for _ in range(15))))
        sage: uweights = [randint(1, 3) for _ in ukeys]
        sage: vweights = [randint(1, 3) for _ in vkeys]
        sage: arcs = {k: [0, 0] for k in ukeys + vkeys}
        sage: for k, u in zip(ukeys, uweights):
        ....:     arcs[k][0] = u
        sage: for k, v in zip(vkeys, vweights):
        ....:     arcs[k][1] = v
        sage: scratch = array('q', [0] * (2 * (n + 1)))
        sage: S = crossing_arcs_sweep_sorted(n, array('q', ukeys), array('q', uweights),
        ....:                                   array('q', vkeys), array('q', vweights), scratch)
        sage: S == crossing_arcs_sweep(n, arcs)
        True
        sage: all(x == 0 for x in scratch)
        True
        sage: S = crossing_arcs_sweep_sorted(n, array('q', ukeys), array('q', uweights), scratch=scratch)
        sage: S == crossing_arcs_sweep(n, {k: [u, u] for k, u in zip(ukeys, uweights)}, True)
        True
        sage: all(x == 0 for x in scratch)
        True

    TESTS::

        sage: e = array('q')
        sage: crossing_arcs_sweep_sorted(1, e, e), crossing_arcs_sweep_sorted(1, e, e, e, e)
        (0, 0)
        sage: crossing_arcs_sweep_sorted(0, e, e)
        Traceback (most recent call last):
        ...
        ValueError: n must be positive
        sage: crossing_arcs_sweep_sorted(4, array('q', [8]), array('q', [1, 1]))
        Traceback (most recent call last):
        ...
        ValueError: ukeys and uweights must have the same length
        sage: crossing_arcs_sweep_sorted(4, array('q', [8]), array('q', [1]), array('q', [8]), None)
        Traceback (most recent call last):
        ...
        ValueError: vkeys and vweights must be both given or both None
        sage: crossing_arcs_sweep_sorted(4, array('q', [13, 8]), array('q', [1, 1]))
        Traceback (most recent call last):
        ...
        ValueError: ukeys must be strictly increasing
        sage: crossing_arcs_sweep_sorted(4, array('q', [5]), array('q', [1]))
        Traceback (most recent call last):
        ...
        ValueError: invalid arc in ukeys
        sage: crossing_arcs_sweep_sorted(4, array('i', [8]), array('q', [1]))
        Traceback (most recent call last):
        ...
        TypeError: ukeys must be an array of typecode 'q'
        sage: crossing_arcs_sweep_sorted(4, e, e, scratch=array('q', [0] * 9))
        Traceback (most recent call last):
        ...
        ValueError: scratch must have length at least 2 * (n + 1)
    """
    if n <= 0:
        raise ValueError("n must be positive")
    cdef bint symmetric = vkeys is None
    if symmetric != (vweights is None):
        raise ValueError("vkeys and vweights must be both given or both None")

    cdef long long *ku = _as_longlongs(ukeys, "ukeys")
    cdef long long *wu = _as_longlongs(uweights, "uweights")
    cdef Py_ssize_t lu = len(ukeys)
    if len(uweights) != lu:
        raise ValueError("ukeys and uweights must have the same length")
    cdef long long *kv = NULL
    cdef long long *wv = NULL
    cdef Py_ssize_t lv = 0
    if not symmetric:
        kv = _as_longlongs(vkeys, "vkeys")
        wv = _as_longlongs(vweights, "vweights")
        lv = len(vkeys)
        if len(vweights) != lv:
            raise ValueError("vkeys and vweights must have the same length")
    if check:
        _check_arcs(n, ku, lu, "ukeys")
        if not symmetric:
            _check_arcs(n, kv, lv, "vkeys")

    cdef long long *fu
    cdef long long *fv
    cdef bint owned = scratch is None
    if owned:
        fu = <long long *> calloc(2 * (n + 1), sizeof(long long))
        if fu == NULL:
            raise MemoryError
    else:
        if len(scratch) < 2 * (n + 1):
            raise ValueError("scratch must have length at least 2 * (n + 1)")
        fu = _as_longlongs(scratch, "scratch")
    # the positions go from 0 to n - 1, stored at 1..n of each tree
    fv = fu + (n + 1)

    cdef long long S = 0
    cdef Py_ssize_t iu, iv, ju, jv, t
    cdef long long last, first, key

    if symmetric:
        iu = 0
        while iu < lu:
            last = ku[iu] // n
            ju = iu
            while ju < lu and ku[ju] // n == last:
                first = ku[ju] - last * n
                S += 2 * wu[ju] * _fenwick_prefix(fu, first + 1)
                ju += 1
            for t in range(iu, ju):
                first = ku[t] - last * n
                _fenwick_add(fu, n, first + 1, wu[t])
                _fenwick_add(fu, n, last, -wu[t])
            iu = ju
        for t in range(lu):
            last = ku[t] // n
            _fenwick_clear(fu, n, ku[t] - last * n + 1)
            _fenwick_clear(fu, n, last)
    else:
        iu = iv = 0
        while iu < lu or iv < lv:
            # the next group is the smallest right endpoint of the two heads
            if iv == lv or (iu < lu and ku[iu] < kv[iv]):
                last = ku[iu] // n
            else:
                last = kv[iv] // n
            ju = iu
            while ju < lu and ku[ju] // n == last:
                S += wu[ju] * _fenwick_prefix(fv, ku[ju] - last * n + 1)
                ju += 1
            jv = iv
            while jv < lv and kv[jv] // n == last:
                S += wv[jv] * _fenwick_prefix(fu, kv[jv] - last * n + 1)
                jv += 1
            for t in range(iu, ju):
                _fenwick_add(fu, n, ku[t] - last * n + 1, wu[t])
                _fenwick_add(fu, n, last, -wu[t])
            for t in range(iv, jv):
                _fenwick_add(fv, n, kv[t] - last * n + 1, wv[t])
                _fenwick_add(fv, n, last, -wv[t])
            iu = ju
            iv = jv
        for t in range(lu):
            last = ku[t] // n
            _fenwick_clear(fu, n, ku[t] - last * n + 1)
            _fenwick_clear(fu, n, last)
        for t in range(lv):
            last = kv[t] // n
            _fenwick_clear(fv, n, kv[t] - last * n + 1)
            _fenwick_clear(fv, n, last)

    if owned:
        free(fu)
    return S
