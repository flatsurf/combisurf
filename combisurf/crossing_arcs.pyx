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

The pairs of arcs sharing a startpoint make the other term of an
intersection number. They are counted by :func:`startpoint_sweep_sorted`,
for the leaves of two curves kept apart as in
:class:`~combisurf.geometric_intersection.GeometricIntersectionMatrix`, and
by :func:`startpoint_sweep_weighted`, for the leaves of weighted multicurves
listed together as in
:meth:`~combisurf.geometric_intersection.GeometricIntersection.geometric_intersection`.

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

from combisurf.partial_sums cimport fenwick_add, fenwick_prefix, fenwick_clear


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
                    S += 2 * u * fenwick_prefix(fu, first + 1)
                else:
                    v = weights[1]
                    S += v * fenwick_prefix(fu, first + 1) + u * fenwick_prefix(fv, first + 1)
                k += 1
            for t in range(k0, k):
                key = keys[t]
                first = key - last * n
                weights = arcs[keys[t]]
                u = weights[0]
                fenwick_add(fu, n, first + 1, u)
                fenwick_add(fu, n, last, -u)
                if not symmetric:
                    v = weights[1]
                    fenwick_add(fv, n, first + 1, v)
                    fenwick_add(fv, n, last, -v)
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
                S += 2 * wu[ju] * fenwick_prefix(fu, first + 1)
                ju += 1
            for t in range(iu, ju):
                first = ku[t] - last * n
                fenwick_add(fu, n, first + 1, wu[t])
                fenwick_add(fu, n, last, -wu[t])
            iu = ju
        for t in range(lu):
            last = ku[t] // n
            fenwick_clear(fu, n, ku[t] - last * n + 1)
            fenwick_clear(fu, n, last)
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
                S += wu[ju] * fenwick_prefix(fv, ku[ju] - last * n + 1)
                ju += 1
            jv = iv
            while jv < lv and kv[jv] // n == last:
                S += wv[jv] * fenwick_prefix(fu, kv[jv] - last * n + 1)
                jv += 1
            for t in range(iu, ju):
                fenwick_add(fu, n, ku[t] - last * n + 1, wu[t])
                fenwick_add(fu, n, last, -wu[t])
            for t in range(iv, jv):
                fenwick_add(fv, n, kv[t] - last * n + 1, wv[t])
                fenwick_add(fv, n, last, -wv[t])
            iu = ju
            iv = jv
        for t in range(lu):
            last = ku[t] // n
            fenwick_clear(fu, n, ku[t] - last * n + 1)
            fenwick_clear(fu, n, last)
        for t in range(lv):
            last = kv[t] // n
            fenwick_clear(fv, n, kv[t] - last * n + 1)
            fenwick_clear(fv, n, last)

    if owned:
        free(fu)
    return S


cdef int _check_leaves(int n, long long *angles, Py_ssize_t num, str name) except -1:
    cdef Py_ssize_t k
    for k in range(num):
        if not 0 <= angles[k] <= n - 2:
            raise ValueError(f"invalid angle in {name}")
    return 0


def startpoint_sweep_sorted(int n, array.array uranks not None, array.array ustarts not None,
                            array.array uangles not None,
                            array.array vranks=None, array.array vstarts=None, array.array vangles=None,
                            array.array scratch=None, bint check=True):
    r"""
    Return the number of pairs of arcs with the same startpoint that cross,
    the arcs being the leaves of two curves listed in increasing rank.

    The leaves of a curve are the conjugates of the curve and of its inverse.
    Each leaf is an arc starting at the angle of its first letter and ending
    at the angle of the inverse of its last letter; the leaves of all the
    curves are ranked by their cyclic order at infinity, so that the leaves
    with a given startpoint have consecutive ranks. Two leaves `A` and `B`
    with the same startpoint and ``rank(A) < rank(B)`` cross when the angle of
    ``A`` is smaller than the one of ``B``. This returns the number of such
    pairs with `A` a leaf of `u` and `B` a leaf of `v` plus the number of
    such pairs with `A` a leaf of `v` and `B` a leaf of `u`.

    This is the term of
    :class:`~combisurf.geometric_intersection.GeometricIntersectionMatrix`
    counting the intersections at the vertex.

    INPUT:

    - ``n`` -- positive integer, the number of half-edges

    - ``uranks``, ``ustarts``, ``uangles`` -- three arrays of typecode
      ``'q'`` of the same length, one entry per leaf of `u` in increasing
      order of rank: its rank, its startpoint and the angle from its
      startpoint to its endpoint minus one, an integer in ``0 .. n - 2``

    - ``vranks``, ``vstarts``, ``vangles`` -- (default: ``None``) the same
      for `v`; when they are ``None``, `v` is taken equal to `u`

    - ``scratch`` -- (default: ``None``) an array of typecode ``'q'``, of
      length at least ``2 * (n + 1)`` and filled with zeros, used as the
      memory of the sweep and filled with zeros again on return; it saves an
      allocation for each call when many of them are made with the same ``n``

    - ``check`` -- boolean (default: ``True``); whether to check that the
      lengths match, that the angles lie in ``0 .. n - 2`` and that the ranks
      are strictly increasing

    OUTPUT: an integer

    ALGORITHM:

    The two lists of leaves are merged by rank and the result is swept group
    by group, a group being the leaves with a given startpoint. Inside a group
    each leaf of `u` is queried against the leaves of `v` of the group
    inserted so far, and conversely, with two Fenwick trees indexed by the
    angles. The cells touched by a group are set back to zero at its end, so
    the cost is `O((|u| + |v|) \log(n))`. When `v` is `u` a single tree is
    enough and each pair is counted twice.

    EXAMPLES:

    Two leaves starting at the same half-edge, the one of smaller rank having
    the smaller angle::

        sage: from array import array
        sage: from combisurf.crossing_arcs import startpoint_sweep_sorted
        sage: n = 8
        sage: q = lambda *x: array('q', x)
        sage: startpoint_sweep_sorted(n, q(3), q(0), q(1), q(4), q(0), q(5))
        1
        sage: startpoint_sweep_sorted(n, q(3), q(0), q(5), q(4), q(0), q(1))
        0
        sage: startpoint_sweep_sorted(n, q(3), q(0), q(1), q(4), q(2), q(5))
        0
        sage: startpoint_sweep_sorted(n, q(3, 4), q(0, 0), q(1, 5))
        2

    It agrees with the pure Python version, and leaves ``scratch`` as it
    found it::

        sage: from combisurf import crossing_arcs_naive
        sage: n = 9
        sage: ranks = sorted(sample(range(60), 40))
        sage: starts = sorted(randrange(n) for _ in ranks)
        sage: angles = [randrange(n - 1) for _ in ranks]
        sage: side = [randrange(2) for _ in ranks]
        sage: u = [[x[k] for k in range(40) if side[k] == 0] for x in (ranks, starts, angles)]
        sage: v = [[x[k] for k in range(40) if side[k] == 1] for x in (ranks, starts, angles)]
        sage: scratch = array('q', [0] * (2 * (n + 1)))
        sage: S = startpoint_sweep_sorted(n, *[array('q', x) for x in u + v], scratch=scratch)
        sage: S == crossing_arcs_naive.startpoint_sweep_sorted(n, *u, *v)
        True
        sage: S = startpoint_sweep_sorted(n, *[array('q', x) for x in u], scratch=scratch)
        sage: S == crossing_arcs_naive.startpoint_sweep_sorted(n, *u)
        True
        sage: all(x == 0 for x in scratch)
        True

    TESTS::

        sage: e = array('q')
        sage: startpoint_sweep_sorted(1, e, e, e), startpoint_sweep_sorted(1, e, e, e, e, e, e)
        (0, 0)
        sage: startpoint_sweep_sorted(0, e, e, e)
        Traceback (most recent call last):
        ...
        ValueError: n must be positive
        sage: startpoint_sweep_sorted(4, q(0), q(0), q(0, 1))
        Traceback (most recent call last):
        ...
        ValueError: uranks, ustarts and uangles must have the same length
        sage: startpoint_sweep_sorted(4, q(0), q(0), q(0), q(1), None, q(0))
        Traceback (most recent call last):
        ...
        ValueError: vranks, vstarts and vangles must be all given or all None
        sage: startpoint_sweep_sorted(4, q(0), q(0), q(3))
        Traceback (most recent call last):
        ...
        ValueError: invalid angle in uangles
        sage: startpoint_sweep_sorted(4, q(1, 0), q(0, 0), q(0, 1))
        Traceback (most recent call last):
        ...
        ValueError: uranks must be strictly increasing
        sage: startpoint_sweep_sorted(4, array('i', [0]), q(0), q(0))
        Traceback (most recent call last):
        ...
        TypeError: uranks must be an array of typecode 'q'
        sage: startpoint_sweep_sorted(4, e, e, e, scratch=array('q', [0] * 9))
        Traceback (most recent call last):
        ...
        ValueError: scratch must have length at least 2 * (n + 1)
    """
    if n <= 0:
        raise ValueError("n must be positive")
    cdef bint symmetric = vranks is None
    if not (symmetric == (vstarts is None) == (vangles is None)):
        raise ValueError("vranks, vstarts and vangles must be all given or all None")

    cdef long long *ru = _as_longlongs(uranks, "uranks")
    cdef long long *su = _as_longlongs(ustarts, "ustarts")
    cdef long long *au = _as_longlongs(uangles, "uangles")
    cdef Py_ssize_t lu = len(uranks)
    if len(ustarts) != lu or len(uangles) != lu:
        raise ValueError("uranks, ustarts and uangles must have the same length")
    cdef long long *rv = NULL
    cdef long long *sv = NULL
    cdef long long *av = NULL
    cdef Py_ssize_t lv = 0
    if not symmetric:
        rv = _as_longlongs(vranks, "vranks")
        sv = _as_longlongs(vstarts, "vstarts")
        av = _as_longlongs(vangles, "vangles")
        lv = len(vranks)
        if len(vstarts) != lv or len(vangles) != lv:
            raise ValueError("vranks, vstarts and vangles must have the same length")
    cdef Py_ssize_t t
    if check:
        _check_leaves(n, au, lu, "uangles")
        for t in range(1, lu):
            if ru[t] <= ru[t - 1]:
                raise ValueError("uranks must be strictly increasing")
        if not symmetric:
            _check_leaves(n, av, lv, "vangles")
            for t in range(1, lv):
                if rv[t] <= rv[t - 1]:
                    raise ValueError("vranks must be strictly increasing")

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
    # the angles go from 0 to n - 2, stored at 1..n - 1 of each tree
    cdef Py_ssize_t size = n - 1
    fv = fu + (n + 1)

    cdef long long S = 0
    cdef long long startpoint
    cdef Py_ssize_t iu, iv, ju, jv, iu0, iv0

    if symmetric:
        iu = 0
        while iu < lu:
            startpoint = su[iu]
            iu0 = iu
            iu += 1
            while iu < lu and su[iu] == startpoint:
                iu += 1
            if iu - iu0 == 1:
                # a single arc through this startpoint crosses nothing
                continue
            fenwick_add(fu, size, au[iu0], 1)
            for t in range(iu0 + 1, iu):
                S += 2 * fenwick_prefix(fu, au[t])
                fenwick_add(fu, size, au[t], 1)
            for t in range(iu0, iu):
                fenwick_clear(fu, size, au[t])
    else:
        iu = iv = 0
        while iu < lu or iv < lv:
            # the leaves sharing a startpoint are consecutive in the cyclic
            # order, so the head of smaller rank opens the group
            if iv == lv or (iu < lu and ru[iu] < rv[iv]):
                startpoint = su[iu]
            else:
                startpoint = sv[iv]
            iu0 = iu
            while iu < lu and su[iu] == startpoint:
                iu += 1
            iv0 = iv
            while iv < lv and sv[iv] == startpoint:
                iv += 1
            if iu0 == iu or iv0 == iv:
                # only one of the two curves goes through this startpoint
                continue
            ju = iu0
            jv = iv0
            while ju < iu or jv < iv:
                if jv == iv or (ju < iu and ru[ju] < rv[jv]):
                    S += fenwick_prefix(fv, au[ju])
                    fenwick_add(fu, size, au[ju], 1)
                    ju += 1
                else:
                    S += fenwick_prefix(fu, av[jv])
                    fenwick_add(fv, size, av[jv], 1)
                    jv += 1
            for t in range(iu0, iu):
                fenwick_clear(fu, size, au[t])
            for t in range(iv0, iv):
                fenwick_clear(fv, size, av[t])

    if owned:
        free(fu)
    return S


def startpoint_sweep_weighted(int n, array.array starts not None, array.array angles not None,
                              array.array uweights not None, array.array vweights=None,
                              array.array scratch=None, bint check=True):
    r"""
    Return the weighted number of pairs of arcs with the same startpoint that
    cross, the arcs being the leaves of a family of curves listed in cyclic
    order at infinity.

    This is :func:`startpoint_sweep_sorted` when the leaves of both sides
    are in a single list and carry weights: the leaf of index `k` has
    `u`-weight ``uweights[k]`` and `v`-weight ``vweights[k]``, and the value
    returned is the sum of `u(A) v(B) + v(A) u(B)` over the pairs of leaves
    `A` before `B` with the same startpoint and the angle of `A` smaller than
    the one of `B`. This is the term of
    :meth:`~combisurf.geometric_intersection.GeometricIntersection.geometric_intersection`
    counting the intersections at the vertex, the weights being the
    multiplicities of the curves in the two multicurves.

    INPUT:

    - ``n`` -- positive integer, the number of half-edges

    - ``starts``, ``angles`` -- two arrays of typecode ``'q'`` of the same
      length, one entry per leaf in cyclic order at infinity: its startpoint
      and the angle from its startpoint to its endpoint minus one, an integer
      in ``0 .. n - 2``; the leaves with the same startpoint must be
      consecutive

    - ``uweights`` -- an array of typecode ``'q'`` of the same length, the
      `u`-weights of the leaves

    - ``vweights`` -- (default: ``None``) the same for the `v`-weights; when
      it is ``None``, the `v`-weights are taken equal to the `u`-weights

    - ``scratch`` -- (default: ``None``) as in :func:`startpoint_sweep_sorted`

    - ``check`` -- boolean (default: ``True``); whether to check that the
      lengths match and that the angles lie in ``0 .. n - 2``

    OUTPUT: an integer, even when ``vweights`` is ``None``

    ALGORITHM:

    The list is swept group by group, a group being the leaves with a given
    startpoint, each leaf being queried against the weights of the other side
    inserted so far in two Fenwick trees indexed by the angles. The cells
    touched by a group are set back to zero at its end, so the cost is
    `O(|starts| \log(n))`. When ``vweights`` is ``None`` a single tree is
    enough.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.crossing_arcs import startpoint_sweep_weighted
        sage: q = lambda *x: array('q', x)
        sage: startpoint_sweep_weighted(8, q(0, 0), q(1, 5), q(2, 0), q(0, 3))
        6
        sage: startpoint_sweep_weighted(8, q(0, 0), q(1, 5), q(2, 3))
        12

    It agrees with the pure Python version, and leaves ``scratch`` as it
    found it::

        sage: from combisurf import crossing_arcs_naive
        sage: n = 9
        sage: starts = sorted(randrange(n) for _ in range(40))
        sage: angles = [randrange(n - 1) for _ in starts]
        sage: uweights = [randrange(3) for _ in starts]
        sage: vweights = [randrange(3) for _ in starts]
        sage: scratch = array('q', [0] * (2 * (n + 1)))
        sage: S = startpoint_sweep_weighted(n, *[array('q', x) for x in (starts, angles, uweights, vweights)], scratch=scratch)
        sage: S == crossing_arcs_naive.startpoint_sweep_weighted(n, starts, angles, uweights, vweights)
        True
        sage: S = startpoint_sweep_weighted(n, *[array('q', x) for x in (starts, angles, uweights)], scratch=scratch)
        sage: S == crossing_arcs_naive.startpoint_sweep_weighted(n, starts, angles, uweights)
        True
        sage: all(x == 0 for x in scratch)
        True

    TESTS::

        sage: e = array('q')
        sage: startpoint_sweep_weighted(1, e, e, e), startpoint_sweep_weighted(1, e, e, e, e)
        (0, 0)
        sage: startpoint_sweep_weighted(0, e, e, e)
        Traceback (most recent call last):
        ...
        ValueError: n must be positive
        sage: startpoint_sweep_weighted(4, q(0), q(0), q(1), q(1, 1))
        Traceback (most recent call last):
        ...
        ValueError: starts, angles and weights must have the same length
        sage: startpoint_sweep_weighted(4, q(0), q(3), q(1))
        Traceback (most recent call last):
        ...
        ValueError: invalid angle in angles
        sage: startpoint_sweep_weighted(4, e, e, e, scratch=array('q', [0] * 9))
        Traceback (most recent call last):
        ...
        ValueError: scratch must have length at least 2 * (n + 1)
    """
    if n <= 0:
        raise ValueError("n must be positive")
    cdef bint symmetric = vweights is None
    cdef long long *s = _as_longlongs(starts, "starts")
    cdef long long *a = _as_longlongs(angles, "angles")
    cdef long long *wu = _as_longlongs(uweights, "uweights")
    cdef long long *wv = NULL
    cdef Py_ssize_t l = len(starts)
    if len(angles) != l or len(uweights) != l:
        raise ValueError("starts, angles and weights must have the same length")
    if not symmetric:
        wv = _as_longlongs(vweights, "vweights")
        if len(vweights) != l:
            raise ValueError("starts, angles and weights must have the same length")
    if check:
        _check_leaves(n, a, l, "angles")

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
    # the angles go from 0 to n - 2, stored at 1..n - 1 of each tree
    cdef Py_ssize_t size = n - 1
    fv = fu + (n + 1)

    cdef long long S = 0
    cdef long long startpoint
    cdef Py_ssize_t i = 0, i0, t
    while i < l:
        startpoint = s[i]
        i0 = i
        i += 1
        while i < l and s[i] == startpoint:
            i += 1
        if i - i0 == 1:
            # a single arc through this startpoint crosses nothing
            continue
        if symmetric:
            for t in range(i0, i):
                S += 2 * wu[t] * fenwick_prefix(fu, a[t])
                fenwick_add(fu, size, a[t], wu[t])
            for t in range(i0, i):
                fenwick_clear(fu, size, a[t])
        else:
            for t in range(i0, i):
                S += wu[t] * fenwick_prefix(fv, a[t]) + wv[t] * fenwick_prefix(fu, a[t])
                fenwick_add(fu, size, a[t], wu[t])
                fenwick_add(fv, size, a[t], wv[t])
            for t in range(i0, i):
                fenwick_clear(fu, size, a[t])
                fenwick_clear(fv, size, a[t])

    if owned:
        free(fu)
    return S
