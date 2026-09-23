# distutils: include_dirs = combisurf/src
r"""
Counting crossing chords of a circle

The chords here are arcs between the positions ``0``, ..., ``n - 1`` of a
circle. Two of them, `A = (i_0, j_0)` and `B = (i_1, j_1)` with `i_0 < j_0`
and `i_1 < j_1`, cross when their endpoints strictly interleave, that is
`i_0 < i_1 < j_0 < j_1` or `i_1 < i_0 < j_1 < j_0`; arcs that share an
endpoint do not cross. Each arc carries a `u`-weight and a `v`-weight, and
:func:`crossing_arcs_sweep_sorted` returns

.. MATH::

    S = \sum_{A, B} u(A) v(B) + v(A) u(B)

over the pairs `(A, B)` of crossing arcs with `i_0 < i_1`.

This is the term, counting the pairs of arcs with four distinct endpoints,
of
:meth:`~combisurf.geometric_intersection.GeometricIntersection.geometric_intersection`,
where the positions are the angles at the vertex and the arcs the
consecutive pairs of letters of the curves, given to
:func:`crossing_arcs_sweep_sorted` as two sorted arrays of arcs, one per
curve, which is also how
:class:`~combisurf.geometric_intersection.GeometricIntersectionMatrix` keeps
them.

The pairs of arcs sharing a startpoint make the other term of an
intersection number. They are counted by :func:`startpoint_sweep_sorted`,
for the leaves of two curves kept apart as in
:class:`~combisurf.geometric_intersection.GeometricIntersectionMatrix`, and
by :func:`startpoint_sweep_weighted`, for the leaves of weighted multicurves
listed together as in
:meth:`~combisurf.geometric_intersection.GeometricIntersection.geometric_intersection`.

Those leaves are the ones of a conjugate tree holding the curves and their
inverses in the free group where ``h ^ 1`` is the inverse of the letter
``h``: :func:`tree_add_with_inverse` adds a curve and its inverse, and
:func:`cyclically_sorted_leaf_arcs` lists the leaves in their cyclic order at
infinity.

The sweeps keep their prefix sums in Fenwick trees (:ref:`fenwick1994`).
Counting crossing chords is the easy case of counting segment intersections
(:ref:`chazelle1986`), the cyclic order of the endpoints giving the sweep
order.

EXAMPLES::

    sage: from array import array
    sage: from combisurf.crossing_arcs import crossing_arcs_sweep_sorted
    sage: n = 4
    sage: crossing_arcs_sweep_sorted(n, array('q', [2 * n + 0]), array('q', [1]),
    ....:                               array('q', [3 * n + 1]), array('q', [1]))
    1

.. SEEALSO::

    The brute forces this is tested against, straight from the definitions
    above, live in ``test/test_geometric_intersection.py``.
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
from libc.limits cimport INT_MAX
from libc.stdlib cimport calloc, free, malloc, qsort

from combisurf.conjugate_tree cimport ConjugateTree, CT_EALPHABET, CT_ENEGATIVE, CT_ETOOLARGE
from combisurf.partial_sums cimport fenwick_add, fenwick_prefix, fenwick_clear


cdef struct _weighted_arc:
    long long key
    long long weight


cdef int _cmp_weighted_arcs(const void *x, const void *y) noexcept nogil:
    cdef long long a = (<const _weighted_arc *> x).key
    cdef long long b = (<const _weighted_arc *> y).key
    return (a > b) - (a < b)


def word_arcs(int n, angles, list words, weights):
    r"""
    Return the arcs of the words ``words[0]``, ``words[2]``, ``words[4]``,
    ... as sorted keys and weights, in the format of
    :func:`crossing_arcs_sweep_sorted`.

    Each pair of cyclically consecutive letters ``w[p - 1], w[p]`` of such a
    word ``w`` is the arc between the positions ``angles[w[p]]`` and
    ``angles[w[p - 1] ^ 1]``: the curve comes in through the half-edge
    ``w[p - 1] ^ 1`` and goes out through ``w[p]``. The arc from ``first`` to
    ``last`` with ``first < last`` has the key ``last * n + first`` and the
    word ``words[2 * j]`` gives each of its arcs the weight ``weights[j]``.

    INPUT:

    - ``n`` -- positive integer, the number of positions on the circle

    - ``angles`` -- a sequence of length ``n`` of integers in ``0 .. n - 1``,
      the position of each letter

    - ``words`` -- a list of words, each an array of typecode ``'i'`` or a
      list, on the letters ``0 .. n - 1``; only the words of even index are
      read, which in a
      :class:`~combisurf.conjugate_tree.ConjugateTree` holding each word
      next to its inverse are the words without their inverses

    - ``weights`` -- a sequence of non-negative integers, of length exactly
      ``(len(words) + 1) // 2``, one entry per word read (``words[0]``,
      ``words[2]``, ...); the weight of the word ``words[2 * j]`` is
      ``weights[j]`` and must fit in a C ``long long``; the words of weight
      ``0`` are skipped

    OUTPUT: a pair ``(keys, key_weights)`` of arrays of typecode ``'q'``:
    the distinct keys by increasing order and, for each of them, the sum of
    the weights of its occurrences. The sum of all the weights read from
    ``words`` is exact, but see :func:`crossing_arcs_sweep_sorted` for the
    bound past which combining two such outputs by a sweep overflows.

    ALGORITHM:

    The arcs are written to a C array and sorted by ``qsort``, and the
    occurrences of a key, now consecutive, are merged. On two curves of
    length 8 in the one-vertex map of genus 32 (``n = 128``), the crossing
    arcs term of
    :meth:`~combisurf.geometric_intersection.GeometricIntersection.geometric_intersection`
    takes 1.9 us with two calls of this function and
    :func:`crossing_arcs_sweep_sorted`.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.crossing_arcs import word_arcs
        sage: word_arcs(4, [0, 2, 1, 3], [array('i', [0, 2]), array('i', [3, 1])], [1])
        (array('q', [9, 12]), array('q', [1, 1]))
        sage: word_arcs(4, [0, 2, 1, 3], [[0, 2], [3, 1], [0], [1]], [3, 2])
        (array('q', [8, 9, 12]), array('q', [2, 3, 3]))
        sage: word_arcs(4, [0, 2, 1, 3], [[0, 2], [3, 1], [0], [1]], [0, 2])
        (array('q', [8]), array('q', [2]))

    Feeding the output of two calls to :func:`crossing_arcs_sweep_sorted`
    counts the crossing arcs of two curves; here the identity angles put the
    letters at their own position on a hexagon, the word ``[0, 3]`` gives the
    arcs ``(0, 2)`` and ``(1, 3)``, the word ``[1, 4]`` gives ``(0, 4)`` and
    ``(1, 5)``, and only ``(0, 2)`` and ``(1, 5)`` cross::

        sage: from combisurf.crossing_arcs import crossing_arcs_sweep_sorted
        sage: n = 6
        sage: angles = list(range(n))
        sage: ukeys, uweights = word_arcs(n, angles, [array('i', [0, 3])], [1])
        sage: ukeys, uweights
        (array('q', [12, 19]), array('q', [1, 1]))
        sage: vkeys, vweights = word_arcs(n, angles, [array('i', [1, 4])], [1])
        sage: vkeys, vweights
        (array('q', [24, 31]), array('q', [1, 1]))
        sage: crossing_arcs_sweep_sorted(n, ukeys, uweights, vkeys, vweights)
        1

    TESTS::

        sage: word_arcs(4, [0, 2, 1, 3], [], [])
        (array('q'), array('q'))
        sage: word_arcs(4, [0, 2, 1, 3], [[0, 1]], [1])
        Traceback (most recent call last):
        ...
        ValueError: degenerate arc in word 0: the letter 1 is followed by its inverse
        sage: word_arcs(4, [0, 2, 1, 3], [[0, 4]], [1])
        Traceback (most recent call last):
        ...
        ValueError: invalid letter in word 0
        sage: word_arcs(4, [0, 2, 1], [], [])
        Traceback (most recent call last):
        ...
        ValueError: angles must have length n
        sage: word_arcs(4, [0, 2, 1, 4], [], [])
        Traceback (most recent call last):
        ...
        ValueError: invalid position in angles
        sage: word_arcs(0, [], [], [])
        Traceback (most recent call last):
        ...
        ValueError: n must be positive

    ``weights`` must have exactly ``(len(words) + 1) // 2`` entries: a
    shorter ``weights`` would be read out of bounds under
    ``boundscheck=False``, and a longer one points to a mismatch with
    ``words``::

        sage: word_arcs(4, [0, 2, 1, 3], [[0, 2], [3, 1]], [1, 2])
        Traceback (most recent call last):
        ...
        ValueError: weights must have (len(words) + 1) // 2 entries
        sage: word_arcs(4, [0, 2, 1, 3], [[0, 2], [3, 1]] * 2, [1])
        Traceback (most recent call last):
        ...
        ValueError: weights must have (len(words) + 1) // 2 entries

    Each weight must fit in a C ``long long``::

        sage: word_arcs(4, [0, 2, 1, 3], [[0, 2], [3, 1]], [2**63])
        Traceback (most recent call last):
        ...
        OverflowError: weight 9223372036854775808 does not fit in a long long
    """
    if n <= 0:
        raise ValueError("n must be positive")
    cdef array.array a_ang = angles if type(angles) is array.array and angles.typecode == 'i' else array.array('i', angles)
    cdef int *ang = a_ang.data.as_ints
    if len(a_ang) != n:
        raise ValueError("angles must have length n")
    cdef int c
    for c in range(n):
        if not 0 <= ang[c] < n:
            raise ValueError("invalid position in angles")

    cdef Py_ssize_t num_words = len(words)
    if len(weights) != (num_words + 1) // 2:
        raise ValueError("weights must have (len(words) + 1) // 2 entries")

    cdef Py_ssize_t i, total = 0
    for i in range(0, num_words, 2):
        wt = weights[i >> 1]
        if not -(1 << 63) <= wt < (1 << 63):
            raise OverflowError(f"weight {wt} does not fit in a long long")
        if wt:
            total += len(words[i])

    cdef array.array a_q = array.array('q', [])
    cdef array.array a_keys = array.clone(a_q, total, False)
    cdef array.array a_weights = array.clone(a_q, total, False)
    if total == 0:
        return a_keys, a_weights

    cdef _weighted_arc *arcs = <_weighted_arc *> malloc(total * sizeof(_weighted_arc))
    if arcs == NULL:
        raise MemoryError
    cdef array.array w
    cdef int *wd
    cdef Py_ssize_t l, p, num = 0, m
    cdef long long weight, first, last
    cdef int letter, previous
    cdef long long *keys = a_keys.data.as_longlongs
    cdef long long *key_weights = a_weights.data.as_longlongs
    try:
        for i in range(0, num_words, 2):
            weight = weights[i >> 1]
            if not weight:
                continue
            obj = words[i]
            w = obj if type(obj) is array.array and obj.typecode == 'i' else array.array('i', obj)
            wd = w.data.as_ints
            l = len(w)
            if l == 0:
                continue
            previous = wd[l - 1]
            for p in range(l):
                letter = wd[p]
                if not 0 <= letter < n or not 0 <= (previous ^ 1) < n:
                    raise ValueError(f"invalid letter in word {i}")
                first = ang[letter]
                last = ang[previous ^ 1]
                if first == last:
                    raise ValueError(f"degenerate arc in word {i}: the letter {previous} is followed by its inverse")
                if last < first:
                    first, last = last, first
                arcs[num].key = last * n + first
                arcs[num].weight = weight
                num += 1
                previous = letter

        qsort(arcs, num, sizeof(_weighted_arc), _cmp_weighted_arcs)

        m = 0
        keys[0] = arcs[0].key
        key_weights[0] = arcs[0].weight
        for p in range(1, num):
            if arcs[p].key == keys[m]:
                key_weights[m] += arcs[p].weight
            else:
                m += 1
                keys[m] = arcs[p].key
                key_weights[m] = arcs[p].weight
    finally:
        free(arcs)

    array.resize(a_keys, m + 1)
    array.resize(a_weights, m + 1)
    return a_keys, a_weights


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

    This computes the sum `S` of the module documentation, where the arc of
    key ``ukeys[k]`` has `u`-weight ``uweights[k]``, the arc of key
    ``vkeys[k]`` has `v`-weight ``vweights[k]``, and all the other weights
    are zero.

    INPUT:

    - ``n`` -- positive integer, the number of positions on the circle

    - ``ukeys``, ``uweights`` -- two arrays of typecode ``'q'`` and of the
      same length; ``ukeys`` holds strictly increasing keys ``last * n +
      first`` (with ``0 <= first < last < n``) and ``uweights`` the
      `u`-weights of these arcs, non-negative

    - ``vkeys``, ``vweights`` -- (default: ``None``) the same for the
      `v`-weights; when they are ``None``, the `v`-weights are taken equal to
      the `u`-weights, which counts the crossings of the `u`-arcs among
      themselves

    - ``scratch`` -- (default: ``None``) an array of typecode ``'q'``, of
      length at least ``2 * (n + 1)`` and filled with zeros, used as the
      memory of the sweep and filled with zeros again on return; it saves an
      allocation for each call when many of them are made with the same ``n``

    - ``check`` -- boolean (default: ``True``); whether to check that the
      keys are sorted and valid

    OUTPUT: the integer `S` of the module documentation, held together with
    the running sweep in a C ``long long``. Writing `U` for the sum of
    ``uweights`` and `V` for the sum of ``vweights`` (or `U` again when
    ``vweights`` is ``None``), `S` is at most `2 U V` and is exact only while
    `2 U V < 2^{63}`; past that bound it silently wraps around, since it is
    a sum of products of weights rather than of the weights themselves.

    ALGORITHM:

    The arcs are visited by merging the two arrays group by group, a group
    being the arcs with a given right endpoint. The `u`-arcs of a group are
    queried against the `v`-weights inserted so far, the `v`-arcs against the
    `u`-weights, and both are inserted afterwards. So there is no dictionary
    and no sort, and the cost is `O((|ukeys| + |vkeys|) \log(n))`. The cells
    of the Fenwick trees touched by the insertions are set back to zero at
    the end, in the same time.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.crossing_arcs import crossing_arcs_sweep_sorted
        sage: n = 4
        sage: crossing_arcs_sweep_sorted(n, array('q', [2 * n + 0]), array('q', [1]),
        ....:                               array('q', [3 * n + 1]), array('q', [1]))
        1
        sage: crossing_arcs_sweep_sorted(n, array('q', [2 * n + 0, 3 * n + 1]), array('q', [1, 1]))
        2

    On the hexagon of :func:`word_arcs`, ``ukeys`` holds the two arcs
    ``(0, 2)`` and ``(1, 3)`` of the word ``[0, 3]``, which already cross
    each other, and ``vkeys`` the two arcs ``(0, 4)`` and ``(1, 5)`` of the
    word ``[1, 4]``, of which only ``(1, 5)`` crosses an arc of ``ukeys``;
    this also leaves ``scratch`` as it found it::

        sage: n = 6
        sage: ukeys, uweights = array('q', [12, 19]), array('q', [1, 1])
        sage: vkeys, vweights = array('q', [24, 31]), array('q', [1, 1])
        sage: scratch = array('q', [0] * (2 * (n + 1)))
        sage: crossing_arcs_sweep_sorted(n, ukeys, uweights, vkeys, vweights, scratch)
        1
        sage: all(x == 0 for x in scratch)
        True
        sage: crossing_arcs_sweep_sorted(n, ukeys, uweights, scratch=scratch)
        2
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

    Reusing ``scratch`` across the two calls above leaves it as it found it::

        sage: scratch = array('q', [0] * (2 * (n + 1)))
        sage: startpoint_sweep_sorted(n, q(3), q(0), q(1), q(4), q(0), q(5), scratch)
        1
        sage: all(x == 0 for x in scratch)
        True
        sage: startpoint_sweep_sorted(n, q(3, 4), q(0, 0), q(1, 5), scratch=scratch)
        2
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


cdef array.array _int_array = array.array('i', [])


cdef array.array _as_int_array(w):
    r"""
    Return ``w`` if it is an array of typecode ``'i'`` and a copy of it as
    one otherwise.
    """
    if type(w) is array.array and (<array.array> w).ob_descr.typecode == b'i':
        return <array.array> w
    return array.array('i', w)


def tree_add_with_inverse(ConjugateTree T not None, w):
    r"""
    Add the free group word ``w`` and, if it is new, its inverse to the
    conjugate tree ``T``.

    The letter ``h ^ 1`` is the inverse of the letter ``h``, so the alphabet
    of ``T``, if it is known, must have even size.

    INPUT:

    - ``T`` -- a :class:`~combisurf.conjugate_tree.ConjugateTree` every word
      of which has been added by this function, so that the words are closed
      under inverse. The only cheap evidence of the contrary is an odd number
      of words, on which this function refuses to run: in particular once
      :meth:`~combisurf.conjugate_tree.ConjugateTree.process` has added a
      single word.

    - ``w`` -- a non-empty cyclically reduced word, given as any sequence of
      integers; unlike :meth:`~combisurf.conjugate_tree.ConjugateTree.process`,
      this function does not convert ``w`` with
      :func:`~combisurf.word.word_init`

    OUTPUT: a pair ``(i, exponent)`` where ``w`` is conjugate to the
    ``exponent``-th power of the ``i``-th word of ``T``. If ``w`` is new,
    its primitive root gets the index ``i``, which is even, and its inverse
    the index ``i + 1``.

    EXAMPLES::

        sage: from combisurf.conjugate_tree import ConjugateTree
        sage: from combisurf.crossing_arcs import tree_add_with_inverse
        sage: T = ConjugateTree(4)
        sage: tree_add_with_inverse(T, [0, 2, 0, 3])
        (0, 1)
        sage: T.words()
        [array('i', [0, 2, 0, 3]), array('i', [2, 1, 3, 1])]
        sage: tree_add_with_inverse(T, [2, 0, 2, 0])
        (2, 2)
        sage: T.words()
        [array('i', [0, 2, 0, 3]), array('i', [2, 1, 3, 1]), array('i', [2, 0]), array('i', [1, 3])]
        sage: tree_add_with_inverse(T, [0, 3, 0, 2])
        (0, 1)
        sage: tree_add_with_inverse(T, [1, 2, 1, 3])
        (1, 1)
        sage: tree_add_with_inverse(T, [0, 2, 0, 2, 0, 2])
        (2, 3)

    TESTS:

    A word already present, then a power of a present word::

        sage: T = ConjugateTree(4)
        sage: tree_add_with_inverse(T, [0, 2])
        (0, 1)
        sage: tree_add_with_inverse(T, [2, 0])
        (0, 1)
        sage: tree_add_with_inverse(T, [3, 1])
        (1, 1)
        sage: tree_add_with_inverse(T, [1, 3, 1, 3])
        (1, 2)
        sage: T.num_words()
        2

    A word of length one, which is its own conjugate only::

        sage: tree_add_with_inverse(T, [1])
        (2, 1)
        sage: tree_add_with_inverse(T, [0, 0])
        (3, 2)

    Invalid input leaves the tree as it was::

        sage: T = ConjugateTree(4)
        sage: tree_add_with_inverse(T, [0, 2, 3])
        Traceback (most recent call last):
        ...
        ValueError: w must be cyclically reduced
        sage: tree_add_with_inverse(T, [2, 1, 0])
        Traceback (most recent call last):
        ...
        ValueError: w must be cyclically reduced
        sage: tree_add_with_inverse(T, [])
        Traceback (most recent call last):
        ...
        ValueError: empty word in input
        sage: tree_add_with_inverse(T, [0, 4])
        Traceback (most recent call last):
        ...
        ValueError: invalid word: letter 4 not in the alphabet {0, 1, ..., 3}
        sage: tree_add_with_inverse(T, [0, -1])
        Traceback (most recent call last):
        ...
        ValueError: invalid word: must be made of non-negative integers
        sage: T.num_words(), T.num_states()
        (0, 1)

    A tree whose words are not closed under inverse::

        sage: T = ConjugateTree()
        sage: T.process([1])
        1
        sage: tree_add_with_inverse(T, [0])
        Traceback (most recent call last):
        ...
        ValueError: the words of this tree are not closed under inverse
        sage: T.num_words(), T.num_states()
        (1, 2)

    An alphabet of odd size is not the alphabet of a free group::

        sage: T = ConjugateTree(3)
        sage: tree_add_with_inverse(T, [2])
        Traceback (most recent call last):
        ...
        ValueError: the alphabet size (=3) must be even
        sage: T.num_words(), T.num_states()
        (0, 1)

    A long word, whose inverse does not fit in the buffer on the stack::

        sage: w = [0, 2] * 50 + [0, 3]
        sage: T = ConjugateTree(4)
        sage: tree_add_with_inverse(T, w)
        (0, 1)
        sage: list(T.word(1)) == [h ^^ 1 for h in reversed(w)]
        True
        sage: T._check()
    """
    cdef array.array a = _as_int_array(w)
    cdef int *v = a.data.as_ints
    cdef Py_ssize_t size = len(a)
    cdef int n = T.T.alphabet_size
    cdef int l, j, h, status, i, r, check
    cdef int stack_buf[64]
    cdef int *buf = stack_buf
    cdef const int *src

    if n % 2:
        raise ValueError(f"the alphabet size (={n}) must be even")
    if size == 0:
        raise ValueError("empty word in input")
    if size > INT_MAX // 2:
        T._raise(CT_ETOOLARGE, a)
    l = <int> size

    # NOTE: all the checks come before any insertion, since a node is never
    # removed from the tree; ct_process checks the letters again
    for j in range(l):
        h = v[j]
        if h < 0:
            T._raise(CT_ENEGATIVE, a)
        if n and h >= n:
            T._raise(CT_EALPHABET, a)
    for j in range(l):
        if v[j] ^ 1 == v[j + 1 if j + 1 < l else 0]:
            raise ValueError("w must be cyclically reduced")
    if T.T.nwords % 2:
        raise ValueError("the words of this tree are not closed under inverse")

    # w and the inverse of its primitive root, which is at most as long
    T._reserve(2, 2 * l)

    T._process(v, l, &status)
    if status <= 0:
        # w is conjugate to a power of a word already present
        i = -status
        if l % T.T.wlen[i]:
            raise RuntimeError("conjugate tree: a word conjugate to a power of a word "
                               "whose length does not divide its own")
        return (i, l // T.T.wlen[i])

    # NOTE: the inverse is built outside of the word buffer of T, which
    # ct_process writes to
    i = T.T.nwords - 1
    r = T.T.wlen[i]
    if r > 64:
        buf = <int *> malloc(r * sizeof(int))
        if buf == NULL:
            raise MemoryError
    src = T.T.wbuf + T.T.wstart[i]
    for j in range(r):
        buf[j] = src[r - 1 - j] ^ 1
    try:
        T._process(buf, r, &check)
    finally:
        if buf != stack_buf:
            free(buf)
    if check != 1:
        # NOTE: the inverse of a cyclically reduced word is not conjugate to a
        # power of it in a free group, and the words were closed under
        # inverse, so the inverse of the new word w is new and primitive
        raise RuntimeError("conjugate tree: the inverse of a new word is not new and primitive")
    return (i, status)


def cyclically_sorted_leaf_arcs(ConjugateTree T not None, angles):
    r"""
    Return the leaves of the conjugate tree ``T`` in their cyclic order at
    infinity, each one described by its word, its first letter and the angle
    it turns.

    The letter ``h ^ 1`` is the inverse of the letter ``h`` and ``angles[h]``
    is the position of the half-edge ``h`` around the vertex. The order is
    the one of
    :meth:`~combisurf.conjugate_tree.ConjugateTree.cyclically_sorted_leaves`
    with ``order = angles`` and ``pivot[b] = angles[b ^ 1]``: below a node,
    the angles are measured from the half-edge through which the curve came
    in.

    The leaf of the conjugate ``(i, k)`` (see
    :meth:`~combisurf.conjugate_tree.ConjugateTree.leaf_as_conjugate`) of the
    word ``w = T.word(i)`` starts with the letter ``w[k]`` and ends with the
    reverse ``w[k - 1] ^ 1`` of the letter before it. Its angle is
    ``(angles[w[k - 1] ^ 1] - angles[w[k]]) % n - 1`` where ``n`` is the
    length of ``angles``.

    INPUT:

    - ``T`` -- a conjugate tree

    - ``angles`` -- a permutation of ``{0, 1, ..., n - 1}``, where ``n`` is
      even and larger than every letter of ``T``

    OUTPUT: a triple of arrays of typecode ``'q'`` with one entry per leaf:
    the index ``i`` of the word of the leaf, its first letter and its angle

    EXAMPLES::

        sage: from combisurf.conjugate_tree import ConjugateTree
        sage: from combisurf.crossing_arcs import cyclically_sorted_leaf_arcs
        sage: T = ConjugateTree()
        sage: T.process([0, 2, 1, 3])
        1
        sage: T.process([2, 0, 3, 1])
        1
        sage: angles = [0, 2, 1, 3]
        sage: cyclically_sorted_leaf_arcs(T, angles)
        (array('q', [1, 0, 1, 0, 1, 0, 1, 0]),
         array('q', [0, 0, 2, 2, 1, 1, 3, 3]),
         array('q', [2, 0, 2, 0, 2, 0, 2, 0]))

    It describes the leaves of
    :meth:`~combisurf.conjugate_tree.ConjugateTree.cyclically_sorted_leaves`::

        sage: n = len(angles)
        sage: pivot = [angles[b ^^ 1] for b in range(n)]
        sage: ans = []
        sage: for s in T.cyclically_sorted_leaves(angles, pivot):
        ....:     i, k = T.leaf_as_conjugate(s)
        ....:     w = T.word(i)
        ....:     ans.append((i, w[k], (angles[w[k - 1] ^^ 1] - angles[w[k]]) % n - 1))
        sage: list(zip(*cyclically_sorted_leaf_arcs(T, angles))) == ans
        True

    TESTS::

        sage: cyclically_sorted_leaf_arcs(T, [0, 1])
        Traceback (most recent call last):
        ...
        ValueError: the letters of this tree do not fit in an alphabet of size 2
        sage: cyclically_sorted_leaf_arcs(T, [0, 2, 1])
        Traceback (most recent call last):
        ...
        ValueError: the length of angles (=3) must be even
    """
    cdef array.array a_ang = T._order_array(angles)
    cdef int *ang = a_ang.data.as_ints
    cdef int n = len(a_ang)
    if n % 2:
        raise ValueError(f"the length of angles (={n}) must be even")

    cdef array.array a_pivot = array.clone(_int_array, n, False)
    cdef int *pivot = a_pivot.data.as_ints
    cdef int c
    for c in range(n):
        pivot[c] = ang[c ^ 1]

    cdef array.array a_leaves = array.clone(_int_array, T.T.nstates, False)
    cdef int *out = a_leaves.data.as_ints
    cdef int num = T._sorted_leaves(a_ang, a_pivot, out)

    cdef array.array a_q = array.array('q', [])
    cdef array.array a_word = array.clone(a_q, num, False)
    cdef array.array a_first = array.clone(a_q, num, False)
    cdef array.array a_turn = array.clone(a_q, num, False)
    cdef long long *word = a_word.data.as_longlongs
    cdef long long *first = a_first.data.as_longlongs
    cdef long long *turn = a_turn.data.as_longlongs
    cdef int j, s, i, l, k, a, b, t
    for j in range(num):
        s = out[j]
        # the conjugate (i, k) of the leaf, as in leaf_as_conjugate
        i = T.T.tword[s]
        l = T.T.wlen[i]
        k = (T.T.tstart[s] - T.T.dep[T.T.parent[s]]) % l
        if k < 0:
            k += l
        a = T.T.wbuf[T.T.wstart[i] + k]
        b = T.T.wbuf[T.T.wstart[i] + (k - 1 if k else l - 1)] ^ 1
        t = ang[b] - ang[a]
        if t < 0:
            t += n
        word[j] = i
        first[j] = a
        turn[j] = t - 1
    return (a_word, a_first, a_turn)


def _leaf_weights(array.array word_index not None, weights):
    r"""
    Return the weights of the leaves of a family of words, from the index of
    the word of each leaf.

    This builds the ``uweights`` and ``vweights`` of
    :func:`startpoint_sweep_weighted` from the first array returned by
    :func:`cyclically_sorted_leaf_arcs`.

    INPUT:

    - ``word_index`` -- an array of typecode ``'q'``, the index ``i`` of the
      word of each leaf

    - ``weights`` -- a sequence of integers, the weight of the words of
      indices ``2 * j`` and ``2 * j + 1`` being ``weights[j]``

    OUTPUT: an array of typecode ``'q'``, the entry ``k`` being
    ``weights[word_index[k] >> 1]``

    ALGORITHM:

    A loop in C. On the 32 leaves of two curves of length 8 it takes 0.2 us,
    against 1.5 us for ``array('q', map(table.__getitem__, word_index))``
    with a table holding each weight twice.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.crossing_arcs import _leaf_weights
        sage: _leaf_weights(array('q', [0, 3, 1, 2, 3]), [5, 7])
        array('q', [5, 7, 5, 7, 7])

    TESTS::

        sage: _leaf_weights(array('q', [4]), [5, 7])
        Traceback (most recent call last):
        ...
        ValueError: invalid word index 4
        sage: _leaf_weights(array('q'), [])
        array('q')
    """
    cdef long long *index = _as_longlongs(word_index, "word_index")
    cdef Py_ssize_t num = len(word_index)
    cdef array.array a_table = array.array('q', weights)
    cdef long long *table = a_table.data.as_longlongs
    cdef Py_ssize_t size = len(a_table)
    cdef array.array a_ans = array.clone(a_table, num, False)
    cdef long long *ans = a_ans.data.as_longlongs
    cdef Py_ssize_t k
    cdef long long j
    for k in range(num):
        j = index[k] >> 1
        if not 0 <= j < size:
            raise ValueError(f"invalid word index {index[k]}")
        ans[k] = table[j]
    return a_ans


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

    OUTPUT: an integer, even when ``vweights`` is ``None``; exact only while
    twice the product of the total `u`-weight and the total `v`-weight given
    is below `2^{63}`, for the same reason as the bound given in the OUTPUT
    of :func:`crossing_arcs_sweep_sorted`

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

    Reusing ``scratch`` across the two calls above leaves it as it found it::

        sage: scratch = array('q', [0] * (2 * (8 + 1)))
        sage: startpoint_sweep_weighted(8, q(0, 0), q(1, 5), q(2, 0), q(0, 3), scratch)
        6
        sage: all(x == 0 for x in scratch)
        True
        sage: startpoint_sweep_weighted(8, q(0, 0), q(1, 5), q(2, 3), scratch=scratch)
        12
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
