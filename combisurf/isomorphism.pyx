# distutils: sources = combisurf/src/isomorphism.c
# distutils: include_dirs = combisurf/src
r"""
Isomorphisms, automorphisms and canonical labels of maps

This module computes canonical labellings, canonical forms, automorphism
groups and isomorphisms of maps given by their face permutation, possibly
decorated by integer values on half-edges. It computes exactly what
:mod:`combisurf.isomorphism_naive` computes, which is its reference: the
functions of the same names take the same arguments and return the same
values (the same canonical labelling, the same canonical form, the
automorphisms in the same order, the same isomorphism), and raise the same
exceptions on malformed input. The conventions (face permutations with
``-1`` on inactive half-edges, the decorations ``data`` kept exactly and
``classes`` kept up to a relabelling of their values, the breadth-first
walk and the key of a root) are described in
:mod:`combisurf.isomorphism_naive`. The methods of
:class:`~combisurf.oriented_map.OrientedMap` (such as
:meth:`~combisurf.oriented_map.OrientedMap.canonical_labelling`) call this
module.

On a connected map, every active half-edge is tried as a root: the
relabelling from it and its key are computed in C (``combisurf/src/isomorphism.c``)
and the minimal key is kept, which costs `O(n^2)` for `n` half-edges. A
disconnected map is currently handed to the function of the same name in
:mod:`combisurf.isomorphism_naive`.

The decorations are stored in C ``int`` arrays: a value of ``data``
outside the range of a C ``int`` raises an ``OverflowError`` (the values
on inactive half-edges are ignored and never converted). The values of
``classes`` are renumbered before the conversion, so they may be any
hashable values.

Besides the functions of :mod:`combisurf.isomorphism_naive`,
:func:`canonical_key` returns the canonical form flattened as bytes, a
compact hashable key for the exploration of large collections of maps up
to isomorphism, and on request the canonical labelling and the number of
automorphisms.

EXAMPLES::

    sage: from combisurf.isomorphism import canonical_labelling, automorphisms, isomorphism, canonical_key
    sage: fp = [2, 5, 4, 1, 0, 3]  # the faces (0,1,2)(~0,~2,~1)
    sage: canonical_labelling(fp)
    array('i', [0, 1, 2, 3, 4, 5])
    sage: len(automorphisms(fp))
    6
    sage: isomorphism(fp, (), [4, 3, 0, 5, 2, 1], ())
    array('i', [0, 1, 4, 5, 2, 3])
    sage: canonical_key(fp) == canonical_key([4, 3, 0, 5, 2, 1])
    True

The results agree with the reference implementation::

    sage: from combisurf import isomorphism_naive
    sage: fp = [2, 3, 4, 5, 0, 1]
    sage: c = [0, 1, 0, 1, 0, 1]
    sage: automorphisms(fp, classes=(c,)) == isomorphism_naive.automorphisms(fp, classes=(c,))
    True
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
from libc.string cimport memcpy, memcmp

import operator

from combisurf import isomorphism_naive


cdef array.array _int_array = array.array('i', [])


cdef inline array.array _new(Py_ssize_t size):
    return array.clone(_int_array, size, False)


cdef array.array _as_int_array(a):
    r"""
    Return ``a`` if it is an array of typecode ``'i'`` and a copy of it as
    such an array otherwise.
    """
    if type(a) is array.array and (<array.array> a).ob_descr.typecode == b'i':
        return <array.array> a
    return array.array('i', a)


cdef class _Map:
    r"""
    A face permutation and its decorations converted to C ``int`` arrays.

    ``dec`` holds the ``nd`` arrays of ``data`` then the ``nc`` renumbered
    arrays of ``classes``, each of ``n`` entries (``0`` on the inactive
    half-edges). ``ne`` is the number of active edges, ``na`` the number of
    active half-edges.
    """
    cdef array.array fp
    cdef array.array dec
    cdef int n
    cdef int nd
    cdef int nc
    cdef int ne
    cdef int na


cdef _check_input(fp, data, classes):
    r"""
    Check the lengths of ``fp`` and of the arrays of ``data`` and
    ``classes``, with the messages of ``isomorphism_naive._check_input``.
    """
    n = len(fp)
    if n % 2:
        raise ValueError("fp must have even length")
    if not isinstance(data, tuple):
        raise TypeError("data must be a tuple")
    for d in data:
        if len(d) != n:
            raise ValueError("data arrays must have length {}".format(n))
    if not isinstance(classes, tuple):
        raise TypeError("classes must be a tuple")
    for c in classes:
        if len(c) != n:
            raise ValueError("classes arrays must have length {}".format(n))


cdef _Map _convert(fp, tuple data, tuple classes):
    r"""
    Return the :class:`_Map` of ``fp`` decorated by ``data`` and
    ``classes``, whose lengths are checked already.

    Raise ``ValueError`` if ``fp`` is not a permutation of its active
    half-edges and ``OverflowError`` if a value of ``data`` on an active
    half-edge is outside the range of a C ``int``.
    """
    cdef _Map m = _Map.__new__(_Map)
    cdef int n, h, x, j, off, i
    cdef int *f
    cdef int *dec
    cdef int *a
    cdef array.array seen
    m.fp = _as_int_array(fp)
    n = m.n = len(m.fp)
    m.nd = len(data)
    m.nc = len(classes)
    f = m.fp.data.as_ints

    # a map from the active half-edges to themselves that is one-to-one is
    # a permutation of them
    seen = _new(n)
    for h in range(n):
        seen.data.as_ints[h] = 0
    m.na = 0
    for h in range(n):
        x = f[h]
        if x == -1:
            continue
        if x < 0 or x >= n or f[x] == -1 or seen.data.as_ints[x]:
            raise ValueError("fp must be a permutation of its active half-edges")
        seen.data.as_ints[x] = 1
        m.na += 1
    m.ne = iso_num_edges(n, f)

    m.dec = _new((m.nd + m.nc) * n)
    dec = m.dec.data.as_ints
    for j in range(m.nd):
        d = data[j]
        off = j * n
        if type(d) is array.array and (<array.array> d).ob_descr.typecode == b'i':
            a = (<array.array> d).data.as_ints
            memcpy(dec + off, a, n * sizeof(int))
            for h in range(n):
                if f[h] == -1:
                    dec[off + h] = 0
        else:
            for h in range(n):
                if f[h] == -1:
                    dec[off + h] = 0
                else:
                    dec[off + h] = d[h]
    for j in range(m.nc):
        c = classes[j]
        off = (m.nd + j) * n
        # renumbered by first appearance in half-edge order, which leaves
        # the key unchanged (it renumbers them again in the order of the
        # new labels) and puts the values in range(n) as iso_key needs
        names = {}
        for h in range(n):
            if f[h] == -1:
                dec[off + h] = 0
            else:
                v = c[h]
                i = names.setdefault(v, len(names))
                dec[off + h] = i
    return m


cdef _Map _prepare(fp, data, classes):
    r"""
    Check the input as ``isomorphism_naive._check_input`` does and convert
    it.
    """
    _check_input(fp, data, classes)
    return _convert(fp, data, classes)


cdef int _run(_Map m, array.array key, int *relabelling, int *roots, int *nroots) except? -2:
    r"""
    Run ``iso_canonical`` on ``m``: write the minimal key into ``key`` (of
    ``ISO_KEY_SIZE`` entries) and return its length, or ``ISO_DISCONNECTED``.
    """
    cdef array.array work = _new(ISO_WORK_SIZE(m.n, m.nd, m.nc))
    cdef int *dec = m.dec.data.as_ints
    return iso_canonical(m.n, m.fp.data.as_ints, m.nd, dec, m.nc, dec + m.nd * m.n,
                         key.data.as_ints, relabelling, roots, nroots, work.data.as_ints)


cdef array.array _complete(array.array r, int n, int ne):
    r"""
    Complete in place the relabelling ``r`` (``-1`` on the half-edges of the
    inactive edges) by the labels ``2*ne, 2*ne + 1, ...`` given to the
    inactive edges in increasing order, as ``isomorphism_naive._complete``.
    """
    cdef int e, k = 2 * ne
    cdef int *a = r.data.as_ints
    for e in range(n // 2):
        if a[2 * e] == -1:
            a[2 * e] = k
            a[2 * e + 1] = k + 1
            k += 2
    return r


def canonical_labelling(fp, data=(), classes=()):
    r"""
    Return the canonical labelling of the map with face permutation ``fp``
    decorated by ``data`` and ``classes``.

    The result is the one of
    :func:`combisurf.isomorphism_naive.canonical_labelling`.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``data`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept exactly (default: no decoration); the values on the
      active half-edges are integers in the range of a C ``int``

    - ``classes`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept up to a relabelling of their values (default: no
      decoration)

    OUTPUT: an ``array('i')`` ``r`` that maps each half-edge to its
    canonical label. It is a permutation of ``range(len(fp))`` that
    preserves the pairs ``{2e, 2e+1}``. The active edges get the labels
    ``0, ..., 2*ne - 1`` (``ne`` the number of active edges) and the
    inactive edges the following ones, in increasing order.

    EXAMPLES::

        sage: from combisurf.isomorphism import canonical_labelling
        sage: canonical_labelling([2, 5, 4, 1, 0, 3])
        array('i', [0, 1, 2, 3, 4, 5])
        sage: canonical_labelling([2, 5, 4, 1, 0, 3], ([0, 0, 1, 0, 0, 0],))
        array('i', [1, 0, 5, 4, 3, 2])
        sage: canonical_labelling([2, 5, 4, 1, 0, 3], classes=([5, 5, 2, 5, 5, 5],))
        array('i', [1, 0, 5, 4, 3, 2])

    TESTS:

    A folded edge, a label gap, a single half-edge and the empty map::

        sage: canonical_labelling([2, -1, 0, 3])
        array('i', [2, 3, 1, 0])
        sage: canonical_labelling([4, -1, -1, -1, 5, 0])
        array('i', [2, 3, 4, 5, 0, 1])
        sage: canonical_labelling([0, -1])
        array('i', [0, 1])
        sage: canonical_labelling([])
        array('i')

    A disconnected map, handed to :mod:`combisurf.isomorphism_naive`::

        sage: canonical_labelling([2, 5, 4, 1, 0, 3, 6, -1])
        array('i', [2, 3, 4, 5, 6, 7, 0, 1])

    Malformed input::

        sage: canonical_labelling([0, -1, 2])
        Traceback (most recent call last):
        ...
        ValueError: fp must have even length
        sage: canonical_labelling([0, -1], [[0, 0]])
        Traceback (most recent call last):
        ...
        TypeError: data must be a tuple
        sage: canonical_labelling([0, -1], ([0],))
        Traceback (most recent call last):
        ...
        ValueError: data arrays must have length 2
        sage: canonical_labelling([1, 1])
        Traceback (most recent call last):
        ...
        ValueError: fp must be a permutation of its active half-edges

    A value of ``data`` outside the range of a C ``int``, and one on an
    inactive half-edge, which is ignored::

        sage: canonical_labelling([1, 0], ([2**40, 0],))
        Traceback (most recent call last):
        ...
        OverflowError: ...
        sage: canonical_labelling([0, -1], ([0, 2**40],))
        array('i', [0, 1])
    """
    cdef _Map m = _prepare(fp, data, classes)
    cdef array.array key = _new(ISO_KEY_SIZE(m.n, m.nd, m.nc))
    cdef array.array r = _new(m.n)
    if _run(m, key, r.data.as_ints, NULL, NULL) == ISO_DISCONNECTED:
        return isomorphism_naive.canonical_labelling(fp, data, classes)
    return _complete(r, m.n, m.ne)


cdef tuple _form(_Map m, array.array key, array.array r, tuple data):
    r"""
    Return the canonical form of ``m`` from its minimal key and the
    relabelling ``r`` that reaches it; the values of ``data`` are read from
    the original arrays, as in ``isomorphism_naive._key``.
    """
    cdef int k = 2 * m.ne
    cdef int h, i, pos
    cdef int *kd = key.data.as_ints
    rinv = [-1] * k
    for h in range(m.n):
        if m.fp.data.as_ints[h] != -1:
            rinv[r.data.as_ints[h]] = h
    out = [tuple([kd[i] for i in range(k)])]
    for d in data:
        out.append(tuple([d[x] for x in rinv if x != -1]))
    pos = k + m.nd * m.na
    while len(out) < 1 + m.nd + m.nc:
        out.append(tuple([kd[i] for i in range(pos, pos + m.na)]))
        pos += m.na
    return tuple(out)


def canonical_form(fp, data=(), classes=()):
    r"""
    Return the canonical form of the map with face permutation ``fp``
    decorated by ``data`` and ``classes``.

    The result is the one of
    :func:`combisurf.isomorphism_naive.canonical_form`.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``data`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept exactly (default: no decoration); the values on the
      active half-edges are integers in the range of a C ``int``

    - ``classes`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept up to a relabelling of their values (default: no
      decoration)

    OUTPUT: a tuple ``(fp_c, d_1, ..., d_k, c_1, ..., c_m)`` of tuples:
    ``fp_c`` is the face permutation in the canonical labels (of length
    ``2*ne``, ``-1`` on the inactive half-edges of folded edges), ``d_i``
    the values of the ``i``-th array of ``data`` read in the increasing
    order of the canonical labels of the active half-edges, and ``c_i`` the
    values of the ``i``-th array of ``classes`` read in the same order and
    renumbered ``0, 1, 2, ...`` by first appearance. Two decorated maps are
    isomorphic if and only if their canonical forms are equal.

    EXAMPLES::

        sage: from combisurf.isomorphism import canonical_form
        sage: canonical_form([2, 5, 4, 1, 0, 3])
        ((2, 5, 4, 1, 0, 3),)
        sage: canonical_form([4, 3, 0, 5, 2, 1])
        ((2, 5, 4, 1, 0, 3),)
        sage: canonical_form([2, 5, 4, 1, 0, 3], ([0, 0, 1, 0, 0, 0],))
        ((2, 5, 4, 1, 0, 3), (0, 0, 0, 0, 0, 1))
        sage: canonical_form([2, 3, 4, 5, 0, 1], ([1, 0, 1, 0, 1, 0],), ([7, 3, 7, 3, 7, 3],))
        ((2, 3, 4, 5, 0, 1), (0, 1, 0, 1, 0, 1), (0, 1, 0, 1, 0, 1))

    TESTS:

    A folded edge, a label gap, the empty map and a disconnected map::

        sage: canonical_form([2, -1, 0, 3])
        ((0, 2, 1, -1),)
        sage: canonical_form([4, -1, -1, -1, 5, 0])
        ((1, 2, 0, -1),)
        sage: canonical_form([], ([],), ([], []))
        ((), (), (), ())
        sage: canonical_form([0, 1, 2, 3], classes=([0, 1, 1, 0],))
        ((0, 1, 2, 3), (0, 1, 0, 1))

    A value of ``data`` outside the range of a C ``int``::

        sage: canonical_form([1, 0], ([0, -2**40],))
        Traceback (most recent call last):
        ...
        OverflowError: ...
    """
    cdef _Map m = _prepare(fp, data, classes)
    cdef array.array key = _new(ISO_KEY_SIZE(m.n, m.nd, m.nc))
    cdef array.array r = _new(m.n)
    if _run(m, key, r.data.as_ints, NULL, NULL) == ISO_DISCONNECTED:
        return isomorphism_naive.canonical_form(fp, data, classes)
    return _form(m, key, r, data)


def automorphisms(fp, data=(), classes=()):
    r"""
    Return the automorphisms of the map with face permutation ``fp``
    decorated by ``data`` and ``classes``.

    The result is the one of
    :func:`combisurf.isomorphism_naive.automorphisms`, in the same order.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``data`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept exactly (default: no decoration); the values on the
      active half-edges are integers in the range of a C ``int``

    - ``classes`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept up to a relabelling of their values (default: no
      decoration)

    OUTPUT: the list of the automorphisms, each an ``array('i')``
    permutation of ``range(len(fp))`` that commutes with ``fp`` on the
    active half-edges, preserves the pairs ``{2e, 2e+1}``, the arrays of
    ``data`` and the partitions defined by the arrays of ``classes``, and is
    the identity on the inactive edges. The identity comes first. On a
    connected map, with ``b_0 < b_1 < ...`` the roots of minimal key, the
    automorphism ``j`` sends the half-edge of label ``i`` from ``b_0`` to
    the half-edge of label ``i`` from ``b_j``.

    EXAMPLES::

        sage: from combisurf.isomorphism import automorphisms
        sage: A = automorphisms([2, 5, 4, 1, 0, 3])
        sage: len(A)
        6
        sage: A[:2]
        [array('i', [0, 1, 2, 3, 4, 5]), array('i', [1, 0, 5, 4, 3, 2])]
        sage: automorphisms([2, 5, 4, 1, 0, 3], ([0, 0, 1, 0, 0, 0],))
        [array('i', [0, 1, 2, 3, 4, 5])]
        sage: fp = [2, 3, 4, 5, 0, 1]
        sage: c = [0, 1, 0, 1, 0, 1]
        sage: len(automorphisms(fp)), len(automorphisms(fp, (c,))), len(automorphisms(fp, classes=(c,)))
        (6, 3, 6)

    TESTS:

    A folded edge, a label gap, the empty map and a disconnected map::

        sage: automorphisms([2, -1, 0, 3])
        [array('i', [0, 1, 2, 3])]
        sage: automorphisms([4, -1, -1, -1, 5, 0])
        [array('i', [0, 1, 2, 3, 4, 5])]
        sage: automorphisms([])
        [array('i')]
        sage: automorphisms([-1, -1])
        [array('i', [0, 1])]
        sage: automorphisms([0, 1, 2, 3], classes=([0, 1, 1, 2],))
        [array('i', [0, 1, 2, 3]), array('i', [3, 2, 1, 0])]
    """
    cdef _Map m = _prepare(fp, data, classes)
    cdef int n = m.n
    cdef int count, j, h
    cdef array.array key = _new(ISO_KEY_SIZE(n, m.nd, m.nc))
    cdef array.array ra = _new(n)
    cdef array.array roots = _new(n)
    cdef array.array rb, rbinv, queue, g
    if _run(m, key, ra.data.as_ints, roots.data.as_ints, &count) == ISO_DISCONNECTED:
        return isomorphism_naive.automorphisms(fp, data, classes)
    if count == 0:
        return [array.array('i', range(n))]
    rb = _new(n)
    rbinv = _new(n)
    queue = _new(n)
    result = []
    for j in range(count):
        # iso_relabelling_from writes -1 into rbinv at the labels of the
        # inactive halves of folded edges; the automorphism maps those
        # halves too, so the inverse is rebuilt on every labelled half-edge
        iso_relabelling_from(n, m.fp.data.as_ints, roots.data.as_ints[j],
                             rb.data.as_ints, rbinv.data.as_ints, queue.data.as_ints)
        for h in range(n):
            if rb.data.as_ints[h] != -1:
                rbinv.data.as_ints[rb.data.as_ints[h]] = h
        g = _new(n)
        for h in range(n):
            if ra.data.as_ints[h] == -1:
                g.data.as_ints[h] = h
            else:
                g.data.as_ints[h] = rbinv.data.as_ints[ra.data.as_ints[h]]
        result.append(g)
    return result


def automorphism_generators(fp, data=(), classes=()):
    r"""
    Return generators of the automorphism group of the map with face
    permutation ``fp`` decorated by ``data`` and ``classes``.

    The result is the one of
    :func:`combisurf.isomorphism_naive.automorphism_generators`.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``data`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept exactly (default: no decoration); the values on the
      active half-edges are integers in the range of a C ``int``

    - ``classes`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept up to a relabelling of their values (default: no
      decoration)

    OUTPUT: a list of automorphisms (as in :func:`automorphisms`) that
    generate the group; empty when the group is trivial. The automorphisms
    are scanned in the order of :func:`automorphisms` and one is kept when
    it is not in the group generated by those already kept.

    EXAMPLES::

        sage: from combisurf.isomorphism import automorphism_generators
        sage: automorphism_generators([2, 5, 4, 1, 0, 3])
        [array('i', [1, 0, 5, 4, 3, 2]), array('i', [2, 3, 4, 5, 0, 1])]
        sage: c = [0, 1, 0, 1, 0, 1]
        sage: automorphism_generators([2, 3, 4, 5, 0, 1], (c,))
        [array('i', [2, 3, 4, 5, 0, 1])]
        sage: automorphism_generators([2, 3, 4, 5, 0, 1], classes=(c,))
        [array('i', [1, 0, 3, 2, 5, 4]), array('i', [2, 3, 4, 5, 0, 1])]

    TESTS::

        sage: automorphism_generators([2, -1, 0, 3])
        []
        sage: automorphism_generators([0, 1, 2, 3])
        [array('i', [0, 1, 3, 2]), array('i', [1, 0, 2, 3]), array('i', [2, 3, 0, 1])]
    """
    n = len(fp)
    gens = []
    group = {tuple(range(n))}
    for g in automorphisms(fp, data, classes):
        if tuple(g) not in group:
            gens.append(g)
            group = isomorphism_naive._group_closure(gens, n)
    return gens


def isomorphism(fp1, data1, fp2, data2, classes1=(), classes2=()):
    r"""
    Return an isomorphism from the first decorated map to the second, or
    ``None`` if they are not isomorphic.

    The result is the one of :func:`combisurf.isomorphism_naive.isomorphism`.

    INPUT:

    - ``fp1``, ``fp2`` -- face permutations (lists or ``array('i')``,
      ``-1`` on inactive half-edges)

    - ``data1``, ``data2`` -- tuples with the same number of arrays indexed
      by half-edges, of lengths ``len(fp1)`` and ``len(fp2)``, kept
      exactly; the values on the active half-edges are integers in the
      range of a C ``int``

    - ``classes1``, ``classes2`` -- tuples with the same number of arrays
      indexed by half-edges, of lengths ``len(fp1)`` and ``len(fp2)``, kept
      up to a relabelling of their values (default: no decoration)

    OUTPUT: ``None`` or an ``array('i')`` ``r`` of length ``len(fp1)`` with
    ``r[fp1[h]] == fp2[r[h]]`` and ``data2[i][r[h]] == data1[i][h]`` for
    every active ``h``, such that ``classes1[i][h] -> classes2[i][r[h]]``
    on the active ``h`` is a well-defined bijection of values for each
    ``i``, that preserves the pairs ``{2e, 2e+1}``. When
    ``len(fp1) == len(fp2)``, ``r`` is a permutation of ``range(len(fp1))``
    that sends the inactive edges of the first map to those of the second in
    increasing order; otherwise ``r`` is ``-1`` on the half-edges of the
    inactive edges of the first map.

    EXAMPLES::

        sage: from combisurf.isomorphism import isomorphism
        sage: fp1 = [2, 5, 4, 1, 0, 3]  # the faces (0,1,2)(~0,~2,~1)
        sage: fp2 = [4, 3, 0, 5, 2, 1]  # the faces (0,2,1)(~0,~1,~2)
        sage: fp3 = [2, 3, 4, 5, 0, 1]  # the faces (0,1,2)(~0,~1,~2)
        sage: isomorphism(fp1, (), fp2, ())
        array('i', [0, 1, 4, 5, 2, 3])
        sage: isomorphism(fp1, (), fp3, ()) is None
        True
        sage: isomorphism(fp1, ([0, 0, 1, 0, 0, 0],), fp2, ([0, 0, 0, 0, 0, 1],))
        array('i', [3, 2, 5, 4, 1, 0])
        sage: c = [0, 1, 0, 1, 0, 1]
        sage: isomorphism(fp3, (), fp3, (), (c,), ([5, 2, 5, 2, 5, 2],))
        array('i', [0, 1, 2, 3, 4, 5])

    TESTS::

        sage: isomorphism(fp1, ((0,) * 6,), fp2, ())
        Traceback (most recent call last):
        ...
        ValueError: data1 and data2 must have the same number of arrays
        sage: isomorphism(fp1, (), fp2, (), ((0,) * 6,), ())
        Traceback (most recent call last):
        ...
        ValueError: classes1 and classes2 must have the same number of arrays
        sage: isomorphism([2, -1, 0, 3, -1, -1], (), [-1, -1, 4, -1, 2, 5], ())
        array('i', [2, 3, 4, 5, 0, 1])
        sage: isomorphism([-1, -1, 4, -1, 2, 5, -1, -1], (), [2, -1, 0, 3], ())
        array('i', [-1, -1, 0, 1, 2, 3, -1, -1])
        sage: isomorphism([2, 5, 4, 1, 0, 3, 6, -1], (), [0, -1, 4, 7, 6, 3, 2, 5], ())
        array('i', [2, 3, 4, 5, 6, 7, 0, 1])
        sage: isomorphism([2, 5, 4, 1, 0, 3], (), [2, 5, 4, 1, 0, 3, 6, -1], ()) is None
        True
    """
    _check_input(fp1, data1, classes1)
    _check_input(fp2, data2, classes2)
    if len(data1) != len(data2):
        raise ValueError("data1 and data2 must have the same number of arrays")
    if len(classes1) != len(classes2):
        raise ValueError("classes1 and classes2 must have the same number of arrays")
    cdef _Map m1 = _convert(fp1, data1, classes1)
    cdef _Map m2 = _convert(fp2, data2, classes2)
    cdef int n1 = m1.n, n2 = m2.n
    cdef array.array key1 = _new(ISO_KEY_SIZE(n1, m1.nd, m1.nc))
    cdef array.array key2 = _new(ISO_KEY_SIZE(n2, m2.nd, m2.nc))
    cdef array.array ra = _new(n1)
    cdef array.array rb = _new(n2)
    cdef int len1, len2, h, e, i
    len1 = _run(m1, key1, ra.data.as_ints, NULL, NULL)
    len2 = _run(m2, key2, rb.data.as_ints, NULL, NULL)
    if len1 == ISO_DISCONNECTED or len2 == ISO_DISCONNECTED:
        return isomorphism_naive.isomorphism(fp1, data1, fp2, data2, classes1, classes2)
    if (m1.ne != m2.ne or len1 != len2
            or memcmp(key1.data.as_ints, key2.data.as_ints, len1 * sizeof(int)) != 0):
        return None
    r = array.array('i', [-1] * n1)
    rbinv = [-1] * (2 * m2.ne)
    for h in range(n2):
        if rb.data.as_ints[h] != -1:
            rbinv[rb.data.as_ints[h]] = h
    for h in range(n1):
        if ra.data.as_ints[h] != -1:
            r[h] = rbinv[ra.data.as_ints[h]]
    if n1 == n2:
        inactive2 = [e for e in range(n2 // 2) if m2.fp.data.as_ints[2 * e] == -1 and m2.fp.data.as_ints[2 * e + 1] == -1]
        i = 0
        for e in range(n1 // 2):
            if r[2 * e] == -1:
                r[2 * e] = 2 * inactive2[i]
                r[2 * e + 1] = 2 * inactive2[i] + 1
                i += 1
    return r


def canonical_key(fp, data=(), classes=(), relabelling=False, automorphism_count=False):
    r"""
    Return the canonical form of the map with face permutation ``fp``
    decorated by ``data`` and ``classes``, flattened as bytes, and on
    request its canonical labelling and its number of automorphisms.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``data`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept exactly (default: no decoration); the values on the
      active half-edges are integers in the range of a C ``int``

    - ``classes`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept up to a relabelling of their values (default: no
      decoration)

    - ``relabelling`` -- boolean (default ``False``); whether to return the
      canonical labelling as well

    - ``automorphism_count`` -- boolean (default ``False``); whether to
      return the number of automorphisms as well

    OUTPUT: the key, a ``bytes`` object: the entries of the canonical form
    (see :func:`canonical_form`) as C ``int`` in the machine byte order,
    the face permutation first, then each array of ``data``, then each
    array of ``classes``. Two maps with the same numbers of arrays of
    ``data`` and of ``classes`` have equal keys if and only if their
    canonical forms are equal, that is, if and only if they are isomorphic.
    With ``relabelling=True`` or ``automorphism_count=True``, a tuple
    ``(key, relabelling)``, ``(key, count)`` or ``(key, relabelling,
    count)``, where the relabelling is the one of
    :func:`canonical_labelling` and ``count`` the number of automorphisms
    (the length of the list of :func:`automorphisms`).

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.isomorphism import canonical_key, canonical_form
        sage: fp = [2, 5, 4, 1, 0, 3]
        sage: key = canonical_key(fp)
        sage: array('i', key)
        array('i', [2, 5, 4, 1, 0, 3])
        sage: canonical_key([4, 3, 0, 5, 2, 1]) == key
        True
        sage: canonical_key([2, 3, 4, 5, 0, 1]) == key
        False
        sage: canonical_key(fp, relabelling=True, automorphism_count=True)
        (b'...', array('i', [0, 1, 2, 3, 4, 5]), 6)

    With decorations::

        sage: key, r = canonical_key(fp, ([0, 0, 1, 0, 0, 0],), relabelling=True)
        sage: array('i', key)
        array('i', [2, 5, 4, 1, 0, 3, 0, 0, 0, 0, 0, 1])
        sage: r
        array('i', [1, 0, 5, 4, 3, 2])
        sage: key, count = canonical_key([2, 3, 4, 5, 0, 1], classes=([7, 3, 7, 3, 7, 3],), automorphism_count=True)
        sage: array('i', key), count
        (array('i', [2, 3, 4, 5, 0, 1, 0, 1, 0, 1, 0, 1]), 6)

    The key is the canonical form flattened::

        sage: fp = [4, -1, -1, -1, 5, 0]
        sage: canonical_form(fp, ([1, 0, 0, 0, 2, 3],))
        ((1, 2, 0, -1), (2, 3, 1))
        sage: array('i', canonical_key(fp, ([1, 0, 0, 0, 2, 3],)))
        array('i', [1, 2, 0, -1, 2, 3, 1])

    TESTS:

    The empty map::

        sage: canonical_key([], relabelling=True, automorphism_count=True)
        (b'', array('i'), 1)
        sage: canonical_key([-1, -1], relabelling=True, automorphism_count=True)
        (b'', array('i', [0, 1]), 1)

    A disconnected map, handed to :mod:`combisurf.isomorphism_naive`: the
    two edges ``(0)(~0)(1)(~1)``, with their classes::

        sage: key, r, count = canonical_key([0, 1, 2, 3], classes=([0, 1, 1, 2],), relabelling=True, automorphism_count=True)
        sage: array('i', key), r, count
        (array('i', [0, 1, 2, 3, 0, 1, 0, 2]), array('i', [1, 0, 2, 3]), 2)

    A value of ``data`` outside the range of a C ``int``, and one on an
    inactive half-edge, which is ignored::

        sage: canonical_key([1, 0], ([2**40, 0],))
        Traceback (most recent call last):
        ...
        OverflowError: ...
        sage: array('i', canonical_key([0, -1], ([3, None],)))
        array('i', [0, -1, 3])
    """
    cdef _Map m = _prepare(fp, data, classes)
    cdef array.array key = _new(ISO_KEY_SIZE(m.n, m.nd, m.nc))
    cdef array.array r = _new(m.n)
    cdef int count, length
    length = _run(m, key, r.data.as_ints, NULL, &count)
    if length == ISO_DISCONNECTED:
        form = isomorphism_naive.canonical_form(fp, data, classes)
        k = array.array('i', [operator.index(x) for t in form for x in t]).tobytes()
        out = [k]
        if relabelling:
            out.append(isomorphism_naive.canonical_labelling(fp, data, classes))
        if automorphism_count:
            out.append(len(isomorphism_naive.automorphisms(fp, data, classes)))
    else:
        out = [(<char *> key.data.as_ints)[:length * sizeof(int)]]
        if relabelling:
            out.append(_complete(r, m.n, m.ne))
        if automorphism_count:
            out.append(count if count else 1)
    if len(out) == 1:
        return out[0]
    return tuple(out)
