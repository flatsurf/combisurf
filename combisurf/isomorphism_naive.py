r"""
Isomorphisms, automorphisms and canonical labels of maps, reference implementation

This module computes canonical labellings, canonical forms, automorphism
groups and isomorphisms of connected maps given by their face permutation,
possibly decorated by values on half-edges. It works on plain arrays and is
written for readability: every active half-edge is tried as a root, which
costs `O(n^2)` for `n` half-edges, with no pruning. It is the reference
against which a faster implementation is to be tested; the methods of
:class:`~combisurf.oriented_map.OrientedMap` (such as
:meth:`~combisurf.oriented_map.OrientedMap.canonical_labelling`) call it.

A map is given by its face permutation ``fp``, a list or ``array('i')`` of
length ``n`` with ``-1`` on the inactive half-edges. The twin of the
half-edge ``h`` is ``h ^ 1``; a folded edge has one active half-edge whose
twin is inactive. An edge may be entirely inactive (a gap in the labels).

Decorations come in two tuples of arrays of length ``n`` indexed by
half-edges (the values on inactive half-edges are ignored):

- ``data``: an isomorphism must preserve each array exactly;

- ``classes``: an isomorphism must preserve each array up to a relabelling
  of its values, that is, an isomorphism ``r`` from a map with the array
  ``c`` to a map with the array ``c'`` must come with a bijection ``s`` of
  values such that ``c'[r[h]] == s(c[h])`` for every active half-edge
  ``h``. Only the partition of the half-edges into classes of equal values
  matters.

The canonical labelling is obtained from a breadth-first walk on
half-edges started at a root: the face of the half-edge taken from the
queue is walked, and each half-edge met on an edge not labelled yet gets
the next even label, its twin the next odd label. The *key* of a root is
the face permutation in the new labels, followed by each array of
``data`` read in the order of the new labels, followed by each array of
``classes`` read in the same order and renumbered ``0, 1, 2, ...`` by first
appearance in that order. The renumbering makes the key independent of the
names of the values, and costs `O(n)` per root like the rest of the key.
The canonical labelling is the one of a root of minimal key.

EXAMPLES::

    sage: from combisurf.isomorphism_naive import canonical_labelling, automorphisms, isomorphism
    sage: fp = [2, 5, 4, 1, 0, 3]  # the faces (0,1,2)(~0,~2,~1)
    sage: canonical_labelling(fp)
    array('i', [0, 1, 2, 3, 4, 5])
    sage: len(automorphisms(fp))
    6
    sage: isomorphism(fp, (), [4, 3, 0, 5, 2, 1], ())
    array('i', [0, 1, 4, 5, 2, 3])

A decoration of the two faces of ``(0,1,2)(~0,~1,~2)``, kept exactly or up
to relabelling::

    sage: fp = [2, 3, 4, 5, 0, 1]  # the faces (0,1,2)(~0,~1,~2)
    sage: c = [0, 1, 0, 1, 0, 1]
    sage: len(automorphisms(fp, (c,)))
    3
    sage: len(automorphisms(fp, classes=(c,)))
    6
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

from array import array
from collections import deque


def _check_input(fp, data, classes=()):
    r"""
    Return ``(n, ne)``, the length of ``fp`` and its number of active edges,
    after checking that ``fp`` has even length and each array of ``data``
    and of ``classes`` has length ``n``.

    TESTS::

        sage: from combisurf.isomorphism_naive import _check_input
        sage: _check_input([0, -1, 2], ())
        Traceback (most recent call last):
        ...
        ValueError: fp must have even length
        sage: _check_input([0, -1], ([0],))
        Traceback (most recent call last):
        ...
        ValueError: data arrays must have length 2
        sage: _check_input([0, -1], (), [[0, 0]])
        Traceback (most recent call last):
        ...
        TypeError: classes must be a tuple
        sage: _check_input([0, -1], (), ([0],))
        Traceback (most recent call last):
        ...
        ValueError: classes arrays must have length 2
        sage: _check_input([-1, -1, 2, -1], ())
        (4, 1)
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
    ne = sum(1 for e in range(n // 2) if fp[2 * e] != -1 or fp[2 * e + 1] != -1)
    return n, ne


def relabelling_from(fp, root):
    r"""
    Return ``(r, k)`` where ``r`` is the relabelling of the half-edges of
    ``fp`` obtained by the breadth-first walk from ``root`` and ``k`` the
    number of labels used.

    The half-edges are taken from a queue that initially holds ``root``.
    For a half-edge ``h`` taken from the queue, its face ``h, fp[h],
    fp[fp[h]], ...`` is walked; each half-edge ``x`` met whose edge has no
    label yet gets the next even label ``2i``, its twin ``x ^ 1`` the label
    ``2i + 1``, and the twin is pushed on the queue when it is active. The
    active half-edge of a folded edge therefore gets an even label and its
    inactive twin the odd one.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``root`` -- an active half-edge

    OUTPUT: a pair ``(r, k)`` with ``r`` an ``array('i')`` of length
    ``len(fp)`` (``-1`` on the half-edges of the edges that are not
    reached) and ``k`` the number of labels used. The map is connected if
    and only if ``k`` is twice its number of active edges.

    EXAMPLES::

        sage: from combisurf.isomorphism_naive import relabelling_from
        sage: fp = [2, 5, 4, 1, 0, 3]  # the faces (0,1,2)(~0,~2,~1)
        sage: relabelling_from(fp, 0)
        (array('i', [0, 1, 2, 3, 4, 5]), 6)
        sage: relabelling_from(fp, 3)
        (array('i', [3, 2, 1, 0, 5, 4]), 6)

    A folded edge and a label gap (the edge ``1`` is inactive)::

        sage: fp = [4, -1, -1, -1, 5, 0]  # the face (0,2,~2)
        sage: relabelling_from(fp, 5)
        (array('i', [2, 3, -1, -1, 1, 0]), 4)

    TESTS:

    A disconnected map::

        sage: relabelling_from([0, 1, 2, 3], 0)
        (array('i', [0, 1, -1, -1]), 2)
    """
    n = len(fp)
    if not (0 <= root < n) or fp[root] == -1:
        raise ValueError("root must be an active half-edge")
    r = array('i', [-1] * n)
    k = 0
    queue = deque([root])
    while queue:
        h = queue.popleft()
        x = h
        while True:
            if r[x] == -1:
                r[x] = k
                r[x ^ 1] = k + 1
                k += 2
                if fp[x ^ 1] != -1:
                    queue.append(x ^ 1)
            x = fp[x]
            if x == h:
                break
    return r, k


def _key(fp, data, r, k, classes=()):
    r"""
    Return the key of the relabelling ``r`` using ``k`` labels: the face
    permutation in the new labels, then each array of ``data`` read in the
    increasing order of the new labels of the active half-edges, then each
    array of ``classes`` read in the same order and renumbered ``0, 1, 2,
    ...`` by first appearance.

    EXAMPLES::

        sage: from combisurf.isomorphism_naive import relabelling_from, _key
        sage: fp = [4, -1, -1, -1, 5, 0]
        sage: r, k = relabelling_from(fp, 5)
        sage: _key(fp, (['a', 'b', 'c', 'd', 'e', 'f'],), r, k)
        ((2, 0, 1, -1), ('f', 'e', 'a'))

    The half-edges ``5, 4, 0`` get the new labels ``0, 1, 2``; the values
    ``7, 3, 7`` of a class array on them become ``0, 1, 0``::

        sage: _key(fp, (), r, k, ([7, 9, 9, 9, 3, 7],))
        ((2, 0, 1, -1), (0, 1, 0))
    """
    fp_new = [-1] * k
    rinv = [-1] * k
    for h, x in enumerate(fp):
        if x != -1:
            fp_new[r[h]] = r[x]
            rinv[r[h]] = h
    key = [tuple(fp_new)]
    for d in data:
        key.append(tuple(d[h] for h in rinv if h != -1))
    for c in classes:
        names = {}
        renumbered = []
        for h in rinv:
            if h != -1:
                v = c[h]
                i = names.get(v)
                if i is None:
                    i = names[v] = len(names)
                renumbered.append(i)
        key.append(tuple(renumbered))
    return tuple(key)


def _minimal_roots(fp, data, classes=()):
    r"""
    Return ``(key, roots)``: the minimal key over the active roots and the
    list of the roots that reach it, each as a pair ``(root, r)`` with
    ``r`` its relabelling, in increasing order of roots.

    Raise ``NotImplementedError`` if ``fp`` is not connected. For a map
    without active half-edge the list of roots is empty.
    """
    n, ne = _check_input(fp, data, classes)
    best = None
    roots = []
    for root in range(n):
        if fp[root] == -1:
            continue
        r, k = relabelling_from(fp, root)
        if k != 2 * ne:
            raise NotImplementedError("disconnected maps are not supported")
        key = _key(fp, data, r, k, classes)
        if best is None or key < best:
            best = key
            roots = [(root, r)]
        elif key == best:
            roots.append((root, r))
    if best is None:
        best = ((),) + tuple(() for _ in data) + tuple(() for _ in classes)
    return best, roots


def _complete(r, n, ne):
    r"""
    Complete in place the relabelling ``r`` (``-1`` on the half-edges of the
    inactive edges) by the labels ``2*ne, 2*ne + 1, ...`` given to the
    inactive edges in increasing order.
    """
    k = 2 * ne
    for e in range(n // 2):
        if r[2 * e] == -1:
            r[2 * e] = k
            r[2 * e + 1] = k + 1
            k += 2
    return r


def canonical_labelling(fp, data=(), classes=()):
    r"""
    Return the canonical labelling of the connected map with face
    permutation ``fp`` decorated by ``data`` and ``classes``.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``data`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept exactly (default: no decoration)

    - ``classes`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept up to a relabelling of their values (default: no
      decoration)

    OUTPUT: an ``array('i')`` ``r`` that maps each half-edge to its
    canonical label. It is a permutation of ``range(len(fp))`` that
    preserves the pairs ``{2e, 2e+1}``. The active edges get the labels
    ``0, ..., 2*ne - 1`` (``ne`` the number of active edges) and the
    inactive edges the following ones, in increasing order.

    Two decorated maps are isomorphic if and only if they are equal once
    relabelled by their canonical labellings.

    EXAMPLES::

        sage: from combisurf.isomorphism_naive import canonical_labelling
        sage: canonical_labelling([2, 5, 4, 1, 0, 3])
        array('i', [0, 1, 2, 3, 4, 5])
        sage: canonical_labelling([2, 5, 4, 1, 0, 3], ([0, 0, 1, 0, 0, 0],))
        array('i', [1, 0, 5, 4, 3, 2])

    The same labelling when the marked half-edge is given as a class::

        sage: canonical_labelling([2, 5, 4, 1, 0, 3], classes=([5, 5, 2, 5, 5, 5],))
        array('i', [1, 0, 5, 4, 3, 2])

    TESTS:

    A folded edge (the faces ``(0,1)(~1)``)::

        sage: canonical_labelling([2, -1, 0, 3])
        array('i', [2, 3, 1, 0])

    A label gap (the edge ``1`` is inactive)::

        sage: canonical_labelling([4, -1, -1, -1, 5, 0])
        array('i', [2, 3, 4, 5, 0, 1])

    A single half-edge and the empty map::

        sage: canonical_labelling([0, -1])
        array('i', [0, 1])
        sage: canonical_labelling([])
        array('i')

    A disconnected map::

        sage: canonical_labelling([0, 1, 2, 3])
        Traceback (most recent call last):
        ...
        NotImplementedError: disconnected maps are not supported
    """
    n, ne = _check_input(fp, data, classes)
    _, roots = _minimal_roots(fp, data, classes)
    r = array('i', roots[0][1]) if roots else array('i', [-1] * n)
    return _complete(r, n, ne)


def canonical_form(fp, data=(), classes=()):
    r"""
    Return the canonical form of the connected map with face permutation
    ``fp`` decorated by ``data`` and ``classes``.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``data`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept exactly (default: no decoration)

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

        sage: from combisurf.isomorphism_naive import canonical_form
        sage: canonical_form([2, 5, 4, 1, 0, 3])
        ((2, 5, 4, 1, 0, 3),)
        sage: canonical_form([4, 3, 0, 5, 2, 1])
        ((2, 5, 4, 1, 0, 3),)
        sage: canonical_form([2, 3, 4, 5, 0, 1])
        ((2, 3, 4, 5, 0, 1),)
        sage: canonical_form([2, 5, 4, 1, 0, 3], ([0, 0, 1, 0, 0, 0],))
        ((2, 5, 4, 1, 0, 3), (0, 0, 0, 0, 0, 1))

    A class array is renumbered, so renaming its values leaves the canonical
    form unchanged::

        sage: canonical_form([2, 3, 4, 5, 0, 1], classes=([0, 1, 0, 1, 0, 1],))
        ((2, 3, 4, 5, 0, 1), (0, 1, 0, 1, 0, 1))
        sage: canonical_form([2, 3, 4, 5, 0, 1], classes=([7, 3, 7, 3, 7, 3],))
        ((2, 3, 4, 5, 0, 1), (0, 1, 0, 1, 0, 1))
        sage: canonical_form([2, 3, 4, 5, 0, 1], ([1, 0, 1, 0, 1, 0],), ([7, 3, 7, 3, 7, 3],))
        ((2, 3, 4, 5, 0, 1), (0, 1, 0, 1, 0, 1), (0, 1, 0, 1, 0, 1))

    TESTS:

    A folded edge, a label gap and the empty map::

        sage: canonical_form([2, -1, 0, 3])
        ((0, 2, 1, -1),)
        sage: canonical_form([4, -1, -1, -1, 5, 0])
        ((1, 2, 0, -1),)
        sage: canonical_form([], ([],))
        ((), ())
        sage: canonical_form([], ([],), ([], []))
        ((), (), (), ())

    A disconnected map::

        sage: canonical_form([0, 1, 2, 3])
        Traceback (most recent call last):
        ...
        NotImplementedError: disconnected maps are not supported
    """
    key, _ = _minimal_roots(fp, data, classes)
    return key


def automorphisms(fp, data=(), classes=()):
    r"""
    Return the automorphisms of the connected map with face permutation
    ``fp`` decorated by ``data`` and ``classes``.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``data`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept exactly (default: no decoration)

    - ``classes`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept up to a relabelling of their values (default: no
      decoration)

    OUTPUT: the list of the automorphisms, each an ``array('i')``
    permutation ``g`` of ``range(len(fp))`` with ``g[fp[h]] == fp[g[h]]``
    and ``d[g[h]] == d[h]`` for every active ``h`` and every ``d`` of
    ``data``, such that for every ``c`` of ``classes`` the map ``c[h] ->
    c[g[h]]`` on the active ``h`` is a well-defined bijection of values,
    that preserves the pairs ``{2e, 2e+1}`` and is the identity
    on the inactive edges. The identity comes first. As the map is
    connected, an automorphism is determined by the image of one active
    half-edge, and there are at most as many automorphisms as active
    half-edges.

    EXAMPLES::

        sage: from combisurf.isomorphism_naive import automorphisms
        sage: A = automorphisms([2, 5, 4, 1, 0, 3])
        sage: len(A)
        6
        sage: A[:2]
        [array('i', [0, 1, 2, 3, 4, 5]), array('i', [1, 0, 5, 4, 3, 2])]

    A decoration that marks one half-edge kills every automorphism::

        sage: automorphisms([2, 5, 4, 1, 0, 3], ([0, 0, 1, 0, 0, 0],))
        [array('i', [0, 1, 2, 3, 4, 5])]

    The two faces of ``(0,1,2)(~0,~1,~2)`` marked by different values: kept
    exactly, only the rotations remain; up to relabelling, the automorphisms
    that swap the two faces (and the two values) remain as well::

        sage: fp = [2, 3, 4, 5, 0, 1]
        sage: c = [0, 1, 0, 1, 0, 1]
        sage: len(automorphisms(fp)), len(automorphisms(fp, (c,))), len(automorphisms(fp, classes=(c,)))
        (6, 3, 6)

    TESTS:

    A folded edge, a label gap, a single half-edge and the empty map::

        sage: automorphisms([2, -1, 0, 3])
        [array('i', [0, 1, 2, 3])]
        sage: automorphisms([4, -1, -1, -1, 5, 0])
        [array('i', [0, 1, 2, 3, 4, 5])]
        sage: automorphisms([0, -1])
        [array('i', [0, 1])]
        sage: automorphisms([])
        [array('i')]

    A disconnected map::

        sage: automorphisms([0, 1, 2, 3])
        Traceback (most recent call last):
        ...
        NotImplementedError: disconnected maps are not supported
    """
    n, ne = _check_input(fp, data, classes)
    _, roots = _minimal_roots(fp, data, classes)
    if not roots:
        return [array('i', range(n))]
    ra = roots[0][1]
    result = []
    for _, rb in roots:
        rbinv = [-1] * (2 * ne)
        for h, x in enumerate(rb):
            if x != -1:
                rbinv[x] = h
        g = array('i', range(n))
        for h, x in enumerate(ra):
            if x != -1:
                g[h] = rbinv[x]
        result.append(g)
    return result


def _group_closure(gens, n):
    r"""
    Return the set of the elements (as tuples) of the group generated by the
    permutations ``gens`` of ``range(n)``.
    """
    identity = tuple(range(n))
    elements = {identity}
    queue = deque([identity])
    while queue:
        x = queue.popleft()
        for g in gens:
            y = tuple(g[i] for i in x)
            if y not in elements:
                elements.add(y)
                queue.append(y)
    return elements


def automorphism_generators(fp, data=(), classes=()):
    r"""
    Return generators of the automorphism group of the connected map with
    face permutation ``fp`` decorated by ``data`` and ``classes``.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``data`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept exactly (default: no decoration)

    - ``classes`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept up to a relabelling of their values (default: no
      decoration)

    OUTPUT: a list of automorphisms (as in :func:`automorphisms`) that
    generate the group; empty when the group is trivial. The automorphisms
    are scanned in the order of :func:`automorphisms` and one is kept when
    it is not in the group generated by those already kept.

    EXAMPLES::

        sage: from combisurf.isomorphism_naive import automorphism_generators
        sage: automorphism_generators([2, 5, 4, 1, 0, 3])
        [array('i', [1, 0, 5, 4, 3, 2]), array('i', [2, 3, 4, 5, 0, 1])]
        sage: automorphism_generators([2, 3, 4, 5, 0, 1])
        [array('i', [1, 0, 3, 2, 5, 4]), array('i', [2, 3, 4, 5, 0, 1])]

    The faces of the last map marked by different values, kept exactly or up
    to relabelling::

        sage: c = [0, 1, 0, 1, 0, 1]
        sage: automorphism_generators([2, 3, 4, 5, 0, 1], (c,))
        [array('i', [2, 3, 4, 5, 0, 1])]
        sage: automorphism_generators([2, 3, 4, 5, 0, 1], classes=(c,))
        [array('i', [1, 0, 3, 2, 5, 4]), array('i', [2, 3, 4, 5, 0, 1])]

    TESTS::

        sage: automorphism_generators([2, -1, 0, 3])
        []
        sage: automorphism_generators([0, -1])
        []
    """
    n = len(fp)
    gens = []
    group = {tuple(range(n))}
    for g in automorphisms(fp, data, classes):
        if tuple(g) not in group:
            gens.append(g)
            group = _group_closure(gens, n)
    return gens


def isomorphism(fp1, data1, fp2, data2, classes1=(), classes2=()):
    r"""
    Return an isomorphism from the first decorated connected map to the
    second, or ``None`` if they are not isomorphic.

    INPUT:

    - ``fp1``, ``fp2`` -- face permutations (lists or ``array('i')``,
      ``-1`` on inactive half-edges)

    - ``data1``, ``data2`` -- tuples with the same number of arrays indexed
      by half-edges, of lengths ``len(fp1)`` and ``len(fp2)``, kept exactly

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

        sage: from combisurf.isomorphism_naive import isomorphism
        sage: fp1 = [2, 5, 4, 1, 0, 3]  # the faces (0,1,2)(~0,~2,~1)
        sage: fp2 = [4, 3, 0, 5, 2, 1]  # the faces (0,2,1)(~0,~1,~2)
        sage: fp3 = [2, 3, 4, 5, 0, 1]  # the faces (0,1,2)(~0,~1,~2)
        sage: r = isomorphism(fp1, (), fp2, ())
        sage: r
        array('i', [0, 1, 4, 5, 2, 3])
        sage: all(r[fp1[h]] == fp2[r[h]] for h in range(6))
        True
        sage: isomorphism(fp1, (), fp3, ()) is None
        True

    With decorations::

        sage: isomorphism(fp1, ([0, 0, 1, 0, 0, 0],), fp2, ([0, 0, 0, 0, 1, 0],))
        array('i', [0, 1, 4, 5, 2, 3])
        sage: isomorphism(fp1, ([0, 0, 1, 0, 0, 0],), fp2, ([0, 0, 0, 0, 0, 1],))
        array('i', [3, 2, 5, 4, 1, 0])
        sage: isomorphism(fp1, ([0, 0, 1, 0, 0, 0],), fp2, ([0, 0, 0, 0, 0, 2],)) is None
        True

    With classes, the values need not match, only the partitions they
    define::

        sage: c = [0, 1, 0, 1, 0, 1]
        sage: isomorphism(fp3, (c,), fp3, ([5, 2, 5, 2, 5, 2],)) is None
        True
        sage: isomorphism(fp3, (), fp3, (), (c,), ([5, 2, 5, 2, 5, 2],))
        array('i', [0, 1, 2, 3, 4, 5])
        sage: isomorphism(fp3, (), fp3, (), (c,), ([0, 0, 1, 1, 0, 1],)) is None
        True

    TESTS:

    Different numbers of arrays::

        sage: isomorphism(fp1, ((0,) * 6,), fp2, ())
        Traceback (most recent call last):
        ...
        ValueError: data1 and data2 must have the same number of arrays
        sage: isomorphism(fp1, (), fp2, (), ((0,) * 6,), ())
        Traceback (most recent call last):
        ...
        ValueError: classes1 and classes2 must have the same number of arrays

    A folded edge and label gaps, with arrays of the same length or not::

        sage: isomorphism([2, -1, 0, 3], (), [-1, -1, 4, -1, 2, 5], ())
        array('i', [2, 3, 4, 5])
        sage: isomorphism([2, -1, 0, 3, -1, -1], (), [-1, -1, 4, -1, 2, 5], ())
        array('i', [2, 3, 4, 5, 0, 1])
        sage: isomorphism([-1, -1, 4, -1, 2, 5, -1, -1], (), [2, -1, 0, 3], ())
        array('i', [-1, -1, 0, 1, 2, 3, -1, -1])

    A disconnected map::

        sage: isomorphism([0, 1, 2, 3], (), [0, 1, 2, 3], ())
        Traceback (most recent call last):
        ...
        NotImplementedError: disconnected maps are not supported
    """
    n1, ne1 = _check_input(fp1, data1, classes1)
    n2, ne2 = _check_input(fp2, data2, classes2)
    if len(data1) != len(data2):
        raise ValueError("data1 and data2 must have the same number of arrays")
    if len(classes1) != len(classes2):
        raise ValueError("classes1 and classes2 must have the same number of arrays")
    key1, roots1 = _minimal_roots(fp1, data1, classes1)
    key2, roots2 = _minimal_roots(fp2, data2, classes2)
    if key1 != key2:
        return None
    r = array('i', [-1] * n1)
    if roots1:
        r1 = roots1[0][1]
        r2 = roots2[0][1]
        r2inv = [-1] * (2 * ne2)
        for h, x in enumerate(r2):
            if x != -1:
                r2inv[x] = h
        for h, x in enumerate(r1):
            if x != -1:
                r[h] = r2inv[x]
    if n1 == n2:
        inactive2 = [e for e in range(n2 // 2) if fp2[2 * e] == -1 and fp2[2 * e + 1] == -1]
        i = 0
        for e in range(n1 // 2):
            if r[2 * e] == -1:
                r[2 * e] = 2 * inactive2[i]
                r[2 * e + 1] = 2 * inactive2[i] + 1
                i += 1
    return r
