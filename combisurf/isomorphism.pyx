# distutils: sources = combisurf/src/isomorphism.c
# distutils: include_dirs = combisurf/src
r"""
Isomorphisms, automorphisms and canonical labels of maps

This module computes canonical labellings, canonical forms, automorphism
groups and isomorphisms of maps given by their face permutation, possibly
decorated by integer values on half-edges. It computes what
:mod:`combisurf.isomorphism_naive` computes, which is its reference: the
functions of the same names take the same arguments and return the same
values (the same canonical labelling, the same canonical form, the
automorphisms in the same order, the same isomorphism), except
:func:`automorphism_generators` on a disconnected map, whose generators
come from the structure of the group and differ from those of the
reference while generating the same group. They raise the same
exceptions on malformed input, with one addition: a face permutation that
is not a permutation of its active half-edges raises a ``ValueError``,
since the C code would read out of bounds on it (the reference loops
forever or raises an ``IndexError`` there). The conventions (face permutations with
``-1`` on inactive half-edges, the decorations ``data`` kept exactly and
``classes`` kept up to a relabelling of their values, the breadth-first
walk and the key of a root) are described in
:mod:`combisurf.isomorphism_naive`. The methods of
:class:`~combisurf.oriented_map.OrientedMap` (such as
:meth:`~combisurf.oriented_map.OrientedMap.canonical_labelling`) call this
module.

On a connected map, every active half-edge of the class of the rule (see
below) is tried as a root: the relabelling from it and its key are
computed in C (``combisurf/src/isomorphism.c``) and the minimal key is
kept, which costs `O(n^2)` for `n` half-edges. A disconnected map is
searched in C as well, component after component, by the search of
:mod:`combisurf.isomorphism_naive` over the sequences of blocks. Without
``classes``, the automorphism group of a disconnected map
is a product of wreath products `\mathrm{Aut}(C) \wr S_k` over the kinds
of `k` isomorphic components `C`: :func:`automorphism_generators` and the
count of :func:`canonical_key` come from this structure, component by
component, without listing the automorphisms (their number is at least
`k!`). With ``classes`` the group is the subgroup of this product that
preserves the partitions of the half-edges by the values of ``classes``:
its generators and its order come from a backtrack search of this
subgroup (the subgroup search of Sims, one base point per component),
exponential in the number of components in the worst case. The search
of the canonical key itself follows every root tied in a component when
there are ``classes``, and returns from a leaf equal to the best one to
the first level where their roots differ, since the subtree left is then
the image of one already searched by an automorphism: `k` tied copies of
a component cost about `k^2` paths of `k` components. The buffers of the C
code are counted in C ``int``: a map
too large for them (about `2^{31} / (14 + 2 nd + 3.5 nc)` half-edges with
``nd`` arrays of ``data`` and ``nc`` of ``classes``) raises an
``OverflowError``.

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

The public functions take the keyword ``root_rule``, one of the integers
:data:`ROOTS_ALL`, :data:`ROOTS_HALF_EDGE_MAX`,
:data:`ROOTS_HALF_EDGE_FEWEST`, :data:`ROOTS_EDGE_MAX` and
:data:`ROOTS_EDGE_FEWEST`, equal to the constants of the same names of
:mod:`combisurf.isomorphism_naive`, which defines the rules: the
half-edges tried as roots. The default is ``ROOTS_HALF_EDGE_FEWEST |
ROOTS_DATA`` (the integer ``10``): decorations narrow the roots tried, and
on an undecorated map it is :data:`ROOTS_HALF_EDGE_FEWEST`, which tries
fewer roots: against
:data:`ROOTS_ALL`, one run each, a GL(2, Z)-orbit exploration of 50
squares took 1.44 s instead of 2.81 s and the canonical labelling of an
origami of 1000 squares 307 us instead of 446 us, while a random
triangulation of 4092 half-edges took 118 us instead of 91.1 us. The
flag :data:`ROOTS_DATA`, combined by ``|`` with one of the four rules
other than :data:`ROOTS_ALL`, appends the values of ``data`` to the
invariants, as in :mod:`combisurf.isomorphism_naive`, so that fewer roots
are tried on a decorated map; the C code sorts each data array once per
call to rank its values, `O(n \log n)`. Any other value, a string or
``ROOTS_ALL | ROOTS_DATA`` included, raises a
``ValueError``. The C code computes the invariants of the half-edges once
per call and tries only the roots of the class of the rule, on a
disconnected map the class chosen at each step among the half-edges not
labelled yet, as the reference does. The automorphism group does not
depend on the rule: on a disconnected map, :func:`automorphism_generators`
and the count of :func:`canonical_key` come from the structure of the
group above under every rule.

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

Under every rule, the results are those of the reference for that rule,
computed in C::

    sage: from combisurf.isomorphism import canonical_form, ROOTS_ALL, ROOTS_HALF_EDGE_FEWEST, ROOTS_HALF_EDGE_MAX
    sage: fp = [2, 5, 4, 1, 0, 3, 6, -1]
    sage: canonical_form(fp)
    ((0, -1, 4, 7, 6, 3, 2, 5),)
    sage: canonical_form(fp, root_rule=ROOTS_HALF_EDGE_FEWEST) == canonical_form(fp)
    True
    sage: canonical_form(fp, root_rule=ROOTS_ALL) == isomorphism_naive.canonical_form(fp, root_rule=isomorphism_naive.ROOTS_ALL)
    True
    sage: canonical_form(fp, root_rule=ROOTS_HALF_EDGE_MAX)
    ((2, 5, 4, 1, 0, 3, 6, -1),)
    sage: canonical_form(fp, root_rule=ROOTS_HALF_EDGE_MAX) == isomorphism_naive.canonical_form(fp, root_rule=isomorphism_naive.ROOTS_HALF_EDGE_MAX)
    True

With the faces of the triangle decorated by ``0`` and ``1``, the rule
``ROOTS_HALF_EDGE_MAX | ROOTS_DATA`` tries three roots instead of six, as
counted by the reference, and gives its results::

    sage: from combisurf.isomorphism import ROOTS_DATA
    sage: d = [0, 1, 0, 1, 0, 1, 0, 0]
    sage: rule = ROOTS_HALF_EDGE_MAX | ROOTS_DATA
    sage: r = [-1] * 8
    sage: [len(isomorphism_naive._roots(fp, r, isomorphism_naive._root_invariants(fp, x, (d,)), x)) for x in (ROOTS_HALF_EDGE_MAX, rule)]
    [6, 3]
    sage: canonical_form(fp, (d,), root_rule=rule)
    ((2, 5, 4, 1, 0, 3, 6, -1), (1, 0, 1, 0, 1, 0, 0))
    sage: canonical_form(fp, (d,), root_rule=rule) == isomorphism_naive.canonical_form(fp, (d,), root_rule=rule)
    True
    sage: canonical_form(fp, (d,), root_rule=ROOTS_ALL | ROOTS_DATA)
    Traceback (most recent call last):
    ...
    ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got 8
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
from libc.string cimport memcpy, memcmp
from math import factorial
from numbers import Integral


cdef array.array _int_array = array.array('i', [])


cdef enum:
    # the value of ROOTS_ALL as a C constant, passed to the C code without
    # the lookup of a module global
    _ROOTS_ALL = 0


ROOTS_ALL = _ROOTS_ALL
ROOTS_HALF_EDGE_MAX = 1
ROOTS_HALF_EDGE_FEWEST = 2
ROOTS_EDGE_MAX = 3
ROOTS_EDGE_FEWEST = 4
ROOTS_DATA = 8


cdef int _check_root_rule(root_rule) except -1:
    r"""
    Return ``root_rule`` as a C ``int`` after checking that it is one of the
    integers ``ROOTS_ALL``, ..., ``ROOTS_EDGE_FEWEST``, or one of the last
    four combined with ``ROOTS_DATA`` by ``|`` (9 to 12), with the message
    of ``isomorphism_naive._check_root_rule``.

    A Python ``int`` is checked first, without the lookup of
    :class:`numbers.Integral`: :func:`canonical_key` may be called millions
    of times in an exploration.
    """
    if type(root_rule) is int:
        if 0 <= root_rule <= 4 or 9 <= root_rule <= 12:
            return root_rule
    elif (not isinstance(root_rule, bool) and isinstance(root_rule, Integral)
            and (0 <= root_rule <= 4 or 9 <= root_rule <= 12)):
        return int(root_rule)
    raise ValueError("root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, "
                     "ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 "
                     "(the last four | ROOTS_DATA), got {!r}".format(root_rule))


def _check_sizes(n, nd, nc):
    r"""
    Raise an ``OverflowError`` if a buffer of the C code for a map of ``n``
    half-edges decorated by ``nd`` arrays of ``data`` and ``nc`` arrays of
    ``classes`` has more entries than a C ``int`` counts.

    The sizes are those of ``combisurf/src/isomorphism.h``, computed here in
    Python integers.

    TESTS::

        sage: from combisurf.isomorphism import _check_sizes
        sage: _check_sizes(10**6, 2, 1)
        sage: _check_sizes(2**28, 2, 1)
        Traceback (most recent call last):
        ...
        OverflowError: map too large: 268435456 half-edges with 2 arrays of data and 1 of classes
        sage: _check_sizes(2**31, 0, 0)
        Traceback (most recent call last):
        ...
        OverflowError: map too large: 2147483648 half-edges with 0 arrays of data and 0 of classes
    """
    key = n * (1 + nd + nc)
    work = 9 * n + 1 + key
    sequence = key + n // 2
    search = (5 + nc) * n + (5 + nc) * (n // 2 + 1) + nc + sequence + key + 1 + 4 * n + 1
    if max(n, key, work, sequence, search) > INT_MAX:
        raise OverflowError("map too large: {} half-edges with {} arrays of data and {} of classes".format(n, nd, nc))


cdef inline array.array _new(Py_ssize_t size):
    return array.clone(_int_array, size, False)


ctypedef fused _sequence:
    list
    tuple


cdef bint _ints_from(_sequence s, const int *f, int *out, int n, bint edges) except -1:
    r"""
    Write the values of ``s`` as C ``int`` into ``out``, of ``n`` entries,
    and return whether every value converted is a Python ``int`` in the
    range of a C ``int``; on ``False``, ``out`` is partly written.

    Without ``edges``, ``s`` has ``n`` entries and ``out[h]`` is ``s[h]``.
    With ``edges``, ``s`` has ``n // 2`` entries, one per edge, and
    ``out[2*e]`` and ``out[2*e + 1]`` are ``s[e]``. If ``f`` is not ``NULL``,
    the values on the half-edges ``h`` with ``f[h] == -1`` (on the edges
    whose two halves have it, with ``edges``) are neither read nor checked
    and ``out`` holds ``0`` there.

    Only the exact type ``int`` is converted here: the conversion of other
    types (``bool``, the integers of Sage, ``float``, the objects with
    ``__index__``) differs between ``array.array`` and Cython, so they are
    left to the generic path of the caller, which then raises the same
    exception as before or none.
    """
    cdef int h, e
    cdef object v
    try:
        if edges:
            for e in range(n // 2):
                if f != NULL and f[2 * e] == -1 and f[2 * e + 1] == -1:
                    out[2 * e] = out[2 * e + 1] = 0
                    continue
                v = s[e]
                if type(v) is not int:
                    return False
                out[2 * e] = out[2 * e + 1] = v
        else:
            for h in range(n):
                if f != NULL and f[h] == -1:
                    out[h] = 0
                    continue
                v = s[h]
                if type(v) is not int:
                    return False
                out[h] = v
    except OverflowError:
        return False
    return True


cdef array.array _as_int_array(a):
    r"""
    Return ``a`` if it is an array of typecode ``'i'`` and a copy of it as
    such an array otherwise.

    A list of Python ``int`` is copied by a typed loop; any other input, or
    a list holding another value, goes through ``array.array('i', a)``,
    which raises the exceptions of :mod:`array` on invalid values.
    """
    cdef array.array out
    cdef Py_ssize_t size
    if type(a) is array.array and (<array.array> a).ob_descr.typecode == b'i':
        return <array.array> a
    if type(a) is list:
        size = len(<list> a)
        if size <= INT_MAX:
            out = _new(size)
            if _ints_from(<list> a, NULL, out.data.as_ints, <int> size, False):
                return out
    return array.array('i', a)


def _decoration_array(d, array.array fp, bint edges=False):
    r"""
    Return the decoration ``d`` of the map of face permutation ``fp`` as an
    array of typecode ``'i'`` indexed by half-edges, or ``None`` if this
    fast conversion does not apply.

    ``d`` is indexed by half-edges, or by edges with ``edges`` (its value
    is then copied to both halves of each edge). The values on the inactive
    half-edges (``fp[h] == -1``; with ``edges``, on the edges whose two
    halves are inactive) are neither read nor checked and replaced by
    ``0``, except that an array of typecode ``'i'`` indexed by half-edges is
    returned as it is.

    The conversion applies to a list, a tuple or an array of typecode
    ``'i'`` of the right length whose values read are Python ``int`` in the
    range of a C ``int``. On any other input it returns ``None``, and
    :meth:`~combisurf.oriented_map.OrientedMap._decoration_data` checks
    and passes ``d`` as before.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.isomorphism import _decoration_array
        sage: fp = array('i', [2, 7, 6, 1, -1, -1, 8, 3, 0, -1])
        sage: d = list(range(10))
        sage: d[4] = d[5] = d[9] = None
        sage: _decoration_array(d, fp)
        array('i', [0, 1, 2, 3, 0, 0, 6, 7, 8, 0])
        sage: _decoration_array((0r, 1r, 'a', 3r, 4r), fp, edges=True)
        array('i', [0, 0, 1, 1, 0, 0, 3, 3, 4, 4])

    The integers of Sage, and values outside the range of a C ``int``, are
    left to the generic path::

        sage: _decoration_array([0] * 10, fp) is None
        True
        sage: d = list(range(10))
        sage: d[8] = int(2**31)
        sage: _decoration_array(d, fp) is None
        True

    TESTS::

        sage: d[8] = 0.5r
        sage: _decoration_array(d, fp) is None
        True
        sage: _decoration_array(list(range(3)), fp) is None
        True
        sage: _decoration_array(array('l', [0] * 10), fp) is None
        True
        sage: a = array('i', [1] * 10)
        sage: _decoration_array(a, fp) is a
        True
        sage: _decoration_array(array('i', [5, 6, 7, 8, 9]), fp, edges=True)
        array('i', [5, 5, 6, 6, 0, 0, 8, 8, 9, 9])
    """
    cdef int n = len(fp)
    cdef int m = n // 2 if edges else n
    cdef int e
    cdef int *f = fp.data.as_ints
    cdef int *a
    cdef array.array out
    if fp.ob_descr.typecode != b'i':
        return None
    if type(d) is array.array and (<array.array> d).ob_descr.typecode == b'i':
        if len(<array.array> d) != m:
            return None
        if not edges:
            return d
        out = _new(n)
        a = (<array.array> d).data.as_ints
        for e in range(m):
            if f[2 * e] == -1 and f[2 * e + 1] == -1:
                out.data.as_ints[2 * e] = out.data.as_ints[2 * e + 1] = 0
            else:
                out.data.as_ints[2 * e] = out.data.as_ints[2 * e + 1] = a[e]
        return out
    out = _new(n)
    if type(d) is list:
        if len(<list> d) == m and _ints_from(<list> d, f, out.data.as_ints, n, edges):
            return out
    elif type(d) is tuple:
        if len(<tuple> d) == m and _ints_from(<tuple> d, f, out.data.as_ints, n, edges):
            return out
    return None


def _constant_on_cycles(array.array d, array.array p):
    r"""
    Return whether the array ``d`` of typecode ``'i'`` is constant on the
    cycles of the partial permutation ``p``, comparing ``d[h]`` and
    ``d[p[h]]`` for every half-edge ``h`` with ``p[h] != -1``.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.isomorphism import _constant_on_cycles
        sage: p = array('i', [1, 0, -1, 3])
        sage: _constant_on_cycles(array('i', [5, 5, 7, 6]), p)
        True
        sage: _constant_on_cycles(array('i', [5, 4, 7, 6]), p)
        False

    TESTS::

        sage: _constant_on_cycles(array('i', [5, 5, 7]), p)
        Traceback (most recent call last):
        ...
        ValueError: d and p must have the same length
        sage: _constant_on_cycles(array('i', [5, 5, 7, 6]), array('i', [1, 0, -1, 4]))
        Traceback (most recent call last):
        ...
        ValueError: p must have values in range(4) or -1
    """
    cdef int n = len(p)
    cdef int h, x
    cdef int *a = d.data.as_ints
    cdef int *q = p.data.as_ints
    if d.ob_descr.typecode != b'i' or p.ob_descr.typecode != b'i':
        raise TypeError("d and p must be arrays of typecode 'i'")
    if len(d) != n:
        raise ValueError("d and p must have the same length")
    for h in range(n):
        x = q[h]
        if x == -1:
            continue
        if x < 0 or x >= n:
            raise ValueError("p must have values in range({}) or -1".format(n))
        if a[x] != a[h]:
            return False
    return True


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
    half-edge is outside the range of a C ``int`` or if the map is too
    large for the buffers of the C code.
    """
    cdef _Map m = _Map.__new__(_Map)
    cdef int n, h, x, j, off, i
    cdef int *f
    cdef int *dec
    cdef int *a
    cdef array.array seen
    cdef Py_ssize_t size = len(fp)
    # every buffer has at most 16 + 4 * (nd + nc) ints per half-edge (plus
    # a constant): below that bound the exact check in Python integers is
    # skipped, which keeps it off the path of small maps
    if size * (16 + 4 * (len(data) + len(classes))) + 64 > INT_MAX:
        _check_sizes(size, len(data), len(classes))
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
        elif not (type(d) is list and len(<list> d) == n and _ints_from(<list> d, f, dec + off, n, False)):
            # a list of Python int is converted by the typed loop above;
            # anything else by the generic conversion of Cython, value by
            # value, which raises on invalid values
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


cdef int _run(_Map m, int rule, array.array key, int *relabelling, int *roots, int *nroots) except? -2:
    r"""
    Run ``iso_canonical`` on ``m`` under the rule ``rule`` (checked by
    ``_check_root_rule``): write the minimal key into ``key`` (of
    ``ISO_KEY_SIZE`` entries) and return its length, or ``ISO_DISCONNECTED``.
    """
    cdef array.array work = _new(ISO_WORK_SIZE(m.n, m.nd, m.nc))
    cdef int *dec = m.dec.data.as_ints
    return iso_canonical(m.n, m.fp.data.as_ints, m.nd, dec, m.nc, dec + m.nd * m.n,
                         <iso_root_rule> rule, key.data.as_ints, relabelling, roots, nroots,
                         work.data.as_ints)


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


cdef int _search(_Map m, int rule, bint all_leaves, array.array sequence, int *relabelling,
                 iso_leaf_callback callback, void *arg) except -1:
    r"""
    Run ``iso_search`` on ``m`` under the rule ``rule``: write the minimal
    sequence of blocks into ``sequence`` (of ``ISO_SEQUENCE_SIZE`` entries)
    and return its length. A Python exception raised by ``callback``
    propagates.
    """
    cdef array.array work = _new(ISO_SEARCH_WORK_SIZE(m.n, m.nd, m.nc))
    cdef int *dec = m.dec.data.as_ints
    cdef int length
    iso_search(m.n, m.fp.data.as_ints, m.nd, dec, m.nc, dec + m.nd * m.n, <iso_root_rule> rule,
               all_leaves, sequence.data.as_ints, &length, relabelling, callback, arg,
               work.data.as_ints)
    return length


cdef int _key_of(_Map m, array.array r, array.array key) except -1:
    r"""
    Write into ``key`` (of ``ISO_KEY_SIZE`` entries) the key of the
    relabelling ``r`` of ``m`` by the ``2*ne`` labels of its active edges
    (``-1`` on the half-edges of the inactive edges) and return its length.
    """
    cdef int n = m.n, k = 2 * m.ne, h
    cdef int *f = m.fp.data.as_ints
    cdef int *a = r.data.as_ints
    cdef int *dec = m.dec.data.as_ints
    cdef array.array rinv = _new(n)
    cdef array.array names = _new(n)
    # the labels of the inactive halves of folded edges keep -1
    for h in range(n):
        rinv.data.as_ints[h] = -1
        names.data.as_ints[h] = -1
    for h in range(n):
        if f[h] != -1:
            rinv.data.as_ints[a[h]] = h
    return iso_key(n, f, m.nd, dec, m.nc, dec + m.nd * n, a, rinv.data.as_ints, k,
                   key.data.as_ints, names.data.as_ints)


cdef class _Automorphisms:
    r"""
    The automorphisms collected from the leaves of ``iso_search``: ``ra``
    is the first leaf reaching the best sequence so far, ``result`` the
    automorphisms ``g[h] = rb^{-1}[ra[h]]`` over the leaves ``rb`` reaching
    it, and ``rbinv`` a work array of ``n`` entries.
    """
    cdef int n
    cdef array.array ra
    cdef array.array rbinv
    cdef list result


cdef int _automorphism_leaf(const int *leaf, int first, void *arg) except -1:
    r"""
    Leaf callback of ``iso_search`` for :func:`automorphisms`: append the
    automorphism of the leaf to the list, restarted at a smaller leaf.
    """
    cdef _Automorphisms s = <_Automorphisms> arg
    cdef int n = s.n, h
    cdef int *ra = s.ra.data.as_ints
    cdef int *rbinv = s.rbinv.data.as_ints
    cdef array.array g
    if first:
        memcpy(ra, leaf, n * sizeof(int))
        s.result = []
    for h in range(n):
        if leaf[h] != -1:
            rbinv[leaf[h]] = h
    g = _new(n)
    for h in range(n):
        if ra[h] == -1:
            g.data.as_ints[h] = h
        else:
            g.data.as_ints[h] = rbinv[ra[h]]
    s.result.append(g)
    return 0


def canonical_labelling(fp, data=(), classes=(), root_rule=ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA):
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

    - ``root_rule`` -- one of the integers :data:`ROOTS_ALL`,
      :data:`ROOTS_HALF_EDGE_MAX`, :data:`ROOTS_HALF_EDGE_FEWEST`,
      :data:`ROOTS_EDGE_MAX` and :data:`ROOTS_EDGE_FEWEST`, or one of the
      last four combined with :data:`ROOTS_DATA` by ``|``: the half-edges
      tried as roots (see :mod:`combisurf.isomorphism_naive`); any other
      value raises a ``ValueError``. The default is
      ``ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA``: decorations narrow the roots
      tried, and on an undecorated map it is :data:`ROOTS_HALF_EDGE_FEWEST`

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
        array('i', [4, 5, 0, 1, 2, 3])
        sage: canonical_labelling([2, 5, 4, 1, 0, 3], classes=([5, 5, 2, 5, 5, 5],))
        array('i', [1, 0, 5, 4, 3, 2])

    TESTS:

    A folded edge, a label gap, a single half-edge and the empty map::

        sage: canonical_labelling([2, -1, 0, 3])
        array('i', [2, 3, 1, 0])
        sage: canonical_labelling([4, -1, -1, -1, 5, 0])
        array('i', [2, 3, 4, 5, 1, 0])
        sage: canonical_labelling([0, -1])
        array('i', [0, 1])
        sage: canonical_labelling([])
        array('i')

    A disconnected map, labelled component after component::

        sage: canonical_labelling([2, 5, 4, 1, 0, 3, 6, -1])
        array('i', [2, 3, 4, 5, 6, 7, 0, 1])
        sage: canonical_labelling([0, 1, 2, 3], classes=([0, 1, 1, 0],))
        array('i', [0, 1, 3, 2])

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

    Under a rule, the triangle ``(0,1,2)(~0,~2,~1)`` beside the folded edge
    ``(3)`` is labelled from the folded edge or from the triangle::

        sage: from combisurf.isomorphism import ROOTS_HALF_EDGE_MAX, ROOTS_EDGE_FEWEST
        sage: canonical_labelling([2, 5, 4, 1, 0, 3, 6, -1], root_rule=ROOTS_EDGE_FEWEST)
        array('i', [2, 3, 4, 5, 6, 7, 0, 1])
        sage: canonical_labelling([2, 5, 4, 1, 0, 3, 6, -1], root_rule=ROOTS_HALF_EDGE_MAX)
        array('i', [0, 1, 2, 3, 4, 5, 6, 7])

    A rule that is not one of the nine integers::

        sage: canonical_labelling([0, -1], root_rule=5)
        Traceback (most recent call last):
        ...
        ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got 5
        sage: canonical_labelling([0, -1], root_rule="edge-max")
        Traceback (most recent call last):
        ...
        ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got 'edge-max'
        sage: canonical_labelling([0, -1], root_rule=int(3))
        array('i', [0, 1])
        sage: from combisurf.isomorphism import ROOTS_DATA
        sage: canonical_labelling([0, -1], root_rule=ROOTS_HALF_EDGE_MAX | ROOTS_DATA)
        array('i', [0, 1])
        sage: canonical_labelling([0, -1], root_rule=13)
        Traceback (most recent call last):
        ...
        ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got 13
    """
    cdef int rule = _check_root_rule(root_rule)
    cdef _Map m = _prepare(fp, data, classes)
    cdef array.array key = _new(ISO_KEY_SIZE(m.n, m.nd, m.nc))
    cdef array.array r = _new(m.n)
    if _run(m, rule, key, r.data.as_ints, NULL, NULL) == ISO_DISCONNECTED:
        _search(m, rule, False, _new(ISO_SEQUENCE_SIZE(m.n, m.nd, m.nc)), r.data.as_ints, NULL, NULL)
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


def canonical_form(fp, data=(), classes=(), root_rule=ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA):
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

    - ``root_rule`` -- the half-edges tried as roots, as in
      :func:`canonical_labelling` (default: ``ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA``)

    OUTPUT: a tuple ``(fp_c, d_1, ..., d_k, c_1, ..., c_m)`` of tuples:
    ``fp_c`` is the face permutation in the canonical labels (of length
    ``2*ne``, ``-1`` on the inactive half-edges of folded edges), ``d_i``
    the values of the ``i``-th array of ``data`` read in the increasing
    order of the canonical labels of the active half-edges, and ``c_i`` the
    values of the ``i``-th array of ``classes`` read in the same order and
    renumbered ``0, 1, 2, ...`` by first appearance. Two decorated maps are
    isomorphic if and only if their canonical forms under the same rule are
    equal.

    EXAMPLES::

        sage: from combisurf.isomorphism import canonical_form
        sage: canonical_form([2, 5, 4, 1, 0, 3])
        ((2, 5, 4, 1, 0, 3),)
        sage: canonical_form([4, 3, 0, 5, 2, 1])
        ((2, 5, 4, 1, 0, 3),)
        sage: canonical_form([2, 5, 4, 1, 0, 3], ([0, 0, 1, 0, 0, 0],))
        ((2, 5, 4, 1, 0, 3), (1, 0, 0, 0, 0, 0))
        sage: canonical_form([2, 3, 4, 5, 0, 1], ([1, 0, 1, 0, 1, 0],), ([7, 3, 7, 3, 7, 3],))
        ((2, 3, 4, 5, 0, 1), (1, 0, 1, 0, 1, 0), (0, 1, 0, 1, 0, 1))

    TESTS:

    A folded edge, a label gap, the empty map and a disconnected map::

        sage: canonical_form([2, -1, 0, 3])
        ((0, 2, 1, -1),)
        sage: canonical_form([4, -1, -1, -1, 5, 0])
        ((2, 0, 1, -1),)
        sage: canonical_form([], ([],), ([], []))
        ((), (), (), ())
        sage: canonical_form([0, 1, 2, 3], classes=([0, 1, 1, 0],))
        ((0, 1, 2, 3), (0, 1, 0, 1))

    A value of ``data`` outside the range of a C ``int``::

        sage: canonical_form([1, 0], ([0, -2**40],))
        Traceback (most recent call last):
        ...
        OverflowError: ...

    Under a rule, the input is checked as under the default one::

        sage: from combisurf.isomorphism import ROOTS_EDGE_MAX
        sage: canonical_form([1, 0], ([0, -2**40],), root_rule=ROOTS_EDGE_MAX)
        Traceback (most recent call last):
        ...
        OverflowError: ...
        sage: canonical_form([2, 5, 4, 1, 0, 3, 6, -1], root_rule=ROOTS_EDGE_MAX)
        ((2, 5, 4, 1, 0, 3, 6, -1),)
    """
    cdef int rule = _check_root_rule(root_rule)
    cdef _Map m = _prepare(fp, data, classes)
    cdef array.array key = _new(ISO_KEY_SIZE(m.n, m.nd, m.nc))
    cdef array.array r = _new(m.n)
    if _run(m, rule, key, r.data.as_ints, NULL, NULL) == ISO_DISCONNECTED:
        _search(m, rule, False, _new(ISO_SEQUENCE_SIZE(m.n, m.nd, m.nc)), r.data.as_ints, NULL, NULL)
        _key_of(m, r, key)
    return _form(m, key, r, data)


cdef array.array _automorphism_from(_Map m, array.array ra, int root, array.array rb,
                                    array.array rbinv, array.array queue):
    r"""
    Return the automorphism sending the half-edge of label ``i`` from the
    relabelling ``ra`` to the half-edge of label ``i`` from ``root``;
    ``rb``, ``rbinv`` and ``queue`` are buffers of ``m.n`` entries.
    """
    cdef int n = m.n, h
    cdef array.array g = _new(n)
    # iso_relabelling_from writes -1 into rbinv at the labels of the
    # inactive halves of folded edges; the automorphism maps those
    # halves too, so the inverse is rebuilt on every labelled half-edge
    iso_relabelling_from(n, m.fp.data.as_ints, root,
                         rb.data.as_ints, rbinv.data.as_ints, queue.data.as_ints)
    for h in range(n):
        if rb.data.as_ints[h] != -1:
            rbinv.data.as_ints[rb.data.as_ints[h]] = h
    for h in range(n):
        if ra.data.as_ints[h] == -1:
            g.data.as_ints[h] = h
        else:
            g.data.as_ints[h] = rbinv.data.as_ints[ra.data.as_ints[h]]
    return g


cdef inline int _find(int *parent, int x):
    r"""
    Return the representative of ``x`` in the union-find ``parent``, with
    path halving.
    """
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


cdef list _leaf_automorphisms(_Map m, int rule, array.array ra):
    r"""
    Return the list of the automorphisms of the disconnected map ``m``, from
    the leaves of ``iso_search`` under the rule ``rule``, in their order;
    ``ra`` is a buffer of ``m.n`` entries.
    """
    cdef _Automorphisms s = _Automorphisms.__new__(_Automorphisms)
    s.n = m.n
    s.ra = ra
    s.rbinv = _new(m.n)
    s.result = []
    _search(m, rule, True, _new(ISO_SEQUENCE_SIZE(m.n, m.nd, m.nc)), NULL, _automorphism_leaf,
            <void *> s)
    return s.result


cdef list _connected_generators(_Map m, array.array ra, array.array roots, int count):
    r"""
    Return the generators of the automorphism group of the connected map
    ``m`` from its ``count`` roots of minimal key ``roots`` and the
    relabelling ``ra`` from the first of them: the automorphism of root
    ``b_j`` is built only when ``b_j`` is outside the orbit of ``b_0`` under
    the generators kept (a union-find over the half-edges).
    """
    cdef int n = m.n, j, h, x, y, b0
    cdef array.array parent, rb, rbinv, queue, g
    cdef int *p
    cdef int *r
    gens = []
    if count <= 1:
        return gens
    parent = _new(n)
    rb = _new(n)
    rbinv = _new(n)
    queue = _new(n)
    p = parent.data.as_ints
    r = roots.data.as_ints
    for h in range(n):
        p[h] = h
    b0 = r[0]
    for j in range(1, count):
        # the orbit of b0 under the generators kept is its class in the
        # union-find; an automorphism is determined by the image of b0
        if _find(p, r[j]) == _find(p, b0):
            continue
        g = _automorphism_from(m, ra, r[j], rb, rbinv, queue)
        for h in range(n):
            x = _find(p, h)
            y = _find(p, g.data.as_ints[h])
            if x != y:
                p[x] = y
        gens.append(g)
    return gens


cdef array.array _identity(int n):
    cdef array.array g = _new(n)
    cdef int h
    for h in range(n):
        g.data.as_ints[h] = h
    return g


cdef list _wreath(_Map m, bint lift, bint labels=False):
    r"""
    Return the structure of the automorphism group of the disconnected map
    ``m`` with its arrays of ``data`` alone, a product of wreath products.

    The components are taken as connected maps of their own (their
    half-edges in increasing order, ``data`` restricted, ``classes``
    ignored) and grouped by canonical key into kinds, a kind being the set
    of the components isomorphic to a given one. The result lists the kinds,
    ordered by the smallest half-edge of their first copy, as tuples
    ``(copies, generators, count, local, roots)``:

    - ``copies`` -- the list of the copies of the kind, ordered by their
      smallest half-edge, each a pair ``(halves, r)`` of an ``array('i')``
      of the half-edges of the component in increasing order and of the
      canonical relabelling ``r`` of the component as a map of its own
      (``r[i]`` is the label of ``halves[i]``); the isomorphism from a copy
      to another sends the half-edge of label ``i`` to the one of label
      ``i`` (see ``_place``)
    - ``generators`` -- the generators of the automorphism group of the
      first copy, lifted to ``m`` (the identity elsewhere); empty unless
      ``lift``
    - ``count`` -- the order of the automorphism group of a copy
    - ``local`` -- the first copy as a map of its own (a ``_Map``, in the
      local indices of ``halves``) if ``labels``, ``None`` otherwise
    - ``roots`` -- the canonical labels of the roots of minimal key of the
      first copy, in increasing order (an ``array('i')`` starting with
      ``0``; its automorphism sending the half-edge of label ``0`` to the
      one of label ``l`` exists exactly for these ``l``) if ``labels``,
      ``None`` otherwise

    The group is the product over the kinds of the wreath product of the
    group of a copy by the symmetric group on the copies.
    """
    cdef int n = m.n, c, h, x, y, i, j, size, cnt, length, top
    cdef int ncomp = 0
    cdef int *f = m.fp.data.as_ints
    cdef int *dec = m.dec.data.as_ints
    cdef array.array comp = _new(n)
    cdef array.array stack = _new(n)
    cdef array.array idx = _new(n)
    cdef array.array halves_c, lf, key, ra, roots, g, gl
    cdef _Map mc
    halves = []
    for h in range(n):
        comp.data.as_ints[h] = -1
    for h in range(n):
        if f[h] == -1 or comp.data.as_ints[h] != -1:
            continue
        comp.data.as_ints[h] = ncomp
        halves.append([])
        stack.data.as_ints[0] = h
        top = 1
        while top:
            top -= 1
            x = stack.data.as_ints[top]
            # an inactive twin belongs to the component but has no face;
            # each half-edge is pushed once, when it gets its component
            for y in (f[x], x ^ 1):
                if comp.data.as_ints[y] == -1:
                    comp.data.as_ints[y] = ncomp
                    if f[y] != -1:
                        stack.data.as_ints[top] = y
                        top += 1
        ncomp += 1
    for h in range(n):
        if comp.data.as_ints[h] != -1:
            halves[comp.data.as_ints[h]].append(h)

    position = {}
    kinds = []
    for c in range(ncomp):
        halves_c = array.array('i', halves[c])
        size = len(halves_c)
        for i in range(size):
            idx.data.as_ints[halves_c.data.as_ints[i]] = i
        lf = _new(size)
        for i in range(size):
            y = f[halves_c.data.as_ints[i]]
            lf.data.as_ints[i] = -1 if y == -1 else idx.data.as_ints[y]
        ld = tuple([array.array('i', [dec[j * n + halves_c.data.as_ints[i]] for i in range(size)])
                    for j in range(m.nd)])
        mc = _convert(lf, ld, ())
        key = _new(ISO_KEY_SIZE(size, mc.nd, 0))
        ra = _new(size)
        roots = _new(size)
        length = _run(mc, _ROOTS_ALL, key, ra.data.as_ints, roots.data.as_ints, &cnt)
        if length == ISO_DISCONNECTED:
            raise RuntimeError("a component of a map is not connected")
        kb = (<char *> key.data.as_ints)[:length * sizeof(int)]
        if kb in position:
            kinds[position[kb]][0].append((halves_c, ra))
            continue
        position[kb] = len(kinds)
        gens = []
        if lift:
            for gl in _connected_generators(mc, ra, roots, cnt):
                g = _identity(n)
                for i in range(size):
                    g.data.as_ints[halves_c.data.as_ints[i]] = halves_c.data.as_ints[gl.data.as_ints[i]]
                gens.append(g)
        if labels:
            kinds.append(([(halves_c, ra)], gens, cnt, mc,
                          array.array('i', sorted([ra.data.as_ints[roots.data.as_ints[i]] for i in range(cnt)]))))
        else:
            kinds.append(([(halves_c, ra)], gens, cnt, None, None))
    return kinds


cdef int _place(array.array g, tuple a, tuple b) except -1:
    r"""
    Write into ``g`` the isomorphism from the copy ``a`` to the copy ``b``
    (see ``_wreath``): the half-edge of label ``i`` in ``a`` goes to the
    one of label ``i`` in ``b``.
    """
    cdef array.array ha = a[0], ra = a[1], hb = b[0], rb = b[1]
    cdef int size = len(ha), i
    cdef array.array rinv = _new(size)
    for i in range(size):
        rinv.data.as_ints[rb.data.as_ints[i]] = i
    for i in range(size):
        g.data.as_ints[ha.data.as_ints[i]] = hb.data.as_ints[rinv.data.as_ints[ra.data.as_ints[i]]]
    return 0


cdef list _wreath_generators(_Map m):
    r"""
    Return the generators of the automorphism group of the disconnected map
    ``m`` without arrays of ``classes``: for each kind of isomorphic
    components (see ``_wreath``), the generators of the group of its first
    copy, then the transposition of its first two copies when it has at
    least two, then the cycle of its copies when it has at least three.
    """
    cdef array.array g
    cdef int t, k
    gens = []
    for copies, aut, _count, _local, _roots in _wreath(m, True):
        gens.extend(aut)
        k = len(copies)
        if k >= 2:
            g = _identity(m.n)
            _place(g, copies[0], copies[1])
            _place(g, copies[1], copies[0])
            gens.append(g)
        if k >= 3:
            g = _identity(m.n)
            for t in range(k):
                _place(g, copies[t], copies[(t + 1) % k])
            gens.append(g)
    return gens


cdef _wreath_order(_Map m):
    r"""
    Return the number of automorphisms of the disconnected map ``m`` without
    arrays of ``classes``, as a Python integer: the product over the kinds
    of `|Aut(C)|^k k!` (see ``_wreath``).
    """
    order = 1
    for copies, _aut, count, _local, _roots in _wreath(m, False):
        order *= count ** len(copies) * factorial(len(copies))
    return order


cdef class _Stabiliser:
    r"""
    The state of the search of :func:`_stabiliser` on a disconnected map.

    The levels are the components, in the order of ``_wreath`` (kinds in
    order, copies in order). Per level ``j``: ``kind[j]`` its kind and
    ``off[j]`` the offset of its half-edges in ``lab``, by canonical label
    (``lab[off[j] + i]`` is its half-edge of label ``i``). Per kind ``t``:
    its copies are the levels ``start[t], ..., start[t + 1] - 1``; its
    first copy as a map of its own has the face permutation ``cfp``, the
    canonical relabelling ``cra`` (local index to label) and its inverse
    ``crainv``, all of ``csize[t]`` entries at the offset ``coff[t]``; the
    labels of its roots of minimal key are ``labels[loff[t]], ...,
    labels[loff[t + 1] - 1]``, in increasing order.

    The search state: ``fw`` and ``bw`` the forward and backward maps of the
    values of each array of ``classes`` (``-1`` where unset, the array
    ``q`` at the offset ``q * n``), ``trail`` the entries of ``fw`` set, in
    order, ``used`` the copies already images, and per level its image
    ``target``, the index ``choice`` of its label in ``labels`` and the
    length ``mark`` of the trail before its step. ``rb``, ``rbinv``,
    ``queue``, ``perm`` and ``ident`` are buffers of the size of the
    largest component.
    """
    cdef int n
    cdef int nd
    cdef int nc
    cdef int ncomp
    cdef int ntrail
    cdef array.array fp
    cdef array.array dec
    cdef array.array kind, off, lab, start, coff, csize, cfp, cra, crainv, loff, labels
    cdef array.array fw, bw, trail, used, target, choice, mark
    cdef array.array rb, rbinv, queue, perm, ident

    cdef int relabel(self, int t, int l) except -1:
        r"""
        Write into ``perm`` the automorphism of the first copy of the kind
        ``t`` sending the half-edge of label ``0`` to the one of label
        ``l``, on the labels: ``perm[i]`` is the label of the image of the
        half-edge of label ``i``.
        """
        cdef int o = self.coff.data.as_ints[t], size = self.csize.data.as_ints[t], h
        cdef int *cra = self.cra.data.as_ints + o
        cdef int *rb = self.rb.data.as_ints
        cdef int *rbinv = self.rbinv.data.as_ints
        cdef int *perm = self.perm.data.as_ints
        # iso_relabelling_from writes -1 into rbinv at the labels of the
        # inactive halves of folded edges, which are rebuilt from rb
        iso_relabelling_from(size, self.cfp.data.as_ints + o, self.crainv.data.as_ints[o + l],
                             rb, rbinv, self.queue.data.as_ints)
        for h in range(size):
            rbinv[rb[h]] = h
        for h in range(size):
            perm[h] = cra[rbinv[h]]
        return 0

    cdef int step(self, int a, int b, int *perm) except -1:
        r"""
        Extend the maps of values by the isomorphism from the level ``a`` to
        the level ``b`` sending the half-edge of label ``i`` to the one of
        label ``perm[i]``; return ``0`` on a conflict (the entries set
        before it stay on the trail) and ``1`` otherwise.
        """
        cdef int n = self.n, size = self.csize.data.as_ints[self.kind.data.as_ints[a]]
        cdef int *f = self.fp.data.as_ints
        cdef int *cls = self.dec.data.as_ints + self.nd * n
        cdef int *la = self.lab.data.as_ints + self.off.data.as_ints[a]
        cdef int *lb = self.lab.data.as_ints + self.off.data.as_ints[b]
        cdef int *fw = self.fw.data.as_ints
        cdef int *bw = self.bw.data.as_ints
        cdef int *trail = self.trail.data.as_ints
        cdef int i, q, h, y, v, w, base
        for i in range(size):
            h = la[i]
            if f[h] == -1:
                continue
            y = lb[perm[i]]
            for q in range(self.nc):
                base = q * n
                v = cls[base + h]
                w = cls[base + y]
                if fw[base + v] == -1:
                    if bw[base + w] != -1:
                        return 0
                    fw[base + v] = w
                    bw[base + w] = v
                    trail[self.ntrail] = base + v
                    self.ntrail += 1
                elif fw[base + v] != w:
                    return 0
        return 1

    cdef void undo(self, int length):
        r"""
        Unset the entries of the maps of values set after the first
        ``length`` entries of the trail.
        """
        cdef int *fw = self.fw.data.as_ints
        cdef int *bw = self.bw.data.as_ints
        cdef int e
        while self.ntrail > length:
            self.ntrail -= 1
            e = self.trail.data.as_ints[self.ntrail]
            bw[(e // self.n) * self.n + fw[e]] = -1
            fw[e] = -1

    cdef int extend(self, int i) except -1:
        r"""
        Search images for the levels ``i + 1, ..., ncomp - 1``, each a copy
        of its kind not used yet, through any label of its roots, in
        order; return ``1`` at the first leaf, the choices left in
        ``target`` and ``choice``, and ``0`` if there is none.
        """
        cdef int *kind = self.kind.data.as_ints
        cdef int *start = self.start.data.as_ints
        cdef int *loff = self.loff.data.as_ints
        cdef int *labels = self.labels.data.as_ints
        cdef int *used = self.used.data.as_ints
        cdef int *target = self.target.data.as_ints
        cdef int *choice = self.choice.data.as_ints
        cdef int *mark = self.mark.data.as_ints
        cdef int j = i + 1, t, b, placed
        if j < self.ncomp:
            target[j] = start[kind[j]]
            choice[j] = 0
        while True:
            if j == self.ncomp:
                return 1
            if j <= i:
                return 0
            t = kind[j]
            placed = 0
            while target[j] < start[t + 1]:
                b = target[j]
                if used[b] or choice[j] >= loff[t + 1] - loff[t]:
                    target[j] += 1
                    choice[j] = 0
                    continue
                mark[j] = self.ntrail
                self.relabel(t, labels[loff[t] + choice[j]])
                if self.step(j, b, self.perm.data.as_ints):
                    used[b] = 1
                    placed = 1
                    break
                self.undo(mark[j])
                choice[j] += 1
            if placed:
                j += 1
                if j < self.ncomp:
                    target[j] = start[kind[j]]
                    choice[j] = 0
            else:
                j -= 1
                if j > i:
                    self.undo(mark[j])
                    used[target[j]] = 0
                    choice[j] += 1

    cdef array.array leaf(self, int i):
        r"""
        Return the automorphism of the map given by the choices of the
        levels ``i, ..., ncomp - 1`` (the identity on the other levels and
        on the inactive edges).
        """
        cdef array.array g = _identity(self.n)
        cdef int *lab = self.lab.data.as_ints
        cdef int *perm = self.perm.data.as_ints
        cdef int j, t, k, oa, ob
        for j in range(i, self.ncomp):
            t = self.kind.data.as_ints[j]
            self.relabel(t, self.labels.data.as_ints[self.loff.data.as_ints[t] + self.choice.data.as_ints[j]])
            oa = self.off.data.as_ints[j]
            ob = self.off.data.as_ints[self.target.data.as_ints[j]]
            for k in range(self.csize.data.as_ints[t]):
                g.data.as_ints[lab[oa + k]] = lab[ob + perm[k]]
        return g


cdef _Stabiliser _stabiliser_state(_Map m):
    r"""
    Return the :class:`_Stabiliser` of the disconnected map ``m``, its
    search state cleared.
    """
    cdef _Stabiliser s = _Stabiliser.__new__(_Stabiliser)
    cdef int n = m.n, size, i, largest = 0
    cdef array.array halves_c, r
    cdef _Map local
    kind = []
    off = []
    lab = []
    start = [0]
    coff = []
    csize = []
    cfp = []
    cra = []
    crainv = []
    loff = [0]
    labels = []
    for t, (copies, _aut, _count, local, roots) in enumerate(_wreath(m, False, True)):
        size = len(copies[0][0])
        largest = max(largest, size)
        coff.append(len(cfp))
        csize.append(size)
        cfp.extend(local.fp)
        r = copies[0][1]
        cra.extend(r)
        inv = [0] * size
        for i in range(size):
            inv[r.data.as_ints[i]] = i
        crainv.extend(inv)
        labels.extend(roots)
        loff.append(len(labels))
        for halves_c, r in copies:
            kind.append(t)
            off.append(len(lab))
            for i in range(size):
                inv[r.data.as_ints[i]] = halves_c.data.as_ints[i]
            lab.extend(inv)
        start.append(len(kind))
    s.n = n
    s.nd = m.nd
    s.nc = m.nc
    s.ncomp = len(kind)
    s.ntrail = 0
    s.fp = m.fp
    s.dec = m.dec
    s.kind = array.array('i', kind)
    s.off = array.array('i', off)
    s.lab = array.array('i', lab)
    s.start = array.array('i', start)
    s.coff = array.array('i', coff)
    s.csize = array.array('i', csize)
    s.cfp = array.array('i', cfp)
    s.cra = array.array('i', cra)
    s.crainv = array.array('i', crainv)
    s.loff = array.array('i', loff)
    s.labels = array.array('i', labels)
    s.fw = array.array('i', [-1]) * (m.nc * n)
    s.bw = array.array('i', [-1]) * (m.nc * n)
    s.trail = _new(m.nc * n)
    s.used = array.array('i', [0]) * s.ncomp
    s.target = _new(s.ncomp)
    s.choice = _new(s.ncomp)
    s.mark = _new(s.ncomp)
    s.rb = _new(largest)
    s.rbinv = _new(largest)
    s.queue = _new(largest)
    s.perm = _new(largest)
    s.ident = array.array('i', range(largest))
    return s


cdef tuple _stabiliser(_Map m):
    r"""
    Return the generators and the order (a Python integer) of the group of
    the automorphisms of the disconnected map ``m`` that preserve the
    partitions defined by its arrays of ``classes``.

    This group `H` is the stabiliser of these partitions in the group `G`
    of ``m`` with its ``data`` alone (see ``_wreath``), computed by a
    subgroup search of Sims. The base is ``b_1, ..., b_c``, the half-edge of
    canonical label ``0`` of each component, in the order of ``_wreath``;
    the pointwise stabiliser of ``b_1, ..., b_{i-1}`` in `G` fixes their
    components and permutes the others as the wreath products do, so the
    orbit of ``b_i`` in it is made of the half-edges whose label is the one
    of a root of minimal key in the copies of the kind of ``b_i`` at the
    levels ``>= i``. The levels are processed from the last to the first:
    each point ``x`` of that orbit outside the orbit of ``b_i`` under the
    generators found so far (a union-find over the half-edges) is tried by
    a depth-first search for an element of `H`, the identity on the levels
    ``< i``, sending ``b_i`` to ``x``, which maps the following levels one
    after the other to unused copies of their kinds and rejects a step
    whose values contradict the bijections of values built so far. The
    first element found is a generator. The generators found at the levels
    ``>= i`` generate the stabiliser of ``b_1, ..., b_{i-1}`` in `H`, so
    the order of `H` is the product over the levels of the size of the
    orbit of ``b_i`` right after its level.

    The search is exponential in the number of components in the worst
    case, when many copies tie and the values of ``classes`` reject the
    candidates only deep in the search.
    """
    cdef _Stabiliser s = _stabiliser_state(m)
    cdef int n = m.n, i, t, b, q, h, x, y, bi, nl
    cdef array.array parent = _new(n)
    cdef array.array sizes = _new(n)
    cdef array.array g
    cdef int *p = parent.data.as_ints
    cdef int *sz = sizes.data.as_ints
    cdef int *used = s.used.data.as_ints
    cdef int *start = s.start.data.as_ints
    cdef int *loff = s.loff.data.as_ints
    cdef int *labels = s.labels.data.as_ints
    cdef int *lab = s.lab.data.as_ints
    cdef int *off = s.off.data.as_ints
    for h in range(n):
        p[h] = h
        sz[h] = 1
    gens = []
    order = 1
    for i in range(s.ncomp - 1, -1, -1):
        t = s.kind.data.as_ints[i]
        nl = loff[t + 1] - loff[t]
        bi = lab[off[i]]
        if start[t + 1] - start[t] == 1 and nl == 1:
            continue
        # the levels < i are fixed: their values are mapped to themselves
        for b in range(i):
            used[b] = 1
            s.step(b, b, s.ident.data.as_ints)
        for b in range(max(start[t], i), start[t + 1]):
            for q in range(nl):
                x = lab[off[b] + labels[loff[t] + q]]
                if x == bi or _find(p, x) == _find(p, bi):
                    continue
                s.target.data.as_ints[i] = b
                s.choice.data.as_ints[i] = q
                s.mark.data.as_ints[i] = s.ntrail
                used[b] = 1
                s.relabel(t, labels[loff[t] + q])
                if s.step(i, b, s.perm.data.as_ints) and s.extend(i):
                    g = s.leaf(i)
                    gens.append(g)
                    for h in range(n):
                        x = _find(p, h)
                        y = _find(p, g.data.as_ints[h])
                        if x != y:
                            p[x] = y
                            sz[y] += sz[x]
                s.undo(s.mark.data.as_ints[i])
                for h in range(i, s.ncomp):
                    used[h] = 0
        s.undo(0)
        for b in range(i):
            used[b] = 0
        order *= sz[_find(p, bi)]
    return gens, order


def automorphisms(fp, data=(), classes=(), root_rule=ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA):
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

    - ``root_rule`` -- the half-edges tried as roots, as in
      :func:`canonical_labelling` (default: ``ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA``); the group
      is the same under every rule, listed in the order of the search of
      the rule

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

    The four edges ``(0)(~0)(1)(~1)(2)(~2)(3)(~3)``, which may be flipped
    and permuted, and with classes that pair them::

        sage: len(automorphisms(list(range(8))))
        384
        sage: len(automorphisms(list(range(8)), classes=([0, 1, 0, 1, 2, 3, 2, 3],)))
        32

    Under a rule::

        sage: from combisurf.isomorphism import ROOTS_HALF_EDGE_FEWEST
        sage: A = automorphisms([2, 3, 4, 5, 0, 1, 6, -1], root_rule=ROOTS_HALF_EDGE_FEWEST)
        sage: len(A), sorted(A) == sorted(automorphisms([2, 3, 4, 5, 0, 1, 6, -1]))
        (6, True)
    """
    cdef int rule = _check_root_rule(root_rule)
    cdef _Map m = _prepare(fp, data, classes)
    cdef int n = m.n
    cdef int count, j
    cdef array.array key = _new(ISO_KEY_SIZE(n, m.nd, m.nc))
    cdef array.array ra = _new(n)
    cdef array.array roots = _new(n)
    cdef array.array rb, rbinv, queue
    if _run(m, rule, key, ra.data.as_ints, roots.data.as_ints, &count) == ISO_DISCONNECTED:
        return _leaf_automorphisms(m, rule, ra)
    if count == 0:
        return [array.array('i', range(n))]
    rb = _new(n)
    rbinv = _new(n)
    queue = _new(n)
    result = []
    for j in range(count):
        result.append(_automorphism_from(m, ra, roots.data.as_ints[j], rb, rbinv, queue))
    return result


def automorphism_generators(fp, data=(), classes=(), root_rule=ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA):
    r"""
    Return generators of the automorphism group of the map with face
    permutation ``fp`` decorated by ``data`` and ``classes``.

    The result is the one of
    :func:`combisurf.isomorphism_naive.automorphism_generators` on a
    connected map; on a disconnected map the generators come from the
    structure of the group and generate the same group.

    INPUT:

    - ``fp`` -- a face permutation (list or ``array('i')``, ``-1`` on
      inactive half-edges)

    - ``data`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept exactly (default: no decoration); the values on the
      active half-edges are integers in the range of a C ``int``

    - ``classes`` -- a tuple of arrays of length ``len(fp)`` indexed by
      half-edges, kept up to a relabelling of their values (default: no
      decoration)

    - ``root_rule`` -- the half-edges tried as roots, as in
      :func:`canonical_labelling` (default: ``ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA``): on a
      connected map the generators may differ from a rule to another, the
      group generated does not; on a disconnected map the generators do not
      depend on the rule

    OUTPUT: a list of automorphisms (as in :func:`automorphisms`) that
    generate the group; empty when the group is trivial. There are three
    cases.

    - On a connected map the automorphisms are scanned in the order of
      :func:`automorphisms` under the rule, and one is kept when it is not
      in the group generated by those already kept. The scan runs over the
      roots of minimal key: the automorphism of root ``b_j`` is built only when
      ``b_j`` is outside the orbit of ``b_0`` under the generators kept (a
      union-find over the half-edges), so the list of automorphisms and the
      group closure are never computed.

    - On a disconnected map without ``classes`` the group is the product,
      over the kinds of isomorphic components, of `\mathrm{Aut}(C) \wr
      S_k` (``data`` links no two components). The kinds are ordered by
      the smallest half-edge of their first copy `C`, and the copies of a
      kind likewise; for each kind the generators are those of
      `\mathrm{Aut}(C)` (by the connected scan), then the transposition of
      the first two copies if `k \geq 2`, then the cycle of the copies if
      `k \geq 3`. The copies are exchanged by the isomorphism sending the
      half-edge of canonical label `i` of one to the one of label `i` of
      the other. These generators differ from those of
      :func:`combisurf.isomorphism_naive.automorphism_generators` (fewer
      of them) and generate the same group. The components are labelled
      under :data:`ROOTS_ALL` whatever the rule.

    - On a disconnected map with ``classes`` the group `H` is the subgroup
      of the previous one (for ``data`` alone) of the automorphisms that
      preserve the partitions defined by the arrays of ``classes``. The
      half-edges of canonical label `0` of the components, in the order
      above, form a base `b_1, \ldots, b_c`; the levels are processed from
      the last to the first, and at the level `i` each image of `b_i` under
      the automorphisms fixing the components of `b_1, \ldots, b_{i-1}`
      that is outside its orbit under the generators found so far is
      tried by a depth-first search for an element of `H` sending `b_i`
      to it, which maps the following components one after the other to
      copies of their kinds. The first element found is a generator; the
      generators come in the order found. The search is exponential in the
      number of components in the worst case, when many components are
      isomorphic and the values of ``classes`` reject the candidates only
      deep in the search.

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
        [array('i', [1, 0, 2, 3]), array('i', [2, 3, 0, 1])]

    Three triangles: the two generators of the group of one triangle, the
    transposition of the first two and the cycle of the three::

        sage: fp = [x + 6*i for i in range(3) for x in [2, 5, 4, 1, 0, 3]]
        sage: automorphism_generators(fp)
        [array('i', [1, 0, 5, 4, 3, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]),
         array('i', [2, 3, 4, 5, 0, 1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]),
         array('i', [6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17]),
         array('i', [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 0, 1, 2, 3, 4, 5])]

    The four edges ``(0)(~0)(1)(~1)(2)(~2)(3)(~3)`` with classes that pair
    them: the generators generate the 32 automorphisms::

        sage: fp = list(range(8))
        sage: c = [0, 1, 0, 1, 2, 3, 2, 3]
        sage: gens = automorphism_generators(fp, classes=(c,))
        sage: gens
        [array('i', [0, 1, 2, 3, 5, 4, 7, 6]),
         array('i', [0, 1, 2, 3, 6, 7, 4, 5]),
         array('i', [1, 0, 3, 2, 4, 5, 6, 7]),
         array('i', [2, 3, 0, 1, 4, 5, 6, 7]),
         array('i', [4, 5, 6, 7, 0, 1, 2, 3])]
        sage: PermutationGroup([[x + 1 for x in g] for g in gens], domain=list(range(1, 9))).order()
        32

    Under a rule, the same group; on a disconnected map, the same
    generators::

        sage: from combisurf.isomorphism import ROOTS_EDGE_MAX
        sage: automorphism_generators(fp, classes=(c,), root_rule=ROOTS_EDGE_MAX) == gens
        True
        sage: from combisurf.isomorphism import ROOTS_HALF_EDGE_FEWEST
        sage: gens = automorphism_generators([2, 3, 4, 5, 0, 1], classes=([0, 1, 0, 1, 0, 1],), root_rule=ROOTS_HALF_EDGE_FEWEST)
        sage: PermutationGroup([[x + 1 for x in g] for g in gens], domain=list(range(1, 7))).order()
        6
    """
    cdef int rule = _check_root_rule(root_rule)
    cdef _Map m = _prepare(fp, data, classes)
    cdef int n = m.n
    cdef int count
    cdef array.array key = _new(ISO_KEY_SIZE(n, m.nd, m.nc))
    cdef array.array ra = _new(n)
    cdef array.array roots = _new(n)
    if _run(m, rule, key, ra.data.as_ints, roots.data.as_ints, &count) != ISO_DISCONNECTED:
        return _connected_generators(m, ra, roots, count)
    if m.nc == 0:
        return _wreath_generators(m)
    return _stabiliser(m)[0]


def isomorphism(fp1, data1, fp2, data2, classes1=(), classes2=(), root_rule=ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA):
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

    - ``root_rule`` -- the half-edges tried as roots in both maps, as in
      :func:`canonical_labelling` (default: ``ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA``)

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

    Under a rule::

        sage: from combisurf.isomorphism import ROOTS_HALF_EDGE_FEWEST
        sage: isomorphism([2, 5, 4, 1, 0, 3, 6, -1], (), [0, -1, 4, 7, 6, 3, 2, 5], (), root_rule=ROOTS_HALF_EDGE_FEWEST)
        array('i', [2, 3, 4, 5, 6, 7, 0, 1])
        sage: isomorphism(fp1, (), fp3, (), root_rule=ROOTS_HALF_EDGE_FEWEST) is None
        True
    """
    cdef int rule = _check_root_rule(root_rule)
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
    cdef array.array seq1, seq2
    len1 = _run(m1, rule, key1, ra.data.as_ints, NULL, NULL)
    len2 = _run(m2, rule, key2, rb.data.as_ints, NULL, NULL)
    if len1 == ISO_DISCONNECTED or len2 == ISO_DISCONNECTED:
        # a connected map has one block, a disconnected one several
        if len1 != len2:
            return None
        seq1 = _new(ISO_SEQUENCE_SIZE(n1, m1.nd, m1.nc))
        seq2 = _new(ISO_SEQUENCE_SIZE(n2, m2.nd, m2.nc))
        len1 = _search(m1, rule, False, seq1, ra.data.as_ints, NULL, NULL)
        len2 = _search(m2, rule, False, seq2, rb.data.as_ints, NULL, NULL)
        if (len1 != len2
                or memcmp(seq1.data.as_ints, seq2.data.as_ints, len1 * sizeof(int)) != 0):
            return None
    elif (m1.ne != m2.ne or len1 != len2
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


def canonical_key(fp, data=(), classes=(), relabelling=False, automorphism_count=False, root_rule=ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA):
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

    - ``root_rule`` -- the half-edges tried as roots, as in
      :func:`canonical_labelling` (default: ``ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA``); the key and
      the relabelling are those of the rule, the count does not depend on
      it

    OUTPUT: the key, a ``bytes`` object: the entries of the canonical form
    (see :func:`canonical_form`) as C ``int`` in the machine byte order,
    the face permutation first, then each array of ``data``, then each
    array of ``classes``. Two maps with the same numbers of arrays of
    ``data`` and of ``classes`` have equal keys under the same rule if and
    only if their canonical forms are equal, that is, if and only if they
    are isomorphic.
    With ``relabelling=True`` or ``automorphism_count=True``, a tuple
    ``(key, relabelling)``, ``(key, count)`` or ``(key, relabelling,
    count)``, where the relabelling is the one of
    :func:`canonical_labelling` and ``count`` the number of automorphisms
    (the length of the list of :func:`automorphisms`), a Python integer
    that may exceed 64 bits. On a disconnected map without ``classes`` it
    is the order `\prod |\mathrm{Aut}(C)|^k k!` of the product of wreath
    products over the kinds of `k` isomorphic components `C`, computed
    without listing the automorphisms; with ``classes`` it is the order
    of the subgroup that preserves the partitions of ``classes``, computed
    by the search of :func:`automorphism_generators`. With ``classes`` the
    search of the key on a disconnected map follows every root tied in a
    component; a leaf equal to the best one so far differs from it by an
    automorphism fixing the components where their roots agree, so the
    search returns at once to the first level where their roots differ
    and tries its next tied root. On `k` tied copies of a component each
    tied root then costs about one path of `k` components, about `k^2`
    paths in all instead of a number of leaves growing as the
    automorphisms (`2^k k!` for `k` isolated edges with one class value).

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
        array('i', [2, 5, 4, 1, 0, 3, 1, 0, 0, 0, 0, 0])
        sage: r
        array('i', [4, 5, 0, 1, 2, 3])
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

    A disconnected map, searched component after component: the two edges
    ``(0)(~0)(1)(~1)``, with their classes::

        sage: key, r, count = canonical_key([0, 1, 2, 3], classes=([0, 1, 1, 2],), relabelling=True, automorphism_count=True)
        sage: array('i', key), r, count
        (array('i', [0, 1, 2, 3, 0, 1, 0, 2]), array('i', [1, 0, 2, 3]), 2)

    Twenty isolated edges have `2^{20} \cdot 20!` automorphisms::

        sage: canonical_key(list(range(40)), automorphism_count=True)[1]
        2551082656125828464640000

    With one class value, six isolated edges keep their `2^6 \cdot 6!`
    automorphisms, and two values on the halves of each edge leave the
    `6!` permutations of the edges, alone or followed by the flip of every
    edge (which exchanges the two values)::

        sage: canonical_key(list(range(12)), classes=([0] * 12,), automorphism_count=True)[1]
        46080
        sage: canonical_key(list(range(12)), classes=([0, 1] * 6,), automorphism_count=True)[1]
        1440

    A value of ``data`` outside the range of a C ``int``, and one on an
    inactive half-edge, which is ignored::

        sage: canonical_key([1, 0], ([2**40, 0],))
        Traceback (most recent call last):
        ...
        OverflowError: ...
        sage: array('i', canonical_key([0, -1], ([3, None],)))
        array('i', [0, -1, 3])

    Under a rule, the key is the canonical form of that rule flattened::

        sage: from combisurf.isomorphism import ROOTS_HALF_EDGE_MAX
        sage: fp = [2, 3, 4, 5, 0, 1, 6, -1]
        sage: key, r, count = canonical_key(fp, relabelling=True, automorphism_count=True, root_rule=ROOTS_HALF_EDGE_MAX)
        sage: array('i', key), r, count
        (array('i', [2, 3, 4, 5, 0, 1, 6, -1]), array('i', [0, 1, 2, 3, 4, 5, 6, 7]), 6)
        sage: array('i', canonical_key(fp))
        array('i', [0, -1, 4, 5, 6, 7, 2, 3])
        sage: canonical_key(fp, root_rule="all")
        Traceback (most recent call last):
        ...
        ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got 'all'
    """
    cdef int rule = _check_root_rule(root_rule)
    cdef _Map m = _prepare(fp, data, classes)
    cdef array.array key = _new(ISO_KEY_SIZE(m.n, m.nd, m.nc))
    cdef array.array r = _new(m.n)
    cdef int count, length
    order = 0
    length = _run(m, rule, key, r.data.as_ints, NULL, &count)
    if length == ISO_DISCONNECTED:
        _search(m, rule, False, _new(ISO_SEQUENCE_SIZE(m.n, m.nd, m.nc)), r.data.as_ints, NULL, NULL)
        if automorphism_count:
            order = _stabiliser(m)[1] if m.nc else _wreath_order(m)
        length = _key_of(m, r, key)
    out = [(<char *> key.data.as_ints)[:length * sizeof(int)]]
    if relabelling:
        out.append(_complete(r, m.n, m.ne))
    if automorphism_count:
        if order:
            out.append(order)
        else:
            out.append(count if count else 1)
    if len(out) == 1:
        return out[0]
    return tuple(out)
