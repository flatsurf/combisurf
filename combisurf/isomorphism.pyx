# distutils: sources = combisurf/src/omap.c combisurf/src/perm.c
# distutils: include_dirs = combisurf/src
r"""
Isomorphisms, automorphisms and canonical labels of maps

This module computes canonical labellings, canonical forms, automorphism
groups and isomorphisms of maps given by their face permutation, possibly
decorated by integer values on half-edges. The methods of
:class:`~combisurf.oriented_map.OrientedMap` (such as
:meth:`~combisurf.oriented_map.OrientedMap.canonical_labelling`) call this
module.

A map on `n` half-edges is given by its face permutation ``fp``, a list
or ``array('i')`` of length `n` (even) with ``-1`` on the inactive
half-edges; the twin of the half-edge ``h`` is ``h ^ 1``, and a folded
edge, whose twin is inactive, is a fixed point of ``fp``. ``fp`` must be
a permutation of its active half-edges: otherwise a ``ValueError`` is
raised, since the C code would read out of bounds on it; other malformed
input raises a ``TypeError`` or a ``ValueError``. A map may be decorated
by two tuples of arrays indexed by the half-edges: the values of
``data`` are kept exactly by an isomorphism, while those of ``classes``
are kept up to a relabelling of their values, so that only the partition
of the half-edges by the values of each array matters. The canonical
labelling of a decorated map is a relabelling of its half-edges, and its
canonical form is the map read in these labels: the face permutation and
the decorations in the canonical labels, the values of ``classes``
renumbered by first appearance. Two decorated maps are isomorphic if and
only if their canonical forms under the same rule are equal. On a
connected map, the canonical labelling is the relabelling by the
breadth-first walk from a root (:ref:`omap-relabelling`) that gives the
minimal key (:ref:`omap-keys`), the first such root in increasing order.

On a connected map, every active half-edge of the class of the rule (see
below) is tried as a root: the relabelling from it and its key are
computed in C (``combisurf/src/omap.c``) and the minimal key is
kept, which costs `O(n^2)` for `n` half-edges. A disconnected map is
searched in C as well, component after component, for the
lexicographically minimal sequence of blocks, one block per component
(:ref:`omap-blocks`). Without
``classes``, the automorphism group of a disconnected map
is a product of wreath products `\mathrm{Aut}(C) \wr S_k` over the kinds
of `k` isomorphic components `C`: :func:`automorphism_generators` and the
count of :func:`canonical_key` come from this structure, component by
component, without listing the automorphisms (their number is at least
`k!`). With ``classes`` the group is the subgroup of this product that
preserves the partitions of the half-edges by the values of ``classes``:
its generators and its order come from a backtrack search of this
subgroup (the subgroup search of Sims, one base point per component),
which maps the components one after the other to copies of their kinds.
Once a placement of a component that maps no new value of ``classes`` has
failed, the later placements of that component mapping no new value are
skipped: the automorphism of the copy relating the two placements carries
the values, so it exchanges the two images. Tied copies whose values are
all mapped after the first components, as `k` isolated edges with one or
two class values in any pattern, then take a time polynomial in `k`. The
search stays exponential in the worst case, where it is a search for the
automorphisms of a structure on the class values: on a graph given by
its edges as isolated edges (the two values of an edge being its
endpoints), it is a backtrack search for the automorphisms of the graph,
which can take far longer than the search of the canonical key, as on a
random cubic graph or a grid. The search
of the canonical key itself, when there are ``classes``, follows among
the roots tied in a component the first one of each way of extending the
numbering of the class values, since two tied roots extending it alike
are exchanged by an automorphism fixing the components labelled before;
it also returns from a leaf equal to the best one to the first level
where their roots differ, since the subtree left is then the image of
one already searched. `k` tied copies of a component then cost one path
per distinct numbering at each level (one for isolated edges all with
the same value, two for `b + 1` edges with the value 0 beside `b` with
the value 1) instead of a number of leaves growing as the automorphisms.
The search stays exponential in the worst case: tied roots with
different numberings are all followed, and their subtrees may part only
deep, as for a graph on the class values given by its edges as isolated
edges, such as a random cubic graph. The buffers of the C code are
counted in C ``int``: a map
too large for them (about `2^{31} / (14 + 2 nd + 3.5 nc)` half-edges with
``nd`` arrays of ``data`` and ``nc`` of ``classes``) raises an
``OverflowError``.

The searches can be interrupted: Ctrl-C raises a ``KeyboardInterrupt``,
and the ``alarm`` of ``cysignals`` an ``AlarmInterrupt``, in the C
searches as in the search of the automorphism group of a disconnected
map. On the `8 \times 8` grid given by its edges as isolated edges, whose
automorphism count takes far longer than the alarm::

    sage: from cysignals.alarm import alarm
    sage: from combisurf.isomorphism import canonical_key
    sage: g = graphs.Grid2dGraph(8, 8).relabel(inplace=False)
    sage: edges = list(g.edges(labels=False))
    sage: fp = list(range(2 * len(edges)))
    sage: classes = ([x for e in edges for x in e],)
    sage: alarm(0.5); canonical_key(fp, (), classes, automorphism_count=True)
    Traceback (most recent call last):
    ...
    AlarmInterrupt

The decorations are stored in C ``int`` arrays: a value of ``data``
outside the range of a C ``int`` raises an ``OverflowError`` (the values
on inactive half-edges are ignored and never converted). The values of
``classes`` are renumbered before the conversion, so they may be any
hashable values.

Beside :func:`canonical_labelling`, :func:`canonical_form`,
:func:`automorphisms`, :func:`automorphism_generators` and
:func:`isomorphism`, :func:`canonical_key` returns the canonical form
flattened as bytes, a compact hashable key for the exploration of large
collections of maps up to isomorphism, and on request the canonical
labelling and the number of automorphisms.

The public functions take the keyword ``root_rule``, one of the integers
:data:`ROOTS_ALL`, :data:`ROOTS_HALF_EDGE_MAX`,
:data:`ROOTS_HALF_EDGE_FEWEST`, :data:`ROOTS_EDGE_MAX` and
:data:`ROOTS_EDGE_FEWEST`: the rule of the half-edges tried as roots.
Under :data:`ROOTS_ALL` every active half-edge is tried; under the other
rules each active half-edge gets an invariant (the degrees of its vertex
and of its face, and under the edge rules those of its twin as well), the
half-edges of equal invariant form a class, and only the half-edges of
one class are tried: the class of the largest invariant under ``*_MAX``,
the one with the fewest half-edges under ``*_FEWEST``, ties broken by the
largest invariant (:ref:`omap-root-rules`). The default is
``ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA`` (the integer ``10``): decorations
narrow the roots tried, and on an undecorated map it is
:data:`ROOTS_HALF_EDGE_FEWEST`, which tries fewer roots: against
:data:`ROOTS_ALL`, this pays on maps with many tied roots, at the price of
the invariants computed to choose the roots, which costs on maps where few
roots tie anyway. The flag :data:`ROOTS_DATA`, combined by ``|`` with one
of the four rules other than :data:`ROOTS_ALL`, appends the values of
``data`` to the invariants, so that fewer roots are tried on a decorated
map; the C code sorts each data array once per call to rank its values,
`O(n \log n)`. Any other value, a string or ``ROOTS_ALL | ROOTS_DATA``
included, raises a ``ValueError``. The C code computes the invariants of
the half-edges once per call and tries only the roots of the class of the
rule, on a disconnected map the class chosen at each step among the
half-edges not labelled yet. The automorphism group does not depend on the
rule: on a disconnected map, :func:`automorphism_generators` and the count
of :func:`canonical_key` come from the structure of the group above under
every rule.

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

The results agree with those of :mod:`combisurf.isomorphism_naive`, a pure
Python implementation they are compared with::

    sage: from combisurf import isomorphism_naive
    sage: fp = [2, 3, 4, 5, 0, 1]
    sage: c = [0, 1, 0, 1, 0, 1]
    sage: automorphisms(fp, classes=(c,)) == isomorphism_naive.automorphisms(fp, classes=(c,))
    True

Under every rule, the results, computed in C, agree with those of
``isomorphism_naive`` under the same rule::

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
``ROOTS_HALF_EDGE_MAX | ROOTS_DATA`` tries three roots instead of six, counted
by the private helpers of ``isomorphism_naive``, and its results agree
with those of ``isomorphism_naive``::

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
from cpython.number cimport PyNumber_Index
from cpython.ref cimport PyObject
from libc.limits cimport INT_MAX
from libc.string cimport memcpy, memcmp
from cysignals.signals cimport sig_on, sig_off, sig_block, sig_unblock
from combisurf.perm cimport _perm_check, _perm_id
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
    Return ``root_rule`` as a C ``int`` after checking that it is a valid
    rule of the roots.

    The valid rules are the integers ``0`` to ``4`` (:data:`ROOTS_ALL`,
    :data:`ROOTS_HALF_EDGE_MAX`, :data:`ROOTS_HALF_EDGE_FEWEST`,
    :data:`ROOTS_EDGE_MAX`, :data:`ROOTS_EDGE_FEWEST`) and ``9`` to ``12``
    (the last four combined with :data:`ROOTS_DATA` by ``|``), given as a
    Python ``int`` or as any other instance of :class:`numbers.Integral`
    except ``bool``. Any other value raises ``ValueError("root_rule must be
    one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX,
    ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12
    (the last four | ROOTS_DATA), got X")``, with ``X`` the ``repr`` of
    ``root_rule``.

    A Python ``int`` is checked first, without the lookup of
    :class:`numbers.Integral`: :func:`canonical_key` may be called millions
    of times in an exploration.

    The examples are those of :func:`check_root_rule`, which calls it.
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


def check_root_rule(root_rule):
    r"""
    Return ``root_rule`` as a Python ``int`` after checking that it is a
    valid rule of the roots.

    INPUT:

    - ``root_rule`` -- an object

    OUTPUT: ``root_rule`` as a Python ``int`` if it is one of the integers
    :data:`ROOTS_ALL`, :data:`ROOTS_HALF_EDGE_MAX`,
    :data:`ROOTS_HALF_EDGE_FEWEST`, :data:`ROOTS_EDGE_MAX` and
    :data:`ROOTS_EDGE_FEWEST` (``0`` to ``4``), or one of the last four
    combined with :data:`ROOTS_DATA` by ``|`` (``9`` to ``12``); any
    integer type is accepted except ``bool``. Any other value raises a
    ``ValueError``.

    EXAMPLES::

        sage: from combisurf.isomorphism import check_root_rule, ROOTS_ALL, ROOTS_EDGE_FEWEST, ROOTS_HALF_EDGE_MAX, ROOTS_DATA
        sage: check_root_rule(ROOTS_EDGE_FEWEST)
        4
        sage: type(check_root_rule(3)), check_root_rule(3)
        (<class 'int'>, 3)
        sage: type(check_root_rule(int(3))), check_root_rule(int(3))
        (<class 'int'>, 3)
        sage: type(check_root_rule(ROOTS_HALF_EDGE_MAX | ROOTS_DATA)), check_root_rule(ROOTS_HALF_EDGE_MAX | ROOTS_DATA)
        (<class 'int'>, 9)

    TESTS::

        sage: check_root_rule(5)
        Traceback (most recent call last):
        ...
        ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got 5
        sage: check_root_rule(ROOTS_ALL | ROOTS_DATA)
        Traceback (most recent call last):
        ...
        ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got 8
        sage: check_root_rule(13)
        Traceback (most recent call last):
        ...
        ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got 13
        sage: check_root_rule("all")
        Traceback (most recent call last):
        ...
        ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got 'all'
        sage: check_root_rule(1.0r)
        Traceback (most recent call last):
        ...
        ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got 1.0
        sage: check_root_rule(True)
        Traceback (most recent call last):
        ...
        ValueError: root_rule must be one of the integers 0 to 4 (ROOTS_ALL, ROOTS_HALF_EDGE_MAX, ROOTS_HALF_EDGE_FEWEST, ROOTS_EDGE_MAX, ROOTS_EDGE_FEWEST) or 9 to 12 (the last four | ROOTS_DATA), got True
    """
    return _check_root_rule(root_rule)


def _check_sizes(n, nd, nc):
    r"""
    Raise an ``OverflowError`` if a buffer of the C code for a map of ``n``
    half-edges decorated by ``nd`` arrays of ``data`` and ``nc`` arrays of
    ``classes`` has more entries than a C ``int`` counts.

    The sizes are those of the size macros of the C code
    (``DOMAP_ISO_KEY_SIZE``, ``DOMAP_CANONICAL_WORK_SIZE``,
    ``DOMAP_ISO_SEQUENCE_SIZE``, ``DOMAP_ISO_SEARCH_WORK_SIZE``,
    ``DOMAP_INVARIANTS_WORK_SIZE`` and ``DOMAP_STABILISER_WORK_SIZE``,
    documented in ``doc/c/omap.rst``), computed here in Python integers.

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
    invariants = 2 * n + 1
    stabiliser = 3 * nc * n + 5 * (n // 2) + 8 * n
    if max(n, key, work, sequence, search, invariants, stabiliser) > INT_MAX:
        raise OverflowError("map too large: {} half-edges with {} arrays of data and {} of classes".format(n, nd, nc))


cdef inline array.array _new(Py_ssize_t size):
    r"""
    Return a new array of typecode ``'i'`` with ``size`` entries, left
    uninitialised.

    The array is cloned from an empty module-level array of typecode
    ``'i'``, which skips the parsing of the typecode of ``array.array``; the
    caller writes every entry it reads.
    """
    return array.clone(_int_array, size, False)


ctypedef fused list_or_tuple:
    list
    tuple


cdef bint _ints_from(list_or_tuple s, const int *f, int *out, int n, bint edges) except -1:
    r"""
    Write the values of ``s`` as C ``int`` into ``out``, of ``n`` entries,
    and return whether every value converted is an integer (a Python ``int``
    or an object with ``__index__``) in the range of a C ``int``; on
    ``False``, ``out`` is partly written.

    Without ``edges``, ``s`` has ``n`` entries and ``out[h]`` is ``s[h]``.
    With ``edges``, ``s`` has ``n // 2`` entries, one per edge, and
    ``out[2*e]`` and ``out[2*e + 1]`` are ``s[e]``. If ``f`` is not ``NULL``,
    the values on the half-edges ``h`` with ``f[h] == -1`` (on the edges
    whose two halves have it, with ``edges``) are neither read nor checked
    and ``out`` holds ``0`` there.

    A value of exact type ``int`` is converted directly, any other one by
    ``PyNumber_Index`` (``operator.index``), so that ``float`` is refused
    and the integers of Sage are accepted. The ``TypeError`` or
    ``OverflowError`` of a value that is not an integer or is out of range
    is caught and turned into the return value ``False``: the caller then
    converts ``s`` by a slower path that raises its own exception. ``s``
    is read only; ``out`` is owned by the caller.
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
                    v = PyNumber_Index(v)
                out[2 * e] = out[2 * e + 1] = v
        else:
            for h in range(n):
                if f != NULL and f[h] == -1:
                    out[h] = 0
                    continue
                v = s[h]
                if type(v) is not int:
                    v = PyNumber_Index(v)
                out[h] = v
    except (OverflowError, TypeError):
        return False
    return True


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
    ``'i'`` of the right length whose values read are integers (with
    ``__index__``) in the range of a C ``int``. On any other input it
    returns ``None``, and
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

    The integers of Sage and the objects with ``__index__`` are converted
    like Python ``int``; values outside the range of a C ``int`` are left
    to the generic path::

        sage: _decoration_array([0] * 10, fp)
        array('i', [0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        sage: _decoration_array([1] * 10, fp)
        array('i', [1, 1, 1, 1, 0, 0, 1, 1, 1, 0])
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


def _constant_on_cycles(array.array d not None, array.array p not None):
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
        sage: _constant_on_cycles(None, p)
        Traceback (most recent call last):
        ...
        TypeError: Argument 'd' has incorrect type (expected array.array, got NoneType)
        sage: _constant_on_cycles(array('i', [5, 5, 7, 6]), None)
        Traceback (most recent call last):
        ...
        TypeError: Argument 'p' has incorrect type (expected array.array, got NoneType)
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


cdef inline int _num_active(const omap_struct *m):
    r"""
    Return the number of active half-edges of ``m``, those ``h`` with
    ``m.fp[h] != -1``.
    """
    cdef int h, k = 0
    for h in range(m.n):
        if m.fp[h] != -1:
            k += 1
    return k


cdef _check_decorated_oriented_map(fp, data, classes):
    r"""
    Check the shape of the map of face permutation ``fp`` decorated by
    ``data`` and ``classes``: the length of ``fp`` and the types and lengths
    of the decorations.

    INPUT:

    - ``fp`` -- a sequence with a length (the face permutation)

    - ``data``, ``classes`` -- the decorations, each expected to be an
      instance of ``tuple`` of sequences of length ``len(fp)``

    OUTPUT: ``None``. With ``n = len(fp)``, the checks, in this order, and
    the exceptions they raise:

    - ``n`` is even, else ``ValueError("fp must have even length")``;
    - ``data`` is an instance of ``tuple``, else ``TypeError("data must be
      a tuple")``;
    - every array of ``data`` has length ``n``, else ``ValueError("data
      arrays must have length n")``;
    - ``classes`` is an instance of ``tuple``, else ``TypeError("classes
      must be a tuple")``;
    - every array of ``classes`` has length ``n``, else
      ``ValueError("classes arrays must have length n")``.

    The values are not read: that ``fp`` is a permutation of its active
    half-edges and that the values of ``data`` are in the range of a C
    ``int`` are checked by ``_convert``, when the map is converted for the
    C code.

    The examples are those of :func:`check_decorated_oriented_map`, which
    calls it.
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


def check_decorated_oriented_map(fp, data, classes):
    r"""
    Check the length of the face permutation ``fp`` and the types and
    lengths of the decorations ``data`` and ``classes``.

    The values are not read: a face permutation that is not a permutation
    of its active half-edges, or a value of ``data`` outside the range of a
    C ``int``, passes this check. The public functions of this module, such
    as :func:`canonical_labelling`, raise those errors when they convert the
    map for the C code.

    INPUT:

    - ``fp`` -- a face permutation (a sequence of even length)

    - ``data`` -- a tuple of arrays of length ``len(fp)``

    - ``classes`` -- a tuple of arrays of length ``len(fp)``

    OUTPUT: ``None``, or a ``ValueError`` or ``TypeError`` on the first
    check that fails.

    EXAMPLES::

        sage: from combisurf.isomorphism import check_decorated_oriented_map
        sage: check_decorated_oriented_map([2, 5, 4, 1, 0, 3], ([0] * 6,), ([1, 2, 1, 2, 1, 2],))

    TESTS::

        sage: check_decorated_oriented_map([0, -1, 2], (), ())
        Traceback (most recent call last):
        ...
        ValueError: fp must have even length
        sage: check_decorated_oriented_map([0, -1], [[0, 0]], ())
        Traceback (most recent call last):
        ...
        TypeError: data must be a tuple
        sage: check_decorated_oriented_map([0, -1], ([0],), ())
        Traceback (most recent call last):
        ...
        ValueError: data arrays must have length 2
        sage: check_decorated_oriented_map([0, -1], (), [[0, 0]])
        Traceback (most recent call last):
        ...
        TypeError: classes must be a tuple
        sage: check_decorated_oriented_map([0, -1], (), ([0, 0, 0],))
        Traceback (most recent call last):
        ...
        ValueError: classes arrays must have length 2

    A face permutation that is not a permutation passes::

        sage: check_decorated_oriented_map([1, 1], (), ())
    """
    _check_decorated_oriented_map(fp, data, classes)


cdef array.array _convert(domap_struct *m, fp, tuple data, tuple classes):
    r"""
    Fill ``m`` with the map of face permutation ``fp`` decorated by
    ``data`` and ``classes`` and return the array that owns the buffers
    ``m`` points into.

    INPUT:

    - ``m`` -- the decorated map to fill, owned by the caller; its previous
      content is ignored

    - ``fp`` -- the face permutation, ``-1`` on the inactive half-edges: an
      array of typecode ``'i'``, a list, or any object that has a length
      (``len(fp)``) and that ``array.array('i', fp)`` accepts

    - ``data`` -- a tuple of ``nd`` arrays of ``n = len(fp)`` integers each,
      indexed by half-edges

    - ``classes`` -- a tuple of ``nc`` arrays of ``n`` hashable values
      each, indexed by half-edges

    The lengths are not checked here: ``n`` must be even and every array of
    ``data`` and ``classes`` must have ``n`` entries, as
    ``_check_decorated_oriented_map`` checks; a shorter array raises an
    ``IndexError`` or is read out of bounds.

    OUTPUT: an array of typecode ``'i'`` of ``(2 + nd + nc) * n`` entries:
    a copy of ``fp``, the vertex permutation computed from it by the C
    function ``_omap_vp_from_fp`` (``-1`` on the inactive half-edges), the
    ``nd`` arrays of ``data``, then the ``nc`` arrays of ``classes``
    renumbered ``0, 1, 2, ...`` by first appearance in half-edge order, each
    of ``n`` entries with ``0`` on the inactive half-edges. ``m.map.fp``,
    ``m.map.vp``, ``m.data`` and ``m.classes`` point into it and
    ``m.map.n``, ``m.nd`` and ``m.nc`` are set. ``fp`` is copied so that the
    buffers stay valid if the caller resizes its own array. The array is
    never resized, and ``m`` is valid as long as it is alive: the caller
    keeps a reference to it.

    Raise ``ValueError("fp must be a permutation of its active
    half-edges")`` if ``fp`` is not one, ``ValueError("fp has inconsistent
    length")`` if iterating ``fp`` does not give ``len(fp)`` values,
    ``TypeError`` if a value of ``data`` on an active half-edge is not an
    integer, ``OverflowError`` if such a value is outside the range of a C
    ``int``, and ``OverflowError("map too large: n half-edges with nd
    arrays of data and nc of classes")`` if a buffer of the C code for this
    map has more entries than a C ``int`` counts. The conversion of ``fp``
    by ``array.array`` raises its own ``TypeError`` or ``OverflowError``.
    """
    cdef int n, h, j, off, i, nd = len(data), nc = len(classes)
    cdef int *f
    cdef int *dec
    cdef int *a
    cdef array.array buf, seen, tmp
    cdef Py_ssize_t size = len(fp)
    # every buffer has at most 16 + 4 * (nd + nc) ints per half-edge (plus
    # a constant): below that bound the exact check in Python integers is
    # skipped, which keeps it off the path of small maps
    if size * (16 + 4 * (len(data) + len(classes))) + 64 > INT_MAX:
        _check_sizes(size, len(data), len(classes))
    n = <int> size
    buf = _new((2 + nd + nc) * size)
    f = buf.data.as_ints
    if type(fp) is array.array and (<array.array> fp).ob_descr.typecode == b'i':
        memcpy(f, (<array.array> fp).data.as_ints, n * sizeof(int))
    elif not (type(fp) is list and size <= INT_MAX
              and _ints_from(<list> fp, NULL, f, n, False)):
        # anything else, or a list holding a value that is not a C int,
        # goes through array.array, which raises its own exceptions
        tmp = array.array('i', fp)
        if len(tmp) != size:
            raise ValueError("fp has inconsistent length")
        memcpy(f, tmp.data.as_ints, n * sizeof(int))

    seen = _new(n)
    if _perm_check(f, n, seen.data.as_ints) != 1:
        raise ValueError("fp must be a permutation of its active half-edges")
    _omap_vp_from_fp(f + n, f, n)

    dec = f + 2 * n
    for j in range(nd):
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
            # anything else value by value through PyNumber_Index, which
            # raises TypeError on a value that is not an integer (the plain
            # conversion of Cython goes through __int__ and truncates a
            # float), and OverflowError on one out of the range of a C int
            for h in range(n):
                if f[h] == -1:
                    dec[off + h] = 0
                else:
                    dec[off + h] = PyNumber_Index(d[h])
    for j in range(nc):
        c = classes[j]
        off = (nd + j) * n
        # renumbered by first appearance in half-edge order, which leaves
        # the key unchanged (it renumbers them again in the order of the
        # new labels) and puts the values in range(n) as _domap_iso_key needs
        names = {}
        for h in range(n):
            if f[h] == -1:
                dec[off + h] = 0
            else:
                v = c[h]
                i = names.setdefault(v, len(names))
                dec[off + h] = i
    m.map.vp = f + n
    m.map.fp = f
    m.map.n = n
    m.nd = nd
    m.data = dec
    m.nc = nc
    m.classes = dec + nd * n
    return buf


cdef array.array _prepare(domap_struct *m, fp, data, classes):
    r"""
    Check and convert the map of face permutation ``fp`` decorated by
    ``data`` and ``classes`` into ``m``, and return the array that owns the
    buffers ``m`` points into.

    INPUT:

    - ``m`` -- the decorated map to fill, owned by the caller

    - ``fp``, ``data``, ``classes`` -- a face permutation and its
      decorations, as taken by the public functions of this module

    OUTPUT: the array of typecode ``'i'`` holding the face permutation, the
    vertex permutation and the decorations in C ``int``, which the caller
    keeps alive as long as it uses ``m``.

    With ``n = len(fp)``, raise first the exceptions of
    ``_check_decorated_oriented_map``, in this order:
    ``ValueError("fp must have even length")``, ``TypeError("data must be a
    tuple")``, ``ValueError("data arrays must have length n")``,
    ``TypeError("classes must be a tuple")`` and ``ValueError("classes
    arrays must have length n")``. Then those of ``_convert``:
    ``ValueError("fp must be a permutation of its active half-edges")`` if
    ``fp`` is not one, ``ValueError("fp has inconsistent length")`` if
    iterating ``fp`` does not give ``len(fp)`` values, ``TypeError`` if a
    value of ``data`` on an active half-edge is not an integer,
    ``OverflowError`` if such a value is outside the range of a C ``int``,
    and ``OverflowError("map too large: n half-edges with nd arrays of data
    and nc of classes")`` if a buffer of the C code for this map has more
    entries than a C ``int`` counts; the conversion of ``fp`` by
    ``array.array`` raises its own ``TypeError`` or ``OverflowError``.
    """
    _check_decorated_oriented_map(fp, data, classes)
    # _check_decorated_oriented_map accepts any instance of tuple, _convert
    # only the exact type
    return _convert(m, fp, tuple(data), tuple(classes))


cdef int _run(const domap_struct *m, int rule, array.array key, int *relabelling, int *roots, int *nroots) except? -2:
    r"""
    Compute the minimal key of the connected map ``m`` over the roots of the
    rule ``rule``, by the C function ``_domap_canonical``.

    INPUT:

    - ``m`` -- a decorated map of ``n = m.map.n`` half-edges

    - ``rule`` -- a rule of the roots (``0`` to ``4`` or ``9`` to ``12``,
      as returned by ``_check_root_rule``), which must be valid; another
      value is not checked

    - ``key`` -- an array of typecode ``'i'`` of at least
      ``DOMAP_ISO_KEY_SIZE(n, m.nd, m.nc)`` entries, owned by the caller

    - ``relabelling`` -- ``NULL`` or a buffer of ``n`` ints owned by the
      caller

    - ``roots`` -- ``NULL`` or a buffer of ``n`` ints owned by the caller

    - ``nroots`` -- ``NULL`` or a pointer to an ``int`` of the caller

    OUTPUT: the length of the minimal key, written into ``key`` (``0`` when
    ``m`` has no active half-edge), or ``DOMAP_ISO_DISCONNECTED`` when
    ``m`` is not connected, in which case no buffer is written. On a
    connected map, ``relabelling`` receives the relabelling from the first
    root of minimal key (``-1`` on the half-edges of the inactive edges),
    ``roots`` the roots of minimal key in increasing order and ``*nroots``
    their number, which is the number of automorphisms of ``m``.

    The work buffer of the C code is allocated here. Ctrl-C and the
    ``alarm`` of ``cysignals`` interrupt the C code and raise their
    exception (``KeyboardInterrupt``, ``AlarmInterrupt``).
    """
    cdef array.array work = _new(DOMAP_CANONICAL_WORK_SIZE(m.map.n, m.nd, m.nc))
    cdef int length
    # an interrupt leaves _domap_canonical by a longjmp back to sig_on: every
    # buffer it writes is owned here or by the caller, and it allocates none
    sig_on()
    length = _domap_canonical(key.data.as_ints, relabelling, roots, nroots, m,
                              <domap_iso_root_rule> rule, work.data.as_ints)
    sig_off()
    return length


cdef array.array _complete_relabelling(array.array r, int n, int ne):
    r"""
    Complete in place the relabelling ``r`` of a map with ``ne`` active
    edges into a permutation, by labelling its inactive edges in increasing
    order.

    INPUT:

    - ``r`` -- an array of typecode ``'i'`` of at least ``n`` entries: a
      relabelling that gives the labels ``0, ..., 2*ne - 1`` to the
      half-edges of the active edges and ``-1`` to both half-edges of each
      inactive edge; not checked

    - ``n`` -- the number of half-edges, even

    - ``ne`` -- the number of active edges

    OUTPUT: ``r`` itself, where the edges ``e`` with ``r[2*e] == -1``, in
    increasing order, get the labels ``2*ne, 2*ne + 1``, then ``2*ne + 2,
    2*ne + 3``, and so on.

    The examples are those of :func:`complete_relabelling`, which calls it.
    """
    cdef int e, k = 2 * ne
    cdef int *a = r.data.as_ints
    for e in range(n // 2):
        if a[2 * e] == -1:
            a[2 * e] = k
            a[2 * e + 1] = k + 1
            k += 2
    return r


def complete_relabelling(array.array r not None, int ne):
    r"""
    Complete in place the relabelling ``r`` of a map with ``ne`` active
    edges into a permutation, by labelling its inactive edges in increasing
    order.

    INPUT:

    - ``r`` -- an array of typecode ``'i'`` of even length: a relabelling
      that gives the labels ``0, ..., 2*ne - 1`` to the half-edges of the
      active edges and ``-1`` to both half-edges of each inactive edge (not
      checked)

    - ``ne`` -- the number of active edges

    OUTPUT: ``r`` itself, where the edges ``e`` with ``r[2*e] == -1``, in
    increasing order, get the labels ``2*ne, 2*ne + 1``, then ``2*ne + 2,
    2*ne + 3``, and so on.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.isomorphism import complete_relabelling
        sage: r = array('i', [2, 3, -1, -1, 0, 1])
        sage: complete_relabelling(r, 2)
        array('i', [2, 3, 4, 5, 0, 1])
        sage: complete_relabelling(r, 2) is r
        True

    TESTS::

        sage: complete_relabelling(array('l', [2, 3, -1, -1, 0, 1]), 2)
        Traceback (most recent call last):
        ...
        TypeError: r must be an array of typecode 'i'
        sage: complete_relabelling(None, 2)
        Traceback (most recent call last):
        ...
        TypeError: Argument 'r' has incorrect type (expected array.array, got NoneType)
    """
    if r.ob_descr.typecode != b'i':
        raise TypeError("r must be an array of typecode 'i'")
    return _complete_relabelling(r, len(r), ne)


cdef int _search(const domap_struct *m, int rule, bint all_leaves, array.array sequence, int *relabelling,
                 domap_iso_leaf_callback callback, void *arg) except -1:
    r"""
    Search the minimal sequence of blocks of the map ``m``, connected or
    not, over its labellings component after component, by the C function
    ``_domap_iso_search``.

    The block of a component labelled from a root is its number of labels,
    its face permutation in these local labels, then its decorations read
    in the same order (the values of ``classes`` renumbered by first
    appearance, continuing the numbering of the components labelled
    before). At each step the roots of the rule ``rule`` among the
    half-edges not labelled yet are tried, and only those of minimal block
    are followed, depth first.

    INPUT:

    - ``m`` -- a decorated map of ``n = m.map.n`` half-edges

    - ``rule`` -- a valid rule of the roots, as returned by
      ``_check_root_rule``

    - ``all_leaves`` -- whether to follow every tied root, so that
      ``callback`` sees every leaf reaching the minimum (one per
      automorphism); otherwise tied roots exchanged by an automorphism are
      skipped

    - ``sequence`` -- an array of typecode ``'i'`` of at least
      ``DOMAP_ISO_SEQUENCE_SIZE(n, m.nd, m.nc)`` entries, owned by the
      caller

    - ``relabelling`` -- ``NULL`` or a buffer of ``n`` ints owned by the
      caller, which receives the relabelling of the first leaf reaching the
      minimal sequence (``-1`` on the half-edges of the inactive edges)

    - ``callback``, ``arg`` -- ``NULL``, or a function called with ``arg``
      on every leaf smaller than the best one so far (with ``first`` set)
      and, with ``all_leaves``, on every leaf equal to it (with ``first``
      unset); the leaf buffer is valid during the call only

    OUTPUT: the length of the minimal sequence, written into ``sequence``.

    A Python exception raised by ``callback`` (its return value ``-1``)
    stops the search and propagates; the outputs are then undefined. The
    work buffer of the C code is allocated here. Ctrl-C and the ``alarm``
    of ``cysignals`` interrupt the search and raise their exception.
    """
    cdef array.array work = _new(DOMAP_ISO_SEARCH_WORK_SIZE(m.map.n, m.nd, m.nc))
    cdef int length
    # an interrupt leaves _domap_iso_search by a longjmp back to sig_on: every
    # buffer it writes is owned here or by the caller, and it allocates none;
    # the longjmp also resets the count of sig_on. An exception of the
    # callback returns from _domap_iso_search normally, hence the sig_off in
    # the finally
    sig_on()
    try:
        _domap_iso_search(sequence.data.as_ints, &length, relabelling, m, <domap_iso_root_rule> rule,
                          all_leaves, callback, arg, work.data.as_ints)
    finally:
        sig_off()
    return length


cdef int _key_of(const domap_struct *m, array.array r, array.array key) except -1:
    r"""
    Write into ``key`` the key of the relabelling ``r`` of ``m`` and return
    its length.

    INPUT:

    - ``m`` -- a decorated map of ``n = m.map.n`` half-edges and ``ne``
      active edges

    - ``r`` -- an array of typecode ``'i'`` of ``n`` entries giving the
      labels ``0, ..., 2*ne - 1`` to the half-edges of the active edges and
      ``-1`` to those of the inactive edges, such as the relabelling of a
      leaf of ``_domap_iso_search``; read only

    - ``key`` -- an array of typecode ``'i'`` of at least
      ``DOMAP_ISO_KEY_SIZE(n, m.nd, m.nc)`` entries, owned by the caller

    OUTPUT: the length of the key, ``2*ne + (m.nd + m.nc) * na`` for ``na``
    active half-edges: the face permutation in the new labels (``-1`` at
    the labels of the inactive halves of folded edges), then each array of
    ``data``, then each array of ``classes`` renumbered by first
    appearance, read in the order of the new labels of the active
    half-edges.
    """
    cdef int n = m.map.n, k = 2 * omap_num_edges(&m.map), h
    cdef const int *f = m.map.fp
    cdef int *a = r.data.as_ints
    cdef array.array rinv = _new(n)
    cdef array.array names = _new(n)
    # the labels of the inactive halves of folded edges keep -1
    for h in range(n):
        rinv.data.as_ints[h] = -1
        names.data.as_ints[h] = -1
    for h in range(n):
        if f[h] != -1:
            rinv.data.as_ints[a[h]] = h
    return _domap_iso_key(key.data.as_ints, m, a, rinv.data.as_ints, k, names.data.as_ints)


# State of the leaf callback of ``_domap_iso_search`` for ``automorphisms``: ``ra``
# is the first leaf reaching the best sequence so far, ``result`` the list of
# the automorphisms ``g[h] = rb^{-1}[ra[h]]`` over the leaves ``rb`` reaching
# it, and ``rbinv`` a work array of ``n`` entries. ``ra``, ``rbinv`` and the
# list are owned by the caller of ``_domap_iso_search``, which outlives the search;
# ``result`` is borrowed and the list is emptied in place, never rebound.
cdef struct _leaf_state:
    int n
    int *ra
    int *rbinv
    PyObject *result


ctypedef int (*_unchecked_leaf_callback)(const int *leaf, int first, void *arg) noexcept


cdef int _automorphism_leaf(const int *leaf, int first, void *arg) except -1:
    r"""
    Append the automorphism of the leaf ``leaf`` to the list of the
    ``_leaf_state`` ``arg``, emptied first when ``first`` is set: the leaf
    callback of ``_domap_iso_search`` for :func:`automorphisms`.

    The automorphism is ``g[h] = leaf^{-1}[ra[h]]``, with ``ra`` the first
    leaf reaching the best sequence so far (stored into the state when
    ``first`` is set), and ``g[h] = h`` on the half-edges of the inactive
    edges. Return ``0``, or ``-1`` with a Python exception set (a failed
    allocation), which stops the search.

    The callback runs inside the ``sig_on`` of ``_search``, where an
    interrupt jumps out at once, and a jump out of the Python code of
    ``_collect_automorphism`` (an allocation, a list append) could leave the
    interpreter in an inconsistent state. ``sig_block`` delays the interrupt
    to the ``sig_unblock``, reached when ``_collect_automorphism`` has
    returned and released its references, so that the jump skips no
    reference count. A ``sig_check`` in the callback, with no ``sig_on`` on
    this path, would instead leave the search uninterruptible between two
    leaves reaching the best sequence, which can be the whole search. The
    call goes through a pointer without exception check, so that the
    ``sig_unblock`` also runs when ``_collect_automorphism`` raises; its
    ``-1``, with the exception set, is then returned to ``_domap_iso_search``.
    """
    cdef int status
    sig_block()
    status = (<_unchecked_leaf_callback> _collect_automorphism)(leaf, first, arg)
    sig_unblock()
    return status


cdef int _collect_automorphism(const int *leaf, int first, void *arg) except -1:
    r"""
    Append to the list of the ``_leaf_state`` ``arg`` the automorphism of the
    leaf ``leaf``, the list emptied first when ``first`` is set.

    INPUT:

    - ``leaf`` -- the relabelling of a leaf (``s.n`` ints, ``-1`` on the
      half-edges of the inactive edges), read during the call only

    - ``first`` -- whether the leaf is smaller than every leaf before; it
      is then copied into ``s.ra`` and the list is emptied in place

    - ``arg`` -- a pointer to the ``_leaf_state`` ``s``, whose buffers
      ``ra`` and ``rbinv`` of ``s.n`` ints and list are owned by the caller
      of the search

    OUTPUT: ``0``, after appending to the list the new array ``g`` of
    typecode ``'i'`` with ``g[h] = leaf^{-1}[ra[h]]``, and ``g[h] = h``
    where ``ra[h] == -1``; ``rbinv`` is overwritten with the inverse of
    ``leaf``. A failed allocation raises ``MemoryError``.
    """
    cdef _leaf_state *s = <_leaf_state *> arg
    cdef int n = s.n, h
    cdef int *ra = s.ra
    cdef int *rbinv = s.rbinv
    cdef list result = <list> s.result
    cdef array.array g
    if first:
        memcpy(ra, leaf, n * sizeof(int))
        del result[:]
    for h in range(n):
        if leaf[h] != -1:
            rbinv[leaf[h]] = h
    g = _new(n)
    for h in range(n):
        if ra[h] == -1:
            g.data.as_ints[h] = h
        else:
            g.data.as_ints[h] = rbinv[ra[h]]
    result.append(g)
    return 0


def canonical_labelling(fp, data=(), classes=(), root_rule=ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA):
    r"""
    Return the canonical labelling of the map with face permutation ``fp``
    decorated by ``data`` and ``classes``.

    On a connected map, the canonical labelling is the relabelling by the
    breadth-first walk (:ref:`omap-relabelling`) from the first root, in
    increasing order among the roots of the rule, whose key is minimal
    (:ref:`omap-keys`). On a disconnected map, it is the labelling of the
    first leaf of the search, component after component, that reaches the
    minimal sequence of blocks (:ref:`omap-blocks`). The inactive edges then
    get the labels left, in increasing order.

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
      tried as roots (see :ref:`omap-root-rules`); any other
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
    cdef domap_t m
    cdef array.array _buf = _prepare(m, fp, data, classes)
    cdef array.array key = _new(DOMAP_ISO_KEY_SIZE(m.map.n, m.nd, m.nc))
    cdef array.array r = _new(m.map.n)
    if _run(m, rule, key, r.data.as_ints, NULL, NULL) == DOMAP_ISO_DISCONNECTED:
        _search(m, rule, False, _new(DOMAP_ISO_SEQUENCE_SIZE(m.map.n, m.nd, m.nc)), r.data.as_ints, NULL, NULL)
    return _complete_relabelling(r, m.map.n, omap_num_edges(&m.map))


cdef tuple _form(const domap_struct *m, array.array key, array.array r, tuple data):
    r"""
    Return the canonical form of ``m`` from its minimal key ``key`` and the
    relabelling ``r`` that reaches it.

    INPUT:

    - ``m`` -- a decorated map of ``n = m.map.n`` half-edges, ``ne`` active
      edges and ``na`` active half-edges

    - ``key`` -- the minimal key of ``m``, an array of typecode ``'i'``:
      the face permutation in the canonical labels (``2*ne`` entries), the
      ``m.nd`` arrays of ``data`` then the ``m.nc`` renumbered arrays of
      ``classes`` (``na`` entries each); read only

    - ``r`` -- the canonical relabelling, an array of typecode ``'i'`` of
      ``n`` entries (``-1`` on the half-edges of the inactive edges); read
      only

    - ``data`` -- the arrays of ``data`` as given by the caller

    OUTPUT: the tuple ``(fp_c, d_1, ..., d_nd, c_1, ..., c_nc)`` of tuples
    returned by :func:`canonical_form`. ``fp_c`` and the ``c_i`` are read
    from ``key``; each ``d_i`` is read from the original array of ``data``
    at the half-edges in the order of their canonical labels, so that its
    values keep their Python type (an integer of Sage stays one) instead of
    the C ``int`` of the key.
    """
    cdef int k = 2 * omap_num_edges(&m.map)
    cdef int na = _num_active(&m.map)
    cdef int h, i, pos
    cdef int *kd = key.data.as_ints
    rinv = [-1] * k
    for h in range(m.map.n):
        if m.map.fp[h] != -1:
            rinv[r.data.as_ints[h]] = h
    out = [tuple([kd[i] for i in range(k)])]
    for d in data:
        out.append(tuple([d[x] for x in rinv if x != -1]))
    pos = k + m.nd * na
    while len(out) < 1 + m.nd + m.nc:
        out.append(tuple([kd[i] for i in range(pos, pos + na)]))
        pos += na
    return tuple(out)


def canonical_form(fp, data=(), classes=(), root_rule=ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA):
    r"""
    Return the canonical form of the map with face permutation ``fp``
    decorated by ``data`` and ``classes``.

    The canonical form is the map read in its canonical labels (see
    :func:`canonical_labelling`): the key of the canonical labelling
    (:ref:`omap-keys`), split into the face permutation and the
    decorations. On a connected map this key is the minimal key over the
    roots of the rule; on a disconnected map it is the key of the labelling
    of the first leaf reaching the minimal sequence of blocks
    (:ref:`omap-blocks`). The values of ``data`` are those of the arrays
    given, so that they keep their Python type.

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
    cdef domap_t m
    cdef array.array _buf = _prepare(m, fp, data, classes)
    cdef array.array key = _new(DOMAP_ISO_KEY_SIZE(m.map.n, m.nd, m.nc))
    cdef array.array r = _new(m.map.n)
    if _run(m, rule, key, r.data.as_ints, NULL, NULL) == DOMAP_ISO_DISCONNECTED:
        _search(m, rule, False, _new(DOMAP_ISO_SEQUENCE_SIZE(m.map.n, m.nd, m.nc)), r.data.as_ints, NULL, NULL)
        _key_of(m, r, key)
    return _form(m, key, r, data)


cdef array.array _automorphism_from(const domap_struct *m, array.array ra, int root, array.array rb,
                                    array.array rbinv, array.array queue):
    r"""
    Return the automorphism of the connected map ``m`` sending the half-edge
    of label ``i`` in the relabelling ``ra`` to the half-edge of label ``i``
    in the relabelling from ``root``.

    INPUT:

    - ``m`` -- a connected decorated map of ``n = m.map.n`` half-edges

    - ``ra`` -- the relabelling from a root of minimal key, an array of
      typecode ``'i'`` of ``n`` entries (``-1`` on the half-edges of the
      inactive edges); read only

    - ``root`` -- an active half-edge whose relabelling gives the same key
      as ``ra``; not checked

    - ``rb``, ``rbinv``, ``queue`` -- work arrays of typecode ``'i'`` of
      ``n`` entries each, owned by the caller and overwritten

    OUTPUT: a new array ``g`` of typecode ``'i'`` of ``n`` entries,
    ``g[h] = rb^{-1}[ra[h]]`` with ``rb`` the relabelling from ``root``,
    and ``g[h] = h`` on the half-edges of the inactive edges.
    """
    cdef int n = m.map.n, h
    cdef array.array g = _new(n)
    # _omap_relabelling_from writes -1 into rbinv at the labels of the
    # inactive halves of folded edges; the automorphism maps those
    # halves too, so the inverse is rebuilt on every labelled half-edge
    _omap_relabelling_from(rb.data.as_ints, rbinv.data.as_ints, m.map.fp,
                           m.map.n, root, queue.data.as_ints)
    for h in range(n):
        if rb.data.as_ints[h] != -1:
            rbinv.data.as_ints[rb.data.as_ints[h]] = h
    for h in range(n):
        if ra.data.as_ints[h] == -1:
            g.data.as_ints[h] = h
        else:
            g.data.as_ints[h] = rbinv.data.as_ints[ra.data.as_ints[h]]
    return g


cdef list _leaf_automorphisms(const domap_struct *m, int rule, array.array ra):
    r"""
    Return the list of the automorphisms of the map ``m``, connected or not,
    from the leaves of the search of ``_search`` under the rule ``rule``.

    Every tied root is followed, so the leaves reaching the minimal
    sequence of blocks are as many as the automorphisms; with ``ra`` the
    first of them, the leaf ``rb`` gives the automorphism ``g[h] =
    rb^{-1}[ra[h]]`` (``g[h] = h`` on the half-edges of the inactive edges).

    INPUT:

    - ``m`` -- a decorated map of ``n = m.map.n`` half-edges

    - ``rule`` -- a valid rule of the roots, as returned by
      ``_check_root_rule``

    - ``ra`` -- a work array of typecode ``'i'`` of ``n`` entries, owned by
      the caller and overwritten with the first leaf

    OUTPUT: a list of new arrays of typecode ``'i'``, the automorphisms in
    the order of their leaves, the identity first. Ctrl-C and the ``alarm``
    of ``cysignals`` interrupt the search and raise their exception.
    """
    cdef _leaf_state s
    cdef array.array rbinv = _new(m.map.n)
    cdef list result = []
    s.n = m.map.n
    s.ra = ra.data.as_ints
    s.rbinv = rbinv.data.as_ints
    s.result = <PyObject *> result
    _search(m, rule, True, _new(DOMAP_ISO_SEQUENCE_SIZE(m.map.n, m.nd, m.nc)), NULL, _automorphism_leaf,
            &s)
    return result


cdef list _connected_generators(const domap_struct *m, array.array ra, array.array roots, int count):
    r"""
    Return the generators of the automorphism group of the connected map
    ``m`` from its roots of minimal key.

    The automorphism of the root ``b_j`` sends the half-edge of label ``i``
    from ``b_0`` to the one of label ``i`` from ``b_j``; the roots are
    scanned in increasing order, and the automorphism of ``b_j`` is built
    and kept only when ``b_j`` is outside the orbit of ``b_0`` under the
    generators kept (a union-find over the half-edges). An automorphism of a
    connected map is determined by the image of ``b_0``, so the generators
    kept generate the group.

    INPUT:

    - ``m`` -- a connected decorated map of ``n = m.map.n`` half-edges

    - ``ra`` -- the relabelling from the first root ``b_0``, an array of
      typecode ``'i'`` of ``n`` entries; read only

    - ``roots`` -- the roots of minimal key ``b_0 < b_1 < ...`` in its first
      ``count`` entries, an array of typecode ``'i'``; read only

    - ``count`` -- the number of roots of minimal key

    OUTPUT: a list of new arrays of typecode ``'i'`` of ``n`` entries, empty
    when ``count <= 1``.
    """
    cdef int n = m.map.n, j, h, b0
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
        if uf_find_halving(p, r[j]) == uf_find_halving(p, b0):
            continue
        g = _automorphism_from(m, ra, r[j], rb, rbinv, queue)
        for h in range(n):
            uf_link(p, uf_find_halving(p, h),
                    uf_find_halving(p, g.data.as_ints[h]))
        gens.append(g)
    return gens


cdef array.array _identity(int n):
    r"""
    Return the identity permutation of ``range(n)`` as a new array of
    typecode ``'i'``.
    """
    cdef array.array g = _new(n)
    _perm_id(g.data.as_ints, n)
    return g


cdef list _wreath(const domap_struct *m, bint lift, bint labels=False):
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
    - ``local`` -- the face permutation of the first copy as a map of its
      own (an ``array('i')`` in the local indices of ``halves``) if
      ``labels``, ``None`` otherwise
    - ``roots`` -- the canonical labels of the roots of minimal key of the
      first copy, in increasing order (an ``array('i')`` starting with
      ``0``; its automorphism sending the half-edge of label ``0`` to the
      one of label ``l`` exists exactly for these ``l``) if ``labels``,
      ``None`` otherwise

    The group is the product over the kinds of the wreath product of the
    group of a copy by the symmetric group on the copies.

    Each component is labelled under :data:`ROOTS_ALL` by ``_run``, so the
    structure does not depend on the rule of the caller. A component found
    disconnected raises ``RuntimeError("a component of a map is not
    connected")``, which signals a bug: the components are the connected
    components of ``m``. Ctrl-C and the ``alarm`` of ``cysignals``
    interrupt the C code and raise their exception.
    """
    cdef int n = m.map.n, c, h, x, y, i, j, size, cnt, length, top
    cdef int ncomp = 0
    cdef const int *f = m.map.fp
    cdef const int *dec = m.data
    cdef array.array comp = _new(n)
    cdef array.array stack = _new(n)
    cdef array.array idx = _new(n)
    cdef array.array halves_c, lf, key, ra, roots, g, gl, _bc
    cdef domap_t mc
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
        _bc = _convert(mc, lf, ld, ())
        key = _new(DOMAP_ISO_KEY_SIZE(size, m.nd, 0))
        ra = _new(size)
        roots = _new(size)
        length = _run(mc, _ROOTS_ALL, key, ra.data.as_ints, roots.data.as_ints, &cnt)
        if length == DOMAP_ISO_DISCONNECTED:
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
            kinds.append(([(halves_c, ra)], gens, cnt, lf,
                          array.array('i', sorted([ra.data.as_ints[roots.data.as_ints[i]] for i in range(cnt)]))))
        else:
            kinds.append(([(halves_c, ra)], gens, cnt, None, None))
    return kinds


cdef int _place(array.array g, tuple a, tuple b) except -1:
    r"""
    Write into ``g`` the isomorphism from the component ``a`` to the
    isomorphic component ``b``, which sends the half-edge of canonical label
    ``i`` in ``a`` to the one of canonical label ``i`` in ``b``.

    INPUT:

    - ``g`` -- an array of typecode ``'i'`` indexed by the half-edges of
      the map, owned by the caller; only its entries at the half-edges of
      ``a`` are written

    - ``a``, ``b`` -- two isomorphic components, each a pair ``(halves,
      r)`` of arrays of typecode ``'i'`` of the same size: ``halves`` the
      half-edges of the component in increasing order and ``r`` its
      canonical relabelling as a map of its own (``r[i]`` the label of
      ``halves[i]``), as listed by ``_wreath``

    OUTPUT: ``0``; the input is not checked.
    """
    cdef array.array ha = a[0], ra = a[1], hb = b[0], rb = b[1]
    cdef int size = len(ha), i
    cdef array.array rinv = _new(size)
    for i in range(size):
        rinv.data.as_ints[rb.data.as_ints[i]] = i
    for i in range(size):
        g.data.as_ints[ha.data.as_ints[i]] = hb.data.as_ints[rinv.data.as_ints[ra.data.as_ints[i]]]
    return 0


cdef list _wreath_generators(const domap_struct *m):
    r"""
    Return the generators of the automorphism group of the disconnected map
    ``m`` with its arrays of ``data`` alone, its arrays of ``classes`` being
    ignored.

    The group is the product, over the kinds of isomorphic components, of
    the wreath product `\mathrm{Aut}(C) \wr S_k` of the group of a copy
    `C` by the symmetric group on its `k` copies. The kinds are ordered by
    the smallest half-edge of their first copy, the copies of a kind
    likewise. For each kind the generators are those of the group of its
    first copy (identity elsewhere), then the transposition of its first
    two copies when `k \geq 2`, then the cycle of its copies when `k \geq
    3`; a copy goes to another by the isomorphism sending the half-edge of
    canonical label ``i`` to the one of label ``i``, the identity on the
    other components.

    OUTPUT: a list of new arrays of typecode ``'i'`` of ``m.map.n`` entries,
    the identity on the half-edges of the inactive edges.

    .. SEEALSO::

        ``_wreath``, which computes the kinds and their copies.
    """
    cdef array.array g
    cdef int t, k
    gens = []
    for copies, aut, _count, _local, _roots in _wreath(m, True):
        gens.extend(aut)
        k = len(copies)
        if k >= 2:
            g = _identity(m.map.n)
            _place(g, copies[0], copies[1])
            _place(g, copies[1], copies[0])
            gens.append(g)
        if k >= 3:
            g = _identity(m.map.n)
            for t in range(k):
                _place(g, copies[t], copies[(t + 1) % k])
            gens.append(g)
    return gens


cdef _wreath_order(const domap_struct *m):
    r"""
    Return the number of automorphisms of the disconnected map ``m`` with
    its arrays of ``data`` alone, its arrays of ``classes`` being ignored,
    as a Python integer.

    The group is the product, over the kinds of `k` isomorphic components
    `C`, of the wreath products `\mathrm{Aut}(C) \wr S_k`, so its order is
    the product of `|\mathrm{Aut}(C)|^k k!` over the kinds; it is computed
    from the order of the group of one copy of each kind, without listing
    the automorphisms, and may exceed every C integer type.

    .. SEEALSO::

        ``_wreath``, which computes the kinds and the orders of their
        groups.
    """
    order = 1
    for copies, _aut, count, _local, _roots in _wreath(m, False):
        order *= count ** len(copies) * factorial(len(copies))
    return order


# State of the generator callback of ``_domap_stabiliser`` for
# ``automorphism_generators``: ``n`` the number of half-edges and ``result``
# the list of the generators, owned by the caller of ``_domap_stabiliser``,
# which outlives the search; ``result`` is borrowed and the list is appended
# in place, never rebound.
cdef struct _generator_state:
    int n
    PyObject *result


ctypedef int (*_unchecked_generator_callback)(const int *g, void *arg) noexcept


cdef int _stabiliser_generator(const int *g, void *arg) except -1:
    r"""
    Append a copy of the generator ``g`` to the list of the
    ``_generator_state`` ``arg``: the generator callback of
    ``_domap_stabiliser`` for :func:`automorphism_generators`.

    ``g`` holds ``n`` ints (``g[h]`` the image of ``h``, with ``n`` the
    field of the state) and is valid during the call only, hence the copy.
    Return ``0``, or ``-1`` with a Python exception set (a failed
    allocation), which stops the search.

    The callback runs inside the ``sig_on`` of ``_stabiliser``, where an
    interrupt jumps out at once, and a jump out of the Python code of
    ``_collect_generator`` (an allocation, a list append) could leave the
    interpreter in an inconsistent state. ``sig_block`` delays the interrupt
    to the ``sig_unblock``, reached when ``_collect_generator`` has returned
    and released its references, so that the jump skips no reference count.
    A ``sig_check`` in the callback, with no ``sig_on`` on this path, would
    instead leave the search uninterruptible between two generators, which
    can be the whole search. The call goes through a pointer without
    exception check, so that the ``sig_unblock`` also runs when
    ``_collect_generator`` raises; its ``-1``, with the exception set, is
    then returned to ``_domap_stabiliser``.
    """
    cdef int status
    sig_block()
    status = (<_unchecked_generator_callback> _collect_generator)(g, arg)
    sig_unblock()
    return status


cdef int _collect_generator(const int *g, void *arg) except -1:
    r"""
    Append to the list of the ``_generator_state`` ``arg`` a copy of the
    generator ``g``.

    INPUT:

    - ``g`` -- a permutation of the half-edges (``s.n`` ints, ``g[h]`` the
      image of ``h``), read during the call only

    - ``arg`` -- a pointer to the ``_generator_state`` ``s``, whose list is
      owned by the caller of the search and appended in place

    OUTPUT: ``0``, after appending a new array of typecode ``'i'`` holding
    the ``s.n`` entries of ``g``. A failed allocation raises
    ``MemoryError``.
    """
    cdef _generator_state *s = <_generator_state *> arg
    cdef array.array c = _new(s.n)
    memcpy(c.data.as_ints, g, s.n * sizeof(int))
    (<list> s.result).append(c)
    return 0


cdef tuple _stabiliser(const domap_struct *m, bint generators):
    r"""
    Return the generators (an empty list unless ``generators``) and the
    order (a Python integer) of the group of the automorphisms of the
    disconnected map ``m`` that preserve the partitions defined by its
    arrays of ``classes``.

    The group `H` is the subgroup, preserving the partitions, of the group
    of ``m`` with its arrays of ``data`` alone, a product over the kinds of
    isomorphic components of wreath products. The components are laid out
    for the C function ``_domap_stabiliser`` (their kinds, their half-edges
    by canonical label, the face permutation, canonical relabelling and
    labels of the roots of minimal key of the first copy of each kind), in
    arrays owned here and kept unresized while the C code reads them.
    ``_domap_stabiliser`` runs a subgroup search of Sims with one base point
    per component, the half-edge of canonical label ``0``, and writes the
    size of the orbit of each base point under the stabiliser of the
    previous ones in `H`: the order of `H` is their product, computed in
    Python integers since it may exceed every C integer type. With
    ``generators``, each generator found is copied into the list by the
    callback ``_stabiliser_generator``.

    OUTPUT: a pair ``(gens, order)``: ``gens`` a list of new arrays of
    typecode ``'i'`` of ``m.map.n`` entries, the identity on the half-edges
    of the inactive edges, and ``order`` a Python integer.

    The work buffer of the C code is allocated here. Ctrl-C and the
    ``alarm`` of ``cysignals`` interrupt the search and raise their
    exception.
    """
    cdef int n = m.map.n, size, i
    cdef array.array halves_c, r
    cdef array.array kind = _new(0), off = _new(0), lab = _new(0)
    cdef array.array start = array.array('i', [0]), coff = _new(0), csize = _new(0)
    cdef array.array cfp = _new(0), cra = _new(0), crainv = _new(0)
    cdef array.array loff = array.array('i', [0]), labels = _new(0)
    cdef array.array orbits, work
    cdef domap_wreath_t w
    cdef _generator_state state
    cdef domap_stabiliser_callback callback = NULL
    for t, (copies, _aut, _count, local, roots) in enumerate(_wreath(m, False, True)):
        size = len(copies[0][0])
        coff.append(len(cfp))
        csize.append(size)
        cfp.extend(local)
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
    # the arrays are complete: they are not resized while w points into them
    w.ncomp = len(kind)
    w.kind = kind.data.as_ints
    w.off = off.data.as_ints
    w.lab = lab.data.as_ints
    w.start = start.data.as_ints
    w.coff = coff.data.as_ints
    w.csize = csize.data.as_ints
    w.cfp = cfp.data.as_ints
    w.cra = cra.data.as_ints
    w.crainv = crainv.data.as_ints
    w.loff = loff.data.as_ints
    w.labels = labels.data.as_ints
    orbits = _new(w.ncomp)
    work = _new(DOMAP_STABILISER_WORK_SIZE(n, m.nc))
    gens = []
    state.n = n
    state.result = <PyObject *> gens
    if generators:
        callback = _stabiliser_generator
    # as in _search: an interrupt leaves _domap_stabiliser by a longjmp back
    # to sig_on, and an exception of the callback returns from it normally,
    # hence the sig_off in the finally
    sig_on()
    try:
        _domap_stabiliser(orbits.data.as_ints, m, w, callback, &state, work.data.as_ints)
    finally:
        sig_off()
    order = 1
    for i in range(w.ncomp):
        order *= orbits.data.as_ints[i]
    return gens, order


def automorphisms(fp, data=(), classes=(), root_rule=ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA):
    r"""
    Return the automorphisms of the map with face permutation ``fp``
    decorated by ``data`` and ``classes``.

    The automorphisms come in the order in which the search meets them. On
    a connected map, with ``b_0 < b_1 < ...`` the roots of the rule of
    minimal key (:ref:`omap-keys`), there is one automorphism per root, in
    increasing order of the roots, the identity (from ``b_0``) first. On a
    disconnected map, the search component after component
    (:ref:`omap-blocks`) follows every tied root, and its leaves reaching
    the minimal sequence of blocks are as many as the automorphisms: with
    ``r_0`` the first of them, the leaf ``r`` gives the automorphism
    sending the half-edge of label ``i`` in ``r_0`` to the one of label
    ``i`` in ``r``, in the order of the leaves, the identity first.

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
    cdef domap_t m
    cdef array.array _buf = _prepare(m, fp, data, classes)
    cdef int n = m.map.n
    cdef int count, j
    cdef array.array key = _new(DOMAP_ISO_KEY_SIZE(n, m.nd, m.nc))
    cdef array.array ra = _new(n)
    cdef array.array roots = _new(n)
    cdef array.array rb, rbinv, queue
    if _run(m, rule, key, ra.data.as_ints, roots.data.as_ints, &count) == DOMAP_ISO_DISCONNECTED:
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

    On a connected map, the generators are automorphisms of the list of
    :func:`automorphisms`, kept in its order when they are not in the group
    generated by those kept before. On a disconnected map, they come from
    the structure of the group: the generators of the group of a copy of
    each kind of isomorphic components and the permutations of its copies
    without ``classes``, the generators found by a subgroup search with
    ``classes`` (see below).

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
      the other. There are at most two generators per kind beside those of
      `\mathrm{Aut}(C)`. The components are labelled under
      :data:`ROOTS_ALL` whatever the rule.

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
      generators come in the order found. After a placement of a component
      mapping no new value of ``classes`` failed, the later placements of
      that component mapping no new value are skipped, an automorphism of
      the copy exchanging the two images. Tied copies whose values are all
      mapped after the first components (`k` isolated edges with one or two
      class values in any pattern) then take a time polynomial in `k`. The
      search stays exponential in the worst case, when the placements map
      new values that the values of ``classes`` reject only deep in the
      search, as on a graph on the class values given by its edges as
      isolated edges.

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
    cdef domap_t m
    cdef array.array _buf = _prepare(m, fp, data, classes)
    cdef int n = m.map.n
    cdef int count
    cdef array.array key = _new(DOMAP_ISO_KEY_SIZE(n, m.nd, m.nc))
    cdef array.array ra = _new(n)
    cdef array.array roots = _new(n)
    if _run(m, rule, key, ra.data.as_ints, roots.data.as_ints, &count) != DOMAP_ISO_DISCONNECTED:
        return _connected_generators(m, ra, roots, count)
    if m.nc == 0:
        return _wreath_generators(m)
    return _stabiliser(m, True)[0]


def isomorphism(fp1, data1, fp2, data2, classes1=(), classes2=(), root_rule=ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA):
    r"""
    Return an isomorphism from the first decorated map to the second, or
    ``None`` if they are not isomorphic.

    On the active edges, the isomorphism returned sends the half-edge of
    canonical label ``i`` in the first map to the one of canonical label
    ``i`` in the second, the canonical labellings being those of
    :func:`canonical_labelling` under the same rule. The maps are
    isomorphic when their minimal keys (:ref:`omap-keys`), or on
    disconnected maps their minimal sequences of blocks
    (:ref:`omap-blocks`), are equal; ``None`` is returned otherwise.

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
    _check_decorated_oriented_map(fp1, data1, classes1)
    _check_decorated_oriented_map(fp2, data2, classes2)
    if len(data1) != len(data2):
        raise ValueError("data1 and data2 must have the same number of arrays")
    if len(classes1) != len(classes2):
        raise ValueError("classes1 and classes2 must have the same number of arrays")
    cdef domap_t m1, m2
    cdef array.array _buf1 = _convert(m1, fp1, tuple(data1), tuple(classes1))
    cdef array.array _buf2 = _convert(m2, fp2, tuple(data2), tuple(classes2))
    cdef int n1 = m1.map.n, n2 = m2.map.n
    cdef int ne1 = omap_num_edges(&m1.map), ne2 = omap_num_edges(&m2.map)
    cdef int len1, len2, h, e, i
    cdef array.array key1, key2, ra, rb, seq1, seq2
    cdef array.array work = None
    # an isomorphism is a bijection of the active half-edges preserving the
    # twins, hence the numbers of active edges and half-edges; the invariants
    # of _domap_invariants_differ, O(n), are compared before the two searches
    if ne1 != ne2 or _num_active(&m1.map) != _num_active(&m2.map):
        return None
    if m1.nc:
        work = _new(DOMAP_INVARIANTS_WORK_SIZE(n1 if n1 > n2 else n2))
    if _domap_invariants_differ(m1, m2, NULL if work is None else work.data.as_ints):
        return None
    key1 = _new(DOMAP_ISO_KEY_SIZE(n1, m1.nd, m1.nc))
    key2 = _new(DOMAP_ISO_KEY_SIZE(n2, m2.nd, m2.nc))
    ra = _new(n1)
    rb = _new(n2)
    len1 = _run(m1, rule, key1, ra.data.as_ints, NULL, NULL)
    len2 = _run(m2, rule, key2, rb.data.as_ints, NULL, NULL)
    if len1 == DOMAP_ISO_DISCONNECTED or len2 == DOMAP_ISO_DISCONNECTED:
        # a connected map has one block, a disconnected one several
        if len1 != len2:
            return None
        seq1 = _new(DOMAP_ISO_SEQUENCE_SIZE(n1, m1.nd, m1.nc))
        seq2 = _new(DOMAP_ISO_SEQUENCE_SIZE(n2, m2.nd, m2.nc))
        len1 = _search(m1, rule, False, seq1, ra.data.as_ints, NULL, NULL)
        len2 = _search(m2, rule, False, seq2, rb.data.as_ints, NULL, NULL)
        if (len1 != len2
                or memcmp(seq1.data.as_ints, seq2.data.as_ints, len1 * sizeof(int)) != 0):
            return None
    elif (ne1 != ne2 or len1 != len2
            or memcmp(key1.data.as_ints, key2.data.as_ints, len1 * sizeof(int)) != 0):
        return None
    r = array.array('i', [-1] * n1)
    rbinv = [-1] * (2 * ne2)
    for h in range(n2):
        if rb.data.as_ints[h] != -1:
            rbinv[rb.data.as_ints[h]] = h
    for h in range(n1):
        if ra.data.as_ints[h] != -1:
            r[h] = rbinv[ra.data.as_ints[h]]
    if n1 == n2:
        inactive2 = [e for e in range(n2 // 2) if m2.map.fp[2 * e] == -1 and m2.map.fp[2 * e + 1] == -1]
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
    search of the key on a disconnected map follows, among the roots tied
    in a component, the first one of each way of extending the numbering
    of the class values: two tied roots extending it alike are exchanged
    by an automorphism fixing the components labelled before, so their
    subtrees reach the same sequences. A leaf equal to the best one so far
    differs from it by an automorphism fixing the components where their
    roots agree, so the search returns at once to the first level where
    their roots differ and tries its next tied root. On `k` tied copies of
    a component the search follows one path per distinct numbering at
    each level: `k` isolated edges with one class value take one path
    instead of `2^k k!` leaves. Tied roots with different numberings are
    all followed, so the search remains exponential in the worst case.
    The count goes through the search of :func:`automorphism_generators`,
    which has a cut of its own: polynomial on tied copies whose values are
    all mapped after the first components, such as isolated edges with one
    or two class values, it stays exponential in the worst case.

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

    The key alone of tied copies follows one path per numbering of the class
    values: 21 isolated edges with the class value ``0`` beside 20 with
    ``1``::

        sage: fp = list(range(82))
        sage: c = [0] * 42 + [1] * 40
        sage: array('i', canonical_key(fp, classes=(c,))) == array('i', fp + c)
        True

    Their count `21! \cdot 2^{21} \cdot 20! \cdot 2^{20}`, by the search of
    :func:`automorphism_generators`, whose placements of an edge mapping no
    new value are skipped after the first that failed::

        sage: count = canonical_key(fp, classes=(c,), automorphism_count=True)[1]
        sage: count
        273336954172212502469698893548133708673843200000000
        sage: count == factorial(21) * 2**21 * factorial(20) * 2**20
        True

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
    cdef domap_t m
    cdef array.array _buf = _prepare(m, fp, data, classes)
    cdef array.array key = _new(DOMAP_ISO_KEY_SIZE(m.map.n, m.nd, m.nc))
    cdef array.array r = _new(m.map.n)
    cdef int count, length
    order = 0
    length = _run(m, rule, key, r.data.as_ints, NULL, &count)
    if length == DOMAP_ISO_DISCONNECTED:
        _search(m, rule, False, _new(DOMAP_ISO_SEQUENCE_SIZE(m.map.n, m.nd, m.nc)), r.data.as_ints, NULL, NULL)
        if automorphism_count:
            order = _stabiliser(m, False)[1] if m.nc else _wreath_order(m)
        length = _key_of(m, r, key)
    out = [(<char *> key.data.as_ints)[:length * sizeof(int)]]
    if relabelling:
        out.append(_complete_relabelling(r, m.map.n, omap_num_edges(&m.map)))
    if automorphism_count:
        if order:
            out.append(order)
        else:
            out.append(count if count else 1)
    if len(out) == 1:
        return out[0]
    return tuple(out)
