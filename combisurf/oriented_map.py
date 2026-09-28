r"""
Combinatorial maps on oriented surfaces

The main class in this module is :class:`OrientedMap` which describes a
cell decomposition of an oriented surface.

A *walk* on a map is a word (see :mod:`combisurf.word`) whose letters are
active half-edges of the map and such that each half-edge ends at the vertex
where the next one starts. It can be given as a string such as ``"0,~1"``, a
list of integers or an ``array('i')``. A walk is *closed* if moreover its last
half-edge ends at the vertex where its first one starts.

A *multiwalk* is a finite collection of walks with non-negative integer
multiplicities. It can be given as a single walk (with multiplicity one), a
list of walks ``[w0, w1, ...]`` (each with multiplicity one) or a list of pairs
``[(w0, m0), (w1, m1), ...]``. It is taken as given: backtracking is not
simplified, repeated walks are not merged and conjugate closed walks are not
identified.
"""
# ****************************************************************************
#  This file is part of combisurf
#
#       Copyright (C) 2018 Mark Bell
#                     2018-2026 Vincent Delecroix
#                     2026 Oscar Fontaine
#                     2024 Kai Fu
#                     2026 Juliette Schabanel
#                     2018 Saul Schleimer
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

import collections
import itertools
import numbers
from array import array

from sage.structure.richcmp import op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE, rich_to_bool

from combisurf.misc import array_hash
from combisurf.permutation import (perm_init, perm_check, perm_trim, perm_cycles, perm_on_array, perm_on_edge_array,
                          perm_invert, perm_conjugate, perm_conjugate_transposition_inplace, perm_cycle_string, perm_dense_cycles, perm_cycles_lengths,
                          perm_cycles_to_string, perm_on_list, perm_on_edge_list, perm_cycle_type,
                          perm_num_cycles, str_to_cycles, str_to_cycles_and_data, perm_compose, perm_from_base64_str,
                          uint_base64_str, uint_from_base64_str, perm_base64_str,
                          perm_orbit, perm_orbit_size,
                          perms_are_transitive, perms_orbits, perm_edge_orbits, edge_relabelling_from)
from combisurf.word import word_init, word_string

def check_relabelling(arg, ne):
    r"""
    EXAMPLES::

        sage: from combisurf.oriented_map import check_relabelling
        sage: from combisurf.permutation import perm_cycle_string
        sage: p = check_relabelling("(0,1,2)", 3)
        sage: perm_cycle_string(p, edge_like=True)
        '(0,1,2)(~0,~1,~2)'
        sage: p = check_relabelling("(0,1,2)(~0,~1,~2)", 3)
        sage: perm_cycle_string(p, edge_like=True)
        '(0,1,2)(~0,~1,~2)'
        sage: p = check_relabelling("(0,~1,2)", 3)
        sage: perm_cycle_string(p, edge_like=True)
        '(0,~1,2)(~0,1,~2)'
        sage: p = check_relabelling("(0,1,2,~0,~1,~2)", 3)
        sage: perm_cycle_string(p, edge_like=True)
        '(0,1,2,~0,~1,~2)'
    """
    n = 2 * ne
    if isinstance(arg, str):
        p = perm_init(arg, n, edge_like=True, partial=True)
    else:
        p = perm_init(arg, n, partial=True)

    if len(p) != n:
        raise ValueError("len(p) = {} while n = {}".format(len(p), n))

    for h in range(n):
        if p[h] == -1 and p[h ^ 1] == -1:
            p[h] = h
            p[h ^ 1] = h ^ 1
        elif p[h] == -1:
            p[h] = p[h ^ 1] ^ 1
        elif p[h ^ 1] == -1:
            p[h ^ 1] = p[h] ^ 1
        elif p[h ^ 1] != p[h] ^ 1:
            raise ValueError("invalid input for relabelling (arg={})".format(arg))

    if not perm_check(p, n):
        raise ValueError("invalid input for relabelling (arg={})".format(arg))

    return p


# half-edge versus dart

# TODO: make sure to remove trailing -1 to make equality consistent
# TODO: should we maintain num_active_darts, num_active_edges

# TODO: add support for half-edges and edges data
#  have a python data protocol for map
#    realloc
#    hash
#    relabel
#    TODO: keep best relabelling generic stuff with reasonable cost
# TODO: ambiguity oriented could refer to the orientation of edges
class OrientedMap:
    r"""
    Map on oriented surface.

    An ``OrientedMap`` is encoded by three permutations called the *vertex permutation*,
    the *edge permutation* and the *face permutation*. The edge permutation is always
    implicit and the vertex and face permutations are always abbreviated as ``vp`` and
    ``fp``. The cycles in the cycle decomposition of ``vp`` and ``fp`` encode
    respectively the vertices and the faces of the map. The domain of the permutations
    is the set of *half-edges* of the map. Each *half-edge* could either be
    paired to make an edge, or be alone and form a *folded* edge.

    EXAMPLES:

        sage: from combisurf import OrientedMap

    The most user friendly way to initialize an ``OrientedMap`` is to use
    strings that describe the cycle decomposition of the vertex and face
    permutations. Here, edges correspond to pair of elements of the form ``i``
    and ``~i``. Here is a map made of three edges::

        sage: OrientedMap(vp="(0,~2,~1,~3,1,~0,3)(2)", fp="(0,1,~2,2)(~0,3,~1,~3)")
        OrientedMap("(0,~2,~1,~3,1,~0,3)(2)", "(0,1,~2,2)(~0,3,~1,~3)")

    Since vertex and face permutations are redundant, you can provide only one
    of them::

        sage: OrientedMap(vp="(0,~2,~1,~3,1,~0,3)(2)")
        OrientedMap("(0,~2,~1,~3,1,~0,3)(2)", "(0,1,~2,2)(~0,3,~1,~3)")
        sage: OrientedMap(fp="(0,1,~2,2)(~0,3,~1,~3)")
        OrientedMap("(0,~2,~1,~3,1,~0,3)(2)", "(0,1,~2,2)(~0,3,~1,~3)")

    Under the hood, the vertex permutation of the above example
    turns out to be `0 \mapsto 5`, `1 \mapsto 6`, `2 \mapsto 1`,
    `3 \mapsto 7`, `4 \mapsto 4`, `5 \mapsto 3`, `6 \mapsto 0`
    `7 \mapsto 2`. Or in cycle notation  ``(0,5,3,7,2,1,6)(4)``.
    More precisely, half-edges in the map use non-negative integers and an
    even half-edge is always paired to the half-edge with the next integer (so
    that ``i`` and ``~i`` in the string notation correspond to ``2i`` and
    ``2i+1`` in the actual permutation)::

        sage: m = OrientedMap(vp="(0,~2,~1,~3,1,~0,3)(2)", fp="(0,1,~2,2)(~0,3,~1,~3)")
        sage: m.vertex_permutation()
        array('i', [5, 6, 1, 7, 4, 3, 0, 2])
        sage: from combisurf.permutation import perm_cycle_string
        sage: perm_cycle_string(m.vertex_permutation())
        '(0,5,3,7,2,1,6)(4)'

    One can alternatively initialize an oriented map with a permutation either as a
    list of images::

        sage: OrientedMap(vp=[5, 6, 1, 7, 4, 3, 0, 2])
        OrientedMap("(0,~2,~1,~3,1,~0,3)(2)", "(0,1,~2,2)(~0,3,~1,~3)")

    Or as a list of disjoint cycles (including fixed points)::

        sage: OrientedMap(vp=[[0,5,3,7,2,1,6], [4]])
        OrientedMap("(0,~2,~1,~3,1,~0,3)(2)", "(0,1,~2,2)(~0,3,~1,~3)")

    It is not necessary for edges to be consecutive integers. One can for example
    initialize a map on a torus with edges `2` and `5`::

        sage: OrientedMap(vp="(2,5,~2,~5)")
        OrientedMap("(2,5,~2,~5)", "(2,5,~2,~5)")

    The corresponding vertex and face permutation of the above examples are lists of
    length 16 where only the images of `4`, `5`, `10` and `11` are different from
    `-1`::

        sage: OrientedMap(vp="(2,5,~2,~5)").vertex_permutation()
        array('i', [-1, -1, -1, -1, 10, 11, -1, -1, -1, -1, 5, 4])

    In cycle notation, if an half-edge is not mentioned then the corresponding edge is folded::

        sage: OrientedMap(vp="(0,1)(~0)")
        OrientedMap("(0,1)(~0)", "(0,~0,1)")
        sage: OrientedMap(vp="(0,1)")
        OrientedMap("(0,1)", "(0,1)")
        sage: OrientedMap(vp=[[0,2],[1]])
        OrientedMap("(0,1)(~0)", "(0,~0,1)")

    Folded edge could also be initialized from list of images as follows::

        sage: OrientedMap(vp=[2,1,0,-1])
        OrientedMap("(0,1)(~0)", "(0,~0,1)")
    """
    __slots__ = ['_vp', '_fp', '_mutable']

    def __init__(self, vp=None, fp=None, mutable=False, check=True):
        r"""
        INPUT:

        - ``vp`` -- ``None`` or data to initialize the vertex permutation

        - ``fp`` -- ``None`` or data to initialize the face permutation

        - ``mutable`` (boolean, default ``False``) -- whether the resulting map should be mutable
          or immutable

        - ``check`` (boolean, default ``True``) -- whether to perform consistency checks of
          the data
        """
        if vp is None and fp is None:
            self._vp = array("i", [])
            self._fp = array("i", [])
            self._mutable = mutable
            return

        if vp is not None:
            vp = perm_init(vp, partial=True, edge_like=isinstance(vp, str))
            if len(vp) % 2 == 1:
                vp.append(-1)
            ne = len(vp) // 2

        if fp is not None:
            fp = perm_init(fp, partial=True, edge_like=isinstance(fp, str))
            if len(fp) % 2 == 1:
                fp.append(-1)
            ne = len(fp) // 2

        if vp is None:
            vp = array('i', [-1] * (2 * ne))
            for i in range(2 * ne):
                if fp[i] == -1:
                    continue
                ii = fp[i ^ 1 if fp[i ^ 1] != -1 else i]
                vp[ii] = i

        if fp is None:
            fp = array('i', [-1] * (2 * ne))
            for i in range(2 * ne):
                if vp[i] == -1:
                    continue
                ii = vp[i] ^ 1 if (vp[i] != -1 and vp[vp[i] ^ 1] != -1) else vp[i]
                if ii != -1:
                    fp[ii] = i

        self._vp = vp
        self._fp = fp
        self._mutable = mutable

        if check:
            self._check(ValueError)

    def _clear_trailing_edges(self):
        vp = self._vp
        fp = self._fp
        while vp and vp[-2] == -1:
            vp.pop()
            vp.pop()
            fp.pop()
            fp.pop()

    def _half_edge_string(self, e):
        return '~%d' % (e // 2) if e % 2 else '%d' % (e // 2)

    def __str__(self):
        vp_cycles = perm_cycles(self._vp)
        vertices = perm_cycles_to_string(vp_cycles, edge_like=True)
        fp_cycles = perm_cycles(self._fp)
        faces = perm_cycles_to_string(fp_cycles, edge_like=True)
        return f"OrientedMap(\"{vertices}\", \"{faces}\")"

    __repr__ = __str__

    def _ep(self, h):
        r"""
        Return the image of ``h`` under the edge permutation.
        """
        if self._vp[h] == -1:
            return -1
        elif self._vp[h ^ 1] == -1:
            return h
        else:
            return h ^ 1

    def _check(self, error=RuntimeError):
        if not (hasattr(self, '_vp') and hasattr(self, '_fp')):
            raise error("missing attributes: these must be _vp, _ep, _fp, _data")
        if not perm_check(self._vp):
            raise error(f"vp is not a permutation: {self._vp}")
        if not perm_check(self._fp):
            raise error(f"fp is not a permutation: {self._fp}")
        if len(self._vp) != len(self._fp):
            raise error("vp and fp have different lengths")
        if len(self._vp) % 2:
            raise error("vp and fp must have even lengths")

        ne = len(self._vp) // 2
        if ne == 0:
            return

        for h in range(2 * ne):
            if (self._vp[h] == -1) != (self._fp[h] == -1):
                raise error(f"vp (={self._vp}) and fp (={self._fp}) with different domains")

        for e in range(ne):
            if self._vp[2 * e] == -1 and self._vp[2 * e + 1] != -1:
                raise error(f"half-edge {2 * e + 1} is active but its twin {2 * e} is not")

        if self._vp[-2] == -1 or self._fp[-2] == -1:
            raise error("trailing inactive edges")

        for h in range(2 * ne):
            if self._vp[h] != -1 and self._fp[self._ep(self._vp[h])] != h:
                raise error(f"fev relation not satisfied at half-edge h={self._half_edge_string(h)}")

    def __getstate__(self):
        a = [self._mutable]
        a.extend(self._fp)
        return a

    def __setstate__(self, arg):
        # We do not know how many slots we have in data
        self._mutable = arg[0]
        self._fp = array('i', arg[1:])
        n = len(self._fp)

        self._vp = array('i', [-1] * n)
        for i in range(n):
            if self._fp[i] == -1:
                continue
            ii = (i ^ 1) if self._fp[i ^ 1] != -1 else i
            self._vp[self._fp[ii]] = i

    def set_immutable(self):
        r"""
        Make the oriented map immutable.

        EXAMPLES:

        By default, ``OrientedMap`` are immutable::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(fp="(0,1,~1,~0)")
            sage: m.is_mutable()
            False

        But they can be made mutable by passing the ``mutable=True`` option::

            sage: m = OrientedMap(fp="(0,1,~1,~0)", mutable=True)
            sage: m.is_mutable()
            True
            sage: m.set_immutable()
            sage: m.is_mutable()
            False
        """
        self._mutable = False

    def is_mutable(self):
        return self._mutable

    def _assert_mutable(self):
        if not self._mutable:
            raise ValueError("immutable map; use a mutable copy instead")

    def _assert_immutable(self):
        if self._mutable:
            raise ValueError("mutable map; use the method .set_immutable() to make it immutable")

    def _assert_connected(self):
        if not self.is_connected():
            raise ValueError("non-connected map")

    def __hash__(self):
        self._assert_immutable()
        return array_hash(self._vp)

    def _check_half_edge(self, h):
        r"""
        Check whether ``h`` is a valid half-edge and raise a ``TypeError`` or
        a ``ValueError`` otherwise.

        If the half-edge is valid return ``h`` as a Python integer.

        TESTS::

            sage: from combisurf import OrientedMap
            sage: OrientedMap("(0,1,~1)")._check_half_edge(0)
            0
            sage: OrientedMap("(0,1,~1)")._check_half_edge(-4)
            Traceback (most recent call last):
            ...
            ValueError: half-edge number out of range (=-4)
            sage: OrientedMap("(0,1,~1)")._check_half_edge(12)
            Traceback (most recent call last):
            ...
            ValueError: half-edge number out of range (=12)
            sage: OrientedMap("(0,1,~1)")._check_half_edge(1)
            Traceback (most recent call last):
            ...
            ValueError: inactive half-edge (=1)
        """
        if not isinstance(h, numbers.Integral):
            raise TypeError(f"invalid half-edge {h} of type {type(h).__name__}")
        h = int(h)
        if h < 0 or h >= len(self._vp):
            raise ValueError(f"half-edge number out of range (={h})")
        if self._vp[h] == -1:
            raise ValueError(f"inactive half-edge (={h})")
        return h

    def _check_half_edge_or_negative(self, h):
        if not isinstance(h, numbers.Integral):
            raise TypeError(f"invalid half-edge {h} of type {type(h).__name__}")
        h = int(h)
        if h >= 0:
            self._check_half_edge(h)
        return h

    # TODO: this does not make any sense, self._ep(h) == 0 is not testing
    # at all that the edge is folded
    def _check_half_edge_folded(self, h):
        if not isinstance(h, numbers.Integral):
            raise TypeError(f"invalid half-edge {h} of type {type(h).__name__}")
        h = int(h)
        if h < 0 or h >= len(self._vp):
            raise ValueError(f"half-edge number out of range (={h})")
        if self._ep(h) == 0:
            raise ValueError(f"the edge corresponding to {h} is folded")

    def _check_edge(self, e):
        if not isinstance(e, numbers.Integral):
            raise TypeError(f"invalid edge {e}")
        e = int(e)
        if e < 0 or 2 * e >= len(self._vp):
            raise ValueError(f"edge number out of range e={e}")
        if self._vp[2 * e] == -1:
            raise ValueError(f"inactive edge e={e}")
        return e

    def _check_face_index(self, f):
        if not isinstance(f, numbers.Integral):
            raise TypeError(f"invalid face {f}")
        f = int(f)
        if f < 0 or f >= self.num_faces():
            raise ValueError(f"face number out of range f={f}")
        return f

    def _check_walk(self, w, closed=False):
        r"""
        Check that ``w`` is a walk on this map and return it as an ``array('i')``.

        See the module documentation for the definition of a walk.

        INPUT:

        - ``w`` -- a string, a list of integers or an array

        - ``closed`` -- boolean (default: ``False``); whether the walk must be closed

        OUTPUT: the walk as an ``array('i')``. If ``w`` is already an
        ``array('i')`` then ``w`` itself is returned, not a copy.

        Raise a ``ValueError`` if ``w`` is not a walk on this map, or not a
        closed one when ``closed=True``.

        TESTS::

            sage: from array import array
            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(fp="(0,1,~0,~1)")

            sage: w = array("i", [0, 2])
            sage: m._check_walk(w) is w
            True
            sage: m._check_walk([0, 2])
            array('i', [0, 2])
            sage: m._check_walk("0,~1")
            array('i', [0, 3])
            sage: m._check_walk(array("I", [0, 2]))
            array('i', [0, 2])
            sage: m._check_walk([0, 3, 5])
            Traceback (most recent call last):
            ...
            ValueError: invalid walk [0, 3, 5]: half-edge number out of range (=5) at position 2

            sage: m = OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~4,~5)(~2,5,~3)")
            sage: m._check_walk([0, 2])
            array('i', [0, 2])
            sage: m._check_walk([0, 0])
            Traceback (most recent call last):
            ...
            ValueError: invalid walk [0, 0]: 0 at position 0 ends at vertex 1 but 0 at position 1 starts at vertex 0
            sage: m._check_walk([0, 2], closed=True)
            Traceback (most recent call last):
            ...
            ValueError: non-closed walk [0, 2]: 2 at position 1 ends at vertex 2 but 0 at position 0 starts at vertex 0
        """
        w = word_init(w)
        h2v = self.half_edge_to_vertex()

        hprev = None # previous half-edge
        for i, h in enumerate(w):
            try:
                self._check_half_edge(h)
            except ValueError as e:
                raise ValueError(f"invalid walk {word_string(w)}: {e} at position {i}") from None
            if hprev is not None and h2v[hprev ^ 1] != h2v[h]:
                raise ValueError(f"invalid walk {word_string(w)}: {hprev} at position {i - 1} ends at vertex {h2v[hprev ^ 1]} but {h} at position {i} starts at vertex {h2v[h]}")
            hprev = h

        if w and closed and h2v[w[0]] != h2v[w[-1] ^ 1]:
            raise ValueError(f"non-closed walk {word_string(w)}: {w[-1]} at position {len(w) - 1} ends at vertex {h2v[w[-1] ^ 1]} but {w[0]} at position 0 starts at vertex {h2v[w[0]]}")

        return w

    def _check_multiwalk(self, words, multiplicities=None, closed=False):
        r"""
        Return the multiwalk given by ``words`` and ``multiplicities`` in normal form.

        See the module documentation for the definitions of walks and
        multiwalks. The normal form does minimal work: the walks are checked and
        converted with :meth:`_check_walk` and the ones with multiplicity zero
        are dropped, but backtracking is not simplified, repeated walks are not
        merged and conjugate closed walks are not identified.

        INPUT:

        - ``words`` -- a single walk, a list of walks ``[w0, w1, ...]`` or a list of
          pairs ``[(w0, m0), (w1, m1), ...]`` where ``m0``, ``m1``, ... are
          non-negative integers. The items of a list must be either all walks or
          all pairs.

        - ``multiplicities`` -- (optional) a list of non-negative integers with the
          same length as ``words``. Only allowed when ``words`` is a list of walks.
          If not provided, the multiplicities are taken from the pairs, or set to
          `1` for a list of walks.

        - ``closed`` -- (boolean, default ``False``) whether the walks must be closed

        OUTPUT: a pair ``(words, multiplicities)`` of lists of the same length,
        where ``words`` contains arrays ``array('i')`` and ``multiplicities``
        contains positive Python integers. The walks with multiplicity zero in
        the input are dropped, so the output can be shorter than the input.

        Raise a ``TypeError`` or a ``ValueError`` if the input is not a valid
        multiwalk.

        TESTS::

            sage: from array import array
            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(fp="(0,1,~0,~1)")

        A single walk::

            sage: m._check_multiwalk(array("i", [0, 2]))
            ([array('i', [0, 2])], [1])
            sage: m._check_multiwalk([0, 2])
            ([array('i', [0, 2])], [1])
            sage: m._check_multiwalk("0,~0")
            ([array('i', [0, 1])], [1])

        A list of walks, with or without multiplicities::

            sage: m._check_multiwalk([])
            ([], [])
            sage: m._check_multiwalk([[0, 2], [], "0,~1"])
            ([array('i', [0, 2]), array('i'), array('i', [0, 3])], [1, 1, 1])
            sage: m._check_multiwalk([[0, 2], [1, 3], [0]], [2, 0, 1])
            ([array('i', [0, 2]), array('i', [0])], [2, 1])
            sage: m._check_multiwalk(iter([[0, 2], [1]]), (3, 4))
            ([array('i', [0, 2]), array('i', [1])], [3, 4])

        A list of pairs (walk, multiplicity)::

            sage: m._check_multiwalk([([0, 2], 1), ([], 3), ("0,1", 2), ([0, 0], 0)])
            ([array('i', [0, 2]), array('i'), array('i', [0, 2])], [1, 3, 2])

        Closed walks::

            sage: m = OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~4,~5)(~2,5,~3)")
            sage: m._check_multiwalk([[0, 1], [0, 2]])
            ([array('i', [0, 1]), array('i', [0, 2])], [1, 1])
            sage: m._check_multiwalk([[0, 1], [0, 2]], closed=True)
            Traceback (most recent call last):
            ...
            ValueError: non-closed walk [0, 2]: 2 at position 1 ends at vertex 2 but 0 at position 0 starts at vertex 0

        Invalid inputs::

            sage: m = OrientedMap(fp="(0,1,~0,~1)")
            sage: m._check_multiwalk([0, 2], [1])
            Traceback (most recent call last):
            ...
            TypeError: invalid multiwalk: multiplicities can only be given with a list of walks
            sage: m._check_multiwalk([([0, 2], 1)], [1])
            Traceback (most recent call last):
            ...
            TypeError: invalid multiwalk: multiplicities can not be given with a list of pairs (walk, multiplicity)
            sage: m._check_multiwalk([[0, 2], [1]], [1])
            Traceback (most recent call last):
            ...
            ValueError: invalid multiwalk: got 2 walks but 1 multiplicities
            sage: m._check_multiwalk([[0, 2], ([0, 2], 2)])
            Traceback (most recent call last):
            ...
            TypeError: invalid multiwalk: the item ([0, 2], 2) at position 1 is a pair (walk, multiplicity) but the first item is a walk
            sage: m._check_multiwalk([([0, 2], 2), [0, 2]])
            Traceback (most recent call last):
            ...
            TypeError: invalid multiwalk: the item [0, 2] at position 1 is a walk but the first item is a pair (walk, multiplicity)
            sage: m._check_multiwalk(([0, 2], 3))
            Traceback (most recent call last):
            ...
            TypeError: invalid multiwalk: the input should either be a walk, a list of walks or a list of pairs (walk, multiplicity)
            sage: m._check_multiwalk([[0], 3])
            Traceback (most recent call last):
            ...
            TypeError: invalid multiwalk: the input should either be a walk, a list of walks or a list of pairs (walk, multiplicity)
            sage: m._check_multiwalk([([0], 1), 3])
            Traceback (most recent call last):
            ...
            TypeError: invalid multiwalk: the input should either be a walk, a list of walks or a list of pairs (walk, multiplicity)
            sage: m._check_multiwalk([([0, 2], "2")])
            Traceback (most recent call last):
            ...
            TypeError: invalid multiwalk: the multiplicity '2' of type str of the walk [0, 2] at position 0 must be an integer
            sage: m._check_multiwalk([[0, 2], [1]], [1, -1])
            Traceback (most recent call last):
            ...
            ValueError: invalid multiwalk: the multiplicity -1 of the walk [1] at position 1 must be non-negative
            sage: m._check_multiwalk([([0, 3, 1, 7], 2)])
            Traceback (most recent call last):
            ...
            ValueError: invalid walk [0, 3, 1, 7]: half-edge number out of range (=7) at position 3
        """
        # a single walk is a string, an array or a sequence of integers
        if isinstance(words, (str, array)):
            single = True
        else:
            if not isinstance(words, (list, tuple)):
                words = list(words)
            single = bool(words) and isinstance(words[0], numbers.Integral)

        if single:
            if multiplicities is not None:
                raise TypeError("invalid multiwalk: multiplicities can only be given with a list of walks")
            return [self._check_walk(words, closed=closed)], [1]

        def is_pair(item):
            # the letters of a walk are integers, so a sequence of length 2
            # whose first entry is not an integer is a pair (walk, multiplicity)
            return (isinstance(item, (tuple, list)) and len(item) == 2 and
                    not isinstance(item[0], numbers.Integral))

        # the first item decides whether words is a list of walks or of pairs
        pairs = bool(words) and is_pair(words[0])
        expected = "a pair (walk, multiplicity)" if pairs else "a walk"
        unexpected = "a walk" if pairs else "a pair (walk, multiplicity)"

        if multiplicities is None:
            if not pairs:
                multiplicities = [1] * len(words)
        else:
            if pairs:
                raise TypeError("invalid multiwalk: multiplicities can not be given with a list of pairs (walk, multiplicity)")
            if not isinstance(multiplicities, (list, tuple, array)):
                multiplicities = list(multiplicities)
            if len(multiplicities) != len(words):
                raise ValueError(f"invalid multiwalk: got {len(words)} walks but {len(multiplicities)} multiplicities")

        words_clean = []
        multiplicities_clean = []
        for i, item in enumerate(words):
            if isinstance(item, numbers.Integral):
                raise TypeError("invalid multiwalk: the input should either be a walk, a list of walks or a list of pairs (walk, multiplicity)")
            if is_pair(item) != pairs:
                raise TypeError(f"invalid multiwalk: the item {item!r} at position {i} is {unexpected} but the first item is {expected}")

            if pairs:
                w, m = item
            else:
                w, m = item, multiplicities[i]

            w = self._check_walk(w, closed=closed)

            if not isinstance(m, numbers.Integral):
                raise TypeError(f"invalid multiwalk: the multiplicity {m!r} of type {type(m).__name__} of the walk {word_string(w)} at position {i} must be an integer")
            m = int(m)
            if m < 0:
                raise ValueError(f"invalid multiwalk: the multiplicity {m} of the walk {word_string(w)} at position {i} must be non-negative")
            if m == 0:
                continue
            words_clean.append(w)
            multiplicities_clean.append(m)

        return words_clean, multiplicities_clean

    @staticmethod
    def _subdivide_false(u, v, half_edge_list, h):
        return 0

    @staticmethod
    def _subdivide_true(u, v, half_edge_list, h):
        if u == v:
            return 2
        elif len(half_edge_list) > 1:
            return 1
        else:
            return 0

    def graph(self, directed=False, subdivide=True, root=None):
        r"""
        Return the sage graph, the embedding and the root edge associated to this map.

        INPUT:

        - ``directed`` (boolean, default ``False``) -- whether to return a SageMath
          ``Graph`` or ``DiGraph``

        - ``subdivide`` -- boolean, non-negative integer, or callable
          ``(u, v, half_edge_list, h) -> k`` (default ``True``) -- if
          ``True`` then insert one vertex on each edge and two vertices on
          each loop so that the resulting graph is simple; if ``False`` (or
          the integer ``0``) then never subdivide, so ``G`` is returned as a
          multigraph with loops, exactly mirroring this map's own
          combinatorics; a positive integer ``k`` inserts exactly ``k``
          vertices on every edge (always simple, like ``True``); a callable
          chooses the subdivision count per edge.

        - ``root`` (``None`` or a valid half-edge) -- if specified, it should be
          a half-edge that is used to identify the external face in planar drawings

        OUTPUT: a 4-tuple ``(G, embedding, root_edge, edge_vertices)`` made of

        - a SageMath graph ``G``

        - the embedding of ``G`` given as a dictionary of cyclic ordering of neighborhoods

        - the root edge as a pair of vertices in ``G``

        - a list of lists of vertices of ``G`` where the list at position ``e`` is
          the list of vertices on ``G`` that are on the path representing the edge ``e``
          in this oriented map

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(vp="(0,~2)(~0,3,1)(~1,~5,2)(~3,4)(~4,5)")
            sage: G, em, root, edge_list = m.graph()
            sage: pos = G.layout_planar(on_embedding=em, external_face=root)
            sage: G.plot(pos=pos, vertex_labels=False, edge_labels=True)
            Graphics object consisting of ... graphics primitives

        With ``subdivide=False``, ``G`` mirrors the map's own combinatorics
        exactly, including a genuine loop (both half-edges of an edge at the
        same vertex, as opposed to a folded edge where only one is active)::

            sage: m = OrientedMap(vp="(0,~0)")
            sage: G, em, root, edge_list = m.graph(subdivide=False)
            sage: G
            Looped multi-graph on 1 vertex
            sage: list(G.edges(labels=False))
            [(0, 0)]
        """
        # NOTE: sage graphs use *clockwise* order for the neighbors
        if self.has_folded_edge():
            raise NotImplementedError

        if root is None:
            root = len(self._vp) - 2

        if subdivide is True:
            # every loop/multi-edge always gets subdivided away, so G is
            # always simple regardless of the map's own combinatorics
            allow_loops = allow_multiedges = False
            subdivide = self._subdivide_true
        elif subdivide is False:
            # never subdivided, so G must keep any loop/multi-edge as-is
            allow_loops = allow_multiedges = True
            subdivide = self._subdivide_false
        elif isinstance(subdivide, numbers.Integral):
            k = int(subdivide)
            if k < 0:
                raise ValueError("subdivide can not be a negative integer")
            # a constant k>=1 subdivides every edge into a simple path;
            # k=0 never subdivides, same as subdivide=False above
            allow_loops = allow_multiedges = (k == 0)
            subdivide = lambda u, v, half_edge_list, h: k
        elif callable(subdivide):
            # a custom callable's per-edge behavior isn't known statically
            # (it might leave some loop/multi-edge unsubdivided), so allow
            # both rather than risk G.add_edge raising below
            allow_loops = allow_multiedges = True
        else:
            raise TypeError("subdivide must be a boolean or a callable")

        vertices = self.vertices()
        half_edge_to_vertex = [-1] * (2 * len(self._vp))
        for i, v in enumerate(vertices):
            for h in v:
                half_edge_to_vertex[h] = i

        # cyclic ordering of neighbors (to specify embedding for SageMath graphs)
        embedding = {}

        if directed:
            from sage.graphs.digraph import DiGraph
            G = DiGraph(len(vertices), multiedges=allow_multiedges, loops=allow_loops)
        else:
            from sage.graphs.graph import Graph
            G = Graph(len(vertices), multiedges=allow_multiedges, loops=allow_loops)

        # group half edges depending on their endpoints
        half_edges = collections.defaultdict(list)
        for e in range(len(self._vp) // 2):
            u = half_edge_to_vertex[2 * e]
            v = half_edge_to_vertex[2 * e + 1]
            if u <= v:
                half_edges[(u, v)].append(2 * e)
            else:
                half_edges[(v, u)].append(2 * e + 1)

        edge_vertices = [None] * (len(self._vp) // 2)
        for (u, v), half_edge_list in half_edges.items():
            for h in half_edge_list:
                k = subdivide(u, v, half_edge_list, h)
                e = h // 2
                edge_vertices[e] = [u] + [G.add_vertex() for _ in range(k)] + [v]
                if h % 2:
                    edge_vertices[e].reverse()

                for i in range(k + 1):
                    G.add_edge(edge_vertices[e][i], edge_vertices[e][i + 1], e)

                # embedding of subdivisions
                for i in range(1, k + 1):
                    embedding[edge_vertices[e][i]] = [edge_vertices[e][i-1], edge_vertices[e][i+1]]

        # embedding of vertices of the original map
        for i, v in enumerate(vertices):
            neighborhood = []
            for h in v:
                e = h // 2
                if h % 2 == 0:
                    neighborhood.append(edge_vertices[e][1])
                else:
                    neighborhood.append(edge_vertices[e][-2])
            embedding[i] = neighborhood

        if root % 2:
            root_edge = edge_vertices[root // 2][:2]
        else:
            root_edge = edge_vertices[root // 2][-2:]

        return G, embedding, root_edge, edge_vertices

    def plot(self, directed=False, subdivide=True, root=None, edge_labels=True, vertex_colors=None, edge_colors=None):
        r"""
        Plot the map.

        INPUT:

        - ``directed``, ``subdivide`` -- options forwarded to :meth:`graph`
        - ``edge_labels``: boolean specifying whether to plot the labels of the edges.
        - ``edge_colors``: dictionary specifying the color to assign to each edge color.
        - ``vertex_colors``: dictionary specifying the color to assign to each vertex color.
        """
        G, em, r, edge_list = self.graph(directed=directed, subdivide=subdivide, root=root)
        pos = G.layout_planar(on_embedding=em, external_face=r)
        if vertex_colors is None:
            vertex_colors={'red':list(range(self.num_vertices()))}
        vertex_colors['#C0C0C0'] = list(range(self.num_vertices(), G.order()))
        if edge_colors is None:
            edge_cols = None
        else:
            edge_cols = collections.defaultdict(list)
            for c, edges in edge_colors.items():
                for e in edges:
                    edge_cols[c].extend(edge_list[e][i:i+2] for i in range(len(edge_list[e])-1))

        return G.plot(pos=pos, vertex_labels=False, edge_labels=edge_labels, vertex_colors=vertex_colors, edge_colors=edge_cols)

    def __eq__(self, other):
        return self._vp == other._vp

    def __ne__(self, other):
        return self._vp != other._vp

    def _cmp_(self, other):
        r"""
        TESTS::

            sage: import itertools
            sage: from combisurf import OrientedMap
            sage: ts = [OrientedMap(fp="(0,1,2)"), OrientedMap(fp="(0,1,2)(~0,~1,~2)"), OrientedMap(fp="(0,1,2)(~0,~2,~1)"), OrientedMap(fp="(0,~0,1)(~1,2,~2)")]
            sage: for t1, t2 in itertools.product(ts, repeat=2):
            ....:     c1 = t1._cmp_(t2)
            ....:     c2 = t2._cmp_(t1)
            ....:     assert c1 == -c2
            ....:     assert (c1 == 0) == (t1 == t2)
        """
        if type(self) is not type(other):
            raise TypeError("can not compare {} with {}".format(type(self).__name__, type(other).__name__))

        return (self._vp > other._vp) - (self._vp < other._vp)

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` and ``other`` according to the operator ``op``.

        EXAMPLES::

            sage: import itertools
            sage: from combisurf import OrientedMap

            sage: ts = [OrientedMap("(0,1,2)"), OrientedMap("(0,1,2)(~0,~1,~2)"), OrientedMap("(0,~0,1)(~1,2,~2)")]
            sage: for m1, m2 in itertools.product(ts, repeat=2):
            ....:     if m1 == m2:
            ....:         assert (m1 <= m2)
            ....:         assert (m1 >= m2)
            ....:         assert not (m1 < m2)
            ....:         assert not (m1 > m2)
            ....:     else:
            ....:         assert (m1 < m2) + (m2 < m1) == 1
            ....:         assert (m1 > m2) + (m2 > m1) == 1
            ....:         assert (m1 < m2) == (m1 <= m2)
            ....:         assert (m1 > m2) == (m1 >= m2)
        """
        if type(self) is not type(other):
            raise TypeError("can not compare {} with {}".format(type(self).__name__, type(other).__name__))

        return rich_to_bool(op, self._cmp_(other))

    def __lt__(self, other):
        return self._richcmp_(other, op_LT)

    def __le__(self, other):
        return self._richcmp_(other, op_LE)

    def __gt__(self, other):
        return self._richcmp_(other, op_GT)

    def __ge__(self, other):
        return self._richcmp_(other, op_GE)

    def copy(self, mutable=None):
        r"""
        Return a copy of this oriented map.

        If ``mutable`` is set to ``True`` or ``False`` return respectively an
        immutable or mutable version of the map.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap([[0,2,4],[1,3,5]], mutable=True)
            sage: s = m.copy()
            sage: s == m
            True
            sage: s.reverse_orientation(0)
            sage: s == m
            False

            sage: m.set_immutable()
            sage: m.copy() is m
            True

            sage: t = m.copy(mutable=True)
            sage: t.reverse_orientation(0)
            sage: s == t
            True

        """
        if mutable is None:
            mutable = self._mutable

        if not self._mutable and not mutable:
            # avoid copies of immutable objects
            return self

        m = OrientedMap.__new__(OrientedMap)
        m._fp = self._fp[:]
        m._vp = self._vp[:]
        m._mutable = mutable
        return m

    def vertex_permutation(self, copy=True):
        r"""
        Return the permutation encoding the vertices of ``self``.

        Going counterclockwise around vertices makes a permutation of the darts
        whose cycles in its cycle decomposition are in bijection with vertices.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: sphere = OrientedMap(fp="(0,1,2)(~0,~2,~1)")
            sage: sphere.vertex_permutation()
            array('i', [5, 2, 1, 4, 3, 0])

        If the triangulation has folded edges, then the vertex permutation is only partial::

            sage: m = OrientedMap(fp="(0,1,2)(~0,3,4)")
            sage: m.vertex_permutation()
            array('i', [4, 8, 1, -1, 2, -1, 0, -1, 6, -1])
        """
        if copy:
            return self._vp[:]
        else:
            return self._vp

    def next_at_vertex(self, h, check=True):
        r"""
        Return the half-edge after ``h`` around the corresponding vertex.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(fp="(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: m.next_at_vertex(0)
            18
            sage: m.next_at_vertex(9)
            7
            sage: m.next_at_vertex(5)
            13
        """
        if check:
            h = self._check_half_edge(h)
        return self._vp[h]

    def previous_at_vertex(self, h, check=True):
        r"""
        Return the half-edge before ``h`` around the corresponding vertex.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(fp="(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: m.previous_at_vertex(9)
            7
            sage: m.previous_at_vertex(8)
            10
            sage: m.previous_at_vertex(4)
            11
        """
        if check:
            h = self._check_half_edge(h)
        return self._fp[self._ep(h)]

    def edge_permutation(self, copy=True):
        r"""
        Return the edge permutation.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap(fp="(0,1,2)(~0,~1,~2)").edge_permutation()
            array('i', [1, 0, 3, 2, 5, 4])
            sage: OrientedMap(fp="(0,1,2)").edge_permutation()
            array('i', [0, -1, 2, -1, 4, -1])
        """
        return array('i', [self._ep(h) for h in range(len(self._vp))])

    def next_in_edge(self, h, check=True):
        r"""
        Return the half-edge which makes an edge with ``h`` (or ``h`` itself if
        the edge is folded).

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(vp="(0)(~0,1,2)")
            sage: m.next_in_edge(0)
            1
            sage: m.next_in_edge(2)
            2
        """
        if check:
            self._check_half_edge(h)
        return self._ep(h)

    previous_in_edge = next_in_edge

    def face_permutation(self, copy=True):
        r"""
        Return the permutation encoding the faces of ``self``.

        Going counterclockwise around faces makes a permutation of the
        darts whose cycles in its cycle decomposition are in bijection with
        vertices.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: sphere = OrientedMap(fp="(0,1,2)(~0,~2,~1)")
            sage: sphere.face_permutation()
            array('i', [2, 5, 4, 1, 0, 3])

        If the triangulation has folded edges, then the face permutation is only partial::

            sage: m = OrientedMap(fp="(0,1,2)(~0,3,4)")
            sage: m.face_permutation()
            array('i', [2, 6, 4, -1, 0, -1, 8, -1, 1, -1])
        """
        if copy:
            return self._fp[:]
        else:
            return self._fp

    def next_in_face(self, h, check=True):
        r"""
        Return the half-edge after ``h`` in the corresponding face.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(fp="(0,1,2)(~1,3,4)")
            sage: m.next_in_face(0)
            2
            sage: m.next_in_face(2)
            4
            sage: m.next_in_face(4)
            0
            sage: m.next_in_face(1)
            Traceback (most recent call last):
            ...
            ValueError: inactive half-edge (=1)
        """
        if check:
            h = self._check_half_edge(h)
        return self._fp[h]

    def previous_in_face(self, h, check=True):
        r"""
        Return the half-edge before ``h`` in the corresponding face.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(fp="(~11,4,~3)(~10,~0,11)(~9,0,10)(~8,9,1)(~7,8,~1)(~6,7,2)(~5,6,~2)(~4,5,3)")
            sage: m.previous_in_face(10)
            9
            sage: m.previous_in_face(1)
            21
            sage: m.previous_in_face(3)
            16
        """
        if check:
            h = self._check_half_edge(h)
        return self._ep(self._vp[h])

    def half_edges(self):
        r"""
        Iterate through the half-edges of this constellation.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: list(OrientedMap(fp="(0,1,2)(~1,3,4)").half_edges())
            [0, 2, 3, 4, 6, 8]
            sage: list(OrientedMap(fp="(2,3)").half_edges())
            [4, 6]
        """
        for h in range(len(self._vp)):
            if self._vp[h] != -1:
                yield h

    def num_half_edges(self):
        r"""
        Return the number of half edges.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap(fp="(0,1,2)(~1,3,4)").num_half_edges()
            6
            sage: OrientedMap(fp="(2,3)").num_half_edges()
            2
        """
        return sum(h != -1 for h in self._vp)

    def has_folded_edge(self):
        r"""
        Return whether this map has a folded edge.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap(fp="(0,1,2)(~0,~1,~2)").has_folded_edge()
            False
            sage: OrientedMap(fp="(0,1,2)").has_folded_edge()
            True
            sage: OrientedMap(fp="(2,~2)").has_folded_edge()
            False
            sage: OrientedMap(fp="(2,3,~3)").has_folded_edge()
            True
        """
        return any((self._vp[2 * e] != -1 and self._vp[2 * e + 1] == -1) for e in range(len(self._vp) // 2))

    def folded_half_edges(self):
        r"""
        Iterate through half-edges on a folded edge.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: list(OrientedMap(fp="(0,1,2)(~0,~1,~2)").folded_half_edges())
            []
            sage: list(OrientedMap(fp="(0,1,2)").folded_half_edges())
            [0, 2, 4]
            sage: list(OrientedMap(fp="(2,5)").folded_half_edges())
            [4, 10]
        """
        vp = self._vp
        for e in range(len(self._vp) // 2):
            if vp[2 * e] != -1 and vp[2 * e + 1] == -1:
                yield 2 * e

    def num_folded_edges(self):
        r"""
        Return the number of folded edges.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap(fp="(0,1,2)(~0,~1,~2)").num_folded_edges()
            0
            sage: OrientedMap(fp="(0,1,2)").num_folded_edges()
            3
            sage: OrientedMap(fp="(2,~2)").num_folded_edges()
            0
            sage: OrientedMap(fp="(2)(3)").num_folded_edges()
            2
        """
        return sum((self._vp[i] != -1 and self._vp[i + 1] == -1) for i in range(0, len(self._vp), 2))

    def num_edges(self):
        r"""
        Return the number of edges.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap(fp="(0,1,2)(~0,~1,~2)").num_edges()
            3
            sage: OrientedMap(fp="(0,1,2)").num_edges()
            3
            sage: OrientedMap("(2,4,~2,~4)").num_edges()
            2
        """
        return sum(self._vp[i] != -1 or self._vp[i + 1] != -1 for i in range(0, len(self._vp), 2))

    def edge_indices(self):
        r"""
        Return the list of edge indices.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap(fp="(0,1,2)(3,4,5)(~0,~3,6)").edge_indices()
            [0, 1, 2, 3, 4, 5, 6]
            sage: OrientedMap(fp="(1,2)(3,5)(~3,6)").edge_indices()
            [1, 2, 3, 5, 6]
        """
        return [e for e in range(len(self._vp) // 2) if self._vp[2 * e] != -1]

    def vertices(self):
        r"""
        Return the list of vertices as lists of half-edges

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(fp="(0,1,2)(3,4,5)(~0,~3,6)")
            sage: m.vertices()
            [[0, 4, 2, 1, 12, 6, 10, 8, 7]]

            sage: OrientedMap("").vertices()
            [[]]
        """
        if not self._vp:
            return [[]]
        return perm_cycles(self._vp, True)

    def half_edge_to_vertex(self):
        r"""
        Return an array whose element at index ``h`` is the index of the vertex
        incident to the half-edge ``h``.

        The indices are the ones of :meth:`vertices`: the entry at ``h`` is the
        ``i`` for which ``h`` belongs to ``self.vertices()[i]``. Inactive
        half-edges get the value ``-1``, and the array has one entry per
        half-edge, so it is empty on a map without edges even though such a map
        still has one vertex.

        .. SEEALSO::

            :meth:`vertices`, :meth:`half_edge_to_face`

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap("(0,1,2)(~0,3,4)(~1,10,~8)(~2,7,~5)(~3,5,6)(~4,15,~13)(~6,16,~15)(~7,8,9)(~9,~12,~17)(~10,11,12)(~11,13,14)(~14,~16,17)", "(0,4,~13,~11,~10,~1)(~0,2,~5,~3)(1,~8,~7,~2)(3,6,~15,~4)(5,7,9,~17,~16,~6)(8,10,12,~9)(11,14,17,~12)(13,15,16,~14)")
            sage: m.half_edge_to_vertex()
            array('i', [0, 1, 0, 2, 0, 3, 1, 4, 1, 5, 4, 3, 4, 6, 3, 7, 7, 2, 7, 8, 2, 9, 9, 10, 9, 8, 10, 5, 10, 11, 5, 6, 6, 11, 11, 8])

        An example with inactive half-edges::

            sage: m = OrientedMap(vp="(2,1,5)(~1,~2,~5)")
            sage: m.half_edge_to_vertex()
            array('i', [-1, -1, 0, 1, 0, 1, -1, -1, -1, -1, 0, 1])
        """
        return perm_dense_cycles(self._vp)

    # TODO: to follow sage Graph convention, we may want to use
    # def vertex_degree(self, h=None)
    # def face_degree(self, h=None)
    def vertex_profile(self, sort=False, reverse=True):
        r"""
        Return the vertex profile of this map.

        The vertex profile is the list of vertex degrees.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap("(0,3)(~0,1,~3,~1)").vertex_profile()
            [2, 4]
            sage: OrientedMap("","").vertex_profile()
            [0]

            sage: OrientedMap("(0,3)(~0,1,~3,~1)").vertex_profile(sort=True, reverse=True)
            [4, 2]
            sage: OrientedMap("(0,3)(~0,1,~3,~1)").vertex_profile(sort=True, reverse=False)
            [2, 4]
        """
        profile = [len(v) for v in self.vertices()]
        if sort:
            profile.sort(reverse=reverse)
        return profile

    def num_vertices(self):
        r"""
        Return the number of vertices.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap(fp="(0,1,2)(3,4,5)(~0,~3,6)").num_vertices()
            1
            sage: OrientedMap(fp="(3,~3)").num_vertices()
            2

            sage: OrientedMap("").num_vertices()
            1
        """
        if not self._vp:
            return 1
        return perm_num_cycles(self._vp)

    def faces(self):
        r"""
        Return the list of edges as lists of half-edges.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(fp="(0,1,2)(3,4,5)(~3,6)")
            sage: m.faces()
            [[0, 2, 4], [6, 8, 10], [7, 12]]

            sage: OrientedMap("").faces()
            [[]]
        """
        if not self._vp:
            return [[]]
        return perm_cycles(self._fp, True)

    def half_edge_to_face(self):
        r"""
        Return an array whose element at index ``h`` is the index of the face
        incident to the half-edge ``h``.

        The indices are the ones of :meth:`faces`: the entry at ``h`` is the
        ``i`` for which ``h`` belongs to ``self.faces()[i]``. Inactive
        half-edges get the value ``-1``, and the array has one entry per
        half-edge, so it is empty on a map without edges even though such a map
        still has one face.

        .. SEEALSO::

            :meth:`faces`, :meth:`half_edge_to_vertex`

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap("(0,1,2)(~0,3,4)(~1,10,~8)(~2,7,~5)(~3,5,6)(~4,15,~13)(~6,16,~15)(~7,8,9)(~9,~12,~17)(~10,11,12)(~11,13,14)(~14,~16,17)", "(0,4,~13,~11,~10,~1)(~0,2,~5,~3)(1,~8,~7,~2)(3,6,~15,~4)(5,7,9,~17,~16,~6)(8,10,12,~9)(11,14,17,~12)(13,15,16,~14)")
            sage: m.half_edge_to_face()
            array('i', [0, 1, 2, 0, 1, 2, 3, 1, 0, 3, 4, 1, 3, 4, 4, 2, 5, 2, 4, 5, 5, 0, 6, 0, 5, 6, 7, 0, 6, 7, 7, 3, 7, 4, 6, 4])

        An example with inactive half-edges::

            sage: m = OrientedMap(vp="(2,1,5)(~1,~2,~5)")
            sage: m.half_edge_to_face()
            array('i', [-1, -1, 0, 1, 1, 2, -1, -1, -1, -1, 2, 0])
        """
        return perm_dense_cycles(self._fp)

    def face_profile(self, sort=False, reverse=True):
        r"""
        Return the face degrees.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap("(0,3)(~0,1,~3,~1)").face_profile()
            [6]
            sage: OrientedMap("").face_profile()
            [0]

        The output is not sorted unless you set ``sort=True``:

            sage: m = OrientedMap(fp="(0,1,2,~0,3)(~1,4,5)(~2,~3,~4,~5)")
            sage: m.face_profile()
            [5, 3, 4]
            sage: m.face_profile(sort=True, reverse=False)
            [3, 4, 5]
            sage: m.face_profile(sort=True, reverse=True)
            [5, 4, 3]
        """
        profile = [len(f) for f in self.faces()]
        if sort:
            profile.sort(reverse=reverse)
        return profile

    def num_faces(self):
        r"""
        Return the number of faces.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap(fp="(0,1,2)(3,4,5)(~0,~3,6)").num_faces()
            3
            sage: OrientedMap(fp="(2,~2)").num_faces()
            1
            sage: OrientedMap("").num_faces()
            1
        """
        if not self._vp:
            return 1
        return perm_num_cycles(self._fp)

    def vertex_degree(self, h, check=True):
        r"""
        Return the degree of the vertex incident to the half-edge ``h``.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(vp="(0,1,2)(3,4,5)(~0,~3,6)")
            sage: m.vertex_degree(2)
            3
        """
        if len(self._vp) == 0:
            return 0
        if check:
            self._check_half_edge(h)
        return perm_orbit_size(self._vp, h)

    def vertex_turn(self, h0, h1):
        r"""
        Return the number of turns from h0 to h1 around a vertex

        EXAMPLES:

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: m.vertex_turn(3, 7)
            3
        """

        vp = self.vertex_permutation(copy=False)
        if h0 == h1:
            return 0
        current = vp[h0]
        turn = 1
        while current != h0 and current != h1:
            turn += 1
            current = vp[current]
        if current == h0:
            raise ValueError("The half-edge {} does not belong to the same vertex than the half-edge {}.".format(h0, h1))
        return turn

    def face_degree(self, h, check=True):
        r"""
        Return the degree of the face incident to the half-edge ``h``.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(vp="(0,1,2)(3,4,5)(~0,~3,6)")
            sage: m.face_degree(2)
            9
        """
        if len(self._vp) == 0:
            return 0
        if check:
            self._check_half_edge(h)
        return perm_orbit_size(self._fp, h)

    def face_turn(self, h0, h1):
        r"""
        Return the number of turns from h0 to h1 around a face

        EXAMPLES:

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: m.face_turn(0, 8)
            5
        """

        fp = self.face_permutation(copy=False)
        if h0 == h1:
            return 0
        current = fp[h0]
        turn = 1
        while current != h0 and current != h1:
            turn += 1
            current = fp[current]
        if current == h0:
            raise ValueError("The half-edge {} does not belong to the same face than the half-edge {}.".format(h0, h1))
        return turn

    def is_connected(self):
        r"""
        Return whether the constellation is connected.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: OrientedMap(fp="(0,1,2)(3,4,5)(~0,~3,6)").is_connected()
            True
            sage: OrientedMap(fp="(0,1,2)(3,4,5)").is_connected()
            False
            sage: OrientedMap(fp="(2,3)(~2,4)").is_connected()
            True
            sage: OrientedMap(fp="(2,3)(4,~4)").is_connected()
            False
        """
        return perms_are_transitive((self._vp, self._fp))

    def connected_components(self):
        r"""
        Return the connected components as a list of lists of edges.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(fp="(0,1,3)(~0,~1,~3)(2,4,5)(~2,~4,~5)")
            sage: m.connected_components()
            [[0, 1, 3], [2, 4, 5]]

        To construct the triangulation induced on each connected component, one can
        use the method :meth:`submap`::

            sage: c0, c1 = m.connected_components()
            sage: m.submap(c0)
            OrientedMap("(0,~3,1,~0,3,~1)", "(0,1,3)(~0,~1,~3)")
            sage: m.submap(c1)
            OrientedMap("(2,~5,4,~2,5,~4)", "(2,4,5)(~2,~4,~5)")

        These components could be relabelled on the edges 0, 1, ... by providing
        the option (``relabel=True``)::

            sage: m.submap(c0, relabel=True)
            OrientedMap("(0,~2,1,~0,2,~1)", "(0,1,2)(~0,~1,~2)")
            sage: m.submap(c1, relabel=True)
            OrientedMap("(0,~2,1,~0,2,~1)", "(0,1,2)(~0,~1,~2)")

        Some more examples::

            sage: OrientedMap(fp="(2,~2)(4)").connected_components()
            [[2], [4]]
        """
        return perm_edge_orbits(self._vp)

    def submap(self, edges, relabel=False, mutable=False, check=True):
        r"""
        Return the submap of this map induced on ``edges``.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: t = OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~3,5)(~2,~5,6)(~7,~6,7)")
            sage: t.submap([0, 3, 5])
            OrientedMap("(0,5,3)(~0,~3)(~5)", "(0,~3,5,~5)(~0,3)")
            sage: t.submap([5, 0, 3])  # isomorphic graph
            OrientedMap("(0,5,3)(~0,~3)(~5)", "(0,~3,5,~5)(~0,3)")

        Note that the result might not be connected::

            sage: t.submap([1, 6])
            OrientedMap("(1)(~1)(6,~6)", "(1,~1)(6)(~6)")
            sage: t.submap([1, 6], relabel=True)
            OrientedMap("(0)(~0)(1,~1)", "(0,~0)(1)(~1)")

        This also works fine with partial maps::

            sage: m = OrientedMap(fp="(2,3,~5,7,~3,5)")
            sage: m.submap([2, 7])
            OrientedMap("(2,7)", "(2,7)")
        """
        if check:
            edges = [self._check_edge(h) for h in edges]
            if len(set(edges)) != len(edges):
                raise ValueError("redundant edges")

        p = [-1] * len(self._vp)
        for i, j in enumerate(edges):
            p[2 * j] = 2 * i
            p[2 * j + 1] = 2 * i + 1

        if relabel:
            ne = len(edges)
        else:
            ne = len(self._vp) // 2

        vp = array('i', [-1] * (2 * ne))

        for e in edges:
            for h in [2 * e, 2 * e + 1]:
                if self._vp[h] == -1:
                    continue
                h_image = self._vp[h]
                while p[h_image] == -1:
                    h_image = self._vp[h_image]
                if relabel:
                    vp[p[h]] = p[h_image]
                else:
                    vp[h] = h_image

        perm_trim(vp)
        return OrientedMap(vp=vp, mutable=mutable, check=check)

    def connected_components_submaps(self, relabel=False, mutable=False):
        r"""
        Run through the connected components as graphs.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(fp="(0,1,3)(~0,~1,~3)(2,4,5)(~2,~4,~5)")
            sage: list(m.connected_components_submaps())
            [OrientedMap("(0,~3,1,~0,3,~1)", "(0,1,3)(~0,~1,~3)"),
             OrientedMap("(2,~5,4,~2,5,~4)", "(2,4,5)(~2,~4,~5)")]
            sage: list(m.connected_components_submaps(relabel=True))
            [OrientedMap("(0,~2,1,~0,2,~1)", "(0,1,2)(~0,~1,~2)"),
             OrientedMap("(0,~2,1,~0,2,~1)", "(0,1,2)(~0,~1,~2)")]
        """
        for comp in self.connected_components():
            yield self.submap(comp, relabel=relabel, mutable=mutable)

    def euler_characteristic(self):
        r"""
        Return the Euler characteristic of this constellation.

        EXAMPLES::

            sage: from combisurf import OrientedMap

        Spheres::

            sage: OrientedMap(fp="(0,1,2)").euler_characteristic()
            2
            sage: OrientedMap(fp="(0,~0)").euler_characteristic()
            2
            sage: OrientedMap(fp="(0)").euler_characteristic()
            2
            sage: OrientedMap(fp="(0,1,2)(~0)(~1)(~2)").euler_characteristic()
            2
            sage: OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~2)(~3,~4)").euler_characteristic()
            2

        A torus::

            sage: OrientedMap(fp="(0,1,2)(~0,~1,~2)").euler_characteristic()
            0

        A genus 2 surface::

            sage: m = OrientedMap(fp="(0,1,2)(~2,3,4)(~4,5,6)(~6,~0,7)(~7,~1,8)(~8,~3,~5)")
            sage: m.euler_characteristic()
            -2

        The empty map is a sphere::

            sage: OrientedMap("").euler_characteristic()
            2
        """
        return self.num_faces() - self.num_edges() + (self.num_vertices() + self.num_folded_edges())

    def genus(self, connected=True, check=2):
        r"""
        INPUT:

        - ``connected`` -- if ``False`` return a tuple if ``True`` and non-connected raises an error


        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: OrientedMap(fp="(0,1,2)").genus()
            0
            sage: OrientedMap(fp="(0,1,~0,~1)").genus()
            1

        If the map is non-connected you get an error::

            sage: OrientedMap(fp="(0,~0)(1,~1)").genus()
            Traceback (most recent call last):
            ...
            ValueError: non-connected map

        In that case, one can obtain the genus of each connected component by setting the
        argument ``connected`` to ``False``::

            sage: OrientedMap(fp="(0,~0)(1,~1)").genus(connected=False)
            [0, 0]
            sage: OrientedMap(fp="(0,2,~0,~2)(1,3)(~1,~3)").genus(connected=False)
            [1, 0]
        """
        if connected:
            if check >= 2:
                self._assert_connected()

            char = self.euler_characteristic()
            two_g = 2 - char
            assert two_g % 2 == 0
            return two_g // 2
        else:
            # TODO: implement something less costly
            return [cc.genus() for cc in self.connected_components_submaps(relabel=True)]

    def _spanning_forest(self, cycles, h2c, roots, used, name):
        r"""
        Return a spanning forest of the cells ``cycles``, rooted at ``roots``.

        This is the common core of the forest and the coforest of
        :meth:`forest_coforest_decomposition`: run on the vertices it builds
        the forest and run on the faces it builds the coforest.

        INPUT:

        - ``cycles`` -- the cells, as the list of their half-edges indexed by
          cell number, that is the output of :meth:`vertices` or :meth:`faces`

        - ``h2c`` -- the inverse map, sending a half-edge to the index of its
          cell, that is the output of :meth:`half_edge_to_vertex` or
          :meth:`half_edge_to_face`

        - ``roots`` -- ``None`` or the cells the trees are rooted at. A tree
          never leaves its connected component. When ``None`` the cells that no
          tree has reached start a new one, which gives exactly one tree per
          connected component.

        - ``used`` -- an array of flags indexed by the edges, read to skip the
          edges already taken and written for the ones this call takes

        - ``name`` -- the name of the argument ``roots`` came from, used in the
          error messages

        OUTPUT: an array of the length of ``cycles`` whose entry is ``-1`` at a
        root, the half-edge joining a cell to its parent at a cell some tree
        reached, and ``-2`` at a cell no tree reached

        This is exercised through :meth:`forest_coforest_decomposition`.
        """
        vp = self._vp
        nc = len(cycles)
        forest = array('i', [-2] * nc)

        if roots is None:
            todo = []
        else:
            # one pass: validate, mark the roots and collect them. Marking as
            # we go makes forest[c] == -1 the test for a repeat, and reading
            # ``roots`` only once lets it be any iterable.
            todo = []
            for c in roots:
                if not isinstance(c, numbers.Integral):
                    raise TypeError(f"invalid entry {c} of type {type(c).__name__} in {name}")
                c = int(c)
                if c < 0 or c >= nc:
                    raise ValueError(f"{name} must consist of integers in range({nc}), got {c}")
                if forest[c] == -1:
                    raise ValueError(f"{name} lists {c} twice")
                forest[c] = -1
                todo.append(c)

        c0 = 0
        while True:
            while todo:
                c = todo.pop()
                for h in cycles[c]:
                    if used[h // 2]:
                        continue
                    # a folded edge is a loop: ep(h) is h, so the cell on the
                    # other side is the one we come from and it never extends
                    # the forest
                    if vp[h ^ 1] == -1:
                        continue
                    h = h ^ 1
                    cc = h2c[h]
                    if forest[cc] == -2:
                        forest[cc] = h
                        used[h // 2] = 1
                        todo.append(cc)
            if roots is not None:
                break
            while c0 < nc and forest[c0] != -2:
                c0 += 1
            if c0 == nc:
                break
            forest[c0] = -1
            todo.append(c0)

        return forest

    def forest_coforest_decomposition(self, root_vertices=None, root_faces=None):
        r"""
        Return a triple ``(forest, coforest, complementary_edges)`` with the
        given roots.

        This is the analogue of :meth:`tree_cotree_decomposition` with several
        roots: the tree is replaced by a spanning forest with one tree per root
        vertex and the cotree by a spanning coforest with one tree per root
        face.

        INPUT:

        - ``root_vertices`` -- (default: ``None``) the vertices the trees of
          the forest are rooted at. A tree never leaves its connected
          component, so each component must contain one of them and a
          ``ValueError`` is raised otherwise. When ``None`` the smallest
          vertex index of each connected component is taken as its root, which
          gives one tree per component.

        - ``root_faces`` -- (default: ``None``) the faces the trees of the
          coforest are rooted at, with the same convention

        OUTPUT: a triple ``(forest, coforest, complementary_edges)`` of arrays
        of integers, in the format of :meth:`tree_cotree_decomposition`

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(vp="(2,3,1,5,6,~5)(~1,4,~2,~6)(~3,7,~4,8)(~7,9,10,11,~9)(~8,~10,~11)")
            sage: m.forest_coforest_decomposition((0,2,4), (1,))
            (array('i', [-1, 8, -1, 20, -1]),
             array('i', [2, -1, 7, 22]),
             array('i', [2, 5, 6, 7, 8, 9]))

        The empty map has one vertex and one face, both roots, and no edge::

            sage: OrientedMap().forest_coforest_decomposition()
            (array('i', [-1]), array('i', [-1]), array('i'))

        Left to itself on a map that is not connected, it roots one tree per
        component::

            sage: m = OrientedMap(fp="(0,1,3)(~0,~1,~3)(2,4,5)(~2,~4,~5)")
            sage: m.is_connected()
            False
            sage: m.connected_components()
            [[0, 1, 3], [2, 4, 5]]
            sage: forest, coforest, comp_edges = m.forest_coforest_decomposition()
            sage: forest
            array('i', [-1, -1])
            sage: coforest
            array('i', [-1, 1, -1, 5])

        the root of each tree being the smallest index it contains::

            sage: [v for v, x in enumerate(forest) if x == -1]
            [0, 1]
            sage: [f for f, x in enumerate(coforest) if x == -1]
            [0, 2]

        Given roots that miss a component, it says so rather than leaving a
        vertex or a face out::

            sage: m.forest_coforest_decomposition((0,), (0, 2))
            Traceback (most recent call last):
            ...
            ValueError: root_vertices must contain a vertex of each connected component
            sage: m.forest_coforest_decomposition((0, 1), (0,))
            Traceback (most recent call last):
            ...
            ValueError: root_faces must contain a face of each connected component
        """
        used = array('i', [0] * (len(self._vp) // 2))

        # the forest on the vertices, then the coforest the same way on the
        # faces of the dual. The two share ``used``, so that the coforest only
        # gets to pick among the edges the forest left.
        forest = self._spanning_forest(self.vertices(), self.half_edge_to_vertex(),
                                       root_vertices, used, "root_vertices")
        coforest = self._spanning_forest(self.faces(), self.half_edge_to_face(),
                                         root_faces, used, "root_faces")

        # a vertex left at -2 was reached by no tree, that is its component
        # holds no root vertex; likewise for the faces
        if -2 in forest:
            raise ValueError("root_vertices must contain a vertex of each connected component")
        if -2 in coforest:
            raise ValueError("root_faces must contain a face of each connected component")

        return (forest, coforest, array('i', [e for e in self.edge_indices() if not used[e]]))

    def tree_cotree_decomposition(self, root_vertex=0, root_face=0):
        r"""
        Return a tree cotree decomposition as a triple ``(tree, cotree, complementary_edges)``.

        INPUT:

        - ``root_vertex`` -- (default: ``0``) the vertex the tree is rooted at

        - ``root_face`` -- (default: ``0``) the face the cotree is rooted at

        OUTPUT: a triple ``(tree, cotree, complementary_edges)`` of arrays of
        integers. The ``tree`` and ``cotree`` have length respectively the
        number of vertices and the number of faces in the map. We describe
        ``tree`` below and the ``cotree`` is similar.

        - ``tree[root_vertex]`` is ``-1``
        - for a non-root vertex ``v``, ``tree[v]`` is a half-edge adjacent to ``v`` and
          going out from the root.

        The last entry ``complementary_edges`` is the array of edge indices that
        are neither part of the tree nor the cotree.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(vp="(2,3,1,5,6,~5)(~1,4,~2,~6)(~3,~4)")
            sage: tree, cotree, comp_edges = m.tree_cotree_decomposition()
            sage: tree
            array('i', [-1, 3, 7])
            sage: cotree
            array('i', [-1, 9, 4])
            sage: comp_edges
            array('i', [5, 6])

        To obtain the edges used in the tree and cotree respectively, one can
        do it as follows (and check that we indeed obtain a partition of
        edges)::

            sage: tree_edges = [h // 2 for h in tree if h != -1]
            sage: tree_edges
            [1, 3]
            sage: cotree_edges = [h // 2 for h in cotree if h != -1]
            sage: cotree_edges
            [4, 2]

        A single tree can not span a map that is not connected::

            sage: m = OrientedMap(fp="(0,1,3)(~0,~1,~3)(2,4,5)(~2,~4,~5)")
            sage: m.tree_cotree_decomposition()
            Traceback (most recent call last):
            ...
            ValueError: a tree cotree decomposition requires a connected map
            sage: m.forest_coforest_decomposition()[0]
            array('i', [-1, -1])

        .. SEEALSO::

            :meth:`forest_coforest_decomposition`
        """
        if not self.is_connected():
            raise ValueError("a tree cotree decomposition requires a connected map")

        return self.forest_coforest_decomposition((root_vertex,), (root_face,))

    def radial_map(self, mapping=False, mutable=False):
        r"""
        Return the radial map of this map.

        The *radial map* of an oriented map is the bipartite quadrangulation
        obtained by adding a vertex in the center of each face, joining this
        added vertex to every corner in the face and removing the original
        edges. The vertices of the radial map are in bijection with the union
        of vertices and faces of the original map. It has as many quadrilateral
        faces as edges in the original map.

        The convention used for labelling is that the half-edge of the radial
        map to the left of `h` in the original map is labelled `2h`. That way,
        vertices of the radial map are either cycles of positively oriented
        edges or cycles of negatively oriented edges. In particular, the
        bipartition of vertices is visible on the labelling.

        INPUT:

        - ``mapping`` -- boolean (default: ``False``); whether to also return
          the list of the images of the half-edges. The image of a half-edge is
          the walk of length two it becomes in the radial map, and the image of
          ``ep(h)`` is the reverse of the image of ``h``.

        - ``mutable`` -- boolean (default: ``False``); whether the result is
          mutable

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(vp="(0,1,~0,2)(~1,~2)")
            sage: m.radial_map()
            OrientedMap("(0,2,1,4)(~0,~2,~5,~1,~4,~3)(3,5)", "(0,~3,5,~2)(~0,4,~1,2)(1,~5,3,~4)")

        An example with inactive half edges::

            sage: m = OrientedMap(vp="(0,3,6,~3)(~0,1,~6,~1)")
            sage: m.radial_map()
            OrientedMap("(0,6,12,7)(~0,~3,~1,~7)(1,2,13,3)(~2,~13,~6,~12)", "(0,~7,12,~6)(~0,7,~1,3)(1,~3,13,~2)(2,~12,6,~13)")

        With ``mapping``, the images of the half-edges are returned as well::

            sage: m = OrientedMap(vp="(0,1,~0,2)(~1,~2)")
            sage: radial, mor = m.radial_map(mapping=True)
            sage: mor
            [array('i', [0, 5]),
             array('i', [4, 1]),
             array('i', [4, 11]),
             array('i', [10, 5]),
             array('i', [8, 7]),
             array('i', [6, 9])]
            sage: all(list(mor[h ^^ 1]) == [mor[h][1] ^^ 1, mor[h][0] ^^ 1]
            ....:     for h in m.half_edges())
            True

        Inactive half-edges have no image::

            sage: m = OrientedMap(vp="(0,3,6,~3)(~0,1,~6,~1)")
            sage: m.radial_map(mapping=True)[1]
            [array('i', [0, 7]),
             array('i', [6, 1]),
             array('i', [4, 27]),
             array('i', [26, 5]),
             None,
             None,
             array('i', [12, 25]),
             array('i', [24, 13]),
             None,
             None,
             None,
             None,
             array('i', [24, 5]),
             array('i', [4, 25])]

        The result is immutable unless ``mutable`` is set::

            sage: m = OrientedMap(vp="(0,1,~0,2)(~1,~2)")
            sage: m.radial_map().is_mutable()
            False
            sage: m.radial_map(mutable=True).is_mutable()
            True

        Raises a ``NotImplementedError`` on maps with folded edge::

            sage: OrientedMap("(0)").radial_map()
            Traceback (most recent call last):
            ...
            NotImplementedError: radial_map is not implemented on a map with a folded edge
        """
        if self.has_folded_edge():
            raise NotImplementedError("radial_map is not implemented on a map with a folded edge")
        n = len(self._vp)
        rvp = array('i', [-1] * (2 * n))
        rfp = array('i', [-1] * (2 * n))
        for h in range(n):
            if self._vp[h] == -1:
                continue
            rvp[2 * h] = 2 * self._vp[h]
            rvp[2 * h + 1] = 2 * self._fp[h] + 1
            rfp[2 * self._fp[h]] = 2 * h + 1
            rfp[2 * (h ^ 1) + 1] = 2 * self._fp[h]
        radial = OrientedMap(vp=rvp, fp=rfp, mutable=mutable)

        if not mapping:
            return radial

        mor = [None] * n
        for e in range(n // 2):
            h = 2 * e
            if self._vp[h] == -1:
                continue
            mor[h] = array('i', [2 * h, 2 * self._fp[h] + 1])
            mor[h + 1] = array('i', [2 * self._fp[h], 2 * h + 1])
        return radial, mor

    def reduced_map(self, forest=None, coforest=None, relabel=False, mapping=False, mutable=False, check=True):
        r"""
        Return the map obtained by contracting ``forest`` and deleting ``coforest``.

        The vertices of the result are the trees of ``forest`` and its faces
        are the trees of ``coforest``. Its edges are the complementary edges,
        which keep their labels unless ``relabel`` is set. If the pair
        (``forest``, ``coforest``) is a tree-cotree decomposition then the
        resulting map has a single vertex and a single face.

        INPUT:

        - ``forest``, ``coforest`` -- (default: ``None``) the half-edges of the
          edges to contract and to delete, in the format returned by
          :meth:`forest_coforest_decomposition`; entries equal to ``-1`` are
          ignored. When both are ``None`` a decomposition is computed. Listing
          the same edge twice raises a ``ValueError``.

        - ``relabel`` -- boolean (default: ``False``); whether to relabel the
          result on ``0, 1, ..., 2 * ne - 1`` so that it has no inactive
          half-edge. Contracting and deleting never renumber, so without this
          the labels of the complementary edges are kept and are sparse. The
          relabelling goes edge by edge, in increasing order, and preserves
          the parity of each half-edge inside its edge.

        - ``mapping`` -- boolean (default: ``False``); whether to also return
          the projection, a list indexed by the half-edges of this map giving
          for each of them the walk it becomes in the result, as an array of
          half-edges of the result (``None`` for inactive half-edges). A
          complementary half-edge is sent to itself, a half-edge of a
          contracted edge to the empty walk and a half-edge of a deleted edge
          to a walk homotopic to it, see below. The walk of ``ep(h)`` is the
          reverse of the walk of ``h``.

        - ``mutable`` -- boolean (default: ``False``); whether the result is
          mutable

        - ``check`` -- boolean (default: ``True``); whether to check that the
          half-edges of ``forest`` and ``coforest`` are ones of this map, that
          neither contains a cycle, and to check the map built. Listing the
          same edge twice is reported whatever its value.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap("(0,5,~4,7,4,~5)(~0,~1,~2,3,2,~6)(1,~3,~7,6)")
            sage: m
            OrientedMap("(0,5,~4,7,4,~5)(~0,~1,~2,3,2,~6)(1,~3,~7,6)", "(0,~6,~7,~4,7,~3,~2,3,1,~0,~5)(~1,6,2)(4,5)")
            sage: print(m.genus(), m.num_vertices(), m.num_faces())
            2 3 3

            sage: r = m.reduced_map()
            sage: r
            OrientedMap("(1,~3,~5,~1,~2,3,2,5)", "(1,~5,2,~1,5,~3,~2,3)")
            sage: print(r.genus(), r.num_vertices(), r.num_faces())
            2 1 1

        With a forest and a coforest that are not spanning, the result has one
        vertex per tree of the forest and one face per tree of the coforest::

            sage: f, c, _ = m.forest_coforest_decomposition(root_vertices=[0, 2], root_faces=[1])
            sage: r = m.reduced_map(forest=f, coforest=c)
            sage: print(r.genus(), r.num_vertices(), r.num_faces())
            2 2 1

        The labels of the complementary edges are kept, and ``relabel`` makes
        them consecutive::

            sage: f, c, comp = m.forest_coforest_decomposition()
            sage: comp
            array('i', [1, 2, 3, 5])
            sage: m.reduced_map(f, c)
            OrientedMap("(1,~3,~5,~1,~2,3,2,5)", "(1,~5,2,~1,5,~3,~2,3)")
            sage: m.reduced_map(f, c, relabel=True)
            OrientedMap("(0,~2,~3,~0,~1,2,1,3)", "(0,~3,1,~0,3,~2,~1,2)")

        With ``mapping``, the projection is returned as well. A half-edge of a
        deleted edge lies on the boundary of a disk, made of the faces of the
        coforest beyond it, and its walk goes around the other side of that
        disk::

            sage: r, proj = m.reduced_map(f, c, mapping=True)
            sage: proj
            [array('i'),
             array('i'),
             array('i', [2]),
             array('i', [3]),
             array('i', [4]),
             array('i', [5]),
             array('i', [6]),
             array('i', [7]),
             array('i', [11]),
             array('i', [10]),
             array('i', [10]),
             array('i', [11]),
             array('i', [2, 5]),
             array('i', [4, 3]),
             array('i'),
             array('i')]
            sage: all(list(proj[h ^^ 1]) == [r._ep(x) for x in reversed(proj[h])]
            ....:     for h in m.half_edges())
            True

        The radial map of the reduced map is the quad system of :meth:`quad_system`
        for the same forest and coforest. The isomorphism between them can be
        read on the projections: a complementary half-edge ``h`` becomes the
        walk ``proj[h]`` of length two in the quad system and the walk of
        ``h`` in the radial map of the reduced map, which gives the image of
        the half-edges of these walks::

            sage: rad_map, rad = r.radial_map(mapping=True)
            sage: q, qproj = m.quad_system(f, c, mapping=True)
            sage: p = {}
            sage: for e in comp:
            ....:     for h in (2 * e, 2 * e + 1):
            ....:         a, b = rad[h]
            ....:         p[a], p[b] = qproj[h]

        These walks need not go through every edge of the radial map, as the
        walk of ``h`` runs through the corners of ``h`` and of
        ``next_in_face(h)``, which the walks of other edges may use as well.
        But the map being connected, the relabelling extends uniquely along the
        face permutation and the edge permutation. It is consistent with the
        value read on the walks and conjugates the face permutations::

            sage: len(p), rad_map.num_half_edges()
            (14, 16)
            sage: rfp = rad_map.face_permutation()
            sage: qfp = q.face_permutation()
            sage: todo = list(p)
            sage: while todo:
            ....:     x = todo.pop()
            ....:     for y, z in ((x ^^ 1, p[x] ^^ 1), (rfp[x], qfp[p[x]])):
            ....:         if y not in p:
            ....:             p[y] = z
            ....:             todo.append(y)
            ....:         elif p[y] != z:
            ....:             raise AssertionError
            sage: len(p) == q.num_half_edges()
            True
            sage: all(qfp[p[x]] == p[rfp[x]] for x in rad_map.half_edges())
            True

        The relabelling also matches the projections on the contracted edges,
        that both send to the empty walk. On the deleted edges the two
        projections are homotopic walks but differ in general, the quad
        system keeping a walk of length two where the reduced map has to go
        around a disk::

            sage: deleted = [h // 2 for h in c if h != -1]
            sage: all([p[x] for y in proj[h] for x in rad[y]] == list(qproj[h])
            ....:     for h in m.half_edges() if h // 2 not in deleted)
            True
            sage: [(proj[2 * e], qproj[2 * e]) for e in deleted]
            [(array('i', [2, 5]), array('i', [20, 9])),
             (array('i', [11]), array('i', [28, 21]))]

        On a sphere the reduced map is empty::

            sage: m = OrientedMap("(0,1,2)(~0,~2,~1)")
            sage: m.genus()
            0
            sage: m.reduced_map(mapping=True)
            (OrientedMap("", ""),
             [array('i'), array('i'), array('i'), array('i'), array('i'), array('i')])

        Folded edges are allowed, but may be neither contracted nor deleted::

            sage: m = OrientedMap(vp="(0,2,~2,4)(~4)")
            sage: m.reduced_map(mapping=True)
            (OrientedMap("(0)", "(0)"),
             [array('i', [0]),
              None,
              None,
              None,
              array('i'),
              array('i'),
              None,
              None,
              array('i'),
              array('i')])

        TESTS::

            sage: m = OrientedMap("(0,5,~4,7,4,~5)(~0,~1,~2,3,2,~6)(1,~3,~7,6)")
            sage: f, c, _ = m.forest_coforest_decomposition()
            sage: m.reduced_map(f, None)
            Traceback (most recent call last):
            ...
            ValueError: forest and coforest must be given together
            sage: m.reduced_map(f, [2, 12])
            Traceback (most recent call last):
            ...
            ValueError: coforest contains a cycle
            sage: m.reduced_map([1, 15, 5], c)
            Traceback (most recent call last):
            ...
            ValueError: forest contains a cycle
            sage: m.reduced_map(f, [1, 12])
            Traceback (most recent call last):
            ...
            ValueError: the edge of the half-edge 1 is listed twice

        .. SEEALSO::

            :meth:`forest_coforest_decomposition`, :meth:`quad_system`,
            :meth:`radial_map`, :meth:`contract_edge`, :meth:`delete_edge`

        ALGORITHM:

        The face permutation of the result is obtained by following, from each
        complementary half-edge, the face permutation composed with the
        involution that exchanges the two half-edges of each deleted edge,
        until the next complementary half-edge. Each half-edge is visited once.

        For the projection, the faces are processed children first in the
        coforest, so that the walk of a deleted half-edge is obtained by
        concatenating the walks already computed on the boundary of its face.
        The cost is linear in the size of the output.
        """
        if forest is None and coforest is None:
            forest, coforest, _ = self.forest_coforest_decomposition()
        elif forest is None or coforest is None:
            raise ValueError("forest and coforest must be given together")

        vp = self._vp
        fp = self._fp
        n = len(vp)

        # kind[e] is 0 for an edge that stays, 1 for a contracted edge and 2
        # for a deleted one
        kind = array('i', [0] * (n // 2))
        contracted = []
        deleted = []
        for k, edges, listed in ((1, forest, contracted), (2, coforest, deleted)):
            for h in edges:
                if h == -1:
                    continue
                if check:
                    h = self._check_half_edge(h)
                if kind[h // 2]:
                    raise ValueError(f"the edge of the half-edge {h} is listed twice")
                if vp[h ^ 1] == -1:
                    raise ValueError(f"the half-edge {h} is folded")
                kind[h // 2] = k
                listed.append(h)

        if check:
            # a contracted cycle or a deleted cocycle would not give a map of
            # the same surface; union-find on the vertices, then on the faces
            for h2c, listed, name in ((self.half_edge_to_vertex(), contracted, "forest"),
                                      (self.half_edge_to_face(), deleted, "coforest")):
                root = list(range(max(h2c, default=-1) + 1))
                for h in listed:
                    u = h2c[h]
                    while root[u] != u:
                        root[u] = u = root[root[u]]
                    v = h2c[h ^ 1]
                    while root[v] != v:
                        root[v] = v = root[root[v]]
                    if u == v:
                        raise ValueError(f"{name} contains a cycle")
                    root[u] = v

        # contracting an edge removes its two half-edges from their faces and
        # deleting it glues the face of each half-edge to that of the other.
        # So the next half-edge in the face of the reduced map is reached by
        # skipping the contracted half-edges and crossing the deleted ones.
        # That is following fp composed with the involution swapping the two
        # half-edges of each deleted edge, a permutation whose orbit from a
        # kept half-edge gets to the next kept one, and each skipped half-edge
        # is visited once overall.
        rfp = array('i', [-1] * n)
        for h in range(n):
            if vp[h] == -1 or kind[h // 2]:
                continue
            y = fp[h]
            while kind[y // 2]:
                y = fp[y] if kind[y // 2] == 1 else fp[y ^ 1]
            rfp[h] = y

        relabelling = None
        if relabel:
            relabelling = array('i', [-1] * n)
            ne = 0
            for e in range(n // 2):
                if vp[2 * e] != -1 and not kind[e]:
                    relabelling[2 * e] = 2 * ne
                    relabelling[2 * e + 1] = 2 * ne + 1
                    ne += 1
            fp_new = array('i', [-1] * (2 * ne))
            for h in range(n):
                if rfp[h] != -1:
                    fp_new[relabelling[h]] = relabelling[rfp[h]]
            rfp = fp_new
        else:
            perm_trim(rfp)
        reduced = OrientedMap(fp=rfp, mutable=mutable, check=check)

        if not mapping:
            return reduced

        # A deleted half-edge c on the face g has its reverse on the parent
        # face of g, that is it lies on the child side of its edge. Pushing c
        # across the disk made of g and its descendants in the coforest gives
        # the walk of ep(c): the boundary of g read after c, in which a deleted
        # half-edge is itself replaced by the walk around its child. Filling
        # in the walks children first makes the whole of it linear in the
        # size of the output.
        h2f = self.half_edge_to_face()
        nf = max(h2f, default=-1) + 1
        child_edge = array('i', [-1] * nf)
        num_children = array('i', [0] * nf)
        for c in deleted:
            g = h2f[c]
            if child_edge[g] != -1:
                raise ValueError(f"the face of the half-edge {c} is the child of two edges of coforest")
            child_edge[g] = c
            num_children[h2f[c ^ 1]] += 1
        order = [g for g in range(nf) if child_edge[g] != -1 and not num_children[g]]
        i = 0
        while i < len(order):
            p = h2f[child_edge[order[i]] ^ 1]
            i += 1
            num_children[p] -= 1
            if not num_children[p] and child_edge[p] != -1:
                order.append(p)
        if len(order) != len(deleted):
            raise ValueError("coforest contains a cycle")

        proj = [None] * n
        for h in range(n):
            if vp[h] != -1 and kind[h // 2] != 2:
                proj[h] = array('i', [] if kind[h // 2] else [h])
        walk = [None] * nf
        for g in order:
            c = child_edge[g]
            w = array('i')
            y = fp[c]
            while y != c:
                k = kind[y // 2]
                if k == 0:
                    w.append(y)
                elif k == 2:
                    w.extend(walk[h2f[y ^ 1]])
                y = fp[y]
            walk[g] = w
            proj[c ^ 1] = w
            proj[c] = array('i', [self._ep(x) for x in reversed(w)])

        if relabelling is not None:
            proj = [None if p is None else array('i', [relabelling[x] for x in p]) for p in proj]

        return reduced, proj

    def quad_system(self, forest=None, coforest=None, relabel=False, mapping=False, mutable=False, check=True):
        r"""
        Return the quad system of this map for the given forest and coforest.

        The quad system is obtained from the radial map by contracting the
        edges of ``forest`` and deleting the edges of ``coforest``, each of
        which amounts to folding the two corners of the corresponding
        quadrilateral, see :meth:`fold_corner`. When ``forest`` and
        ``coforest`` are spanning, the result is a quadrangulation of the same
        surface with two vertices of degree `4g`, `4g` edges and `2g`
        quadrilateral faces.

        INPUT:

        - ``forest``, ``coforest`` -- (default: ``None``) the half-edges of the
          edges to contract and to delete, in the format returned by
          :meth:`forest_coforest_decomposition`; entries equal to ``-1`` are
          ignored. When both are ``None`` a decomposition is computed. Listing
          the same edge twice raises a ``ValueError``.

        - ``relabel`` -- boolean (default: ``False``); whether to relabel the
          result on ``0, 1, ..., 2 * ne - 1`` so that it has no inactive
          half-edge. Folding never renumbers, so without this the labels of the
          radial map are kept and are sparse. The relabelling goes edge by edge
          and preserves the parity of each half-edge inside its edge, so that
          the bipartition of the vertices stays visible on the labels.

        - ``mapping`` -- boolean (default: ``False``); whether to also return
          the projection, a list indexed by the half-edges of this map giving
          for each of them the walk of length zero or two it becomes in the
          quad system

        - ``mutable`` -- boolean (default: ``False``); whether the result is
          mutable

        - ``check`` -- boolean (default: ``True``); whether to check that the
          half-edges of ``forest`` and ``coforest`` are ones of this map, and
          whether to check the map built when ``relabel`` is set. Listing the
          same edge twice is reported whatever its value, being caught by the
          bookkeeping rather than by a check.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: m.genus()
            2
            sage: q = m.quad_system()
            sage: q.num_vertices(), q.num_edges(), q.num_faces()
            (2, 8, 4)
            sage: q.vertex_profile()
            [8, 8]
            sage: q.face_profile()
            [4, 4, 4, 4]

        The projection is the mapping of :meth:`radial_map` followed by the
        folds. It sends a half-edge to a walk of length zero or two, the
        length being zero exactly when the folds identify the two half-edges
        of its image in the radial map. That is always the case for a
        half-edge of a contracted edge and for a half-edge of a monogon face,
        and it happens for some half-edges of deleted edges as well. The walk
        of ``ep(h)`` is the reverse of the walk of ``h``::

            sage: q, proj = m.quad_system(mapping=True)
            sage: sorted(set(len(p) for p in proj if p is not None))
            [0, 2]
            sage: all(list(proj[h ^^ 1]) == [proj[h][1] ^^ 1, proj[h][0] ^^ 1]
            ....:     for h in m.half_edges() if proj[h])
            True

        A half-edge of a deleted edge may have empty image too::

            sage: mm = OrientedMap("(0,~2,4,~1,~0,3)(1,~4)(2,~3)", "(0,~1,~4,~2,~3,~0,3,2)(1,4)")
            sage: forest, coforest, _ = mm.forest_coforest_decomposition()
            sage: [h // 2 for h in coforest if h != -1]
            [1]
            sage: q, proj = mm.quad_system(forest, coforest, mapping=True)
            sage: proj[2], proj[3]
            (array('i'), array('i'))

        Folding never renumbers, so the labels of the radial map are kept and
        are sparse; ``relabel`` compacts them, the projection included::

            sage: list(m.quad_system().half_edges())
            [0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 20, 21, 22, 23, 26, 27]
            sage: list(m.quad_system(relabel=True).half_edges())
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
            sage: q, proj = m.quad_system(relabel=True, mapping=True)
            sage: all(x in list(q.half_edges()) for p in proj for x in p)
            True

        .. SEEALSO::

            :meth:`radial_map`, :meth:`fold_corner`,
            :meth:`forest_coforest_decomposition`

        ALGORITHM:

        The folds are performed on a copy of the radial map, in which the edge
        ``e`` of this map is a quadrilateral face whose four sides alternate
        parity, the parity of a half-edge of the radial map being the side of
        the bipartition its tail belongs to. Contracting ``e`` folds the two
        even sides of that quadrilateral and deleting ``e`` folds the two odd
        ones.

        Folding relabels: :meth:`fold_corner` at ``a`` reads ``b = fp[a]`` and
        makes ``b`` take over the position of ``a ^ 1``. So the sides computed
        on the initial radial map may be dead by the time they are needed. To
        avoid resolving names on the fly we keep, for each edge ``e``, one live
        half-edge ``handle[e]`` lying on its quadrilateral together with the
        inverse map ``owner``. The four sides are then recovered in constant
        time by walking the face of the handle, and the two of the wanted
        parity can both be read before folding, since two half-edges of equal
        parity are never opposite to each other. A handle dies only as the
        ``a`` or the ``a ^ 1`` of a fold: in the first case its quadrilateral
        is the one being collapsed, that is the edge currently processed, and
        in the second case ``b`` is its replacement.

        The projection, on the other hand, is only read once every fold has
        been performed. The dying half-edges are linked to the ones they are
        identified with and a single pass with path compression resolves the
        whole forest of links at the end, so that the method is linear.
        """
        if self.has_folded_edge():
            raise NotImplementedError("quad_system is not implemented on a map with a folded edge")

        if forest is None and coforest is None:
            forest, coforest, _ = self.forest_coforest_decomposition()
        elif forest is None or coforest is None:
            raise ValueError("forest and coforest must be given together")

        vp = self._vp
        n = len(vp)

        quad = self.radial_map(mapping=mapping, mutable=True)
        if mapping:
            quad, mor = quad
        qfp = quad.face_permutation(copy=False)

        # handle[e] is a live half-edge of quad lying on the quadrilateral of
        # the edge e of this map and owner is its inverse; this is what keeps
        # track of where each quadrilateral went as folding relabels half-edges
        handle = array('i', [-1] * (n // 2))
        owner = array('i', [-1] * (2 * n))
        for e in range(n // 2):
            if vp[2 * e] != -1:
                handle[e] = 4 * e + 1
                owner[4 * e + 1] = e

        # the half-edges of the radial map that die get linked to the ones they
        # are identified with; only the projection needs them
        link = array('i', [-1] * (2 * n)) if mapping else None

        for parity, edges in ((0, forest), (1, coforest)):
            for h in edges:
                if h == -1:
                    continue
                if check:
                    h = self._check_half_edge(h)
                x = handle[h // 2]
                if x == -1:
                    raise ValueError(f"the edge of the half-edge {h} is listed twice")
                # the handle is dropped before folding: the two fold arguments
                # are all that is needed from this quadrilateral, and a handle
                # left behind would turn into a stale owner entry that a later
                # fold could use to overwrite a live one
                handle[h // 2] = owner[x] = -1

                # the four sides alternate parity; contracting folds the two
                # even ones and deleting the two odd ones. Two half-edges of
                # equal parity are never opposite, so the first fold does not
                # invalidate the second argument.
                a0 = a1 = -1
                for _ in range(4):
                    if x & 1 == parity:
                        if a0 == -1:
                            a0 = x
                        else:
                            a1 = x
                    x = qfp[x]

                # the checks of fold_corner cannot fire here, so skip them:
                # quad is mutable, a is live, and no neighbour ever lies on a
                # folded edge. Every face of the radial map has degree four and
                # a fold takes two from a face, so no face ever reaches degree
                # one and the monogon branch, the only one that folds an edge
                # onto itself, is never taken.
                for a in (a0, a1):
                    b = qfp[a]
                    quad.fold_corner(a, check=0)
                    if b == a ^ 1:
                        # the edge of a bounded its face on both sides and got
                        # pruned, so nothing survives it and there is nothing
                        # to link. This ends the last quadrilateral of a map of
                        # genus zero.
                        owner[a] = owner[a ^ 1] = -1
                        continue
                    # b takes over the position of a ^ 1, hence lies on the
                    # same face
                    f = owner[a ^ 1]
                    if f != -1:
                        handle[f] = b
                        owner[b] = f
                    owner[a] = owner[a ^ 1] = -1
                    if mapping:
                        link[a] = b ^ 1
                        link[a ^ 1] = b

        if mapping:
            # every link is final here, so a single pass with path compression
            # resolves them all, each slot being compressed at most once
            stack = []
            for x in range(2 * n):
                if link[x] == -1:
                    continue
                y = x
                while link[y] != -1:
                    stack.append(y)
                    y = link[y]
                while stack:
                    link[stack.pop()] = y

        relabelling = None
        if relabel:
            qvp = quad.vertex_permutation(copy=False)
            qfp = quad.face_permutation(copy=False)
            nq = len(qvp)
            relabelling = array('i', [-1] * nq)
            ne = 0
            for y in range(0, nq, 2):
                if qvp[y] != -1:
                    relabelling[y] = 2 * ne
                    relabelling[y + 1] = 2 * ne + 1
                    ne += 1
            fp_new = array('i', [-1] * (2 * ne))
            for y in range(nq):
                if qfp[y] != -1:
                    fp_new[relabelling[y]] = relabelling[qfp[y]]
            quad = OrientedMap(fp=fp_new, mutable=True, check=check)

        if not mutable:
            quad.set_immutable()

        if not mapping:
            return quad

        proj = []
        for h in range(n):
            if vp[h] == -1:
                proj.append(None)
                continue
            a = mor[h][0]
            b = mor[h][1]
            if link[a] != -1:
                a = link[a]
            if link[b] != -1:
                b = link[b]
            if b == a ^ 1:
                proj.append(array('i', []))
            elif relabelling is None:
                proj.append(array('i', [a, b]))
            else:
                proj.append(array('i', [relabelling[a], relabelling[b]]))
        return quad, proj

    def geometric_intersection(self, u, v=None, punctured_faces=False):
        r"""
        Return the geometric intersection between two multiwalks or the self-intersection
        of a single one.

        INPUT:

        - ``u`` -- a multiwalk of closed walks, given in one of the following forms
          (see the module documentation):

          - a single walk ``w``
          - a list of walks ``[w0, w1, ...]``
          - a list of pairs (walk, multiplicity) ``[(w0, m0), (w1, m1), ...]``

        - ``v`` -- (default ``None``) an optional multiwalk of closed walks, in one
          of the same forms

        - ``punctured_faces`` -- either a boolean (default ``False``) or a list of
          face indices. If set to ``False`` or ``True`` then respectively no face or
          all faces are considered as punctured. If a list of integers is provided,
          then the faces corresponding to them are considered as punctured.

        OUTPUT: integer. Either the self-intersection of ``u`` if ``v`` is not
        provided or the geometric intersection between ``u`` and ``v``.

        The self-intersection counts each crossing once, while the intersection
        of ``u`` with itself counts it twice: for a closed curve ``c`` one has
        `i(c, c) = 2 i(c)`.

        EXAMPLES:

        A sample of curves on a pair of pants::

            sage: from combisurf import OrientedMap
            sage: P = OrientedMap(fp="(0)(1)(~0,~1)")
            sage: P.num_vertices(), P.num_faces()
            (1, 3)
            sage: w0 = [0, 3]
            sage: w1 = [0, 2, 0, 0, 2]
            sage: w2 = [0, 2, 1, 2]
            sage: P.geometric_intersection(w0, punctured_faces=True)
            1
            sage: P.geometric_intersection([w0], [w2], punctured_faces=True)
            2
            sage: P.geometric_intersection([(w0, 1), (w1, 3)], [(w2, 2)], punctured_faces=True)
            28

        The self-intersection against the intersection with itself::

            sage: P.geometric_intersection(w1, punctured_faces=True)
            2
            sage: P.geometric_intersection(w1, w1, punctured_faces=True)
            4

        .. SEEALSO::

            Each call builds a new
            :class:`~combisurf.geometric_intersection.GeometricIntersectionPairing`,
            which computes a reduced map. To compute many intersection numbers on
            the same map, build the pairing once and call it instead, and to
            compute all the intersection numbers among a fixed list of curves use
            :meth:`~combisurf.geometric_intersection.GeometricIntersectionPairing.intersection_matrix`.
        """
        from combisurf.geometric_intersection import GeometricIntersectionPairing
        return GeometricIntersectionPairing(self, punctured_faces)(u, v)

    def geometric_intersection_matrix(self, wlist, punctured_faces=False):
        r"""
        Return the ``GeometricIntersectionMatrix`` of ``wlist``.

        The :class:`~combisurf.geometric_intersection.GeometricIntersectionMatrix` is a
        lazy object that can be queried to answer pairwise intersections among the curves
        represented by the walks in ``wlist``. It is much faster than repeated calls
        to :meth:`geometric_intersection`.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.lyndon_word_family import cyclically_reduced_lyndon_words
            sage: torus = OrientedMap(fp="(0,1,~0,~1)")
            sage: I = torus.geometric_intersection_matrix(cyclically_reduced_lyndon_words(2, 1, 4, up_to_inverse=True), punctured_faces=True)
            sage: I
            GeometricIntersectionMatrix of 8 curves on OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)")
            sage: print(I.row(3))
            [1, 1, 2, 0, 3, 1, 3, 1]
            sage: print(I.matrix())
            [0 1 1 1 1 1 2 2]
            [1 0 1 1 2 2 1 1]
            [1 1 0 2 1 3 1 3]
            [1 1 2 0 3 1 3 1]
            [1 2 1 3 0 4 3 5]
            [1 2 3 1 4 0 5 3]
            [2 1 1 3 3 5 0 4]
            [2 1 3 1 5 3 4 0]
        """
        from combisurf.geometric_intersection import GeometricIntersectionPairing
        return GeometricIntersectionPairing(self, punctured_faces).matrix(wlist)

    #############
    # Mutations #
    #############

    # TODO: delete and use the more general relabel (with list of cycles)
    #def swap(self, h0, h1, check=2):
        r"""
         Modify the map by multiplying the vertex and face permutations
         by the transposition ``(h0, h1)`` respectively on left and right.

         If the vertices adjacent to the half-edges ``h0``` and ``h1`` are
         distinct they get merged. If they are attached common vertex it is
         split into two vertices. Similarly for the adjacent faces.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap("(0,~0,1,~1)", "(0)(1)(~0,~1)", mutable=True)
            sage: m.swap(0, 1)
            sage: m
            OrientedMap("(0,1,~1)(~0)", "(0,~0,~1)(1)")
            sage: m.swap(0, 2)
            sage: m
            OrientedMap("(0,~1)(~0)(1)", "(0,~0,~1,1)")
            sage: m.swap(0, 3)
            sage: m
            OrientedMap("(0)(~0)(1)(~1)", "(0,~0)(1,~1)")

        This operation also works when the map has folded edges::

            sage: m = OrientedMap("(1,3,~3)", mutable=True)
            sage: m.swap(2, 7)
            sage: m
            OrientedMap("(1)(3,~3)", "(1)(3)(~3)")
            sage: m.swap(6, 7)
            sage: m
            OrientedMap("(1)(3)(~3)", "(1)(3,~3)")
        """
        """if check >= 1:
            self._assert_mutable()
            h0 = self._check_half_edge(h0)
            h1 = self._check_half_edge(h1)

        vp = self._vp
        fp = self._fp

        fp0_pre = self._ep(vp[h0])
        fp1_pre = self._ep(vp[h1])

        vp[h0], vp[h1] = vp[h1], vp[h0]
        fp[fp0_pre], fp[fp1_pre] = fp[fp1_pre], fp[fp0_pre]"""

    def fold_corner(self, h, check=2):
        r"""
        Fold the half-edge ``h`` onto the next half-edge in its face.

        This operation consists in merging the head of ``next_in_face(h)`` with
        the tail of ``h``, which identifies ``h`` with ``ep(next_in_face(h))``
        and ``ep(h)`` with ``next_in_face(h)``. Unless the face of ``h`` is a
        monogon (see below), the resulting map has one edge less and the face
        of ``h`` loses two from its degree. By convention it is the edge of
        ``h`` that disappears.

        That identification rule settles the two degenerate corners as well.
        First, when ``next_in_face(h)`` is ``ep(h)``, both identifications read
        ``h`` with ``h`` and nothing is glued: the edge is pruned (via a call
        to :meth:`delete_edge`). Secondly, when ``next_in_face(h)`` is ``h``
        itself, that is when the face of ``h`` is a monogon, ``h`` is
        identified with ``ep(h)`` and the edge becomes a folded edge. This is
        the one case in which the edge of ``h`` survives the fold.

        Neither the edge of ``h`` nor that of ``next_in_face(h)`` may be
        folded on entry: a ``NotImplementedError`` is raised when one of them
        is. The edge of ``h`` may well become folded by the fold itself,
        through the monogon case above.

        INPUT:

        - ``h`` -- a half-edge, the one whose edge is folded away

        - ``check`` -- integer (default: ``2``); the level of checks to
          perform. Level ``1`` checks that the map is mutable, that ``h`` is
          one of its half-edges and that no folded edge is involved; level
          ``0`` performs none of that and assumes the caller has done it. No
          check is specific to level ``2`` here, which elsewhere in this class
          guards the expensive ones, as in :meth:`genus` and :meth:`relabel`.

        EXAMPLES::

            sage: from combisurf import OrientedMap

        The edge of ``h`` disappears and its face loses two from its degree::

            sage: m = OrientedMap(fp="(0,1,~0,~1)", mutable=True)
            sage: m.num_edges(), m.face_profile()
            (2, [4])
            sage: m.fold_corner(0)
            sage: m
            OrientedMap("(1)(~1)", "(1,~1)")
            sage: m.num_edges(), m.face_profile()
            (1, [2])

        When ``next_in_face(h)`` is ``ep(h)``, that is when the head of ``h``
        has degree one, the edge is pruned::

            sage: m = OrientedMap("(0)(~0)", mutable=True)
            sage: m
            OrientedMap("(0)(~0)", "(0,~0)")
            sage: m.fold_corner(0)
            sage: m
            OrientedMap("", "")

        When the face of ``h`` is a monogon the edge is glued to itself and
        stays on as a folded edge::

            sage: m = OrientedMap("(0,~0)", mutable=True)
            sage: m.face_profile()
            [1, 1]
            sage: m.fold_corner(0)
            sage: m
            OrientedMap("(0)", "(0)")
            sage: m.num_folded_edges()
            1

        .. SEEALSO::

            :meth:`fold_half_edge`, :meth:`delete_edge`
        """
        if check >= 1:
            self._assert_mutable()
            h = self._check_half_edge(h)

        vp = self._vp
        fp = self._fp

        a = h
        a1 = a ^ 1
        b = fp[a]
        b1 = b ^ 1
        if check >= 1 and (vp[a1] == -1 or vp[b1] == -1):
            raise NotImplementedError("fold_corner is not implemented when the edge of h "
                                      "or that of next_in_face(h) is folded")
        if b == a:
            # the face of a is a monogon and the rule identifies a with ep(a),
            # gluing the edge to itself rather than removing it
            self.fold_half_edge(a, check=0)
            return
        if b == a1:
            # a is a leaf half-edge: its head is a vertex of degree one and the
            # corner is bounded by the edge of a on both sides. Folding it
            # prunes that edge, which still takes two from the face of a.
            self.delete_edge(a // 2, check=0)
            return

        # a and b leave the boundary of their face; a is identified with b1 and
        # a1 with b, so that the edge of a disappears. Everything is read
        # before anything is written.
        pa = self._ep(vp[a])          # the position before a in its face
        pa1 = self._ep(vp[a1])        # the position before a1 in its face
        na1 = fp[a1]                  # the position after a1 in its face
        nb = fp[b]                    # the position after b in its face

        # b takes over the position of a1
        val = nb if na1 == a else na1
        if val == a1:
            val = b
        fp[b] = val
        vp[val] = b ^ 1

        # what came before a1 now comes before b
        if pa1 != a and pa1 != a1 and pa1 != b:
            fp[pa1] = b
            vp[b] = self._ep(pa1)

        # the face of a closes over the positions of a and b
        if nb != a and pa != a and pa != a1:
            val = b if nb == a1 else nb
            fp[pa] = val
            vp[val] = self._ep(pa)

        vp[a] = vp[a1] = fp[a] = fp[a1] = -1

        self._clear_trailing_edges()

    def fold_half_edge(self, h, check=2):
        r"""
        Fold the half-edge ``h`` onto its reverse.

        The half-edge ``h`` is identified with ``ep(h)`` so that the edge
        becomes a folded edge, of which only ``2 * (h // 2)`` stays active.
        The position of ``h`` disappears from its face, the surviving
        half-edge takes over the position of ``ep(h)``, and the two ends of
        the edge are merged.

        The result depends on ``h`` and not only on its edge: folding ``h``
        and folding ``ep(h)`` give different maps in general, which is why
        this operation takes a half-edge where :meth:`contract_edge` and
        :meth:`delete_edge` take an edge.

        A ``ValueError`` is raised if the edge of ``h`` is already folded.

        INPUT:

        - ``h`` -- a half-edge, whose position in its face disappears

        - ``check`` -- integer (default: ``2``); the level of checks to
          perform. Level ``1`` checks that the map is mutable, that ``h`` is
          one of its half-edges and that its edge is not already folded; level
          ``0`` performs none of that and assumes the caller has done it. No
          check is specific to level ``2`` here, which elsewhere in this class
          guards the expensive ones, as in :meth:`genus` and :meth:`relabel`.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: t = OrientedMap(vp="(0,1,2)(~0,~1,~2)", mutable=True)
            sage: t.num_folded_edges()
            0
            sage: t.fold_half_edge(0)
            sage: t
            OrientedMap("(0,~1,~2,1,2)", "(0,2,~1,~2,1)")
            sage: t.num_folded_edges()
            1

        Folding the other half-edge of the same edge gives a different map::

            sage: t0 = OrientedMap(vp="(0,~0,1,~1)", mutable=True)
            sage: t0.fold_half_edge(0)
            sage: t0
            OrientedMap("(0,1,~1)", "(0,~1)(1)")
            sage: t1 = OrientedMap(vp="(0,~0,1,~1)", mutable=True)
            sage: t1.fold_half_edge(1)
            sage: t1
            OrientedMap("(0)(1,~1)", "(0)(1)(~1)")

        The two ends of the edge are merged, and the Euler characteristic is
        preserved because a folded edge carries a vertex of its own::

            sage: t.num_vertices()
            1
            sage: OrientedMap(vp="(0,1,2)(~0,~1,~2)").euler_characteristic()
            0
            sage: t.euler_characteristic()
            0

        .. SEEALSO::

            :meth:`fold_corner`
        """
        if check >= 1:
            self._assert_mutable()
            h = self._check_half_edge(h)

        vp = self._vp
        fp = self._fp

        if check >= 1 and vp[h ^ 1] == -1:
            raise ValueError(f"the edge of the half-edge {h} is already folded")

        s = h & ~1                    # survives, the even half-edge
        d = h | 1                     # dies

        # h leaves the boundary of its face and is identified with its reverse.
        # Everything is read before anything is written. Note that the edge of
        # s is folded in the result, so that ep(s) is s.
        ph = self._ep(vp[h])          # the position before h in its face
        ph1 = self._ep(vp[h ^ 1])     # the position before h ^ 1 in its face
        nh = fp[h]                    # the position after h in its face
        nh1 = fp[h ^ 1]               # the position after h ^ 1 in its face

        # s takes over the position of h ^ 1
        val = nh if nh1 == h else nh1
        if val == d:
            val = s
        fp[s] = val
        vp[val] = s

        # what came before h ^ 1 now comes before s
        if s == h and ph1 != h and ph1 != d:
            fp[ph1] = s
            vp[s] = self._ep(ph1)

        # the face of h closes over the position of h
        if ph != h and ph != d:
            val = s if nh == d else nh
            fp[ph] = val
            # ep(s) is s in the result, but vp[d] is only cleared below, so
            # _ep does not see the fold yet and ph == s is settled by hand
            vp[val] = s if ph == s else self._ep(ph)

        vp[d] = fp[d] = -1

        self._clear_trailing_edges()

    def contract_edge(self, e, check=2):
        r"""
        Contract the edge ``e``.

        Inverse operation of :meth:`split_vertex` except that here we allow
        vertices of degree one.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: vp = "(0,1,~0,~1)"
            sage: fp = "(0,1,~0,~1)"

            sage: vp02 = "(0,2,~0,~1)(1,~2)"
            sage: fp02 = "(0,2,1,~0,~1,~2)"
            sage: m = OrientedMap(vp02, fp02, mutable=True)
            sage: m.contract_edge(2)
            sage: m == OrientedMap(vp, fp)
            True

            sage: vp01 = "(0,2,~1)(~0,~2,1)"
            sage: fp01 = "(0,1,2,~0,~1,~2)"
            sage: m = OrientedMap(vp01, fp01, mutable=True)
            sage: m.contract_edge(2)
            sage: m == OrientedMap(vp, fp)
            True

            sage: vp03 = "(0,2)(~0,~1,~2,1)"
            sage: fp03 = "(0,1,~0,2,~1,~2)"
            sage: m = OrientedMap(vp03, fp03, mutable=True)
            sage: m.contract_edge(2)
            sage: m == OrientedMap(vp, fp)
            True

        Degree 1 vertices::

            sage: m = OrientedMap(vp="(0,1)(~0)(~1)", mutable=True)
            sage: m.contract_edge(1)
            sage: m
            OrientedMap("(0)(~0)", "(0,~0)")

            sage: m = OrientedMap("(0,1)(~0)(~1)", mutable=True)
            sage: m.contract_edge(0)
            sage: m
            OrientedMap("(1)(~1)", "(1,~1)")

            sage: m = OrientedMap(vp="(0,~1)(~0)(1)", mutable=True)
            sage: m.contract_edge(0)
            sage: m
            OrientedMap("(1)(~1)", "(1,~1)")

            sage: m = OrientedMap(vp="(0,~1)(~0)(1)", mutable=True)
            sage: m.contract_edge(1)
            sage: m
            OrientedMap("(0)(~0)", "(0,~0)")

            sage: m = OrientedMap(vp="(0)", mutable=True)
            sage: m.contract_edge(0)
            sage: m
            OrientedMap("", "")

            sage: m = OrientedMap(vp="(0,1)", mutable=True)
            sage: m.contract_edge(1)
            sage: m
            OrientedMap("(0)", "(0)")

            sage: m = OrientedMap(vp="(0,1)", mutable=True)
            sage: m.contract_edge(0)
            sage: m
            OrientedMap("(1)", "(1)")

            sage: m = OrientedMap(vp="(0,~0)", mutable=True)
            sage: m.contract_edge(0)
            sage: m
            OrientedMap("", "")

            sage: m = OrientedMap(vp="(0,~0,1,~1)", mutable=True)
            sage: m.contract_edge(0)
            sage: m
            OrientedMap("(1,~1)", "(1)(~1)")

            sage: m = OrientedMap(vp="(0,1,~1,~0)", mutable=True)
            sage: m.contract_edge(0)
            sage: m
            OrientedMap("(1,~1)", "(1)(~1)")
        """
        if check >= 1:
            self._assert_mutable()
            e = self._check_edge(e)

        vp = self._vp
        fp = self._fp

        h0 = 2 * e
        h1 = self._ep(h0)

        if h0 == h1:
            # folded edge
            if fp[h0] == h0:
                # the vertex at h0 has degree one
                vp[h0] = fp[h0] = -1
            else:
                h0_fp_prev = self.previous_in_face(h0)
                h0_fp_next = self.next_in_face(h0)
                assert h0 != h0_fp_prev
                assert h0 != h0_fp_next
                fp[h0_fp_prev] = h0_fp_next

                h0_vp_prev = self.previous_at_vertex(h0)
                h0_vp_next = self.next_at_vertex(h0)
                assert h0 != h0_vp_prev
                assert h0 != h0_vp_next
                vp[h0_vp_prev] = h0_vp_next

            vp[h0] = fp[h0] = -1

        else:
            # non-folded edge
            h0_vp_prev = self.previous_at_vertex(h0)
            h0_vp_next = self.next_at_vertex(h0)
            h1_vp_next = self.next_at_vertex(h1)
            h0_fp_prev = self.previous_in_face(h0)
            if h0_fp_prev == h1:
                h0_fp_prev = self.previous_in_face(h0_fp_prev)
            h0_fp_next = self.next_in_face(h0)
            if h0_fp_next == h1:
                h0_fp_next = self.next_in_face(h0_fp_next)

            h1_vp_prev = self.previous_at_vertex(h1)
            h1_vp_next = self.next_at_vertex(h1)
            h1_fp_prev = self.previous_in_face(h1)
            if h1_fp_prev == h0:
                h1_fp_prev = self.previous_in_face(h1_fp_prev)
            h1_fp_next = self.next_in_face(h1)
            if h1_fp_next == h0:
                h1_fp_next = self.next_in_face(h1_fp_next)

            if h0 != h0_fp_next:
                fp[h0_fp_prev] = h0_fp_next
            if h1 != h1_fp_next:
                fp[h1_fp_prev] = h1_fp_next

            if h0 != h0_vp_next:
                if h1 != h1_vp_next:
                    # vertices at h0 and h1 have degree > 1
                    vp[h0_vp_prev] = h1_vp_next
                    vp[h1_vp_prev] = h0_vp_next
                else:
                    # vertex at h1 has degree 1 but not h0
                    vp[h0_vp_prev] = h0_vp_next
            else:
                if h1 != h1_vp_next:
                    # vertex at h0 has degree 1 but not h1
                    vp[h1_vp_prev] = h1_vp_next

            vp[h0] = vp[h1] = fp[h0] = fp[h1] = -1

        self._clear_trailing_edges()

    def delete_edge(self, e, check=2):
        r"""
        Delete the edge ``e``.

        EXAMPLES::

            sage: from combisurf import OrientedMap

        In each of the examples below, the edge ``2`` is attached to the map

            sage: vp = "(0,1,~0,~1)"
            sage: fp = "(0,1,~0,~1)"
            sage: m = OrientedMap(vp, fp)
            sage: m
            OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)")

        in a different configuration, and deleting it gives that map back::

            sage: vp20 = "(0,~2,2,1,~0,~1)"
            sage: fp20 = "(0,1,~0,~1,2)(~2)"
            sage: m = OrientedMap(vp20, fp20, mutable=True)
            sage: m.delete_edge(2)
            sage: m
            OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)")

            sage: vp10 = "(0,2,1,~2,~0,~1)"
            sage: fp10 = "(0,~2)(~0,~1,2,1)"
            sage: m = OrientedMap(vp10, fp10, mutable=True)
            sage: m.delete_edge(2)
            sage: m
            OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)")

            sage: vp30 = "(0,2,1,~0,~2,~1)"
            sage: fp30 = "(0,1,~2)(~0,~1,2)"
            sage: m = OrientedMap(vp30, fp30, mutable=True)
            sage: m.delete_edge(2)
            sage: m
            OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)")

            sage: vp00 = "(0,~2,2,1,~0,~1)"
            sage: fp00 = "(0,1,~0,~1,2)(~2)"
            sage: m = OrientedMap(vp00, fp00, mutable=True)
            sage: m.delete_edge(2)
            sage: m
            OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)")

            sage: vp22 = "(0,1,~2,2,~0,~1)"
            sage: fp22 = "(0,2,1,~0,~1)(~2)"
            sage: m = OrientedMap(vp22, fp22, mutable=True)
            sage: m.delete_edge(2)
            sage: m
            OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)")
        """
        if check >= 1:
            self._assert_mutable()
            e = self._check_edge(e)

        vp = self._vp
        fp = self._fp

        h0 = 2 * e
        h1 = self._ep(h0)

        if h0 == h1:
            # folded edge (same behavior as with contract_edge)
            if fp[h0] == h0:
                # the vertex at h0 has degree one
                vp[h0] = fp[h0] = -1
            else:
                h0_fp_prev = self.previous_in_face(h0)
                h0_fp_next = self.next_in_face(h0)
                assert h0 != h0_fp_prev
                assert h0 != h0_fp_next
                fp[h0_fp_prev] = h0_fp_next

                h0_vp_prev = self.previous_at_vertex(h0)
                h0_vp_next = self.next_at_vertex(h0)
                assert h0 != h0_vp_prev
                assert h0 != h0_vp_next
                vp[h0_vp_prev] = h0_vp_next

            vp[h0] = fp[h0] = -1

        else:
            # non-folded edge
            h0_vp_prev = self.previous_at_vertex(h0)
            if h0_vp_prev == h1:
                h0_vp_prev = self.previous_at_vertex(h0_vp_prev)
            h0_vp_next = self.next_at_vertex(h0)
            if h0_vp_next == h1:
                h0_vp_next = self.next_at_vertex(h0_vp_next)
            h0_fp_prev = self.previous_in_face(h0)
            h0_fp_next = self.next_in_face(h0)

            h1_vp_prev = self.previous_at_vertex(h1)
            if h1_vp_prev == h0:
                h1_vp_prev = self.previous_at_vertex(h1_vp_prev)
            h1_vp_next = self.next_at_vertex(h1)
            if h1_vp_next == h0:
                h1_vp_next = self.next_at_vertex(h1_vp_next)
            h1_fp_prev = self.previous_in_face(h1)
            h1_fp_next = self.next_in_face(h1)

            if h0 != h0_vp_next:
                vp[h0_vp_prev] = h0_vp_next
            if h1 != h1_vp_next:
                vp[h1_vp_prev] = h1_vp_next

            if h0 != h0_fp_next:
                if h1 != h1_fp_next:
                    # faces at h0 and h1 have degree > 1
                    fp[h0_fp_prev] = h1_fp_next
                    fp[h1_fp_prev] = h0_fp_next
                else:
                    # face at h1 has degree 1 but not h0
                    fp[h0_fp_prev] = h0_fp_next
            else:
                if h1 != h1_vp_next:
                    # face at h0 has degree 1 but not h1
                    fp[h1_fp_prev] = h1_fp_next

            vp[h0] = vp[h1] = fp[h0] = fp[h1] = -1

        self._clear_trailing_edges()

    def add_edge(self, h0=-1, h1=-1, e=None, check=2):
        r"""
        Add an edge between the corners adjacent to the half-edges ``h0`` and ``h1``.

        This operation is the inverse of :meth:`delete_edge` and dual to :meth:`insert_edge`.
        See also :meth:`contract_edge`.

        INPUT:

        -- ``h0``: integer -- ``-1`` or valid half-edge index. If set to ``-1``
           then the newly added edge is attached to a new vertex.

        -- ``h1``: integer -- ``-1`` or valid half-edge index. If set to ``-1``
           then the newly added edge is attached to a new vertex.

        -  `e``: integer -- an optional edge index for the newly created edge.

        - ``check`` -- whether to check input consistency

        EXAMPLES:

            sage: from combisurf import OrientedMap

        We start from a torus and triangulate it by adding an edge in the middle of the face::

            sage: m = OrientedMap(fp="(0,1,~0,~1)", mutable=True)
            sage: m.add_edge(0, 1)
            sage: m
            OrientedMap("(0,2,1,~0,~2,~1)", "(0,1,~2)(~0,~1,2)")
        """
        if check >= 1:
            self._assert_mutable()
            h0 = self._check_half_edge_or_negative(h0)
            if h1 is not None:
                h1 = self._check_half_edge_or_negative(h1)

        if e is None:
            e = len(self._vp) // 2
        else:
            if not isinstance(e, numbers.Integral):
                raise ValueError("e must be an integer")
            e = int(e)
            if 2 * e > len(self._vp):
                raise NotImplementedError
            if self._vp[2 * e] != -1:
                raise ValueError(f"invalid edge e (={e}) already in use")

        if 2 * e == len(self._vp):
            # make extra space for the new edge
            self._vp.append(-1)
            self._vp.append(-1)
            self._fp.append(-1)
            self._fp.append(-1)

        vp = self._vp
        fp = self._fp

        if h0 < 0:
            if h1 < 0:
                if h0 == h1:
                    # independent loop
                    vp[2 * e] = 2 * e + 1
                    vp[2 * e + 1] = 2 * e
                    fp[2 * e] = 2 * e
                    fp[2 * e + 1] = 2 * e + 1
                else:
                    # independent edge
                    vp[2 * e] = 2 * e
                    vp[2 * e + 1] = 2 * e + 1
                    fp[2 * e] = 2 * e + 1
                    fp[2 * e + 1] = 2 * e
            else:
                # edge attached only at h1
                h1_fp_prev = self.previous_in_face(h1)

                vp[2 * e] = 2 * e
                vp[2 * e + 1] = vp[h1]
                vp[h1] = 2 * e + 1

                fp[h1_fp_prev] = 2 * e + 1
                fp[2 * e + 1] = 2 * e
                fp[2 * e] = h1

        elif h1 < 0:
            # edge attached only at h0
            h0_fp_prev = self.previous_in_face(h0)

            vp[2 * e + 1] = 2 * e + 1
            vp[2 * e] = vp[h0]
            vp[h0] = 2 * e

            fp[h0_fp_prev] = 2 * e
            fp[2 * e] = 2 * e + 1
            fp[2 * e + 1] = h0

        elif h0 == h1:
            # loop at h0
            h0_fp_prev = self.previous_in_face(h0)

            vp[2 * e] = 2 * e + 1
            vp[2 * e + 1] = vp[h0]
            vp[h0] = 2 * e

            fp[h0_fp_prev] = 2 * e + 1
            fp[2 * e + 1] = h0
            fp[2 * e] = 2 * e

        else:
            # edge attached at h0 and h1
            h0_fp_prev = self.previous_in_face(h0)
            h1_fp_prev = self.previous_in_face(h1)

            vp[2 * e] = vp[h0]
            vp[h0] = 2 * e
            vp[2 * e + 1] = vp[h1]
            vp[h1] = 2 * e + 1

            fp[h1_fp_prev] = 2 * e + 1
            fp[2 * e + 1] = h0
            fp[h0_fp_prev] = 2 * e
            fp[2 * e] = h1

    def insert_edge(self, h0=-1, h1=-1, e=None, check=2):
        r"""
        Add an edge.

        This operation is the inverse of :meth:`contract_edge` and dual to
        :meth:`add_edge`. See also :meth:`delete_edge`.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(fp="(0,1,~0,~1)", mutable=True)
            sage: m.insert_edge(0, 1)
            sage: m
            OrientedMap("(0,1,~2)(~0,~1,2)", "(0,2,1,~0,~2,~1)")

            sage: G = OrientedMap(vp = [[0, 2], [1, 4], [3, 5]], mutable=True)
            sage: G_dual = G.dual()
            sage: G.add_edge(0, 4)
            sage: G_dual.insert_edge(0, 4)
            sage: G_dual == G.dual()
            True

        """

        if check >= 1:
            self._assert_mutable()
            h0 = self._check_half_edge_or_negative(h0)
            if h1 is not None:
                h1 = self._check_half_edge_or_negative(h1)

        if e is None:
            e = len(self._vp) // 2
        else:
            if not isinstance(e, numbers.Integral):
                raise ValueError("e must be an integer")
            e = int(e)
            if 2 * e > len(self._vp):
                raise NotImplementedError
            if self._vp[2 * e] != -1:
                raise ValueError(f"invalid edge e (={e}) already in use")

        if 2 * e == len(self._vp):
            # make extra space for the new edge
            self._vp.append(-1)
            self._vp.append(-1)
            self._fp.append(-1)
            self._fp.append(-1)

        vp = self._vp
        fp = self._fp

        if h0 < 0:
            if h1 < 0:
                if h0 != h1:
                    # independent loop
                    vp[2 * e] = 2 * e + 1
                    vp[2 * e + 1] = 2 * e
                    fp[2 * e] = 2 * e
                    fp[2 * e + 1] = 2 * e + 1
                else:
                    # independent edge
                    vp[2 * e] = 2 * e
                    vp[2 * e + 1] = 2 * e + 1
                    fp[2 * e] = 2 * e + 1
                    fp[2 * e + 1] = 2 * e
            else:
                # loop attached after h1
                fp[2 * e] = 2 * e
                fp[2 * e + 1] = fp[h1]
                fp[h1] = 2 * e + 1

                vp[2 * e + 1] = self._ep(h1)
                vp[2 * e] = 2 * e + 1
                vp[fp[2 * e + 1]] = 2 * e

        elif h1 < 0:
            # loop attached after h0
            fp[2 * e + 1] = 2 * e + 1
            fp[2 * e] = fp[h0]
            fp[h0] = 2 * e

            vp[2 * e] = self._ep(h0)
            vp[2 * e + 1] = 2 * e
            vp[fp[2 * e]] = 2 * e + 1

        elif h0 == h1:
            # insert an edge with a vertex after h0
            fp[2 * e] = 2 * e + 1
            fp[2 * e + 1] = fp[h0]
            fp[h0] = 2 * e

            vp[fp[2 * e + 1]] = 2 * e
            vp[2 * e] = self._ep(h0)
            vp[2 * e + 1] = 2 * e + 1

        else:
            # insert an edge between h0 and h1
            fp[2 * e + 1] = fp[h1]
            fp[2 * e] = fp[h0]
            fp[h0] = 2 * e
            fp[h1] = 2 * e + 1

            vp[2 * e] = self._ep(h0)
            vp[2 * e + 1] = self._ep(h1)
            vp[fp[2 * e + 1]] = 2 * e
            vp[fp[2 * e]] = 2 * e + 1

    def reverse_orientation(self, e, check=True):
        r"""
        Change the orientation of the edge ``e``.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m = OrientedMap(fp="(0,1,2)(~0,~1,~2)", mutable=True)
            sage: m.reverse_orientation(0)
            sage: m
            OrientedMap("(0,2,~1,~0,~2,1)", "(0,~1,~2)(~0,1,2)")
            sage: m.reverse_orientation(1)
            sage: m
            OrientedMap("(0,2,1,~0,~2,~1)", "(0,1,~2)(~0,~1,2)")
            sage: m.reverse_orientation(2)
            sage: m
            OrientedMap("(0,~2,1,~0,2,~1)", "(0,1,2)(~0,~1,~2)")

            sage: m = OrientedMap(fp="(0,~5,4)(3,5,6)(1,2,~6)", mutable=True)
            sage: m.reverse_orientation(0)
            sage: m
            OrientedMap("(0,4,5,3,~6,2,1,6,~5)", "(0,~5,4)(1,2,~6)(3,5,6)")
            sage: m.reverse_orientation(5)
            sage: m
            OrientedMap("(0,4,~5,3,~6,2,1,6,5)", "(0,5,4)(1,2,~6)(3,~5,6)")
            sage: m._check()

        One can alternatively use ``relabel`` (which would be slower in that
        situation)::

            sage: m = OrientedMap(fp="(0,~5,4)(1,2,~6)(3,5,6)", mutable=True)
            sage: m1 = m.copy()
            sage: m1.reverse_orientation(5)
            sage: m1.reverse_orientation(6)
            sage: m2 = m.copy()
            sage: m2.relabel("(5,~5)(6,~6)")
            sage: m1 == m2
            True
        """
        if check:
            self._assert_mutable()
            e = self._check_edge(e)

        ep = self._ep
        h = 2 * e
        H = self._ep(h)
        if h == H:
            return

        h_vp_inv = self.previous_at_vertex(h)
        H_vp_inv = self.previous_at_vertex(H)
        h_fp_inv = self.previous_in_face(h)
        H_fp_inv = self.previous_in_face(H)

        perm_conjugate_transposition_inplace(self._vp, h, H, h_vp_inv, H_vp_inv)
        perm_conjugate_transposition_inplace(self._fp, h, H, h_fp_inv, H_fp_inv)

    # TODO: if given as cycles, can we make it O(input size)?
    def relabel(self, p=None, check=2):
        r"""
        Relabel this triangulation inplace.

        If the permutation ``p`` is provided as argument, then it used for
        relabelling. Otherwise, the map is relabelled so that it uses
        consecutive integers starting at ``0`` for its edges.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: m =OrientedMap("(0,5,3,~3)(~0)", "(0,~0,~3,5)(3)", mutable=True)
            sage: m.relabel("(0,~0)")
            sage: m
            OrientedMap("(0)(~0,5,3,~3)", "(0,~3,5,~0)(3)")
            sage: m.relabel("(0,1,~2)(~0,~1,2)")
            sage: m
            OrientedMap("(1)(~1,5,3,~3)", "(1,~3,5,~1)(3)")

            sage: m.set_immutable()
            sage: m.relabel("(0,~1)")
            Traceback (most recent call last):
            ...
            ValueError: immutable map; use a mutable copy instead

        An example where ``p`` is not provided::

            sage: m = OrientedMap(vp="(0,5,3,~3)(~0)", mutable=True)
            sage: m.relabel()
            sage: m
            OrientedMap("(0,2,1,~1)(~0)", "(0,~0,~1,2)(1)")
        """
        # NOTE: if this operation is used often, it would not be too complicated
        # to make it inplace
        if check >= 1:
            self._assert_mutable()

        if p is None:
            p = array('i', [-1] * len(self._vp))
            ee = 0
            for e in range(len(self._vp) // 2):
                if self._vp[2 * e] != -1:
                    p[2 * e] = 2 * ee
                    p[2 * e + 1] = 2 * ee + 1
                    ee += 1

        elif check >= 2:
            p = check_relabelling(p, len(self._vp) // 2)

        self._vp = perm_conjugate(self._vp, p)
        self._fp = perm_conjugate(self._fp, p)

        self._clear_trailing_edges()

    # TODO: should we make it possible to choose the triangulation? Right now
    # we just pick the dual to a path.
    def triangulate(self, h=None):
        r"""
        Triangulate the faces of degree more than 3.

        If a half-edge ``h`` is specified as input, only triangulate the face
        adjacent to it.

        EXAMPLES::

            sage: from combisurf import OrientedMap

        We start from a 5-gon and triangulate first one of its faces and then
        the other one::

            sage: m = OrientedMap(fp="(0,1,2,3,4)(~0,~4,~3,~2,~1)", mutable=True)
            sage: m.triangulate(0)
            sage: m
            OrientedMap("(0,~4)(~0,1,6)(~1,2,5)(~2,3)(~3,4,~6,~5)", "(0,6,4)(~0,~4,~3,~2,~1)(1,5,~6)(2,3,~5)")
            sage: m.triangulate(1)
            sage: m
            OrientedMap("(0,~4,8)(~0,1,6)(~1,~8,~7,2,5)(~2,3)(~3,7,4,~6,~5)", "(0,6,4)(~0,8,~1)(1,5,~6)(2,3,~5)(~2,~7,~3)(~4,7,~8)")

        Trying to triangulate a monogon or a bigon generates a ValueError::

            sage: m = OrientedMap(fp="(0)(~0)", mutable=True)
            sage: m.triangulate(0)
            Traceback (most recent call last):
            ...
            ValueError: can not triangulate a monogon or a bigon
            sage: m = OrientedMap(fp="(0,1)(~0,~1)", mutable=True)
            sage: m.triangulate(0)
            Traceback (most recent call last):
            ...
            ValueError: can not triangulate a monogon or a bigon

        Immutable maps can not be triangulated::

            sage: m = OrientedMap("(0,1,2)(~0,~2,~1)")
            sage: m.triangulate(0)
            Traceback (most recent call last):
            ...
            ValueError: immutable map; use a mutable copy instead
        """
        self._assert_mutable()

        if h is None:
            for f in self.faces():
                if len(f) > 3:
                    self.triangulate(f[0])

        f = perm_orbit(self._fp, h)
        if len(f) == 1 or len(f) == 2:
            raise ValueError("can not triangulate a monogon or a bigon")

        i = 1
        j = len(f) - 1
        while i < j:
            if i + 2 < j:
                self.add_edge(f[i + 1], f[j])
            if i + 1 < j:
                self.add_edge(f[i], f[j])
            i += 1
            j -= 1

    # TODO: consider listing all quotients by looking at blocks under the monodromy group
    def automorphism_quotient(self, mapping=False, mutable=False, check=True):
        r"""
        Return the quotient under the automorphism group.

        EXAMPLES::

            sage: from veerer import *  # not tested

        Veering triangulation example::

            sage: vt = VeeringTriangulation("(0,1,2)(3,4,~0)(5,6,~1)(7,~2,8)(9,~3,~6)(10,~7,~4)(11,~5,12)(13,14,~8)(15,~9,16)(17,18,~10)(19,~17,~11)(20,~13,~12)(21,~14,~18)(22,~21,~15)(23,24,~16)(25,~23,~19)(26,~20,~25)(~26,~24,~22)", "RBBRBRBRRBRBBRBBRRBRRRBBRRB")  # not tested
            sage: len(vt.automorphisms())  # not tested
            2
            sage: qvt = vt.automorphism_quotient()  # not tested
            sage: qvt  # not tested
            VeeringTriangulation("(0,1,2)(~0,3,4)(~1,5,6)(~2,8,7)(~3,~6,9)(~4,10,~7)(~5,12,11)(~8,13,~12)(~10,14,~11)", "RBBRBRBRRBRBBRR")
            sage: (vt.stratum(), qvt.stratum())  # not tested
            (H_4(2^3), Q_1(1^3, -1^3))

            sage: vt.automorphism_quotient(mapping=True)  # not tested
            (VeeringTriangulation("(0,1,2)(~0,3,4)(~1,5,6)(~2,8,7)(~3,~6,9)(~4,10,~7)(~5,12,11)(~8,13,~12)(~10,14,~11)", "RBBRBRBRRBRBBRR"),
             array('i', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 18, 20, 21, 22, 23, 24, 25, 26, 26, 25, 24, 13, 12, 7, 6, 28, 28, 23, 22, 21, 20, 17, 16, 11, 10, 3, 2, 8, 9, 1, 0, 15, 14, 5, 4]))

        Strebel graph example::

            sage: sg = StrebelGraph("(0,~0,~1)(1,2,~2)")  # not tested
            sage: sg.automorphism_quotient()  # not tested
            StrebelGraph("(0,~0,1)")

        TESTS::

            sage: vt = VeeringTriangulation("(~0,1,2)(~1,3,4)(~2,5,6)(~3,7,8)(~4,~7,9)(~6,10,11)(~8,12,13)(~9,14,15)(~10,16,17)(~11,18,~17)(~12,19,~18)(~14,20,~16)(0:1)(~5:1)(~13:1)(~15:1)(~19:1)(~20:1)", "RBBBRRBRBBRRBRBRBBBRR")  # not tested
            sage: vt.automorphism_quotient()  # not tested
            VeeringTriangulation("(~0,1,2)(~1,3,4)(~2,5,6)(~3,7,8)(~4,~7,~6)(~8,9,10)(0:1)(~5:1)(~10:1)", "RBBBRRBRBBR")
        """
        return self.quotient(perms_orbits(self.automorphisms()), mapping, mutable, check)

    def quotient(self, blocks, mapping=False, mutable=False, check=True):
        if check:
            if not all(blocks):
                raise ValueError("each block must be non empty")
            blocks = [[self._check_half_edge(h) for h in block] for block in blocks]
            half_edges = set().union(*blocks)
            if half_edges != set(self.half_edges()):
                raise ValueError("invalid blocks")
            for block in blocks:
                for l in self._half_edges_data:
                    if len(set(l[h] for h in block)) != 1:
                        raise ValueError("block must be constant on half-edges data")
                for l in self._edges_data:
                    if len(set(l[h // 2] for h in block)) != 1:
                        raise ValueError("block must be constant on edges data")

        half_edge_to_block = [-1] * len(self._vp)
        for i, block in enumerate(blocks):
            for j in block:
                half_edge_to_block[j] = i

        if check:
            for block in blocks:
                if len(set(half_edge_to_block[self._vp[h]] for h in block)) != 1:
                    raise ValueError("invalid blocks")
                if len(set(half_edge_to_block[self._fp[h]] for h in block)) != 1:
                    raise ValueError("invalid blocks")

        ne = 0
        block_relabelling = {-1: -1}
        for e in range(len(self._vp) // 2):
            i = half_edge_to_block[2 * e]
            if i in block_relabelling:
                continue
            block_relabelling[i] = 2 * ne
            if self._vp[2 * e + 1] != -1:
                ii = half_edge_to_block[2 * e + 1]
                if ii in block_relabelling:
                    assert i == ii  # folding
                else:
                    block_relabelling[ii] = 2 * ne + 1
            ne += 1

        vp = array('i', [-1] * (2 * ne))
        fp = array('i', [-1] * (2 * ne))
        half_edges_data = [array('i', [0] * (2 * ne)) for _ in self._half_edges_data]
        edges_data = [array('i', [0] * ne) for _ in self._edges_data]
        for i, block in enumerate(blocks):
            h = block[0]
            ii = half_edge_to_block[self._vp[h]]
            vp[block_relabelling[i]] = block_relabelling[ii]

            ii = half_edge_to_block[self._fp[h]]
            fp[block_relabelling[i]] = block_relabelling[ii]

            for ldest, lsrc in zip(half_edges_data, self._half_edges_data):
                ldest[block_relabelling[i]] = lsrc[h]
            for ldest, lsrc in zip(edges_data, self._edges_data):
                ldest[block_relabelling[i] // 2] = lsrc[h // 2]

        quotient = self.from_permutations(vp, fp, half_edges_data, edges_data, mutable, check)
        return (quotient, array('i', [block_relabelling[half_edge_to_block[h]] for h in range(len(self._vp))])) if mapping else quotient

    #################
    # Automorphisms #
    #################

    def _relabelling_from(self, root):
        r"""
        When connected, return a canonical relabelling map obtained from walking
        along the triangulation starting at ``root``.

        The returned relabelling array maps the current edge to the new
        labelling.

        EXAMPLES::

            sage: from veerer import *  # not tested
            sage: from array import array  # not tested

        The torus example (6 symmetries)::

            sage: fp = array('i', [2, 3, 5, 4, 1, 0])  # not tested
            sage: vp = array('i', [4, 5, 1, 0, 2, 3])  # not tested
            sage: T = Triangulation.from_permutations(vp, fp, (array('i', [0]*6),), mutable=True)  # not tested
            sage: T._relabelling_from(3)  # not tested
            array('i', [5, 4, 1, 0, 2, 3])

            sage: p = T._relabelling_from(0)  # not tested
            sage: T.relabel(p)  # not tested
            sage: for i in range(6):  # not tested
            ....:     p = T._relabelling_from(i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     assert S == T

        The sphere example (3 symmetries)::

            sage: fp = array('i', [2, -1, 4, -1, 0, -1])  # not tested
            sage: vp = array('i', [4, -1, 0, -1, 2, -1])  # not tested
            sage: T = Triangulation.from_permutations(vp, fp, (array('i', [0]*6),), mutable=True)  # not tested
            sage: T._relabelling_from(2)  # not tested
            array('i', [4, 5, 0, 1, 2, 3])
            sage: p = T._relabelling_from(0)  # not tested
            sage: T.relabel(p)  # not tested
            sage: for i in range(3):  # not tested
            ....:     p = T._relabelling_from(2 * i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     assert S == T

        An example with no automorphism::

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)", mutable=True)  # not tested
            sage: p = T._relabelling_from(0)  # not tested
            sage: T.relabel(p)  # not tested
            sage: for i in T.half_edges():  # not tested
            ....:     if i == 0: continue
            ....:     p = T._relabelling_from(i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     S._check()
            ....:     assert S != T
        """
        root = self._check_half_edge(root)
        relabelling = array('i', [-1] * len(self._vp))
        fp_new = array('i', [-1] * len(self._vp))
        last = edge_relabelling_from(relabelling, fp_new, self._fp, len(self._vp), root, 0)
        assert fp_new == perm_conjugate(self._fp, relabelling), (fp_new, perm_conjugate(self._fp, relabelling))
        if last != len(self._vp):
            raise ValueError("non-connected constellation")
        return relabelling

    def automorphisms(self):
        r"""
        Return the list of automorphisms of this constellation.

        The output is a list of arrays that are permutations acting on the set
        of half edges.

        For triangulations with boundaries, we allow automorphism to permute
        boundaries. Though, boundary edge have to be mapped on boundary edge.

        EXAMPLES::

            sage: from veerer import *  # not tested

        An example with 4 symmetries in genus 2::

            sage: T = Triangulation("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)")  # not tested
            sage: A = T.automorphisms()  # not tested
            sage: len(A)  # not tested
            4

        And the "sphere octagon" has 8::

            sage: s  = "(0,8,~7)(1,9,~0)(2,10,~1)(3,11,~2)(4,12,~3)(5,13,~4)(6,14,~5)(7,15,~6)"  # not tested
            sage: len(Triangulation(s).automorphisms())  # not tested
            8

        A veering triangulation with 4 symmetries in genus 2::

            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"  # not tested
            sage: cols = "BRBBBRRBBBBR"  # not tested
            sage: V = VeeringTriangulation(fp, cols)  # not tested
            sage: A = V.automorphisms()  # not tested
            sage: len(A)  # not tested
            4

        Examples with boundaries::

            sage: t = Triangulation("(0,1,2)", boundary="(~0:1)(~1:1)(~2:1)")  # not tested
            sage: len(t.automorphisms())  # not tested
            3
            sage: t = Triangulation("(0,1,2)", boundary="(~0:1,~1:1,~2:1)")  # not tested
            sage: len(t.automorphisms())  # not tested
            3
            sage: t = Triangulation("(0,1,2)", boundary="(~0:1,~1:1,~2:2)")  # not tested
            sage: len(t.automorphisms())  # not tested
            1

        Linear families::

            sage: s = StrebelGraph("(0,3,7,~6,~2,1)(2,5,~4,~3,~1,~0)(4,8,~5)(6,~8,~7)")  # not tested
            sage: f = StrebelGraphLinearFamily(s, [(2, 0, 0, 0, 1, 0, 1, 0, 2), (0, 2, 0, 0, 0, 1, 0, 1, 2), (0, 0, 1, 1, 0, 0, 0, 0, 2)])  # not tested
            sage: len(s.automorphisms())  # not tested
            2
            sage: len(f.automorphisms())  # not tested
            2

        A non-connected example::

            sage: t = Triangulation("(0,1,3)(2,4,~4)(~2,5,~5)(6,7,8)")  # not tested
            sage: len(t.automorphisms())  # not tested
            36

        TESTS::

            sage: examples = []  # not tested
            sage: examples.append(Triangulation("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"))  # not tested
            sage: examples.append(Triangulation("(0,8,~7)(1,9,~0)(2,10,~1)(3,11,~2)(4,12,~3)(5,13,~4)(6,14,~5)(7,15,~6)"))  # not tested
            sage: examples.append(Triangulation("(0,1,2)", boundary="(~0:1)(~1:1)(~2:1)"))  # not tested
            sage: examples.append(Triangulation("(0,1,3)(2,4,~4)(~2,5,~5)(6,7,8)"))  # not tested

            sage: examples.append(StrebelGraph("(0,3,7,~6,~2,1)(2,5,~4,~3,~1,~0)(4,8,~5)(6,~8,~7)"))  # not tested
            sage: for G in examples:  # not tested
            ....:     H = G.copy(mutable=True)
            ....:     for a in G.automorphisms():
            ....:         assert H == G
            ....:         H.relabel(a)
            ....:         assert H == G, (G, H, a)
        """
        best_relabellings = self.best_relabelling(return_all=True)[0]
        p0 = perm_invert(best_relabellings[0])
        return [perm_compose(p, p0) for p in best_relabellings]

    def automorphism_gens(self):
        return self.automorphisms()

    def best_relabelling(self, return_all=False):
        r"""
        Return a pair ``(r, data)`` where ``r`` is a relabelling that
        brings this constellation to the canonical one.

        EXAMPLES::

            sage: from veerer import Triangulation, VeeringTriangulation, StrebelGraph  # not tested
            sage: from veerer.permutation import perm_random_centralizer  # not tested

            sage: examples = []  # not tested
            sage: triangles = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"  # not tested
            sage: examples.append(Triangulation(triangles, mutable=True))  # not tested
            sage: examples.append(Triangulation("(0,1,3)(2,4,~4)(~2,5,~5)(6,7,8)", mutable=True))  # not tested
            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"  # not tested
            sage: cols = "BRBBBRRBBBBR"  # not tested
            sage: examples.append(VeeringTriangulation(fp, cols, mutable=True))  # not tested
            sage: fp = "(0,16,~15)(1,19,~18)(2,22,~21)(3,21,~20)(4,20,~19)(5,23,~22)(6,18,~17)(7,17,~16)(8,~1,~23)(9,~2,~8)(10,~3,~9)(11,~4,~10)(12,~5,~11)(13,~6,~12)(14,~7,~13)(15,~0,~14)"  # not tested
            sage: cols = "RRRRRRRRBBBBBBBBBBBBBBBB"  # not tested
            sage: examples.append(VeeringTriangulation(fp, cols, mutable=True))  # not tested
            sage: examples.append(StrebelGraph("(0,6,~5,~3,~1,4,~4,2,~2)(1)(3,~0)(5)(~6)", mutable=True))  # not tested
            sage: examples.append(StrebelGraph("(0,6,~5,~3,~1,4,~4:3,2,~2:3)(1:2)(3:2,~0)(5:2)(~6)", mutable=True))  # not tested

            sage: for G in examples:  # not tested
            ....:     print(G)
            ....:     r, fp, half_edges_data, edges_data = G.best_relabelling()
            ....:     for _ in range(10):
            ....:         p = perm_random_centralizer(G.edge_permutation())
            ....:         G.relabel(p)
            ....:         r2, fp2, half_edges_data2, edges_data2 = G.best_relabelling()
            ....:         assert fp2 == fp, G
            ....:         assert half_edges_data2 == half_edges_data, (G, half_edges_data2, half_edges_data)
            ....:         assert edges_data2 == edges_data, (G, edges_data2, edges_data)
            Triangulation("(0,~1,2)(~0,1,~3)(~2,~4,6)(3,4,~5)(5,~9,~7)(~6,8,7)(~8,~10,11)(9,10,~11)")
            Triangulation("(0,1,3)(2,4,~4)(~2,5,~5)(6,7,8)")
            VeeringTriangulation("(0,~1,2)(~0,1,~3)(~2,~4,6)(3,4,~5)(5,~9,~7)(~6,8,7)(~8,~10,11)(9,10,~11)", "BRBBBRRBBBBR")
            VeeringTriangulation("(0,16,~15)(~0,~14,15)(1,19,~18)(~1,~23,8)(2,22,~21)(~2,~8,9)(3,21,~20)(~3,~9,10)(4,20,~19)(~4,~10,11)(5,23,~22)(~5,~11,12)(6,18,~17)(~6,~12,13)(7,17,~16)(~7,~13,14)", "RRRRRRRRBBBBBBBBBBBBBBBB")
            StrebelGraph("(0,6,~5,~3,~1,4,~4,2,~2)(~0,3)(1)(5)(~6)")
            StrebelGraph("(0,6,~5,~3,~1,4,~4:3,2,~2:3)(~0,3:2)(1:2)(5:2)(~6)")
        """
        n = len(self._vp)
        ne = n // 2

        if not self.is_connected():
            # each component is labelled with consecutive half-edge labels
            # we use canonical labels for each of them, and then use a total ordering on the components
            components = {}
            for cc in self.connected_components():
                # TODO: set check to False
                comp = self.submap(cc, check=True)
                relabelling_best, fp_best, half_edges_data_best, edges_data_best = comp.best_relabelling(return_all=return_all)

                comp_hashable = [fp_best.tobytes()]
                comp_hashable.extend(data.tobytes() for data in half_edges_data_best)
                comp_hashable.extend(data.tobytes() for data in edges_data_best)
                comp_hashable = tuple(comp_hashable)
                if comp_hashable not in components:
                    components[comp_hashable] = []
                if return_all:
                    data = (cc, relabelling_best[0], relabelling_best, fp_best, half_edges_data_best, edges_data_best)
                else:
                    data = (cc, relabelling_best, None, fp_best, half_edges_data_best, edges_data_best)

                components[comp_hashable].append(data)

            relabelling_best = array('i', [-1] * n)
            fp_best = array('i', [-1] * n)
            half_edges_data_best = array('i', [0] * n)
            edges_data_best = array('i', [0] * n)

            shift = 0
            for comp_hashable in sorted(components):
                value = components[comp_hashable]
                for comp, comp_relabelling_best, _, _, _, _ in components[comp_hashable]:
                    # NOTE: elements in comp are edges, not half-edges
                    for i, j in enumerate(comp):
                        i0 = comp_relabelling_best[2 * i]
                        i1 = comp_relabelling_best[2 * i + 1]
                        relabelling_best[2 * j] = shift + i0
                        relabelling_best[2 * j + 1] = shift + i1
                    shift += 2 * len(comp)

            fp_best = perm_conjugate(self._fp, relabelling_best)
            half_edges_data_best = tuple(l[:] for l in self._half_edges_data)
            for ldest, lsrc in zip(half_edges_data_best, self._half_edges_data):
                perm_on_array(ldest, lsrc, relabelling_best, n)
            edges_data_best = tuple(l[:] for l in self._edges_data)
            for ldest, lsrc in zip(edges_data_best, self._edges_data):
                perm_on_edge_array(ldest, lsrc, relabelling_best, n)

            if not return_all:
                return (relabelling_best, fp_best, half_edges_data_best, edges_data_best)

            relabellings = []
            for oc in itertools.product(*[itertools.permutations(components[comp_hashable]) for comp_hashable in sorted(components)]):
                # run through all permutations of isomorphic components
                comps = [data[0] for isom_comps in oc for data in isom_comps]
                for comp_relabellings in itertools.product(*[data[2] for isom_comps in oc for data in isom_comps]):
                    # run through products available relabellings
                    relabelling = array('i', [-1] * n)
                    shift = 0
                    for comp, comp_relabelling in zip(comps, comp_relabellings):
                        # NOTE: elements in comp are edges, not half-edges
                        for i, j in enumerate(comp):
                            i0 = comp_relabelling[2 * i]
                            i1 = comp_relabelling[2 * i + 1]
                            relabelling[2 * j] = shift + i0
                            relabelling[2 * j + 1] = shift + i1
                        shift += 2 * len(comp)
                    relabellings.append(relabelling)

            return (relabellings, fp_best, half_edges_data_best, edges_data_best)

        else:
            # connected case
            fp = self._fp
            half_edges_data = self._half_edges_data
            edges_data = self._edges_data
            relabellings = []

            relabelling_new = array('i', [-1] * n)
            relabelling_best = array('i', [-1] * n)
            fp_new = array('i', [-1] * n)
            fp_best = array('i', [-1] * n)
            half_edges_data_new = tuple(l[:] for l in half_edges_data)
            half_edges_data_best = tuple(l[:] for l in half_edges_data)
            edges_data_new = tuple(l[:] for l in edges_data)
            edges_data_best = tuple(l[:] for l in edges_data)
            k_half_edges = len(half_edges_data)
            k_edges = len(edges_data)

            half_edges = self.half_edges()
            edge_relabelling_from(relabelling_best, fp_best, self._fp, 2 * ne, next(half_edges), 0)
            for i in range(k_half_edges):
                perm_on_array(half_edges_data_best[i], half_edges_data[i], relabelling_best, 2 * ne)
            for i in range(k_edges):
                perm_on_edge_array(edges_data_best[i], edges_data[i], relabelling_best, 2 * ne)

            if return_all:
                relabellings.append(relabelling_best[:])

            for start_half_edge in half_edges:
                # reinitialize relabelling_new as intended by edge_relabelling_from
                for i in range(n):
                    relabelling_new[i] = fp_new[i] = -1
                end_image = edge_relabelling_from(relabelling_new, fp_new, self._fp, n, start_half_edge, 0)
                assert end_image == 2 * ne, (end_image, ne)
                assert sum(x == -1 for x in fp_new) == sum(x == -1 for x in self._fp)
                assert all(x != -1 for x in relabelling_new)

                c = 0
                if fp_new < fp_best:
                    # no need to compare anything else
                    c = -1
                elif fp_new > fp_best:
                    # no need to go further
                    c = 1
                    continue

                for i in range(k_half_edges):
                    perm_on_array(half_edges_data_new[i], half_edges_data[i], relabelling_new, 2 * ne)
                    if not c:
                        if half_edges_data_new[i] < half_edges_data_best[i]:
                            c = -1
                        elif half_edges_data_new[i] > half_edges_data_best[i]:
                            c = 1
                            break
                if c == 1:
                    continue

                for i in range(k_edges):
                    perm_on_edge_array(edges_data_new[i], edges_data[i], relabelling_new, 2 * ne)
                    if not c:
                        if edges_data_new[i] < edges_data_best[i]:
                            c = -1
                        elif edges_data_new[i] > edges_data_best[i]:
                            c = 1
                            break
                if c == 1:
                    continue

                # at this stage either c=0 and relabelling is identical or c=-1 and we found something better
                if c == -1:
                    fp_best, fp_new = fp_new, fp_best
                    relabelling_best, relabelling_new = relabelling_new, relabelling_best
                    half_edges_data_best, half_edges_data_new = half_edges_data_new, half_edges_data_best
                    edges_data_best, edges_data_new = edges_data_new, edges_data_best
                    if return_all:
                        relabellings.clear()
                        relabellings.append(relabelling_best[:])
                elif return_all:
                    assert c == 0
                    relabellings.append(relabelling_new[:])

            return (relabellings, fp_best, half_edges_data_best, edges_data_best) if return_all else (relabelling_best, fp_best, half_edges_data_best, edges_data_best)

    # TODO: expand and clean documentation
    def set_canonical_labels(self, mapping=False, check=True):
        r"""
        Set labels in a canonical way in its automorphism class.

        EXAMPLES::

            sage: from veerer import *  # not tested
            sage: from veerer.permutation import perm_random, perm_random_centralizer  # not tested

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),  # not tested
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: T = Triangulation(t, mutable=True)  # not tested
            sage: T  # not tested
            Triangulation("(0,10,~9)(~0,11,~10)(1,~8,9)(~1,~7,8)(2,~6,7)(~2,~5,6)(3,~4,5)(~3,~11,4)")
            sage: T._check()  # not tested
            sage: T.set_canonical_labels()  # not tested
            sage: T  # not tested
            Triangulation("(0,1,2)(~0,~2,3)(~1,4,5)(~3,6,7)(~4,8,~5)(~6,9,~7)(~8,10,11)(~9,~11,~10)")
            sage: T._check()  # not tested
        """
        if check:
            self._assert_mutable()

        r, fp_best, half_edges_data_best, edges_data_best = self.best_relabelling()
        self._fp = fp_best
        self._vp = perm_conjugate(self._vp, r)
        self._half_edges_data = half_edges_data_best
        self._edges_data = edges_data_best
        self._set_data_pointers()
        if mapping:
            return r

    def _non_isom_easy(self, other):
        r"""
        A quick certificate of non-isomorphism that does not require relabellings.
        """
        return (len(self._vp) != len(other._vp) or
            perm_cycle_type(self._vp) != perm_cycle_type(other._vp) or
            self.num_folded_edges() != other.num_folded_edges() or
            perm_cycle_type(self._fp) != perm_cycle_type(other._fp) or
            any(sorted(l_self) != sorted(l_other) for l_self, l_other in zip(self._half_edges_data, other._half_edges_data)) or
            any(sorted(l_self) != sorted(l_other) for l_self, l_other in zip(self._edges_data, other._edges_data)))

    def is_isomorphic(self, other, certificate=False):
        r"""
        Return whether ``self`` is isomorphic to ``other``.

        INPUT:

        - ``other`` -- an :class:`OrientedMap`

        - ``certificate`` -- optional boolean (default ``False``), whether to
          additionally return the relabelling when ``self`` and ``other`` are
          isomorphic

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(vp="(0,1,2)(~0,~1,~2)")
            sage: m.is_isomorphic(m)
            Traceback (most recent call last):
            ...
            NotImplementedError

        .. TODO::

            The implementation relied on :meth:`best_relabelling` which reads
            the attribute ``_half_edges_data``. That attribute is never
            assigned on an :class:`OrientedMap`, so the whole canonical
            labelling machinery raises ``AttributeError``. Rewrite it, or drop
            the ``_half_edges_data`` / ``_edges_data`` indirection that was
            inherited from veerer.
        """
        raise NotImplementedError

    def dual(self, mutable=None, check=True):
        r"""
        Return the dual map of self.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: G = OrientedMap(vp=[[0,1,2,3]])
            sage: G
            OrientedMap("(0,~0,1,~1)", "(0)(~0,~1)(1)")
            sage: G.dual()
            OrientedMap("(0,1)(~0)(~1)", "(0,~0,1,~1)")

            sage: H = OrientedMap(vp=[[0,2,1,3]])
            sage: H.dual()
            OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)")
            sage: H == H.dual()
            True

            sage: I = OrientedMap(vp=[2,1,0,-1])
            sage: I
            OrientedMap("(0,1)(~0)", "(0,~0,1)")
            sage: I.dual()
            OrientedMap("(0,1,~0)", "(0,1)(~0)")

        Applying four times the dual function return the same map (with same labels)::

            sage: G0 = OrientedMap(vp=[[0,1,2,3]])
            sage: G1 = G0.dual()
            sage: G2 = G1.dual()
            sage: G3 = G2.dual()
            sage: G4 = G3.dual()
            sage: G4 == G0
            True
            sage: G2 == G0
            False

        """
        if mutable is None:
            mutable = self._mutable
        return OrientedMap(fp=self._vp, mutable=mutable)


    def smoothing(self, h, check=True):
        r"""
        Smooth the 4-degree vertex corresponding to h in the direction of h.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: G = OrientedMap(vp=[[0, 2, 4, 6], [1], [3], [5], [7]], mutable=True)
            sage: G.smoothing(0)
            sage: G
            OrientedMap("(0)(~0)(1)(~1)", "(0,~0)(1,~1)")
            sage: H = OrientedMap(vp=[[0, 2, 4, 6], [1, 8], [3, 9], [5, 10], [7, 11]], mutable=True)
            sage: H.smoothing(0)
            sage: H
            OrientedMap("(0,~5)(~0,4)(1,5)(~1,~4)", "(0,4,~1,5)(~0,~5,1,~4)")
            sage: I = OrientedMap(vp=[[0, 2, 3, 4], [1], [5]], mutable=True)
            sage: I.smoothing(0)
            sage: I
            OrientedMap("(0)(~0)", "(0,~0)")
            sage: J = OrientedMap(vp=[[0, 2, 3, 4], [1], [5]], mutable=True)
            sage: J.smoothing(2)
            sage: J
            OrientedMap("(1)(~1)", "(1,~1)")
        """
        if check:
            self._assert_mutable()
            self._check_half_edge(h)

        h1 = self._vp[h]
        h2 = self._vp[h1]
        h3 = self._vp[h2]
        if h==h1 or h==h2 or h==h3 or h!=self._vp[h3]:
            raise ValueError("the case of vertex of degree other than 4 is not implemented yet")
        self._check_half_edge_folded(h)
        self._check_half_edge_folded(h1)
        self._check_half_edge_folded(h2)
        self._check_half_edge_folded(h3)

        oh = self._ep(h)
        oh1 = self._ep(h1)
        oh2 = self._ep(h2)
        oh3 = self._ep(h3)

        self._vp[h] = self._vp[oh3] if self._vp[oh3] != oh3 else h
        self._vp[h1] = self._vp[oh2] if self._vp[oh2] != oh2 else h1
        if h3 != oh:
            self._vp[self._fp[h3]] = h
        if h2 != oh1:
            self._vp[self._fp[h2]] = h1
        self._vp[oh3] = -1
        self._vp[oh2] = -1
        self._vp[h2] = -1
        self._vp[h3] = -1

        self._fp[oh] = self._fp[h3]
        self._fp[oh1] = self._fp[h2]
        self._fp[self._ep(self._vp[h])] = h
        self._fp[self._ep(self._vp[h1])] = h1
        self._fp[h2] = -1
        self._fp[h3] = -1
        self._fp[oh2] = -1
        self._fp[oh3] = -1

        self._clear_trailing_edges()

    def disjoint_union(self, *others, check=True):
        r"""
        Add a copy of others in self. The labels of the half edges of others will be shifted by the sizes of the previous maps.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: M = OrientedMap(vp=[1, 0], mutable=True)
            sage: edge = OrientedMap(vp=[0, 1])
            sage: triangle = OrientedMap(vp=[3, 4, 5, 0, 1, 2])
            sage: M.disjoint_union(edge, triangle)
            sage: M
            OrientedMap("(0,~0)(1)(~1)(2,~3)(~2,4)(3,~4)", "(0)(~0)(1,~1)(2,4,3)(~2,~3,~4)")
        """

        if check:
            self._assert_mutable()
            self._check()
            for m in others:
                m._check()

        for m in others:
            n = len(self._vp)
            for h in range(len(m._vp)):
                self._vp.append(m._vp[h]+n)
                self._fp.append(m._fp[h]+n)


    def merge_vertices(self, *corners, check=True):
        r"""
        Merges the corners in corners.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: M = OrientedMap(vp=[0, 2, 1, 4, 3, 5], mutable=True)
            sage: M.merge_vertices(0, 2, 5)
            sage: M
            OrientedMap("(0,~0,1,~2)(~1,2)", "(0)(~0,~2,~1)(1,2)")
        """

        if check:
            self._assert_mutable()
            self._check()
            for c in corners:
                self._check_half_edge(c)

        for i in range(len(corners)-1):
            c0 = corners[i]
            c1 = corners[i+1]
            nv0 = self.next_at_vertex(c0)
            nv1 = self.next_at_vertex(c1)
            nf0 = self.previous_in_face(c0)
            nf1 = self.previous_in_face(c1)
            self._vp[c0] = nv1
            self._vp[c1] = nv0
            self._fp[nf1] = c0
            self._fp[nf0] = c1


    def move_half_edge(self, h, c, check=True):
        r"""
        Move the half_edge h to the corner after c. If c is negative then create a new vertex and attach h to it.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: G = OrientedMap(vp = [[0, 4, 2], [1], [3], [5]], mutable=True)
            sage: G.move_half_edge(4, 1)
            sage: G
            OrientedMap("(0,1)(~0,2)(~1)(~2)", "(0,2,~2,~0,1,~1)")
            sage: H = OrientedMap(vp = [[0, 4, 2], [1, 6], [3], [5]], mutable=True)
            sage: H.move_half_edge(6, 0)
            sage: H
            OrientedMap("(0,3,2,1)(~0)(~1)(~2)", "(0,~0,1,~1,2,~2,3)")
            sage: I = OrientedMap(vp = [[0, 4, 2], [1, 6], [3], [5]], mutable=True)
            sage: I.move_half_edge(4, 6)
            sage: I
            OrientedMap("(0,1)(~0,3,2)(~1)(~2)", "(0,2,~2,3,~0,1,~1)")

            sage: M = OrientedMap(vp = "(0,~1)(~0,2,5)(1,~2)(3,4)(~3,~4,~5)", mutable = True)
            sage: M.move_half_edge(10, -2)
            sage: M
            OrientedMap("(0,~1)(~0,2)(1,~2)(3,4)(~3,~4,~5)(5)", "(0,2,1)(~0,~1,~2)(3,~5,5,~4)(~3,4)")
        """

        if check:
            self._check()
            self._assert_mutable()
            h = self._check_half_edge(h)
            c = self._check_half_edge_or_negative(c)

        oh = self._ep(h)
        pre_h = self._fp[oh]

        self._vp[pre_h] = self._vp[h]
        if c >= 0:
            self._vp[h] = self._vp[c]
            self._vp[c] = h
            self._fp[oh] = c
        else:
            self._vp[h] = h
            self._fp[oh] = h
        self._fp[self._ep(self._vp[h])] = h
        self._fp[self._ep(self._vp[pre_h])] = pre_h


    def turn_around_vertex(self, h0, h1):
        r"""
        Compute the number of turn from dart h0 to dart h1 around a vertex.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(vp=[[0,1,2], [3]])
            sage: m.turn_around_vertex(0, 2)
            2
        """
        vp = self.vertex_permutation(copy=False)
        current = h0
        if h1 == h0:
            return 0
        else:
            current = vp[h0]
        turn = 1
        while current != h0 and h1 != current:
            turn += 1
            current = vp[current]
        if current == h0:
            raise ValueError("The half-edge {} is not on the same vertex as the half-edge {}.".format(f,e))
        return turn

# - relabel: keep combinatorics but change labellings
# - slide or half_edge_slide (possibly flip as a shortcut)
# - contract_edge, delete_edge
# - smoothing (for even degree vertices)
# - split_vertex (and adding an edge between the newly created one)
# - split_face
# - add_edge(h1, h2=None, h=None): if h2=None => folded and h1=h2 => loop (h is the new name)
# - glue(h1, h2)
# - union(m1, m2, m3, ...): disjoint union
