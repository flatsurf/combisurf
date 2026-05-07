r"""
Combinatorial maps on oriented surfaces

The main class in this module is :class:`OrientedMap` which describes a
cell decomposition of an oriented surface.
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
from combisurf.permutation import (perm_init, perm_check, perm_cycles, perm_on_array, perm_on_edge_array,
                          perm_invert, perm_conjugate, perm_conjugate_transposition_inplace, perm_cycle_string, perm_cycles_lengths,
                          perm_cycles_to_string, perm_on_list, perm_on_edge_list, perm_cycle_type,
                          perm_num_cycles, str_to_cycles, str_to_cycles_and_data, perm_compose, perm_from_base64_str,
                          uint_base64_str, uint_from_base64_str, perm_base64_str,
                          perm_orbit, perm_orbit_size,
                          perms_are_transitive, perms_orbits, perm_edge_orbits, edge_relabelling_from)


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



def remove_trailing_minus_ones(p):
    while p and p[-2] == -1:
        if p[-1] != -1:
            raise ValueError("invalid permutation")
        p.pop()
        p.pop()

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
    implicit and the vertex and face permutations are always abreviated as ``vp`` and
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

    In cycle notation, if an half-edge is not mentionned then the corresponding edge is folded::

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

        if len(vp) != len(fp):
            raise ValueError(f"inconsistent input: vp has length {len(vp)} while fp has length {len(fp)}")

        remove_trailing_minus_ones(vp)
        remove_trailing_minus_ones(fp)

        self._vp = vp
        self._fp = fp
        self._mutable = mutable

        if check:
            self._check(ValueError)


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
        ne = len(self._vp) // 2

        if not (hasattr(self, '_vp') and hasattr(self, '_fp')):
            raise error("missing attributes: these must be _vp, _ep, _fp, _data")
        if not perm_check(self._vp, 2 * ne):
            raise error(f"vp is not a permutation: {self._vp}")
        if not perm_check(self._fp, 2 * ne):
            raise error(f"fp is not a permutation: {self._fp}")

        if self._vp and (self._vp[-2] == -1 or self._fp[-2] == -1):
            raise error("trailing -1 in vertex or face permutation")

        for h in range(2 * ne):
            if (self._vp[h] == -1) != (self._fp[h] == -1):
                raise ValueError(f"vp (={self._vp}) and fp (={self._fp}) with different domains")

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
            ValueError: invalid half-edge (=1); the underlying edge is folded
        """
        if not isinstance(h, numbers.Integral):
            raise TypeError(f"invalid half-edge {h} of type {type(h).__name__}")
        h = int(h)
        if h < 0 or h >= len(self._vp):
            raise ValueError(f"half-edge number out of range (={h})")
        if self._vp[h] == -1:
            raise ValueError(f"invalid half-edge (={h}); the underlying edge is folded")
        return h

    def _check_half_edge_or_negative(self, h):
        if not isinstance(h, numbers.Integral):
            raise TypeError(f"invalid half-edge {h} of type {type(h).__name__}")
        h = int(h)
        if h >= 0:
            if  h >= len(self._vp):
                raise ValueError(f"half-edge number out of range (={h})")
            if self._vp[h] == -1:
                raise ValueError(f"invalid half-edge (={h}); the underlying edge is folded")
        return h

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

        - ``subdivide`` (boolean, default ``False``) -- if ``True`` then insert
          one vertex on each edge and two vertices on each loop so that the
          resulting graph is simple.

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
        """
        # NOTE: sage graphs use *clockwise* order for the neighbors
        if self.has_folded_edge():
            raise NotImplementedError

        if root is None:
            root = len(self._vp) - 2

        if subdivide is True:
            subdivide = self._subdivide_true
        elif subdivide is False:
            subdivide = self._subdivide_false
        elif isinstance(subdivide, numbers.Integral):
            k = int(subdivide)
            if k < 0:
                raise ValueError("subdivide can not be a negative integer")
            subdivide = lambda u, v, half_edge_list, h: k
        elif not callable(subdivide):
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
            G = DiGraph(len(vertices), multiedges=not subdivide, loops=not subdivide)
        else:
            from sage.graphs.graph import Graph
            G = Graph(len(vertices), multiedges=not subdivide, loops=not subdivide)

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
        - ``edge_colors``: dictionnary specifying the color to assign to each edge color.
        - ``vertex_colors``: dictionnary specifying the color to assign to each vertex color.
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
            sage: s.flip(0)
            sage: s == m
            False

            sage: m.set_immutable()
            sage: m.copy() is m
            True

            sage: t = m.copy(mutable=True)
            sage: t.flip(0)
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
            ValueError: invalid half-edge (=1); the underlying edge is folded
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
            sage: OrientedMap(fp="(2,~3)(~2,4)").is_connected()
            True
            sage: OrientedMap(fp="(2,~3)(4,~4)").is_connected()
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
        Return the submap of this constellation induced on ``edges``.

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

        assert perm_check(vp)
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

        while vp and vp[-2] == -1:
            vp.pop()
            vp.pop()
            fp.pop()
            fp.pop()

    def delete_edge(self, e, check=2):
        r"""
        Delete the edge ``e``.

        EXAMPLES::

            sage: from combisurf import OrientedMap

            sage: vp = "(0,1,~0,~1)"
            sage: fp = "(0,1,~0,~1)"
            sage: m = OrientedMap(vp, fp)

            sage: vp20 = "(0,~2,2,1,~0,~1)"
            sage: fp20 = "(0,1,~0,~1,2)(~2)"
            sage: m = OrientedMap(vp20, fp20, mutable=True)
            sage: m.delete_edge(2)
            sage: m

            sage: vp10 = "(0,2,1,~2,~0,~1)"
            sage: fp10 = "(0,~2)(~0,~1,2,1)"
            sage: m = OrientedMap(vp10, fp10, mutable=True)
            sage: m.delete_edge(2)
            sage: m

            sage: vp30 = "(0,2,1,~0,~2,~1)"
            sage: fp30 = "(0,1,~2)(~0,~1,2)"
            sage: m = OrientedMap(vp30, fp30, mutable=True)
            sage: m.delete_edge(2)
            sage: m

            sage: vp00 = "(0,~2,2,1,~0,~1)"
            sage: fp00 = "(0,1,~0,~1,2)(~2)"
            sage: m = OrientedMap(vp00, fp00, mutable=True)
            sage: m.delete_edge(2)
            sage: m

            sage: vp22 = "(0,1,~2,2,~0,~1)"
            sage: fp22 = "(0,2,1,~0,~1)(~2)"
            sage: m = OrientedMap(vp00, fp00, mutable=True)
            sage: m.delete_edge(2)
            sage: m
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

        while vp and vp[-2] == -1:
            vp.pop()
            vp.pop()
            fp.pop()
            fp.pop()

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
            OrientedMap("(0,1,2,~0,~1,~2)", "(0,2,~1)(~0,~2,1)")

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
            sage: m.relabel("(0,1,~2)(~0,~1,2)")
            sage: m

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
        remove_trailing_minus_ones(self._vp)
        remove_trailing_minus_ones(self._fp)

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
            ValueError: immutable graph; make a mutable copy first
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

            sage: from veerer import *

        Veering triangulation example::

            sage: vt = VeeringTriangulation("(0,1,2)(3,4,~0)(5,6,~1)(7,~2,8)(9,~3,~6)(10,~7,~4)(11,~5,12)(13,14,~8)(15,~9,16)(17,18,~10)(19,~17,~11)(20,~13,~12)(21,~14,~18)(22,~21,~15)(23,24,~16)(25,~23,~19)(26,~20,~25)(~26,~24,~22)", "RBBRBRBRRBRBBRBBRRBRRRBBRRB")
            sage: len(vt.automorphisms())
            2
            sage: qvt = vt.automorphism_quotient()
            sage: qvt
            VeeringTriangulation("(0,1,2)(~0,3,4)(~1,5,6)(~2,8,7)(~3,~6,9)(~4,10,~7)(~5,12,11)(~8,13,~12)(~10,14,~11)", "RBBRBRBRRBRBBRR")
            sage: (vt.stratum(), qvt.stratum())  # optional - surface_dynamics
            (H_4(2^3), Q_1(1^3, -1^3))

            sage: vt.automorphism_quotient(mapping=True)
            (VeeringTriangulation("(0,1,2)(~0,3,4)(~1,5,6)(~2,8,7)(~3,~6,9)(~4,10,~7)(~5,12,11)(~8,13,~12)(~10,14,~11)", "RBBRBRBRRBRBBRR"),
             array('i', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 18, 20, 21, 22, 23, 24, 25, 26, 26, 25, 24, 13, 12, 7, 6, 28, 28, 23, 22, 21, 20, 17, 16, 11, 10, 3, 2, 8, 9, 1, 0, 15, 14, 5, 4]))

        Strebel graph example::

            sage: sg = StrebelGraph("(0,~0,~1)(1,2,~2)")
            sage: sg.automorphism_quotient()
            StrebelGraph("(0,~0,1)")

        TESTS::

            sage: vt = VeeringTriangulation("(~0,1,2)(~1,3,4)(~2,5,6)(~3,7,8)(~4,~7,9)(~6,10,11)(~8,12,13)(~9,14,15)(~10,16,17)(~11,18,~17)(~12,19,~18)(~14,20,~16)(0:1)(~5:1)(~13:1)(~15:1)(~19:1)(~20:1)", "RBBBRRBRBBRRBRBRBBBRR")
            sage: vt.automorphism_quotient()
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

            sage: from veerer import *
            sage: from array import array

        The torus example (6 symmetries)::

            sage: fp = array('i', [2, 3, 5, 4, 1, 0])
            sage: vp = array('i', [4, 5, 1, 0, 2, 3])
            sage: T = Triangulation.from_permutations(vp, fp, (array('i', [0]*6),), mutable=True)
            sage: T._relabelling_from(3)
            array('i', [5, 4, 1, 0, 2, 3])

            sage: p = T._relabelling_from(0)
            sage: T.relabel(p)
            sage: for i in range(6):
            ....:     p = T._relabelling_from(i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     assert S == T

        The sphere example (3 symmetries)::

            sage: fp = array('i', [2, -1, 4, -1, 0, -1])
            sage: vp = array('i', [4, -1, 0, -1, 2, -1])
            sage: T = Triangulation.from_permutations(vp, fp, (array('i', [0]*6),), mutable=True)
            sage: T._relabelling_from(2)
            array('i', [4, 5, 0, 1, 2, 3])
            sage: p = T._relabelling_from(0)
            sage: T.relabel(p)
            sage: for i in range(3):
            ....:     p = T._relabelling_from(2 * i)
            ....:     S = T.copy()
            ....:     S.relabel(p)
            ....:     assert S == T

        An example with no automorphism::

            sage: T = Triangulation("(0,1,2)(3,4,5)(~0,~3,6)", mutable=True)
            sage: p = T._relabelling_from(0)
            sage: T.relabel(p)
            sage: for i in T.half_edges():
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

            sage: from veerer import *

        An example with 4 symmetries in genus 2::

            sage: T = Triangulation("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)")
            sage: A = T.automorphisms()
            sage: len(A)
            4

        And the "sphere octagon" has 8::

            sage: s  = "(0,8,~7)(1,9,~0)(2,10,~1)(3,11,~2)(4,12,~3)(5,13,~4)(6,14,~5)(7,15,~6)"
            sage: len(Triangulation(s).automorphisms())
            8

        A veering triangulation with 4 symmetries in genus 2::

            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: V = VeeringTriangulation(fp, cols)
            sage: A = V.automorphisms()
            sage: len(A)
            4

        Examples with boundaries::

            sage: t = Triangulation("(0,1,2)", boundary="(~0:1)(~1:1)(~2:1)")
            sage: len(t.automorphisms())
            3
            sage: t = Triangulation("(0,1,2)", boundary="(~0:1,~1:1,~2:1)")
            sage: len(t.automorphisms())
            3
            sage: t = Triangulation("(0,1,2)", boundary="(~0:1,~1:1,~2:2)")
            sage: len(t.automorphisms())
            1

        Linear families::

            sage: s = StrebelGraph("(0,3,7,~6,~2,1)(2,5,~4,~3,~1,~0)(4,8,~5)(6,~8,~7)")
            sage: f = StrebelGraphLinearFamily(s, [(2, 0, 0, 0, 1, 0, 1, 0, 2), (0, 2, 0, 0, 0, 1, 0, 1, 2), (0, 0, 1, 1, 0, 0, 0, 0, 2)])
            sage: len(s.automorphisms())
            2
            sage: len(f.automorphisms())
            2

        A non-connected example::

            sage: t = Triangulation("(0,1,3)(2,4,~4)(~2,5,~5)(6,7,8)")
            sage: len(t.automorphisms())
            36

        TESTS::

            sage: examples = []
            sage: examples.append(Triangulation("(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"))
            sage: examples.append(Triangulation("(0,8,~7)(1,9,~0)(2,10,~1)(3,11,~2)(4,12,~3)(5,13,~4)(6,14,~5)(7,15,~6)"))
            sage: examples.append(Triangulation("(0,1,2)", boundary="(~0:1)(~1:1)(~2:1)"))
            sage: examples.append(Triangulation("(0,1,3)(2,4,~4)(~2,5,~5)(6,7,8)"))

            sage: examples.append(StrebelGraph("(0,3,7,~6,~2,1)(2,5,~4,~3,~1,~0)(4,8,~5)(6,~8,~7)"))
            sage: for G in examples:
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

            sage: from veerer import Triangulation, VeeringTriangulation, StrebelGraph
            sage: from veerer.permutation import perm_random_centralizer

            sage: examples = []
            sage: triangles = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: examples.append(Triangulation(triangles, mutable=True))
            sage: examples.append(Triangulation("(0,1,3)(2,4,~4)(~2,5,~5)(6,7,8)", mutable=True))
            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: examples.append(VeeringTriangulation(fp, cols, mutable=True))
            sage: fp = "(0,16,~15)(1,19,~18)(2,22,~21)(3,21,~20)(4,20,~19)(5,23,~22)(6,18,~17)(7,17,~16)(8,~1,~23)(9,~2,~8)(10,~3,~9)(11,~4,~10)(12,~5,~11)(13,~6,~12)(14,~7,~13)(15,~0,~14)"
            sage: cols = "RRRRRRRRBBBBBBBBBBBBBBBB"
            sage: examples.append(VeeringTriangulation(fp, cols, mutable=True))
            sage: examples.append(StrebelGraph("(0,6,~5,~3,~1,4,~4,2,~2)(1)(3,~0)(5)(~6)", mutable=True))
            sage: examples.append(StrebelGraph("(0,6,~5,~3,~1,4,~4:3,2,~2:3)(1:2)(3:2,~0)(5:2)(~6)", mutable=True))

            sage: for G in examples:
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

            sage: from veerer import *
            sage: from veerer.permutation import perm_random, perm_random_centralizer

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: T = Triangulation(t, mutable=True)
            sage: T
            Triangulation("(0,10,~9)(~0,11,~10)(1,~8,9)(~1,~7,8)(2,~6,7)(~2,~5,6)(3,~4,5)(~3,~11,4)")
            sage: T._check()
            sage: T.set_canonical_labels()
            sage: T
            Triangulation("(0,1,2)(~0,~2,3)(~1,4,5)(~3,6,7)(~4,8,~5)(~6,9,~7)(~8,10,11)(~9,~11,~10)")
            sage: T._check()
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

    def iso_sig(self):
        r"""
        Return a canonical signature.

        EXAMPLES::

            sage: from veerer import *
            sage: T = Triangulation("(0,3,1)(~0,4,2)(~1,~2,~4)")
            sage: T.iso_sig()
            '5_1__1_2~46098537_0000000000'
            sage: TT = Triangulation.from_string(T.iso_sig())
            sage: TT
            Triangulation("(0,1,2)(~1,3,4)(~2,~4,~3)")
            sage: TT.iso_sig() == T.iso_sig()
            True

            sage: T = Triangulation("(0,10,~6)(1,12,~2)(2,14,~3)(3,16,~4)(4,~13,~5)(5,~1,~0)(6,~17,~7)(7,~14,~8)(8,13,~9)(9,~11,~10)(11,~15,~12)(15,17,~16)")
            sage: T.iso_sig()
            'i_1__1_264a0e8i1mcj3sgr5tkq7xov9dbupwfzhyln_000000000000000000000000000000000000'
            sage: Triangulation.from_string(T.iso_sig())
            Triangulation("(0,1,2)(~0,3,4)(~1,5,6)(~2,7,8)(~3,9,10)(~4,11,12)(~5,~9,13)(~6,14,~12)(~7,~13,15)(~8,~14,16)(~10,~16,17)(~11,~15,~17)")

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: cols = [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE]
            sage: T = VeeringTriangulation(t, cols, mutable=True)
            sage: T.iso_sig()
            'c_1_1_1_2548061cag39ei7dbkfnmjhl_000000000000000000000000_122212122212'

        If we relabel the triangulation, the isomorphic signature does not change::

            sage: from veerer.permutation import perm_random_centralizer
            sage: p = perm_random_centralizer(T.edge_permutation())
            sage: T.relabel(p)
            sage: T.iso_sig()
            'c_1_1_1_2548061cag39ei7dbkfnmjhl_000000000000000000000000_122212122212'

        An isomorphic triangulation can be reconstructed from the isomorphic
        signature via::

            sage: s = T.iso_sig()
            sage: T2 = VeeringTriangulation.from_string(s)
            sage: T == T2
            False
            sage: T.is_isomorphic(T2)
            True

        TESTS::

            sage: from veerer.veering_triangulation import VeeringTriangulation
            sage: from veerer.permutation import perm_random

            sage: t = [(-12, 4, -4), (-11, -1, 11), (-10, 0, 10), (-9, 9, 1),
            ....:      (-8, 8, -2), (-7, 7, 2), (-6, 6, -3), (-5, 5, 3)]
            sage: cols = [RED, RED, RED, RED, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE, BLUE]
            sage: T = VeeringTriangulation(t, cols, mutable=True)
            sage: iso_sig = T.iso_sig()
            sage: for _ in range(10):
            ....:     p = perm_random_centralizer(T.edge_permutation())
            ....:     T.relabel(p)
            ....:     assert T.iso_sig() == iso_sig

            sage: VeeringTriangulation("(0,1,2)(3,4,~1)(5,6,~4)", "RBGGRBG").iso_sig()
            '7_1_1_1_2~4~068~5ac~9~_00000000000000_8128128'
            sage: VeeringTriangulation.from_string('7_1_1_1_2~4~068~5ac~9~_00000000000000_8128128')
            VeeringTriangulation("(0,1,2)(~2,3,4)(~4,5,6)", "GRBGRBG")
        """
        T = self.copy(mutable=True)
        T.set_canonical_labels()
        return T.to_string()

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

        - ``other`` - a constellation

        - ``certificate`` -- optional boolean (default ``False``), whether to
           additionally return the relabelling when ``self`` and ``other`` are
           isomorphic

        EXAMPLES::

            sage: from veerer import Triangulation
            sage: sphere = Triangulation("(0,1,2)(~0,~2,~1)")
            sage: sphere2 = Triangulation("(0,2,1)(~0,~1,~2)")
            sage: torus = Triangulation("(0,1,2)(~0,~1,~2)")
            sage: sphere.is_isomorphic(sphere2)
            True
            sage: sphere.is_isomorphic(torus)
            False

        TESTS::

            sage: from veerer import Triangulation, VeeringTriangulation
            sage: from veerer.permutation import perm_random_centralizer

            sage: T = Triangulation("(0,5,1)(~0,4,2)(~1,~2,~4)(3,6,~5)", mutable=True)
            sage: TT = T.copy()
            sage: for _ in range(10):
            ....:     rel = perm_random_centralizer(TT.edge_permutation())
            ....:     TT.relabel(rel)
            ....:     assert T.is_isomorphic(TT)

            sage: fp = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
            sage: cols = "BRBBBRRBBBBR"
            sage: V = VeeringTriangulation(fp, cols, mutable=True)
            sage: W = V.copy()
            sage: p = perm_random_centralizer(V.edge_permutation())
            sage: W.relabel(p)
            sage: assert V.is_isomorphic(W) is True
            sage: ans, cert = V.is_isomorphic(W, True)
            sage: V.relabel(cert)
            sage: assert V == W
        """
        if type(self) is not type(other):
            raise TypeError("can only check isomorphisms between identical types")

        if self._non_isom_easy(other):
            return (False, None) if certificate else False

        r1, fp1, half_edges_data1, edges_data1 = self.best_relabelling()
        r2, fp2, half_edges_data2, edges_data2 = other.best_relabelling()

        if fp1 != fp2 or half_edges_data1 != half_edges_data2 or edges_data1 != edges_data2:
            return (False, None) if certificate else False
        elif certificate:
            return (True, perm_compose(r1, perm_invert(r2)))
        else:
            return True

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

        remove_trailing_minus_ones(self._fp)
        remove_trailing_minus_ones(self._vp)

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

    def hom(self, images, codomain=None):
        r"""
        Return a morphism defined by the given ``images``.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: m = OrientedMap(fp="(0,1,~0,~1)")
            sage: m.hom([[0,1,2],[2],[0,0],[3,2,1]])
            OrientedMapMorphism(OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)"), OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)"), [array('i', [0, 1, 2]), array('i', [2]), array('i', [0, 0]), array('i', [3, 2, 1])])
        """
        from combisurf.morphism import OrientedMapMorphism_list
        if codomain is None:
            codomain = self
        return OrientedMapMorphism_list(self, codomain, images)

# - relabel: keep combinatorics but change labellings
# - slide or half_edge_slide (possibly flip as a shortcut)
# - contract_edge, delete_edge
# - smoothing (for even degree vertices)
# - split_vertex (and adding an edge between the newly created one)
# - split_face
# - add_edge(h1, h2=None, h=None): if h2=None => folded and h1=h2 => loop (h is the new name)
# - glue(h1, h2)
# - union(m1, m2, m3, ...): disjoint union

