r"""
Vertex colored combinatorial maps

The main class in this module is :class:`ColoredOrientedMap`.
"""
# ****************************************************************************
#  This file is part of combisurf
#
#       Copyright (C) 2026 Juliette Schabanel
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

from combisurf import OrientedMap
from combisurf.permutation import (perm_cycles, perm_cycles_to_string, perm_orbit)
from collections import defaultdict

class ColoredOrientedMap(OrientedMap):
    r"""
    INPUT:
        - ``ecolors`` -- ``None`` or a list or dictionary of data so that ecolors[i] is the data associated to the edge (2i,2i+1).
        - ``vcolors`` -- a dictionary of data so that vcolors[h] is the data associated to the vertex incident to half edge h, and only one half edge incident to a given vertex is present. If the map is atomic (no edges), the color of the vertex can be specified in entry ``-1``.
    """
    def __init__(self, vp=None, fp=None, vcolors={}, ecolors=None, mutable=False, check=True):
        OrientedMap.__init__(self, vp, fp, mutable, check)

        if ecolors is None:
            self._edge_colors = [None]*(len(self._vp)//2)
        else:
            if isinstance(ecolors, list):
                if check:
                    if len(ecolors) != len(self._vp)//2:
                        raise ValueError("The number of colors doesn't match the number of edges")
                self._edge_colors = ecolors
            elif isinstance(ecolors, dict):
                self._edge_colors = [ecolors.get(e) for e in range(len(self._vp)//2)]
            else:
                raise TypeError("ecolors should be a list or dictionary")

        self._vertex_colors = {}
        if len(self._vp) == 0:
            self._vertex_colors = {-1: vcolors.get(-1)}
        else:
            for v in self.vertices():
                col = None
                for h in v:
                    hcol = vcolors.get(h)
                    if hcol is not None:
                        if col is None:
                            col = hcol
                        else:
                            raise ValueError(f"multiple values for {v} in vcolors")
                self._vertex_colors[v[0]] = col

    def __str__ (self):
        vp_cycles = perm_cycles(self._vp)
        vertices = perm_cycles_to_string(vp_cycles, edge_like=True)
        fp_cycles = perm_cycles(self._fp)
        faces = perm_cycles_to_string(fp_cycles, edge_like=True)
        return f"ColoredOrientedMap(\"{vertices}\", \"{faces}\", \n edge colors: {self._edge_colors}, \n vertex colors: {self._vertex_colors})"

    __repr__=__str__

    def copy(self, mutable=None):
        r"""
        Return a copy of this oriented map.

        If ``mutable`` is set to ``True`` or ``False`` return respectively an immutable or mutable version of the map.
        """

        if mutable is None:
            mutable = self._mutable

        if not self._mutable and not mutable:
            # avoid copies of immutable objects
            return self

        m = ColoredOrientedMap.__new__(ColoredOrientedMap)
        m._fp = self._fp[:]
        m._vp = self._vp[:]
        m._mutable = mutable
        m._edge_colors = self._edge_colors
        m._vertex_colors = self._vertex_colors
        return m

    def edge_color(self, e):
        r"""
        Return the color of the edge e.
        """
        return self._edge_colors[e]

    def _v_id(self, h):
        r"""
        Return the id of the vertex v incident to h (= half edge of minimal label incident to v).
        """
        h = self._check_half_edge(h)
        return min(perm_orbit(self._vp, h))

    def vertex_color(self, h):
        r"""
        Return the color of the vertex incident to h.
        """
        if h == -1 and len(self._vp)==0:
            return self._vertex_colors[-1]
        return self._vertex_colors[self._v_id(h)]


    def submap(self, edges, mutable=False, check=True):
        r"""
        Return the submap of this constellation induced on ``edges``.
        """

        unlab_submap = OrientedMap.submap(self, edges, mutable=mutable, check=check)
        relab_submap = OrientedMap.submap(self, edges, relabel=True, mutable=mutable, check=check)

        sub_edge_colors = [self.edge_color(e) for e in edges]
        sub_vert_colors = {(v[0] - 2*sum([e not in edges for e in range(v[0]//2)])): self.vertex_color(v[0]) for v in unlab_submap.vertices()}

        return ColoredOrientedMap(vp=relab_submap._vp, ecolors=sub_edge_colors, vcolors=sub_vert_colors, mutable=mutable, check=check)


    def plot(self, oriented=False, subdivide=True, root=None, edge_labels=True, vertex_colors=None, edge_colors=None):
        r"""
        Plot the map.

        INPUT:
            - ``oriented``: boolean specifying whether edge should be oriented.
            - ``subdivide``: boolean specifying whether multiple edges and loop should be subdivided for pretty plotting.
            - ``edge_labels``: boolean specifying whether to plot the labels of the edges.
            - ``edge_colors``: dictionary specifying the color to assign to each edge color.
            - ``vertex_colors``: dictionary specifying the color to assign to each vertex color.

        """

        if edge_colors is None:
            edge_cols = None
        else:
            edge_cols = defaultdict(list)
            for e, c in enumerate(self._edge_colors):
                col = edge_colors.get(c)
                if col is not None:
                    edge_cols[col].append(e)

        if vertex_colors is None:
            vert_cols = None
        else:
            vert_cols = defaultdict(list)
            for i, v in enumerate(self.vertices()):
                c = vertex_colors.get(self.vertex_color(v[0]))
                if c is not None:
                    vert_cols[c].append(i)

        return OrientedMap.plot(self, oriented=oriented, subdivide=subdivide, root=root, edge_labels=edge_labels, edge_colors=edge_cols, vertex_colors=vert_cols)


    #############
    # Mutations #
    #############

    def set_edge_color(self, h, col, check=True):
        r"""
        Change the color of the half edge h to col
        """
        if check:
            h = self._check_half_edge_or_negative(h)
            self._assert_mutable()
        self._edge_colors[h//2] = col

    def set_vertex_color(self, h, col, check=True):
        r"""
        Change the color of the vertex incident to half edge h to col
        """
        if check:
            h = self._check_half_edge_or_negative(h)
            self._assert_mutable()
        self._vertex_colors[self._v_id(h)] = col


    def add_edge(self, h0=-1, h1=-1, e=None, e_color=None, v0_color=None, v1_color=None, check=2):
        r"""
        Add an edge between the corners of ``h0`` and ``h1`` with color ``col``.
        If ``h0`` or ``h1`` is negative, then new vertices are created and their color can be specified in ``v0_color`` and ``v1_color``.
        If self is atomic and ``v0_color`` is not provided, the vertex keeps its color.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: M = ColoredOrientedMap(vp=[3, 4, 5, 0, 1, 2], vcolors={0:0, 1:2, 3:None}, ecolors = [5, None, 'a'], mutable=True)
            sage: M.add_edge(0, 2, e_color = 42)
            sage: M
            ColoredOrientedMap("(0,3,~1)(~0,2)(1,~3,~2)", "(0,2,~3)(~0,~1,~2)(1,3)",
             edge colors: [5, None, 'a', 42],
             vertex colors: {0: 0, 1: 2, 2: None})

            sage: E = ColoredOrientedMap(vp=[], vcolors={-1:'a'}, mutable=True)
            sage: E.add_edge(-1, -2, v1_color='b')
            sage: E
            ColoredOrientedMap("(0)(~0)", "(0,~0)", edge colors: [None], vertex colors: {0: 'a', 1: 'b'})
        """

        n = len(self._vp)
        if h0 < 0 :
            if n == 0 and v0_color is None:
                    self._vertex_colors[0] = self.vertex_color(-1)
            else:
                self._vertex_colors[n] = v0_color
        if h1 < 0 and h1 != h0:
            self._vertex_colors[n+1] = v1_color
        if n == 0:
            self._vertex_colors.pop(-1)

        OrientedMap.add_edge(self, h0=h0, h1=h1, e=e, check=check)
        self._edge_colors.append(e_color)


    def insert_edge(self, h0=-1, h1=-1, e=None, e_color=None, v_color=None, check=2):
        r"""
        Add an edge by splitting the vertex of ``h0`` and ``h1`` between them with color col. The two new vertices will have the same color as the original one.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: M = ColoredOrientedMap(vp="(0, 1, 2)(~0, 3, ~3)(~1,~2)", vcolors={0:0, 1:2, 3:None}, ecolors = [5, 7, 11, 13], mutable=True)
            sage: M.insert_edge(0, 2, e_color=42)
            sage: M
            ColoredOrientedMap("(0,1,2)(~0,3,~3,~4,~1,~2,4)", "(0,4,~3,~0,2,~1)(1,~4,~2)(3)",
            edge colors: [5, 7, 11, 13, 42],
            vertex colors: {0: 0, 1: 2, 3: None})

            sage: M.insert_edge(-1, -2, e_color='loop', v_color='n')
            sage: M
            ColoredOrientedMap("(0,1,2)(~0,3,~3,~4,~1,~2,4)(5,~5)", "(0,4,~3,~0,2,~1)(1,~4,~2)(3)(5)(~5)",
             edge colors: [5, 7, 11, 13, 42, 'loop'],
             vertex colors: {0: 0, 1: 2, 3: None, 10: 'n'})
        """

        OrientedMap.insert_edge(self, h0=h0, h1=h1, e=e, check=check)
        self._edge_colors.append(e_color)
        if h0 >= 0 and h1 >= 0:
            v0 = min(self._v_id(h0),self._v_id(h1))
            v1 = max(self._v_id(h0),self._v_id(h1))
            self._vertex_colors[v1] = self.vertex_color(v0)
        if h0 < 0 and h1 < 0:
            self._vertex_colors[len(self._vp)-2]=v_color


    def _clean_suppr_edge(self, e, h0, h1):
        OrientedMap.relabel(self)
        self._edge_colors.pop(e)
        vertex_colors = {}
        for h, c in self._vertex_colors.items():
            if h < 2*e:
                vertex_colors[h] = c
            elif h == 2*e:
                if h0 < 2*e:
                    vertex_colors[self._v_id(h0)] = c
                elif h0 > 2*e+1:
                    vertex_colors[self._v_id(h0-2)] = c
            elif h == 2*e+1:
                if h1 < 2*e:
                    vertex_colors[self._v_id(h1)] = c
                elif h1 > 2*e+1:
                    vertex_colors[self._v_id(h1-2)] = c
            else:
                vertex_colors[h-2] = c
        self._vertex_colors = vertex_colors


    def contract_edge(self, e, check=2):
        r"""
        Contract the edge ``e``.
        If both end point have the same color they keep it, if one is ``None``the new vertex gets the others color, otherwise it becomes ``None``.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: M = ColoredOrientedMap(vp="(0, 1, 2)(~0, 3, ~3)(~1,~2)", vcolors={0:0, 1:2, 3:None}, ecolors = [5, 7, 11, 13], mutable=True)
            sage: M.contract_edge(1)
            sage: M
            ColoredOrientedMap("(0,~1,1)(~0,2,~2)", "(0,~2,~0,1)(~1)(2)",  edge colors: [5, 11, 13],  vertex colors: {0: 0, 1: 2})
            sage: M.contract_edge(0)
            sage: M
            ColoredOrientedMap("(0,1,~1,~0)", "(0,~1)(~0)(1)", edge colors: [11, 13], vertex colors: {0: None})
        """

        v0 = min(self._v_id(2*e),self._v_id(2*e+1))
        v1 = max(self._v_id(2*e),self._v_id(2*e+1))
        c0 = self.vertex_color(v0)
        c1 = self.vertex_color(v1)
        if c1 is None or c1 == c0:
            self._vertex_colors.pop(v1)
        else:
            if c0 is None:
                self.set_vertex_color(v0, c1)
            else:
                self.set_vertex_color(v0, None)
            self._vertex_colors.pop(v1)
        h0 = self.next_at_vertex(2*e)
        h1 = self.next_at_vertex(2*e+1)
        OrientedMap.contract_edge(self, e=e, check=check)
        self._clean_suppr_edge(e, h0, h1)



    def delete_edge(self, e, check=2):
        r"""
        Delete the edge ``e``.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: M = ColoredOrientedMap(vp="(0, 1, ~2)(~0, 3, ~3)(~1,2)", vcolors={0:0, 1:2, 3:None}, ecolors = [5, 7, 11, 13], mutable=True)
            sage: M.delete_edge(2)
            sage: M  # not tested
            ColoredOrientedMap("(0,~1)(~0,2,~2)(1)", "(0,~2,~0,~1,1)(2)", edge colors: [5, 11, 13], vertex colors: {0: 0, 1: 2, 2: None})
        """

        h0 = self.next_at_vertex(2*e)
        h1 = self.next_at_vertex(2*e+1)
        OrientedMap.delete_edge(self, e=e, check=check)
        self._clean_suppr_edge(e, h0, h1)


    def reverse_orientation(self, e, check=True):
        r"""
        Change the orientation of the edge ``e``.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: M = ColoredOrientedMap(vp="(0, 1, 2)(~0, 3, ~3)(~1,~2)", vcolors={0:0, 1:2, 3:None}, ecolors = [5, 7, 11, 13], mutable=True)
            sage: M.reverse_orientation(1)  # not tested
             ColoredOrientedMap("(0,~1,2)(~0,3,~3)(1,~2)", "(0,~3,~0,2,1)(~1,~2)(3)",
             edge colors: [5, 7, 11, 13], vertex colors: {0: 0, 1: 2, 2: None})
        """

        if self._v_id(2*e) == 2*e and self._v_id(2*e+1) == 2*e+1:
            c0 = self.vertex_color(2*e)
            self.set_vertex_color(2*e, self.vertex_color(2*e+1))
            self.set_vertex_color(2*e+1, c0)
        elif self._v_id(2*e) == 2*e:
            self._vertex_colors[2*e+1] = self.vertex_color(2*e)
            self._vertex_colors.pop(2*e)
        elif self._v_id(2*e+1) == 2*e+1:
            self._vertex_colors[2*e] = self.vertex_color(2*e+1)
            self._vertex_colors.pop(2*e+1)
        OrientedMap.reverse_orientation(self, e, check=check)


    def move_half_edge(self, h, c, v_color=None, check=True):
        r"""
        Move the half_edge h to the corner after c.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: M = ColoredOrientedMap(vp="(0,1,2)(~0,3,4,~4)(~1,~3,5)(~2,~5,6)(~7,~6,7)", vcolors={0:0, 1:1, 3:0, 5:1, 13:1}, ecolors=[i%2 for i in range(8)], mutable=True)
            sage: M.move_half_edge(9, 7)
            sage: M
            ColoredOrientedMap("(0,1,2)(~0,3,4)(~1,~3,~4,5)(~2,~5,6)(~6,7,~7)", "(0,4,~3,~0,2,6,~7,~6,~5,~4,3,~1)(1,5,~2)(7)",
             edge colors: [0, 1, 0, 1, 0, 1, 0, 1], vertex colors: {0: 0, 1: 1, 3: 0, 5: 1, 13: 1})
            sage: M.move_half_edge(1, 10)
            sage: M
            ColoredOrientedMap("(0,1,2)(~0,~1,~3,~4,5)(~2,~5,6)(3,4)(~6,7,~7)", "(0,5,~2,1,~0,2,6,~7,~6,~5,~4,3,~1)(~3,4)(7)",
             edge colors: [0, 1, 0, 1, 0, 1, 0, 1], vertex colors: {0: 0, 1: 0, 5: 1, 13: 1, 6: 1})
            sage: M.move_half_edge(12, -1, v_color=42)
            sage: M
             ColoredOrientedMap("(0,1,2)(~0,~1,~3,~4,5)(~2,~5)(3,4)(6)(~6,7,~7)", "(0,5,~2,1,~0,2,~5,~4,3,~1)(~3,4)(6,~7,~6)(7)",
             edge colors: [0, 1, 0, 1, 0, 1, 0, 1], vertex colors: {0: 0, 1: 0, 5: 1, 13: 1, 6: 1, 12: 42})
        """

        if check:
            self._assert_mutable()
            h = self._check_half_edge(h)
            c = self._check_half_edge_or_negative(c)

        if h == c:
            return None

        nh = None
        if h == self._v_id(h):
            if self.next_at_vertex(h) == h:
                self._vertex_colors.pop(h)
            else :
                nh = self.next_at_vertex(h)
                colh = self.vertex_color(h)
        recycle = True
        if c < 0:
            self._vertex_colors[h] = v_color
            recycle = False
        elif h < self._v_id(c):
            self._vertex_colors[h] = self.vertex_color(c)
            if h != self._v_id(c):
                self._vertex_colors.pop(self._v_id(c))
            recycle = False

        OrientedMap.move_half_edge(self, h, c, check=check)

        if nh is not None:
            self._vertex_colors[self._v_id(nh)]=colh
            if recycle:
                self._vertex_colors.pop(h)




    def disjoint_union(self, *others, check=True):
        r"""
        Add a copy of others in self. The labels of the half edges of others will be shifted by the sizes of the previous maps.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: M1 = ColoredOrientedMap(vp="(0, 1, 2)(~0, 3, ~3)(~1,~2)", vcolors={0:0, 1:2, 3:None}, ecolors = [5, 7, 11, 13], mutable=True)
            sage: M2 = ColoredOrientedMap(vp="(0, 1)(~0, ~1)", vcolors={0:True, 1:False}, ecolors = ['a', 'b'])
            sage: M1.disjoint_union(M2)
            sage: M1
            ColoredOrientedMap("(0,1,2)(~0,3,~3)(~1,~2)(4,5)(~4,~5)", "(0,~3,~0,2,~1)(1,~2)(3)(4,~5)(~4,5)",
             edge colors: [5, 7, 11, 13, 'a', 'b'], vertex colors: {0: 0, 1: 2, 3: None, 8: True, 9: False})
        """

        shift = len(self._vp)
        OrientedMap.disjoint_union(self, *others, check=check)
        for m in others:
            self._edge_colors.extend(m._edge_colors)
            for v, c in m._vertex_colors.items():
                self._vertex_colors[v+shift]=c
            shift += len(m._vp)



    def merge_vertices(self, *corners, color=None, check=True):
        r"""
        Merges the corners in corners. If color is given then the new vertex takes this color. Otherwise if exactly one color is present among the corners then it keeps it.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: M = ColoredOrientedMap(vp="(0,1,2)(~0,3,~3)(~1,~2)(4,5)(~4,~5)", ecolors=[5, 7, 11, 13, 'a', 'b'], vcolors={0: 0, 1: 2, 3: True, 8: True, 9: False}, mutable = True)
            sage: M.merge_vertices(3, 8)
            sage: M
            ColoredOrientedMap("(0,1,2)(~0,3,~3)(~1,5,4,~2)(~4,~5)", "(0,~3,~0,2,4,~5,~1)(1,~2)(3)(~4,5)",
             edge colors: [5, 7, 11, 13, 'a', 'b'], vertex colors: {0: 0, 1: 2, 9: False, 3: True})
            sage: M.merge_vertices(0, 1, color=42)
            sage: M
            ColoredOrientedMap("(0,3,~3,~0,1,2)(~1,5,4,~2)(~4,~5)", "(0,~3)(~0,2,4,~5,~1)(1,~2)(3)(~4,5)",
             edge colors: [5, 7, 11, 13, 'a', 'b'], vertex colors: {9: False, 3: True, 0: 42})
        """

        if check:
            self._assert_mutable()
            self._check()
            for c in corners:
                self._check_half_edge(c)

        new_color = color
        if color is None :
            for c in corners:
                if self.vertex_color(c) is not None:
                    if new_color is None or new_color == self.vertex_color(c) :
                        new_color = self.vertex_color(c)
                    else:
                        new_color = None
                        break
        new_id = min([self._v_id(c) for c in corners])
        for c in corners:
            if self._v_id(c) in self._vertex_colors:
                self._vertex_colors.pop(self._v_id(c))
        self._vertex_colors[new_id] = new_color
        OrientedMap.merge_vertices(self, *corners, check=False)

    def swap_color(self, c0=0, c1=1, check=True):
        r"""
        Swap colors c0 and c1 in m.
        """

        if check:
            self._assert_mutable()
            self._check()

        for v, c in self._vertex_colors.items():
            if c == c0:
                self._vertex_colors[v]=c1
            elif c == c1:
                self._vertex_colors[v]=c0



    def relabel(self, p=None, check=2):
        raise NotImplementedError

    def dual(self, mutable=None, check=True):
        raise NotImplementedError

    def smoothing(self, h, check=True):
        raise NotImplementedError

    def canonical_labelling(self, half_edge_data=None, edge_data=None, vertex_data=None, face_data=None, check=True):
        r"""
        Return the canonical labelling of this coloured map.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: m = ColoredOrientedMap(fp="(0,1,2)(~0,~2,~1)", ecolors=[0, 1, 1])
            sage: m.canonical_labelling()
            Traceback (most recent call last):
            ...
            NotImplementedError: colours are not supported yet
        """
        raise NotImplementedError("colours are not supported yet")

    def canonical_form(self, half_edge_data=None, edge_data=None, vertex_data=None, face_data=None, check=True):
        r"""
        Return the canonical form of this coloured map.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: m = ColoredOrientedMap(fp="(0,1,2)(~0,~2,~1)", ecolors=[0, 1, 1])
            sage: m.canonical_form()
            Traceback (most recent call last):
            ...
            NotImplementedError: colours are not supported yet
        """
        raise NotImplementedError("colours are not supported yet")

    def automorphisms(self, half_edge_data=None, edge_data=None, vertex_data=None, face_data=None, check=True):
        r"""
        Return the automorphisms of this coloured map.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: m = ColoredOrientedMap(fp="(0,1,2)(~0,~2,~1)", ecolors=[0, 1, 1])
            sage: m.automorphisms()
            Traceback (most recent call last):
            ...
            NotImplementedError: colours are not supported yet
        """
        raise NotImplementedError("colours are not supported yet")

    def is_isomorphic(self, other, certificate=False, half_edge_data=None, edge_data=None, vertex_data=None, face_data=None, check=True):
        r"""
        Return whether this coloured map is isomorphic to another.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: m = ColoredOrientedMap(fp="(0,1,2)(~0,~2,~1)", ecolors=[0, 1, 1])
            sage: m.is_isomorphic(m)
            Traceback (most recent call last):
            ...
            NotImplementedError: colours are not supported yet
        """
        raise NotImplementedError("colours are not supported yet")

    def num_automorphisms(self, half_edge_data=None, edge_data=None, vertex_data=None, face_data=None, check=True):
        r"""
        Return the order of the automorphism group of this coloured map.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: m = ColoredOrientedMap(fp="(0,1,2)(~0,~2,~1)", ecolors=[0, 1, 1])
            sage: m.num_automorphisms()
            Traceback (most recent call last):
            ...
            NotImplementedError: colours are not supported yet
        """
        raise NotImplementedError("colours are not supported yet")

    def automorphism_vertex_action(self, a, check=True):
        r"""
        Return the permutation of the vertices induced by an automorphism.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: m = ColoredOrientedMap(fp="(0,1,2)(~0,~2,~1)", ecolors=[0, 1, 1])
            sage: m.automorphism_vertex_action([0, 1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            NotImplementedError: colours are not supported yet
        """
        raise NotImplementedError("colours are not supported yet")

    def automorphism_face_action(self, a, check=True):
        r"""
        Return the permutation of the faces induced by an automorphism.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: m = ColoredOrientedMap(fp="(0,1,2)(~0,~2,~1)", ecolors=[0, 1, 1])
            sage: m.automorphism_face_action([0, 1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            NotImplementedError: colours are not supported yet
        """
        raise NotImplementedError("colours are not supported yet")

    def canonical_roots(self, half_edge_data=None, edge_data=None, vertex_data=None, face_data=None, check=True):
        r"""
        Return the half-edges that some minimiser of the canonical key
        labels ``0``.

        EXAMPLES::

            sage: from combisurf import ColoredOrientedMap
            sage: m = ColoredOrientedMap(fp="(0,1,2)(~0,~2,~1)", ecolors=[0, 1, 1])
            sage: m.canonical_roots()
            Traceback (most recent call last):
            ...
            NotImplementedError: colours are not supported yet
        """
        raise NotImplementedError("colours are not supported yet")
