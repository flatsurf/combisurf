r"""
Walks on oriented maps

This module introduce one class: :class:`Walk`.
"""
# ****************************************************************************
#  This file is part of combisurf
#
#       Copyright (C) 2026 Oscar Fontaine
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

# TODO : Mettre les notations cohérentes avec oriented map : possibilité de string et 2i -> i // 2i+1 -> ~i


from combisurf import OrientedMap, QuadSystem, Geodesic, StarShapedSpace


def is_subword(u, v):
    r"""
    Test if u is a subword of v in O(|u|+|v|).
    """
    if len(u) == 0:
        return True
    elif len(v) == 0:
        return False
    cnd = 0
    T = [-1]
    for i in range(1, len(u)):
        if u[i] == u[cnd]:
            T.append(u[cnd])
        else:
            T.append(cnd)
            while cnd>=0 and u[i] == u[cnd]:
                cnd = T[cnd]
            cnd += 1
    j = 0
    k = 0
    res = False
    while j < len(v) and not res:
        if u[k] == v[j]:
            j += 1
            k += 1
            if k == len(u):
                res = True
        else:
            k = T[k]
            if k == -1:
                k += 1
                j += 1
    return res



class Walk:

    def __init__(self, oriented_map, walk, check=True):
        r"""
        Methods:
            _oriented_map: the underlying map
            _walk: the walk (as a list of half-edges) in oriented_map
        """

        self._oriented_map = oriented_map
        self._walk = walk

        if check:
            self._check
            

    def _check(self):
        if type(self._oriented_map) != OrientedMap:
            raise TypeError("The underlying surface must be an OrientedMap.")
        size = self._oriented_map.num_half_edges()
        previous = None
        for elt in self._walk:
            if elt >= size:
                raise ValueError(f"edge {elt} does not belong to the oriented map {self._oriented_map}.")
            elif not previous is None:
                current = self.oriented_map._ep(previous)
                while current != self.oriented_map._ep(previous) or current != elt:
                    current = self.oriented_map._vp[current]
                if current != elt:
                    raise ValueError(f"edge {elt} does not follow edge {previous}.")
            previous = elt

    
    def __eq__(self, other):
        return (self._oriented_map == other._oriented_map) and (self._walk == other._walk)

    def __len__(self):
        return len(self._walk)

    def __iter__(self):
        yield from self._walk

    def geodesic(self, quad_system=None):
        
        if quad_system is None:
            quad_system = QuadSystem(self._oriented_map)
        
        geodesic = Geodesic(quad_system)
        for e in self._walk:
            for f in quad_system._proj[e]:
                geodesic.add_edge(f)
        geodesic.canonical()
        return geodesic
    

    def is_homotopic(self, other, quad_system=None):
        r"""
        Return whether self and other are freely homotopic.
        
        EXAMPLES::
        
            sage: from combisurf import OrientedMap, QuadSystem, Geodesic, Walk
            sage: m = OrientedMap(vp=[[0, 2, 4, 6],[7, 8, 5], [9, 10, 12, 11], [3, 15, 1, 13, 14]])
            sage: w1 = Walk(m, [])
            sage: w2 = Walk(m, [4, 8, 11, 9, 7])
            sage: w3 = Walk(m, [6, 5])
            sage: w1.is_homotopic(w2)
            True
            sage: w1.is_homotopic(w3)
            False
            sage: w4 = Walk(m, [6, 8, 11, 9, 7])
            sage: w5 = Walk(m, [2, 15])
            sage: w3.is_homotopic(w4)
            True
            sage: w3.is_homotopic(w5)
            False
            sage: w6 = Walk(m, [11])
            sage: w3.is_homotopic(w6)
            True

        """
        if self._oriented_map != other._oriented_map:
            raise ValueError("The walk belong to differents maps.")

        geodesic_self = self.geodesic()
        geodesic_other = other.geodesic()
        
        if len(geodesic_self) != len(geodesic_other):
            return False

        c = geodesic_self._geodesic.copy()
        c.extend(c)
        return is_subword(geodesic_other._geodesic, c)


    def simplicity(self, edge_to_vertex=None, quad_system=None):

        r"""
        Return whether self lift into a simple walk in the universal covering.
        """
        
        if not self._walk:
            return True
        elif len(self._walk)==2 and self._oriented_map._ep(self._walk[0]) == self._walk[1]:
            return False
        if quad_system is None:
            quad_system = QuadSystem(self._oriented_map)
        if edge_to_vertex is None:
            edge_to_vertex = self._oriented_map.half_edges_to_vertices()

        star = StarShapedSpace(quad_system, edge_to_vertex[self._walk[0]])
        seen_value = {} 
        previous_vertex = 0
        first_values = (0, edge_to_vertex[self._walk[0]])
        seen_value[first_values] = True

        for i in range(len(self._walk)):
            edge = self._walk[i]
            for elt in quad_system._proj[edge]:
                previous_vertex = star.insert_edge(previous_vertex, elt)
            map_vertex = edge_to_vertex[self._oriented_map._ep(edge)]
            if seen_value.get((previous_vertex, map_vertex)) is None:
                seen_value[(previous_vertex, map_vertex)] = True
            elif i == len(self._walk) - 1 and (previous_vertex, map_vertex) == first_values:
                return True
            else:
                return False
        return True

