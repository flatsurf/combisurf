r"""
Morphism
"""
# ****************************************************************************
#  This file is part of combisurf
#
#       Copyright (C) 2026 Vincent Delecroix
#                     2026 Oscar Fontaine
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

from combisurf.word import word_init, word_reduce
from combisurf.oriented_map import OrientedMap

class OrientedMapMorphism:
    r"""
    Morphism between oriented maps.

    A morphism is a partial maps between walks in the domain to walks in the
    codomain. The domain must be a subset of walks stable under concatenation
    and the map should respect concatenation. More precisely, if walks
    ``u`` and ``v`` belongs to the domain and for a morphism ``f`` both
    ``f(u)`` and ``f(v)`` are defined then we impose that ``f(uv)`` exists
    and assume that ``f(u) f(v) = f(uv)``.
    """
    def __init__(self, domain_map, codomain_map):
        self._domain_map = domain_map
        self._codomain_map = codomain_map

    def domain_map(self):
        return self._domain_map

    def codomain_map(self):
        return self._codomain_map


class OrientedMapMorphism_list(OrientedMapMorphism):
    r"""
    EXAMPLES::

        sage: from combisurf import OrientedMap
        sage: m = OrientedMap(fp="(0,1,~0,~1)")
        sage: f = m.hom([[0,2], [3,1], [0], [1]])
        sage: f
        OrientedMapMorphism(OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)"), OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)"), [[0, 2], [3, 1], [0], [1]])
        sage: f([0, 2, 0])
        array('i', [0, 2, 0, 0, 2])
        sage: f([0, 1], reduce=True)
        array('i')
    """
    def __init__(self, domain_map, codomain_map, images):
        OrientedMapMorphism.__init__(self, domain_map, codomain_map)
        images = list(map(word_init, images))
        if len(images) != len(self._domain_map._vp):
            raise ValueError("wrong data to specify images")
        self._images = images

    def __repr__(self):
        return f"OrientedMapMorphism({self._domain_map}, {self._codomain_map}, {[list(self._half_edge_image(h)) for h in self._domain_map.half_edges()]})"

    def _half_edge_image(self, h, check=True):
        if check:
            h = self._domain_map._check_half_edge(h)
        return self._images[h]

    def _array_image(self, w, reduce=False):
        ans = array("i")
        if reduce:
            w = word_reduce(w)
            for h in w:
                word_free_group_mul_inplace(ans, self._half_edge_image(h))
        else:
            for h in w:
                ans += self._half_edge_image(h)
        return ans

    # TODO: make __call__ take a walk on the domain as input
    def __call__(self, w, reduce=False, check=True):
        if check:
            w = word_init(w)
        return self._array_image(w, reduce)
