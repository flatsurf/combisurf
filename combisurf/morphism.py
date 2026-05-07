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
    def __init__(self, domain, codomain):
        self._domain = domain
        self._codomain = codomain

    def domain(self):
        return self._domain

    def codomain(self):
        return self._codomain

    def _image(self, h, check=True):
        r"""
        Return the image of the half-edge ``h`` as an array of integers.
        """
        raise NotImplementedError

    def __repr__(self):
        return f"OrientedMapMorphism({self._domain}, {self._codomain}, {[self._image(h) for h in self._domain.half_edges()]})"

    # TODO: should we have an option to remove spur in domain/codomain?
    def __call__(self, w, reduce=False, check=True):
        if check:
            w = word_init(w)
        ans = array("i")
        if reduce:
            w = word_reduce(w)
            for h in w:
                word_free_group_mul_inplace(ans, self._image(h))
        else:
            for h in w:
                ans += self._image(h)
        return ans


class OrientedMapMorphism_list(OrientedMapMorphism):
    r"""
    EXAMPLES::

        sage: from combisurf import OrientedMap
        sage: m = OrientedMap(fp="(0,1,~0,~1)")
        sage: f = m.hom([[0,2], [3,1], [0], [1]])
        sage: f([0, 2, 0])
        array('i', [0, 2, 0, 0, 2])
        sage: f([0, 1], reduce=True)
        array('i')
    """
    def __init__(self, domain, codomain, images):
        OrientedMapMorphism.__init__(self, domain, codomain)
        images = list(map(word_init, images))
        if len(images) != len(self._domain._vp):
            raise ValueError("wrong data to specify images")
        self._images = images

    def _image(self, h, check=True):
        if check:
            h = self.domain()._check_half_edge(h)
        return self._images[h]
