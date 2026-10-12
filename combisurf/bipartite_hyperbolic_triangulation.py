r"""
Bipartite hyperbolic triangulations

A *bipartite hyperbolic triangulation* is an oriented map (see
:mod:`~combisurf.oriented_map`) whose faces are either *filled faces*, all
triangles, or *punctured faces*, of any degree, such that:

- the triangles are coloured white and black, two triangles sharing an edge
  having distinct colours;

- every vertex with no punctured corner has degree at least `8`.

The punctured faces are removed from the surface, each of them being a
puncture (or a boundary component).

The colours are encoded in the labels: the three half-edges of a white
triangle are even and the three half-edges of a black triangle are odd. Since
the two half-edges of an edge have distinct parities, two triangles sharing an
edge then have distinct colours. There is no parity constraint on the
punctured faces: an edge between two punctured faces may be oriented either
way. The barycentric subdivisions of
:meth:`~combisurf.oriented_map.OrientedMap.barycentric_subdivision` follow
this convention, their white triangles being those whose corners read
(vertex, edge vertex, face vertex) counterclockwise.

With no punctured face, these are the reducing triangulations of
:ref:`colin-de-verdiere-despre-dubois-2024`, §3.1, on which every
non-contractible closed walk is freely homotopic to a strongly reduced closed
walk (Proposition 3.8 there), which is unique (Proposition 3.5 there). By
:ref:`dubois-2024`, Proposition 1.5, the geometric intersection of strongly
reduced closed walks of such a triangulation of a closed surface is their
intersection on the same map with every face punctured (Dubois calls
"reduced" the closed walks called strongly reduced here).

A walk on the map (see :mod:`~combisurf.oriented_map`) is a sequence of
half-edges where each one starts where the previous one ends. The vertex
permutation ``vp`` turns counterclockwise around the vertices and the face of a
half-edge `h` is the face on the left of `h`. The *turn* from `h_1` to `h_2`,
two half-edges such that `h_1` ends where `h_2` starts, at a vertex of degree
`d`, is the number of steps of ``vp`` from ``h1 ^ 1`` to `h_2`, in
`\{0, 1, \ldots, d - 1\}`. As in :mod:`~combisurf.hyperbolic_quadrangulation`,
the *weighted turn* gives the weight `1` to each corner of a filled face (a
*filled corner*) and the weight `4` to each corner of a punctured face; the
*weighted degree* `wd` of a vertex is its degree plus `3` per punctured corner,
and the weighted turn from `h_1` to `h_2` is the sum of the weights of the
corners swept counterclockwise from ``h1 ^ 1`` to `h_2`, that is the corners
on the right of the walk, taken modulo `wd` in `(-wd/2, wd/2]`. A positive
weighted turn turns right and a negative one turns left. The *colour* of the
turn from `h_1` to `h_2` is the parity of `h_1`: *white* when `h_1` is even and
*black* otherwise. When the face of `h_1` is filled, it is the colour of that
triangle; in general it is the colour opposite to the one of the corner
that follows ``h1 ^ 1`` counterclockwise, the first corner on the right.

Dictionary with :ref:`colin-de-verdiere-despre-dubois-2024`, §3.1: the
article counts the triangles on the left of the walk, so that its `k`-turn is
the weighted turn `-k` here, and its red triangles are the white triangles
here. Its turn `k_r` (resp. `k_b`) is the weighted turn `-k` with colour white
(resp. black). Its bad turns `0, 1, -1, 2_r, -2_r` are here the weighted turns
`0`, `\pm 1` and the white weighted turns `\pm 2`. A closed walk is *reduced*
when none of its weighted turns, read cyclically, is bad, and *strongly
reduced* when moreover its weighted turns are not all white `-3` (the
article's `3_r`) and not all black `3` (the article's `-3_b`). The empty walk
is strongly reduced. Reversing a walk negates its weighted turns, keeps the
colour of the turns `\pm 2` (the parities of `h_1` and of ``h2 ^ 1`` are equal
when two corners separate them) and changes the colour of the turns `\pm 3`,
so that it keeps (strong) reduction.

The weight `4` of a punctured corner is the least for which the moves of the
reduction (see :meth:`BipartiteHyperbolicTriangulation.strongly_reduced`)
never go through a punctured face. A side of a turn that contains a punctured
corner weighs at least `4`, so every weighted turn `0, \pm 1, \pm 2, \pm 3`,
the only values the moves read, counts filled corners on a side free of
punctures and means what it means in the article. At most one representative
of a turn lies in `[-3, 3]`, since two would need both sides to weigh at most
`3`, hence to be free of punctures, and the vertex to have degree at most `6`
with no punctured corner, which is excluded. A move changes the weighted turns
at its two ends by one corner of the strip of triangles it goes across, a
filled one, and gives the weighted turns `\mp 2`, `\mp 3` or `\mp 4` to the
vertices on the far side of the strip, with `2`, `3` or `4` corners in the
strip; their other side weighs at least `4` if it contains a punctured corner
and at least `5` otherwise (degree at least `8`), so these values are the
weighted turns there. With the weight `3`, a side made of a lone punctured
corner would weigh `3` and read as `\pm 3`: the boundary of a white punctured
triangle, with the triangle on its left, would read as white turns `-3` only,
and the last move of the reduction would push it across the puncture; at a
vertex with three filled corners and one punctured corner, the far side of a
strip would weigh `3` and the turn `-3` there would be read as `3`.

EXAMPLES:

The barycentric subdivision of a torus made of two quadrilaterals, with one
of them punctured, is a bipartite hyperbolic triangulation of the punctured
torus::

    sage: from combisurf import OrientedMap
    sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
    sage: m = OrientedMap(fp="(0,1,2,3)(~0,~1,~2,~3)")
    sage: t, punctured = m.barycentric_subdivision(punctured_faces=[1])
    sage: T = BipartiteHyperbolicTriangulation(t, punctured)
    sage: T.punctured_faces()
    [1]
"""

from array import array

from combisurf.hyperbolic_quadrangulation import _mod_normalize
from combisurf.oriented_map import OrientedMap

# the weighted turns that are bad whatever their colour, and the ones that the
# end of a flip must avoid once multiplied by the side s of the flip
_BAD = (0, 1, -1)
_NOT_FLIP_END = (-1, 0, 1, 2, 3)


class BipartiteHyperbolicTriangulation:
    r"""
    Bipartite hyperbolic triangulation.

    See :mod:`~combisurf.bipartite_hyperbolic_triangulation` for the
    definitions. The colour of a triangle is the parity of its half-edges:
    white for even half-edges and black for odd half-edges.

    .. NOTE::

        The attributes are the tables of
        :class:`~combisurf.hyperbolic_quadrangulation.HyperbolicQuadrangulation`,
        which turn half-edges into weighted turns and back in constant time:
        ``_h2v`` maps a half-edge to its vertex, ``_upos[h]`` and ``_wpos[h]``
        are the position and the weighted position of ``h`` around its
        vertex, ``_deg`` and ``_wdeg`` the degrees and the weighted degrees,
        and the entry ``_offset[v] + p`` of ``_at`` is the half-edge at the
        weighted position `p` around the vertex ``v``, or ``-1`` inside a
        punctured corner.

    INPUT:

    - ``m`` -- an :class:`~combisurf.oriented_map.OrientedMap` without folded
      edge; it may have inactive half-edges

    - ``punctured_faces`` -- (default: ``False``) a boolean or a list of face
      indices; the punctured faces, ``True`` standing for every face and
      ``False`` for none

    - ``check`` -- boolean (default: ``True``); whether to check that ``m`` is
      a connected bipartite hyperbolic triangulation and that the surface, once
      the punctured faces are removed, has negative Euler characteristic

    EXAMPLES:

    The barycentric subdivision of a torus made of two quadrilaterals, with
    the face ``1`` punctured, has `2 + 4 + 1 = 7` vertices, sixteen edges,
    eight triangles and the punctured face of degree `8`. Every vertex but the
    face vertex lies on the punctured face, and the face vertex has degree
    `8`. Its Euler characteristic is `0`, that is `-1` once the punctured face
    is removed::

        sage: from combisurf import OrientedMap
        sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
        sage: m = OrientedMap(fp="(0,1,2,3)(~0,~1,~2,~3)")
        sage: t, punctured = m.barycentric_subdivision(punctured_faces=[1])
        sage: t.num_vertices(), t.num_edges(), t.face_profile(), punctured
        (7, 16, [3, 8, 3, 3, 3, 3, 3, 3, 3], [1])
        sage: t.vertex_profile()
        [6, 3, 6, 3, 3, 3, 8]
        sage: t.euler_characteristic() - len(punctured)
        -1
        sage: T = BipartiteHyperbolicTriangulation(t, punctured)
        sage: T
        BipartiteHyperbolicTriangulation(OrientedMap("(0,~16,7,4,~20,3)(~0,~1,8)(1,6,~22,5,2,~18)(~2,~3,10)(~4,~5,12)(~6,~7,14)(~8,18,~10,20,~12,22,~14,16)", "(0,8,16)(~0,3,~2,5,~4,7,~6,1)(~1,~18,~8)(2,10,18)(~3,~20,~10)(4,12,20)(~5,~22,~12)(6,14,22)(~7,~16,~14)"), [1])

    The colours of the triangles are the parities of their half-edges::

        sage: [(f, 'white' if f[0] % 2 == 0 else 'black') for i, f in enumerate(t.faces())
        ....:  if i not in punctured]
        [([0, 16, 32], 'white'),
         ([3, 37, 17], 'black'),
         ([4, 20, 36], 'white'),
         ([7, 41, 21], 'black'),
         ([8, 24, 40], 'white'),
         ([11, 45, 25], 'black'),
         ([12, 28, 44], 'white'),
         ([15, 33, 29], 'black')]

    The barycentric subdivision of a closed surface of genus two made of one
    vertex and one face is not hyperbolic: its edge vertices have degree
    `4`::

        sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
        sage: t, _ = m.barycentric_subdivision()
        sage: t.num_vertices(), t.num_edges(), t.num_faces()
        (6, 24, 16)
        sage: t.vertex_profile()
        [16, 4, 4, 4, 4, 16]
        sage: BipartiteHyperbolicTriangulation(t)
        Traceback (most recent call last):
        ...
        ValueError: the vertex 1 touches no punctured face and has degree 4 instead of at least 8

    A face that is not punctured must be a triangle::

        sage: m = OrientedMap(fp="(0,1,2,3)(~0,~1,~2,~3)")
        sage: BipartiteHyperbolicTriangulation(m, [1])
        Traceback (most recent call last):
        ...
        ValueError: the face 0 is not punctured and has degree 4 instead of 3

    The labels must encode the colours of the triangles. Swapping the two
    half-edges of the edge between two triangles gives triangles with mixed
    parities::

        sage: t, punctured = m.barycentric_subdivision(punctured_faces=[1], mutable=True)
        sage: t.reverse_orientation(8)
        sage: BipartiteHyperbolicTriangulation(t, punctured)
        Traceback (most recent call last):
        ...
        ValueError: the triangle 0 with half-edges [0, 17, 32] has mixed parities; the labels must encode the 2-colouring of the triangles: white triangles even, black triangles odd

    The surface must have negative Euler characteristic once its punctured
    faces are removed::

        sage: t, _ = m.barycentric_subdivision()
        sage: BipartiteHyperbolicTriangulation(t)
        Traceback (most recent call last):
        ...
        ValueError: the surface must have negative Euler characteristic once its punctured faces are removed (got 0)

    With ``check=False``, none of the above is checked::

        sage: BipartiteHyperbolicTriangulation(t, check=False)
        BipartiteHyperbolicTriangulation(OrientedMap("(0,~16,7,~23,4,~20,3,~19)(~0,9,~1,8)(1,~17,6,~22,5,~21,2,~18)(~2,11,~3,10)(~4,13,~5,12)(~6,15,~7,14)(~8,18,~10,20,~12,22,~14,16)(~9,19,~11,21,~13,23,~15,17)", "(0,8,16)(~0,~19,~9)(1,9,17)(~1,~18,~8)(2,10,18)(~2,~21,~11)(3,11,19)(~3,~20,~10)(4,12,20)(~4,~23,~13)(5,13,21)(~5,~22,~12)(6,14,22)(~6,~17,~15)(7,15,23)(~7,~16,~14)"), [])
        sage: BipartiteHyperbolicTriangulation(m, check=False)
        BipartiteHyperbolicTriangulation(OrientedMap("(0,~3,2,~1)(~0,3,~2,1)", "(0,1,2,3)(~0,~1,~2,~3)"), [])

    TESTS:

    The input must be a connected oriented map without folded edge::

        sage: BipartiteHyperbolicTriangulation("(0,1,2)(~0,~1,~2)")
        Traceback (most recent call last):
        ...
        TypeError: m must be an oriented map
        sage: BipartiteHyperbolicTriangulation(OrientedMap(vp="(0,1,~1)"))
        Traceback (most recent call last):
        ...
        NotImplementedError: bipartite hyperbolic triangulations with folded edges are not implemented
        sage: m = OrientedMap(fp="(0,1,2)(~0,~1,~2)(3,4,5)(~3,~4,~5)")
        sage: BipartiteHyperbolicTriangulation(m)
        Traceback (most recent call last):
        ...
        ValueError: the map must be connected
    """
    def __init__(self, m, punctured_faces=False, check=True):
        if not isinstance(m, OrientedMap):
            raise TypeError("m must be an oriented map")
        if m.has_folded_edge():
            raise NotImplementedError("bipartite hyperbolic triangulations with folded edges are not implemented")
        m = m.copy(mutable=False)
        punctured = m._check_face_list(punctured_faces)
        faces = m.faces()
        is_punctured = [False] * len(faces)
        for f in punctured:
            is_punctured[f] = True
        h2f = m.half_edge_to_face()

        if check:
            if not m.is_connected():
                raise ValueError("the map must be connected")
            chi = m.euler_characteristic() - len(punctured)
            if chi >= 0:
                raise ValueError(f"the surface must have negative Euler characteristic once its punctured faces are removed (got {chi})")
            for i, f in enumerate(faces):
                if is_punctured[i]:
                    continue
                if len(f) != 3:
                    raise ValueError(f"the face {i} is not punctured and has degree {len(f)} instead of 3")
                if not f[0] % 2 == f[1] % 2 == f[2] % 2:
                    raise ValueError(f"the triangle {i} with half-edges {f} has mixed parities; the labels must encode the 2-colouring of the triangles: white triangles even, black triangles odd")
            for i, v in enumerate(m.vertices()):
                if len(v) < 8 and not any(is_punctured[h2f[h]] for h in v):
                    raise ValueError(f"the vertex {i} touches no punctured face and has degree {len(v)} instead of at least 8")

        h2v = m.half_edge_to_vertex()
        n = len(h2v)
        upos = array('i', [-1] * n)
        wpos = array('i', [-1] * n)
        deg = array('i')
        wdeg = array('i')
        offset = array('i')
        at = array('i')
        for vertex in m.vertices():
            offset.append(len(at))
            p = 0
            for d, h in enumerate(vertex):
                upos[h] = d
                wpos[h] = p
                at.append(h)
                if is_punctured[h2f[h]]:
                    at.extend((-1, -1, -1))
                    p += 4
                else:
                    p += 1
            deg.append(len(vertex))
            wdeg.append(p)

        self._m = m
        self._punctured_faces = punctured
        self._h2v = h2v
        self._upos = upos
        self._wpos = wpos
        self._deg = deg
        self._wdeg = wdeg
        self._offset = offset
        self._at = at

    def __repr__(self):
        return f"BipartiteHyperbolicTriangulation({self._m}, {self._punctured_faces})"

    def map(self):
        r"""
        Return the underlying oriented map.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
            sage: m = OrientedMap(fp="(0,1,2,3)(~0,~1,~2,~3)")
            sage: t, punctured = m.barycentric_subdivision(punctured_faces=[1])
            sage: T = BipartiteHyperbolicTriangulation(t, punctured)
            sage: T.map() == t
            True
        """
        return self._m

    def punctured_faces(self):
        r"""
        Return the sorted list of the indices of the punctured faces.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
            sage: m = OrientedMap(fp="(0,1,2,3)(~0,~1,~2,~3)")
            sage: t, punctured = m.barycentric_subdivision(punctured_faces=[1])
            sage: BipartiteHyperbolicTriangulation(t, punctured).punctured_faces()
            [1]

        Every face punctured::

            sage: BipartiteHyperbolicTriangulation(m, punctured_faces=True).punctured_faces()
            [0, 1]
        """
        return list(self._punctured_faces)

    def _check_turn(self, h1, h2):
        r"""
        Return ``(h1, h2)`` as Python integers after checking that they are
        half-edges and that ``h1`` ends where ``h2`` starts.
        """
        m = self._m
        h1 = m._check_half_edge(h1)
        h2 = m._check_half_edge(h2)
        if self._h2v[h1 ^ 1] != self._h2v[h2]:
            raise ValueError(f"the half-edge {h1} does not end where the half-edge {h2} starts")
        return h1, h2

    def turn(self, h1, h2):
        r"""
        Return the turn from ``h1`` to ``h2``.

        This is the number of steps of the vertex permutation from ``h1 ^ 1``
        to ``h2``, in `\{0, 1, \ldots, d - 1\}` where `d` is the degree of the
        vertex where ``h1`` ends and ``h2`` starts. It does not depend on the
        punctured faces.

        INPUT:

        - ``h1``, ``h2`` -- half-edges, ``h1`` ending where ``h2`` starts

        OUTPUT: an integer in `\{0, 1, \ldots, d - 1\}`

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
            sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
            sage: T = BipartiteHyperbolicTriangulation(m.triangle_system(relabel=True))
            sage: v = T.map().vertices()[0]
            sage: v
            [0, 27, 6, 25, 8, 29, 4, 21, 2, 23]
            sage: [T.turn(v[0] ^^ 1, h) for h in v]
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            sage: T.turn(0, 0)
            Traceback (most recent call last):
            ...
            ValueError: the half-edge 0 does not end where the half-edge 0 starts
        """
        h1, h2 = self._check_turn(h1, h2)
        g = h1 ^ 1
        return (self._upos[h2] - self._upos[g]) % self._deg[self._h2v[g]]

    def weighted_turn(self, h1, h2, check=True):
        r"""
        Return the weighted turn from ``h1`` to ``h2``.

        This is the sum of the weights of the corners swept counterclockwise
        from ``h1 ^ 1`` to ``h2``, the corners on the right of the walk, a
        filled corner weighing `1` and a punctured corner `4`, taken modulo
        the weighted degree `wd` of the vertex in `(-wd/2, wd/2]`. A positive
        weighted turn turns right and a negative one left. The `k`-turn of
        :ref:`colin-de-verdiere-despre-dubois-2024`, which leaves `k`
        triangles on the left, is the weighted turn `-k`; see
        :mod:`~combisurf.bipartite_hyperbolic_triangulation`.

        INPUT:

        - ``h1``, ``h2`` -- half-edges, ``h1`` ending where ``h2`` starts

        - ``check`` -- boolean (default: ``True``); whether to check that
          ``h1`` and ``h2`` are half-edges and that ``h1`` ends where ``h2``
          starts; with ``False``, the arguments must be valid Python integers
          and the call takes constant time without any validation

        OUTPUT: an integer in `(-wd/2, wd/2]`

        EXAMPLES:

        Around a vertex of the triangle system of a closed surface of genus
        two, of degree `10`, the weighted turns are the turns taken in
        `(-5, 5]`::

            sage: from combisurf import OrientedMap
            sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
            sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
            sage: T = BipartiteHyperbolicTriangulation(m.triangle_system(relabel=True))
            sage: v = T.map().vertices()[0]
            sage: [T.weighted_turn(v[0] ^^ 1, h) for h in v]
            [0, 1, 2, 3, 4, 5, -4, -3, -2, -1]

        The boundary of a triangle, which lies on its left, turns sharp left
        at each corner (the article's `1`-turns), and the reverse walk sharp
        right::

            sage: f = T.map().faces()[0]
            sage: T.weighted_turns(f)
            [-1, -1, -1]
            sage: T.weighted_turns([h ^^ 1 for h in reversed(f)])
            [1, 1, 1]

        On the barycentric subdivision of a torus made of two
        quadrilaterals, with one of them punctured, the vertex ``0`` has
        degree `6` and two punctured corners, the corners of the punctured
        face ``1``, so that its weighted degree is `12`. From ``0 ^^ 1``, the
        turns to ``8`` and to ``41`` sweep a punctured corner and the
        weighted turns read the other side::

            sage: m = OrientedMap(fp="(0,1,2,3)(~0,~1,~2,~3)")
            sage: t, punctured = m.barycentric_subdivision(punctured_faces=[1])
            sage: T = BipartiteHyperbolicTriangulation(t, punctured)
            sage: h2f = t.half_edge_to_face()
            sage: v = t.vertices()[0]
            sage: v, [h2f[h] for h in v]
            ([0, 33, 14, 8, 41, 6], [0, 8, 1, 5, 4, 1])
            sage: [T.turn(v[0] ^^ 1, h) for h in v]
            [0, 1, 2, 3, 4, 5]
            sage: [T.weighted_turn(v[0] ^^ 1, h) for h in v]
            [0, 1, 2, 6, -5, -4]
        """
        if check:
            h1, h2 = self._check_turn(h1, h2)
        g = h1 ^ 1
        wpos = self._wpos
        return _mod_normalize(wpos[h2] - wpos[g], self._wdeg[self._h2v[g]])

    def weighted_turns(self, w, closed=True):
        r"""
        Return the list of the weighted turns of the walk ``w``.

        The entry ``i`` is the weighted turn from ``w[i]`` to ``w[i + 1]``,
        the indices being taken modulo the length of ``w`` when ``w`` is
        closed.

        INPUT:

        - ``w`` -- a walk

        - ``closed`` -- boolean (default: ``True``); whether ``w`` is a closed
          walk, which is then checked; if ``False``, ``w`` is a walk between
          two vertices and there is no weighted turn from the last half-edge
          to the first

        OUTPUT: a list of integers, of the length of ``w`` if ``closed`` and
        one less (zero for the empty walk) otherwise

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
            sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
            sage: T = BipartiteHyperbolicTriangulation(m.triangle_system(relabel=True))
            sage: f = T.map().faces()[0]
            sage: T.weighted_turns(f)
            [-1, -1, -1]
            sage: T.weighted_turns(f[:2], closed=False)
            [-1]
            sage: T.weighted_turns([0, 1])
            [0, 0]
            sage: T.weighted_turns(f[:2])
            Traceback (most recent call last):
            ...
            ValueError: non-closed walk ...
        """
        w = self._m._check_walk(w, closed=closed)
        n = len(w)
        count = n if closed else max(n - 1, 0)
        return [self.weighted_turn(w[i], w[i + 1 if i + 1 < n else 0], check=False) for i in range(count)]

    def is_strongly_reduced(self, w, closed=True):
        r"""
        Return whether the walk ``w`` is strongly reduced.

        See :mod:`~combisurf.bipartite_hyperbolic_triangulation` for the
        definitions. A closed walk is strongly reduced when none of its
        weighted turns, read cyclically, is `0`, `\pm 1` or a white `\pm 2`,
        and when its weighted turns are not all white `-3` and not all black
        `3`; the empty walk is strongly reduced. Strong reduction is only
        defined for closed walks: when ``closed`` is ``False``, this returns
        whether the walk between two vertices is *reduced*, that is whether
        none of its weighted turns, read linearly, is bad. The check takes
        time `O(|w|)`.

        INPUT:

        - ``w`` -- a walk

        - ``closed`` -- boolean (default: ``True``); whether ``w`` is a closed
          walk, which is then checked, or a walk between two vertices

        OUTPUT: a boolean

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
            sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
            sage: T = BipartiteHyperbolicTriangulation(m.triangle_system(relabel=True))

        The boundary of a triangle and a backtrack are not reduced::

            sage: T.is_strongly_reduced(T.map().faces()[0])
            False
            sage: T.is_strongly_reduced([0, 1])
            False
            sage: T.is_strongly_reduced([])
            True

        A strongly reduced closed walk, whose weighted turn `2` is black::

            sage: w = [0, 3]
            sage: T.weighted_turns(w), [h % 2 for h in w]
            ([-4, 2], [0, 1])
            sage: T.is_strongly_reduced(w)
            True

        A closed walk whose weighted turns are all white `-3` is reduced but
        not strongly reduced, and the other side of the cylinder of triangles
        on its left has the white weighted turns `3`::

            sage: w = [0, 18, 24]
            sage: T.weighted_turns(w), [h % 2 for h in w]
            ([-3, -3, -3], [0, 0, 0])
            sage: T.is_strongly_reduced(w)
            False
            sage: r = T.strongly_reduced(w)
            sage: r, T.weighted_turns(r)
            ([24, 4, 12], [3, 3, 3])
            sage: T.is_strongly_reduced(r)
            True

        A walk between two vertices is reduced when none of its weighted turns
        is bad, the first and the last half-edges being free::

            sage: f = T.map().faces()[0]
            sage: T.is_strongly_reduced(f[:2], closed=False)
            False
            sage: T.is_strongly_reduced(f[:1], closed=False)
            True
        """
        w = self._m._check_walk(w, closed=closed)
        n = len(w)
        if not n:
            return True
        count = n if closed else n - 1
        wpos = self._wpos
        wdeg = self._wdeg
        h2v = self._h2v
        all_white_left = all_black_right = closed
        for i in range(count):
            h = w[i]
            g = h ^ 1
            t = _mod_normalize(wpos[w[i + 1 if i + 1 < n else 0]] - wpos[g], wdeg[h2v[g]])
            white = not (h & 1)
            if t in _BAD or (white and (t == 2 or t == -2)):
                return False
            if t != -3 or not white:
                all_white_left = False
            if t != 3 or white:
                all_black_right = False
        return not (all_white_left or all_black_right)

    def strongly_reduced(self, w, closed=True, check=True):
        r"""
        Return the strongly reduced closed walk freely homotopic to ``w``.

        See :mod:`~combisurf.bipartite_hyperbolic_triangulation` for the
        definitions. The result is empty when ``w`` is contractible in the
        surface without its punctured faces. Otherwise, it is strongly reduced
        (see :meth:`is_strongly_reduced`) and, by
        :ref:`colin-de-verdiere-despre-dubois-2024`, Proposition 3.5, unique up
        to cyclic rotation in the free homotopy class of ``w``, so that two
        closed walks are freely homotopic exactly when their results are
        rotations of each other. The rotation of the result is not normalised.
        The time is linear in the length of ``w``.

        INPUT:

        - ``w`` -- a closed walk

        - ``closed`` -- boolean (default: ``True``); only closed walks are
          implemented, ``False`` raises a ``NotImplementedError``

        - ``check`` -- boolean (default: ``True``); whether to check that
          ``w`` is a closed walk

        OUTPUT: a list of half-edges, empty when ``w`` is contractible

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
            sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
            sage: t, mor = m.triangle_system(relabel=True, mapping=True)
            sage: T = BipartiteHyperbolicTriangulation(t)

        The boundary of a triangle and a backtrack are contractible::

            sage: T.strongly_reduced(t.faces()[0])
            []
            sage: T.strongly_reduced([0, 1])
            []

        The image of the boundary of the face of ``m`` is contractible too::

            sage: w = [x for h in m.faces()[0] for x in mor[h]]
            sage: len(w)
            16
            sage: T.strongly_reduced(w)
            []

        The image of the edge ``0`` of ``m``, a loop, is strongly reduced; a
        conjugate of it reduces to a rotation of it::

            sage: w = list(mor[0])
            sage: w, T.weighted_turns(w)
            ([0, 3], [-4, 2])
            sage: T.strongly_reduced(w)
            [0, 3]
            sage: w2 = [x for h in [2, 0, 3] for x in mor[h]]
            sage: w2
            [4, 3, 0, 3, 2, 5]
            sage: T.strongly_reduced(w2)
            [3, 0]

        The reversal of the result is the result of the reversal, up to
        rotation::

            sage: T.strongly_reduced([h ^^ 1 for h in reversed(w2)])
            [1, 2]

        On a punctured surface, the boundary of a punctured face is not
        contractible. Here the barycentric subdivision of a torus made of two
        quadrilaterals, with one of them punctured: the boundary of the
        punctured face, which lies on its left, has the white weighted turns
        `2` at the edge vertices, which have two filled corners on its right,
        and is pushed across them::

            sage: m = OrientedMap(fp="(0,1,2,3)(~0,~1,~2,~3)")
            sage: t, punctured = m.barycentric_subdivision(punctured_faces=[1])
            sage: T = BipartiteHyperbolicTriangulation(t, punctured)
            sage: f = t.faces()[1]
            sage: f, T.weighted_turns(f)
            ([1, 6, 5, 10, 9, 14, 13, 2], [-4, 2, -4, 2, -4, 2, -4, 2])
            sage: r = T.strongly_reduced(f)
            sage: r, T.weighted_turns(r)
            ([45, 40, 33, 44, 37, 32, 41, 36], [-2, 6, -2, 6, -2, 6, -2, 6])
            sage: T.is_strongly_reduced(r)
            True

        TESTS:

        The empty walk, a walk that is not closed, and walks between two
        vertices::

            sage: T.strongly_reduced([])
            []
            sage: T.strongly_reduced([0])
            Traceback (most recent call last):
            ...
            ValueError: non-closed walk ...
            sage: T.strongly_reduced([0], closed=False)
            Traceback (most recent call last):
            ...
            NotImplementedError: strongly reduced walks are only implemented for closed walks

        .. SEEALSO::

            :meth:`is_strongly_reduced`,
            :meth:`~combisurf.hyperbolic_quadrangulation.HyperbolicQuadrangulation.geodesic`

        ALGORITHM:

        :ref:`colin-de-verdiere-despre-dubois-2024`, Proposition 3.8. The
        weighted turns are stored as a circular doubly linked list of runs,
        the *compressed homotopy sequence*: maximal runs of `3` and of `-3`,
        and single other values, each with the half-edge that starts it, see
        :class:`TurnRuns`. With `s = \pm 1` the side of the move (`s = -1`
        for the moves of the article, `s = 1` for the moves of the article on
        the reversed walk), ``x`` and ``y`` the weighted turns before and
        after the pattern, the moves are, every value being taken modulo the
        weighted degree of its vertex:

        - spur: ``x, 0, y`` becomes ``x + y``;

        - spike: ``x, s, y`` becomes ``x - s, y - s``, the two sides of the
          triangle being replaced by the third one;

        - bracket: ``x, 2s, (3s)^k, 2s, y`` becomes ``x - s, (-3s)^(k+1), y -
          s``, the walk going along the other side of the strip of triangles;

        - flip: ``x, (3s)^k, 2s, (3s)^l, y`` with a white `2s` and ``s x``,
          ``s y`` not in `\{-1, 0, 1, 2, 3\}` becomes ``x - s, -2s,
          (-3s)^(k-1), -4s, (-3s)^(l-1), -2s, y - s`` when `k, l \geq 1`,
          ``x - s, (-3s)^l, -2s, y - s`` when `k = 0 < l`, ``x - s, -2s,
          (-3s)^k, y - s`` when `l = 0 < k`, and ``x - s, -2s, y - s`` when
          `k = l = 0`.

        When ``x`` and ``y`` are the same weighted turn, the pattern covering
        the rest of the walk, it becomes the single weighted turn ``x - 2s``
        (``x + y`` and the empty walk for a spur). This covers the heart-like
        case of Fig. 3.6 of the article: a flip of ``4s, (3s)^k, 2s,
        (3s)^l`` gives the white weighted turn `2s` at the vertex of ``x``,
        which the next flip removes. A move touches at most six runs and
        inserts at most seven; the scan then resumes two runs before them.
        The moves decrease `3` times the length plus the number of bad turns,
        except for the heart-like case, so that their number is linear. Last,
        a walk whose weighted turns are all white `-3` or all black `3` goes
        around a cylinder of triangles on the side `s` of the sign of its
        turns; it is replaced by the other side of the cylinder, whose
        weighted turns are all `-3s` (Fig. 3.7 of the article).
        """
        if not closed:
            raise NotImplementedError("strongly reduced walks are only implemented for closed walks")
        if check:
            w = self._m._check_walk(w, closed=True)
        if not len(w):
            return []
        return TurnRuns(self, w).strongly_reduced()


class TurnRuns:
    r"""
    The weighted turns of a closed walk of a
    :class:`BipartiteHyperbolicTriangulation`, run-length encoded in a
    circular doubly linked list, and their reduction to the strongly reduced
    closed walk.

    A run is a node with a value, the weighted turn, a length, the number of
    consecutive weighted turns equal to the value, and a first half-edge, the
    half-edge ``h`` such that the first weighted turn of the run is the one
    from ``h`` to the next half-edge of the walk. Only the runs of `\pm 3` are
    merged and can be longer than one. A weighted turn `\pm 3` keeps the
    parity of the half-edges, so that all the turns of a run have the colour
    of its first half-edge.

    INPUT:

    - ``T`` -- a :class:`BipartiteHyperbolicTriangulation`

    - ``w`` -- a list of half-edges of ``T``, a closed walk of positive length

    EXAMPLES::

        sage: from combisurf import OrientedMap
        sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation, TurnRuns
        sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
        sage: T = BipartiteHyperbolicTriangulation(m.triangle_system(relabel=True))
        sage: TurnRuns(T, [0, 1]).strongly_reduced()
        []
    """
    def __init__(self, T, w):
        h2v = self._h2v = T._h2v
        wpos = self._wpos = T._wpos
        wdeg = self._wdeg = T._wdeg
        self._offset = T._offset
        self._at = T._at

        n = len(w)
        val = []
        length = []
        first = []
        for i in range(n):
            h = w[i]
            g = h ^ 1
            t = _mod_normalize(wpos[w[i + 1 if i + 1 < n else 0]] - wpos[g], wdeg[h2v[g]])
            if val and val[-1] == t and (t == 3 or t == -3):
                length[-1] += 1
            else:
                val.append(t)
                length.append(1)
                first.append(h)
        if len(val) > 1 and val[0] == val[-1] and (val[0] == 3 or val[0] == -3):
            # the run at the end goes on at the start: it absorbs the first
            # one and keeps its own first half-edge, which comes before
            length[-1] += length[0]
            del val[0], length[0], first[0]
        k = len(val)
        self._val = val
        self._len = length
        self._first = first
        self._next = list(range(1, k)) + [0]
        self._prev = [k - 1] + list(range(k - 1))
        self._alive = [True] * k
        self._dirty = [True] * k
        self._size = k
        self._num_dirty = k
        self._num_turns = n

    def _half_edge_at(self, h, t):
        r"""
        Return the half-edge at weighted position ``t`` counterclockwise from
        ``h`` around its vertex.
        """
        v = self._h2v[h]
        x = self._at[self._offset[v] + (self._wpos[h] + t) % self._wdeg[v]]
        if x == -1:
            raise RuntimeError(f"the weighted position {t} from the half-edge {h} is inside a punctured corner")
        return x

    def _norm(self, t, h):
        r"""
        Return ``t`` modulo the weighted degree of the vertex of ``h``.
        """
        return _mod_normalize(t, self._wdeg[self._h2v[h]])

    def _new(self, value, length, first):
        i = len(self._val)
        self._val.append(value)
        self._len.append(length)
        self._first.append(first)
        self._next.append(i)
        self._prev.append(i)
        self._alive.append(True)
        self._dirty.append(True)
        self._num_dirty += 1
        self._size += 1
        return i

    def _remove(self, u):
        p = self._prev[u]
        q = self._next[u]
        self._next[p] = q
        self._prev[q] = p
        self._alive[u] = False
        if self._dirty[u]:
            self._dirty[u] = False
            self._num_dirty -= 1
        self._size -= 1

    def _insert_after(self, a, runs):
        r"""
        Insert new runs after the node ``a`` and return the node of the first.
        """
        b = self._next[a]
        first = None
        for value, length, h in runs:
            u = self._new(value, length, h)
            self._next[a] = u
            self._prev[u] = a
            a = u
            if first is None:
                first = u
        self._next[a] = b
        self._prev[b] = a
        return first

    def _merge(self, start, count):
        r"""
        Merge the equal runs of `\pm 3` among the ``count`` adjacencies that
        follow the node ``start`` and return a live node among them.
        """
        val = self._val
        u = start
        for _ in range(count):
            v = self._next[u]
            if u != v and val[u] == val[v] and (val[u] == 3 or val[u] == -3):
                # u comes first and keeps its first half-edge
                self._len[u] += self._len[v]
                self._remove(v)
            else:
                u = v
        return start if self._alive[start] else u

    def _mark_dirty(self, start, count):
        r"""
        Mark dirty ``count`` nodes from two nodes before ``start`` and return
        the first of them.

        A move is looked for from its *key* run: the run `0` of a spur, `\pm
        1` of a spike, the first `\pm 2` of a bracket, the `\pm 2` of a flip.
        Whether a move applies depends only on the runs at most two runs away
        from its key run (the run after a bracket, three runs after its key,
        is rewritten but its value is not read), so that the moves created by
        a change have their key runs at most two runs before or after the
        changed ones.
        """
        u = self._prev[self._prev[start]]
        x = u
        for _ in range(min(count, self._size)):
            if not self._dirty[x]:
                self._dirty[x] = True
                self._num_dirty += 1
            x = self._next[x]
        return u

    def _splice(self, P, middle, N, first_N, runs):
        r"""
        Replace the last weighted turn of ``P``, the runs ``middle`` and the
        first weighted turn of ``N`` by ``runs``, mark the runs around them
        dirty and return the node to resume the scan from.

        ``P`` and ``N`` may be the same run; ``first_N`` is the first
        half-edge of ``N`` once its first weighted turn is removed.
        """
        length = self._len
        removed = 2
        for u in middle:
            removed += length[u]
            self._remove(u)
        new = self._insert_after(P, runs)
        length[P] -= 1
        length[N] -= 1
        self._first[N] = first_N
        if length[P] == 0:
            self._remove(P)
        if N != P and length[N] == 0:
            self._remove(N)
        self._num_turns += sum(r[1] for r in runs) - removed
        start = self._merge(self._prev[new], len(runs) + 1)
        return self._mark_dirty(start, len(runs) + 6)

    def _replace_all(self, nodes, runs):
        r"""
        Replace the whole cycle, made of ``nodes``, by ``runs``, mark it dirty
        and return a node of it.
        """
        for u in nodes:
            self._remove(u)
        a = self._new(*runs[0])
        for value, length, h in runs[1:]:
            self._insert_after(self._prev[a], [(value, length, h)])
        self._num_turns = sum(r[1] for r in runs)
        a = self._merge(a, len(runs))
        return self._mark_dirty(a, self._size)

    def _frame(self, middle):
        r"""
        Return ``(rest, P, N)``: the number of weighted turns of the walk out
        of the runs ``middle`` and the runs before and after them.
        """
        pattern = 0
        for u in middle:
            pattern += self._len[u]
        return self._num_turns - pattern, self._prev[middle[0]], self._next[middle[-1]]

    def _rewrite(self, middle, s, rest, P, N, inner, g_first, g_last):
        r"""
        Replace the part of the walk made of the runs ``middle`` by a walk
        from ``g_first`` to ``g_last`` with the weighted turns ``inner``
        between them, and return the node to resume the scan from.

        The part goes from the start of ``h_0``, the first half-edge of
        ``middle[0]``, to the end of ``h_m``, the first half-edge of ``N``.
        The new walk leaves the start of ``h_0`` one filled corner away from
        ``h_0`` on the side ``s``, along ``g_first``, and arrives at the end
        of ``h_m`` one filled corner away from ``h_m ^ 1`` on the side ``s``,
        along ``g_last``: the weighted turns ``x`` before the part and ``y``
        after it become ``x - s`` and ``y - s``. When they are the same
        weighted turn, ``rest == 1``, it becomes ``x - 2s``, from ``g_last``
        to ``g_first``.
        """
        h0 = self._first[middle[0]]
        x = self._val[P]
        if rest == 1:
            runs = [(self._norm(x - 2 * s, h0), 1, g_last)] + inner
            return self._replace_all(middle + [P], runs)
        hm = self._first[N]
        y = self._val[N]
        runs = [(self._norm(x - s, h0), 1, self._half_edge_at(h0, -x) ^ 1)]
        runs.extend(inner)
        runs.append((self._norm(y - s, hm ^ 1), 1, g_last))
        return self._splice(P, middle, N, self._half_edge_at(hm ^ 1, y), runs)

    def _spur(self, u):
        r"""
        Remove the spur at the run ``u``, of value `0`, and return the node to
        resume from, ``None`` when the walk vanishes.

        Weighted turns ``x, 0, y`` from ``h_{-1}``, ``h_0 = first[u]``,
        ``h_0 ^ 1`` and ``h_2`` become the single weighted turn ``x + y``
        from ``h_{-1}`` to ``h_2``, at the vertex where ``h_0`` starts.
        """
        rest, P, N = self._frame([u])
        if rest <= 1:
            return None
        h0 = self._first[u]
        x = self._val[P]
        y = self._val[N]
        runs = [(self._norm(x + y, h0), 1, self._half_edge_at(h0, -x) ^ 1)]
        return self._splice(P, [u], N, self._half_edge_at(h0, y), runs)

    def _spike(self, u, s):
        r"""
        Remove the spike at the run ``u``, of value ``s = ±1``, and return the
        node to resume from, ``-1`` if there is none.

        The walk ``h_0, h_1`` goes along two sides of the triangle on its side
        ``s``, and is replaced by its third side.
        """
        rest, P, N = self._frame([u])
        if rest < 1:
            return -1
        f = self._half_edge_at(self._first[u], -s)
        return self._rewrite([u], s, rest, P, N, [], f, f)

    def _bracket(self, u, s):
        r"""
        Remove the bracket ``2s, (3s)^k, 2s`` that starts at the run ``u``,
        if there is one, and return the node to resume from, ``-1`` if there
        is none.

        The walk ``h_0, ..., h_{k+2}`` goes along one side of a strip of
        `2k + 3` triangles on its side ``s``; the other side of the strip is
        a walk of length `k + 2` with the weighted turns `-3s`, whose first
        half-edge is at weighted position ``-s`` from ``h_0``.
        """
        val = self._val
        b = self._next[u]
        if val[b] == 3 * s:
            k = self._len[b]
            c = self._next[b]
            middle = [u, b, c]
        else:
            k = 0
            c = b
            middle = [u, c]
        if val[c] != 2 * s:
            return -1
        rest, P, N = self._frame(middle)
        if rest < 1:
            return -1
        g_first = self._half_edge_at(self._first[u], -s)
        g_last = self._half_edge_at(self._first[N] ^ 1, s) ^ 1
        return self._rewrite(middle, s, rest, P, N, [(-3 * s, k + 1, g_first)], g_first, g_last)

    def _flip(self, u, s):
        r"""
        Flip the run ``u``, of value ``2s``, if it is white and is the middle
        of a flip ``x, (3s)^k, 2s, (3s)^l, y``, and return the node to resume
        from, ``-1`` otherwise.

        The walk ``h_0, ..., h_{k+l+1}`` goes along one side of a band of
        triangles on its side ``s``: the triangle ``U_j`` of each ``h_j`` and,
        between consecutive ones at a weighted turn `3s`, a triangle ``D_j``.
        At the corner of the weighted turn `2s`, the triangles ``U_k`` and
        ``U_{k+1}`` share the edge ``c``, at weighted position ``s`` from
        ``h_k ^ 1``, whose other end ``p`` is a vertex of the other side of
        the band, reached by both ``U_k`` and ``U_{k+1}``. The other side of
        the band is the new walk: it has the same length, its weighted turns
        are `-2s` at the first and the last vertices after and before the
        ends, `-3s` along the band and `-(2 + [k > 0] + [l > 0]) s` at ``p``,
        where it arrives at weighted position ``s (1 + [k > 0])`` from ``c ^
        1`` and leaves at weighted position ``-s (1 + [l > 0])``.
        """
        first = self._first
        if first[u] & 1:
            return -1
        val = self._val
        p = self._prev[u]
        q = self._next[u]
        if val[p] == 3 * s:
            k = self._len[p]
            middle = [p, u]
        else:
            k = 0
            middle = [u]
        if val[q] == 3 * s:
            l = self._len[q]
            middle.append(q)
        else:
            l = 0
        rest, P, N = self._frame(middle)
        if rest < 1 or s * val[P] in _NOT_FLIP_END or s * val[N] in _NOT_FLIP_END:
            return -1
        at = self._half_edge_at
        g_first = at(first[middle[0]], -s)
        g_last = at(first[N] ^ 1, s) ^ 1
        c1 = at(first[u] ^ 1, s) ^ 1
        inner = []
        if k:
            inner.append((-2 * s, 1, g_first))
            if k > 1:
                inner.append((-3 * s, k - 1, at(g_first ^ 1, -2 * s)))
            e_in = at(c1, 2 * s) ^ 1
        else:
            e_in = g_first
        apex = -s * (2 + (k > 0) + (l > 0))
        inner.append((self._norm(apex, c1), 1, e_in))
        if l:
            if l > 1:
                inner.append((-3 * s, l - 1, at(c1, -2 * s)))
            inner.append((-2 * s, 1, at(g_last, 2 * s) ^ 1))
        return self._rewrite(middle, s, rest, P, N, inner, g_first, g_last)

    def _cylinder(self, u):
        r"""
        Replace the walk, whose weighted turns are all ``3s`` and form the
        single run ``u``, by the other side of the cylinder of triangles on
        its side ``s``, and return its node.

        Each half-edge ``h_j`` has a triangle ``U_j`` on its side ``s``, and
        at each vertex a triangle ``D_j`` lies between ``U_j`` and
        ``U_{j+1}``. The other side of the cylinder goes along the sides of
        the ``D_j`` away from the walk and has the weighted turns `-3s`. Its
        half-edge in ``D_0`` is at weighted position ``-s`` from the reverse
        of the side of ``U_0`` at weighted position ``s`` from ``h_0 ^ 1``.
        """
        s = self._val[u] // 3
        h = self._first[u]
        f = self._half_edge_at(self._half_edge_at(h ^ 1, s) ^ 1, -s)
        return self._replace_all([u], [(-3 * s, self._num_turns, f)])

    def _reduce(self):
        r"""
        Apply the moves until there is none and return a node of the walk,
        or ``None`` when the walk vanishes.

        All the runs are dirty at first. At a dirty run that is the key run
        of no move, the run is marked clean and the cursor moves on; a move
        marks dirty the runs around the new ones and moves the cursor back to
        the first of them.
        """
        val = self._val
        dirty = self._dirty
        nxt = self._next
        limit = 8 * self._num_turns + 16
        moves = 0
        cur = 0
        while self._num_dirty:
            if not dirty[cur]:
                cur = nxt[cur]
                continue
            t = val[cur]
            if t == 0:
                r = self._spur(cur)
                if r is None:
                    return None
            elif t == 1 or t == -1:
                r = self._spike(cur, t)
            elif t == 2 or t == -2:
                s = t // 2
                r = self._bracket(cur, s)
                if r == -1:
                    r = self._flip(cur, s)
            else:
                r = -1
            if r == -1:
                dirty[cur] = False
                self._num_dirty -= 1
                cur = nxt[cur]
            else:
                cur = r
                moves += 1
                if moves > limit:
                    raise RuntimeError("the reduction does not terminate; the triangulation is not a bipartite hyperbolic triangulation")
        return cur

    def strongly_reduced(self):
        r"""
        Reduce the walk to the strongly reduced one and return it, the empty
        list when it is contractible.
        """
        u = self._reduce()
        if u is None:
            return []
        if self._size == 1:
            t = self._val[u]
            # all white -3 or all black 3
            if (t == -3 or t == 3) and (self._first[u] & 1) == (t == 3):
                u = self._cylinder(u)
        return self._walk(u)

    def _walk(self, u):
        r"""
        Return the walk, starting from the first half-edge of the run ``u``.
        """
        val = self._val
        length = self._len
        first = self._first
        nxt = self._next
        out = []
        for _ in range(self._size):
            h = first[u]
            t = val[u]
            for _ in range(length[u]):
                out.append(h)
                h = self._half_edge_at(h ^ 1, t)
            u = nxt[u]
        return out
