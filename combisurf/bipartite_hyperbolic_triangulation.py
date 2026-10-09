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
non-contractible closed walk is freely homotopic to a unique strongly reduced
closed walk. By :ref:`dubois-2024`, Proposition 1.5, the geometric
intersection of reduced closed walks of such a triangulation of a closed
surface is their intersection on the same map with every face punctured.

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

from combisurf.oriented_map import OrientedMap


class BipartiteHyperbolicTriangulation:
    r"""
    Bipartite hyperbolic triangulation.

    See :mod:`~combisurf.bipartite_hyperbolic_triangulation` for the
    definitions. The colour of a triangle is the parity of its half-edges:
    white for even half-edges and black for odd half-edges.

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

        if check:
            if not m.is_connected():
                raise ValueError("the map must be connected")
            chi = m.euler_characteristic() - len(punctured)
            if chi >= 0:
                raise ValueError(f"the surface must have negative Euler characteristic once its punctured faces are removed (got {chi})")
            faces = m.faces()
            is_punctured = [False] * len(faces)
            for f in punctured:
                is_punctured[f] = True
            for i, f in enumerate(faces):
                if is_punctured[i]:
                    continue
                if len(f) != 3:
                    raise ValueError(f"the face {i} is not punctured and has degree {len(f)} instead of 3")
                if not f[0] % 2 == f[1] % 2 == f[2] % 2:
                    raise ValueError(f"the triangle {i} with half-edges {f} has mixed parities; the labels must encode the 2-colouring of the triangles: white triangles even, black triangles odd")
            h2f = m.half_edge_to_face()
            for i, v in enumerate(m.vertices()):
                if len(v) < 8 and not any(is_punctured[h2f[h]] for h in v):
                    raise ValueError(f"the vertex {i} touches no punctured face and has degree {len(v)} instead of at least 8")

        self._m = m
        self._punctured_faces = punctured

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
