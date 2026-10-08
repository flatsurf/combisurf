r"""
Hyperbolic quadrangulations with boundary and their geodesics

A *quadrangulation with boundary* is an oriented map (see
:mod:`~combisurf.oriented_map`) whose faces are either *interior faces*, all
of degree `4`, or *punctured faces*, of any degree. The punctured faces are
removed from the surface, each of them being a puncture (or a boundary
component). The quadrangulation is *hyperbolic* when every vertex touches a
punctured face or has degree at least `5`. The quad systems of
:meth:`~combisurf.oriented_map.OrientedMap.quad_system` and the radial maps
of :meth:`~combisurf.oriented_map.OrientedMap.radial_map`, punctured at their
bigons, are the main examples.

A walk on the map (see :mod:`~combisurf.oriented_map`) is a sequence of
half-edges where each one starts where the previous one ends. The vertex
permutation ``vp`` turns counterclockwise around the vertices, and the face
whose ``fp``-cycle contains a half-edge `h` is the face on the left of `h`.
Since ``fp[vp[g] ^ 1] == g``, that face reaches the vertex of `g` along
``vp[g] ^ 1`` and leaves it along `g`: the corner between `g` and ``vp[g]``,
counterclockwise from `g`, is a corner of the face of `g`::

          vp[g]
            \
             \   corner of the face of g
              \
               v ---------- g

The *turn* from `h_1` to `h_2`, two half-edges such that `h_1` ends where
`h_2` starts, at a vertex of degree `d`, is the number `t` in
`\{0, 1, \ldots, d - 1\}` of steps of ``vp`` from ``h1 ^ 1`` to `h_2`. It does
not depend on the punctures. A turn `0` is a *spur*: the walk goes back along
the edge it came from.

The reduction of walks works on the *weighted turn* instead. Each corner of
an interior face weighs `1` and each corner of a punctured face weighs `4`;
the *weighted degree* of a vertex is the sum of the weights of its corners,
that is its degree plus `3` per punctured corner. The weighted turn from `h_1`
to `h_2` is the sum of the weights of the corners swept counterclockwise from
``h1 ^ 1`` to `h_2`, taken modulo the weighted degree `wd` and represented in
`(-wd/2, wd/2]`. At a vertex without punctured corner it is the turn, taken
in `(-d/2, d/2]`.

A positive weighted turn turns right and a negative one turns left: walking
along the boundary of an interior face, which lies on the left, gives the
weighted turn `-1` at each corner, a *sharp left turn*. The articles below
count corners clockwise and call this turn `1`; their turns are the opposite
of the weighted turns here.

When one side of a turn sweeps no punctured corner and at most `3` interior
corners, the weighted turn is the signed number of these corners, positive
counterclockwise, as in the integer turns of :ref:`despre-lazarus-2019`,
§4.2, since the other side weighs at least `4`. When both sides sweep a
punctured corner, its absolute value is at least `4`. So the weighted turns
`0, \pm 1, \pm 2, \pm 3` count interior corners on a side without punctured
corner, and the reductions below never go through a punctured corner.

A *bracket* is a part of a closed walk whose weighted turns are
`1, 2, \ldots, 2, 1` or `-1, -2, \ldots, -2, -1`, with any number of `\pm 2`,
possibly none: the walk goes along one side of a strip of interior faces,
and the other side of the strip is shorter by two. A closed walk is a
*geodesic* when it has no spur and no bracket, the indices of the weighted
turns being taken cyclically. On a hyperbolic quadrangulation with boundary,
every closed walk is freely homotopic to a geodesic, obtained by removing
spurs and brackets, and the homotopic geodesics all have the same length:
the minimal length in the free homotopy class, `0` for a contractible closed
walk. See :ref:`erickson-whittlesey-2013`, Lemma 4.3, and
:ref:`despre-lazarus-2019`, Theorem 9 and Corollary 11.

EXAMPLES:

The quad system of the punctured torus has two vertices of degree `5`, five
edges, two interior faces and the bigon of the puncture::

    sage: from combisurf import OrientedMap
    sage: from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation
    sage: m = OrientedMap(vp=[[0, 2, 1, 3]])
    sage: q, mor, bigons = m.quad_system(punctured_faces=True, relabel=True, mapping=True)
    sage: q
    OrientedMap("(0,2,1,3,4)(~0,~4,~2,~1,~3)", "(0,~3,1,~2)(~0,4)(~1,2,~4,3)")
    sage: q.faces(), bigons
    ([[0, 7, 2, 5], [1, 8], [3, 4, 9, 6]], [1])
    sage: Q = HyperbolicQuadrangulation(q, bigons)

The commutator of the two generators goes around the puncture: its
geodesic is the boundary of the bigon::

    sage: w = [h for x in [0, 2, 1, 3] for h in mor[x]]
    sage: w
    [8, 5, 4, 3, 4, 9, 2, 5]
    sage: Q.weighted_turns(w)
    [1, 0, 1, -1, -1, -2, -1, 3]
    sage: Q.geodesic(w)
    [1, 8]
    sage: Q.is_geodesic([1, 8])
    True
"""

from array import array

from combisurf.oriented_map import OrientedMap


def _normalize(t, d):
    r"""
    Return the representative of ``t`` modulo ``d`` in `(-d/2, d/2]`.

    TESTS::

        sage: from combisurf.hyperbolic_quadrangulation import _normalize
        sage: [_normalize(t, 6) for t in range(-6, 7)]
        [0, 1, 2, 3, -2, -1, 0, 1, 2, 3, -2, -1, 0]
        sage: [_normalize(t, 5) for t in range(-5, 6)]
        [0, 1, 2, -2, -1, 0, 1, 2, -2, -1, 0]
    """
    t %= d
    return t - d if 2 * t > d else t


class HyperbolicQuadrangulation:
    r"""
    Hyperbolic quadrangulation with boundary.

    See :mod:`~combisurf.hyperbolic_quadrangulation` for the definitions.

    INPUT:

    - ``m`` -- an :class:`~combisurf.oriented_map.OrientedMap` without folded
      edge; it may have inactive half-edges

    - ``punctured_faces`` -- (default: ``False``) a boolean or a list of face
      indices; the punctured faces, ``True`` standing for every face and
      ``False`` for none

    - ``check`` -- boolean (default: ``True``); whether to check that ``m`` is
      a connected hyperbolic quadrangulation with boundary and that the
      surface, once the punctured faces are removed, has negative Euler
      characteristic

    EXAMPLES:

    The quad system of a closed surface of genus two has two vertices of
    degree `8` and four quadrilaterals::

        sage: from combisurf import OrientedMap
        sage: from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation
        sage: m = OrientedMap(vp=[[0, 2, 1, 3, 4, 6, 5, 7]])
        sage: q, bigons = m.quad_system()
        sage: q.vertex_profile(), q.face_profile(), bigons
        ([8, 8], [4, 4, 4, 4], [])
        sage: HyperbolicQuadrangulation(q)
        HyperbolicQuadrangulation(OrientedMap("(0,2,1,3,4,6,5,7)(~0,~2,~1,~7,~4,~6,~5,~3)", "(0,~3,1,~2)(~0,7,~1,2)(3,~5,6,~4)(4,~7,5,~6)"), [])

    The bigons of a quad system with punctures are given as the punctured
    faces::

        sage: q, bigons = m.quad_system(punctured_faces=True)
        sage: bigons
        [1]
        sage: HyperbolicQuadrangulation(q, bigons)
        HyperbolicQuadrangulation(OrientedMap("(0,2,1,3,4,6,5,7,8)(~0,~8,~2,~1,~7,~4,~6,~5,~3)", "(0,~3,1,~2)(~0,8)(~1,2,~8,7)(3,~5,6,~4)(4,~7,5,~6)"), [1])

    TESTS:

    The closed torus is not hyperbolic::

        sage: q, _ = OrientedMap(vp=[[0, 2, 1, 3]]).quad_system()
        sage: HyperbolicQuadrangulation(q)
        Traceback (most recent call last):
        ...
        ValueError: the surface must have negative Euler characteristic once its punctured faces are removed (got 0)

    A face that is not punctured must be a quadrilateral::

        sage: m = OrientedMap(vp=[[0, 2, 1, 3, 4, 6, 5, 7]])
        sage: HyperbolicQuadrangulation(m)
        Traceback (most recent call last):
        ...
        ValueError: the face 0 is not punctured and has degree 8 instead of 4
        sage: HyperbolicQuadrangulation(m, punctured_faces=True)
        HyperbolicQuadrangulation(OrientedMap("(0,1,~0,~1,2,3,~2,~3)", "(0,1,~0,~3,2,3,~2,~1)"), [0])

    A vertex that touches no punctured face must have degree at least five::

        sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
        sage: r, _ = m.radial_map()
        sage: r.vertex_profile()
        [4, 11, 6, 3, 4]
        sage: HyperbolicQuadrangulation(r)
        Traceback (most recent call last):
        ...
        ValueError: the vertex 0 touches no punctured face and has degree 4 instead of at least 5

    Folded edges are not supported::

        sage: HyperbolicQuadrangulation(OrientedMap(vp="(0,1,~1)"))
        Traceback (most recent call last):
        ...
        NotImplementedError: hyperbolic quadrangulations with folded edges are not implemented
    """
    def __init__(self, m, punctured_faces=False, check=True):
        if not isinstance(m, OrientedMap):
            raise TypeError("m must be an oriented map")
        if m.has_folded_edge():
            raise NotImplementedError("hyperbolic quadrangulations with folded edges are not implemented")
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
                if not is_punctured[i] and len(f) != 4:
                    raise ValueError(f"the face {i} is not punctured and has degree {len(f)} instead of 4")
            for i, v in enumerate(m.vertices()):
                if len(v) < 5 and not any(is_punctured[h2f[h]] for h in v):
                    raise ValueError(f"the vertex {i} touches no punctured face and has degree {len(v)} instead of at least 5")

        # Positions around each vertex, counterclockwise from its first
        # half-edge: _upos counts the half-edges and _wpos the weights of the
        # corners. The weighted position p of the vertex v is the entry
        # _offset[v] + p of _at, which holds the half-edge at that position,
        # or -1 for the three positions inside a punctured corner. This turns
        # half-edges into weighted turns and back in constant time.
        vp = m.vertex_permutation(copy=False)
        n = len(vp)
        hv = array('i', [-1] * n)
        upos = array('i', [-1] * n)
        wpos = array('i', [-1] * n)
        deg = array('i')
        wdeg = array('i')
        offset = array('i')
        at = array('i')
        for h0 in range(n):
            if vp[h0] == -1 or hv[h0] != -1:
                continue
            v = len(deg)
            offset.append(len(at))
            d = p = 0
            h = h0
            while True:
                hv[h] = v
                upos[h] = d
                wpos[h] = p
                at.append(h)
                d += 1
                if is_punctured[h2f[h]]:
                    at.extend((-1, -1, -1))
                    p += 4
                else:
                    p += 1
                h = vp[h]
                if h == h0:
                    break
            deg.append(d)
            wdeg.append(p)

        self._m = m
        self._punctured_faces = punctured
        self._hv = hv
        self._upos = upos
        self._wpos = wpos
        self._deg = deg
        self._wdeg = wdeg
        self._offset = offset
        self._at = at

    def __repr__(self):
        return f"HyperbolicQuadrangulation({self._m}, {self._punctured_faces})"

    def map(self):
        r"""
        Return the underlying oriented map.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation
            sage: q, bigons = OrientedMap(vp=[[0, 2, 1, 3]]).quad_system(punctured_faces=True, relabel=True)
            sage: Q = HyperbolicQuadrangulation(q, bigons)
            sage: Q.map()
            OrientedMap("(0,2,1,3,4)(~0,~4,~2,~1,~3)", "(0,~3,1,~2)(~0,4)(~1,2,~4,3)")
            sage: Q.map() == q
            True
        """
        return self._m

    def punctured_faces(self):
        r"""
        Return the sorted list of the indices of the punctured faces.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: q, bigons = m.quad_system(punctured_faces=True)
            sage: HyperbolicQuadrangulation(q, bigons).punctured_faces()
            [1, 5]
            sage: HyperbolicQuadrangulation(q, [5, 1, 5]).punctured_faces()
            [1, 5]
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
        if self._hv[h1 ^ 1] != self._hv[h2]:
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
            sage: from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation
            sage: q, bigons = OrientedMap(vp=[[0, 2, 1, 3]]).quad_system(punctured_faces=True, relabel=True)
            sage: Q = HyperbolicQuadrangulation(q, bigons)
            sage: q.vertices()
            [[0, 4, 2, 6, 8], [1, 9, 5, 3, 7]]
            sage: [Q.turn(0, h) for h in [1, 9, 5, 3, 7]]
            [0, 1, 2, 3, 4]
            sage: Q.turn(0, 2)
            Traceback (most recent call last):
            ...
            ValueError: the half-edge 0 does not end where the half-edge 2 starts
        """
        h1, h2 = self._check_turn(h1, h2)
        g = h1 ^ 1
        return (self._upos[h2] - self._upos[g]) % self._deg[self._hv[g]]

    def _weighted_turn(self, h1, h2):
        g = h1 ^ 1
        wpos = self._wpos
        return _normalize(wpos[h2] - wpos[g], self._wdeg[self._hv[g]])

    def weighted_turn(self, h1, h2):
        r"""
        Return the weighted turn from ``h1`` to ``h2``.

        This is the sum of the weights of the corners swept counterclockwise
        from ``h1 ^ 1`` to ``h2``, an interior corner weighing `1` and a
        punctured corner `4`, taken modulo the weighted degree `wd` of the
        vertex in `(-wd/2, wd/2]`. A positive weighted turn turns right and a
        negative one left.

        When one side sweeps no punctured corner and at most three interior
        corners, the weighted turn is their signed number, positive
        counterclockwise. When both sides sweep a punctured corner, its
        absolute value is at least `4`.

        INPUT:

        - ``h1``, ``h2`` -- half-edges, ``h1`` ending where ``h2`` starts

        OUTPUT: an integer in `(-wd/2, wd/2]`

        EXAMPLES:

        On the quad system of the punctured torus, the second vertex is
        ``vp``-cycle ``(1, 9, 5, 3, 7)``. The corner between ``1`` and
        ``vp[1] = 9`` belongs to the bigon, the face of ``1``, and weighs `4`;
        the other corners of the vertex are interior, so that its weighted
        degree is `8`::

            sage: from combisurf import OrientedMap
            sage: from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation
            sage: q, bigons = OrientedMap(vp=[[0, 2, 1, 3]]).quad_system(punctured_faces=True, relabel=True)
            sage: Q = HyperbolicQuadrangulation(q, bigons)
            sage: q.vertices()
            [[0, 4, 2, 6, 8], [1, 9, 5, 3, 7]]
            sage: q.faces(), bigons
            ([[0, 7, 2, 5], [1, 8], [3, 4, 9, 6]], [1])

        From ``0``, which ends at that vertex, the weighted turns to its
        half-edges are positive to the right, negative to the left; ``7``
        follows ``0`` in the interior face on the left of ``0``, so that
        going from ``0`` to ``7`` is a sharp left turn::

            sage: [Q.weighted_turn(0, h) for h in [1, 9, 5, 3, 7]]
            [0, 4, -3, -2, -1]
            sage: [Q.turn(0, h) for h in [1, 9, 5, 3, 7]]
            [0, 1, 2, 3, 4]

        The boundary of an interior face, which lies on its left, turns
        sharp left at each corner, and the reverse walk sharp right::

            sage: w = q.faces()[0]
            sage: Q.weighted_turns(w)
            [-1, -1, -1, -1]
            sage: Q.weighted_turns([h ^^ 1 for h in reversed(w)])
            [1, 1, 1, 1]

        With two bigons at a vertex, a turn can sweep a punctured corner on
        each side. Here the weighted degree of the vertex ``(0, 4, 8, 2, 6,
        10)`` is `12`, the corners at ``8`` and ``10`` being punctured::

            sage: q, bigons = OrientedMap(vp=[[0, 2, 4, 1, 3, 5]]).quad_system(punctured_faces=True, relabel=True)
            sage: Q = HyperbolicQuadrangulation(q, bigons)
            sage: q.vertices()[0], q.faces(), bigons
            ([0, 4, 8, 2, 6, 10], [[0, 7, 2, 5], [1, 10], [3, 8], [4, 11, 6, 9]], [1, 2])
            sage: [Q.weighted_turn(1, h) for h in [0, 4, 8, 2, 6, 10]]
            [0, 1, 2, 6, -5, -4]
        """
        h1, h2 = self._check_turn(h1, h2)
        return self._weighted_turn(h1, h2)

    def weighted_turns(self, w):
        r"""
        Return the list of the weighted turns of the closed walk ``w``.

        The entry ``i`` is the weighted turn from ``w[i]`` to ``w[i + 1]``,
        the indices being taken modulo the length of ``w``.

        INPUT:

        - ``w`` -- a closed walk

        OUTPUT: a list of integers, of the length of ``w``

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation
            sage: q, bigons = OrientedMap(vp=[[0, 2, 1, 3]]).quad_system(punctured_faces=True, relabel=True)
            sage: Q = HyperbolicQuadrangulation(q, bigons)
            sage: Q.weighted_turns([1, 8])
            [4, 4]
            sage: Q.weighted_turns([0, 1])
            [0, 0]
            sage: Q.weighted_turns([0, 2])
            Traceback (most recent call last):
            ...
            ValueError: invalid walk [0, 2]: 0 at position 0 ends at vertex 1 but 2 at position 1 starts at vertex 0
        """
        w = self._m._check_walk(w, closed=True)
        n = len(w)
        return [self._weighted_turn(w[i], w[i + 1 if i + 1 < n else 0]) for i in range(n)]

    def is_geodesic(self, w):
        r"""
        Return whether the closed walk ``w`` has no spur and no bracket.

        The weighted turns are read cyclically: ``w`` is a geodesic when none
        of them is `0` and no consecutive weighted turns read
        `1, 2, \ldots, 2, 1` or `-1, -2, \ldots, -2, -1`, with any number of
        `\pm 2`, the run of `\pm 2` possibly going once around the whole walk
        (a *cyclic bracket*). The empty walk is a geodesic.

        INPUT:

        - ``w`` -- a closed walk

        OUTPUT: a boolean

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation
            sage: q, bigons = OrientedMap(vp=[[0, 2, 1, 3]]).quad_system(punctured_faces=True, relabel=True)
            sage: Q = HyperbolicQuadrangulation(q, bigons)

        The boundary of the bigon, the boundary of an interior face and a
        backtrack::

            sage: Q.is_geodesic([1, 8])
            True
            sage: Q.is_geodesic([0, 7, 2, 5])
            False
            sage: Q.is_geodesic([0, 1])
            False
            sage: Q.is_geodesic([])
            True
        """
        t = self.weighted_turns(w)
        n = len(t)
        if 0 in t:
            return False
        for s in (1, -1):
            # last index of a turn s; a bracket is two turns s, possibly the
            # same one when it is the only turn s, with only turns 2 s between
            i0 = -1
            for i in range(n):
                if t[i] == s:
                    i0 = i
            if i0 == -1:
                continue
            straight = True
            for k in range(1, n + 1):
                x = t[(i0 + k) % n]
                if x == s:
                    if straight:
                        return False
                    straight = True
                elif x != 2 * s:
                    straight = False
        return True

    def geodesic(self, w, check=True):
        r"""
        Return a geodesic freely homotopic to the closed walk ``w``.

        The geodesic has no spur and no bracket (see :meth:`is_geodesic`). It
        is obtained by the reduction of :ref:`erickson-whittlesey-2013`, Lemma
        4.3, in time linear in the length of ``w``. Any cyclic rotation of the
        result is a valid answer; the geodesics freely homotopic to ``w`` all
        have the same length, :ref:`despre-lazarus-2019`, Corollary 11.

        INPUT:

        - ``w`` -- a closed walk

        - ``check`` -- boolean (default: ``True``); whether to check that
          ``w`` is a closed walk

        OUTPUT: a list of half-edges, empty when ``w`` is contractible

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.hyperbolic_quadrangulation import HyperbolicQuadrangulation
            sage: m = OrientedMap(vp=[[0, 2, 1, 3, 4, 6, 5, 7]])
            sage: q, mor, bigons = m.quad_system(mapping=True, relabel=True)
            sage: Q = HyperbolicQuadrangulation(q, bigons)

        The image of the relator of the genus two surface is contractible::

            sage: relator = [3, 0, 2, 1, 7, 4, 6, 5]
            sage: w = [h for x in relator for h in mor[x]]
            sage: len(w)
            16
            sage: Q.geodesic(w)
            []

        The image of a generator is a geodesic; the image of a conjugate of it
        is not, and its geodesic has the same length::

            sage: Q.geodesic(mor[0])
            [0, 5]
            sage: w = [h for x in [2, 0, 3] for h in mor[x]]
            sage: w
            [4, 3, 0, 5, 2, 5]
            sage: Q.is_geodesic(w)
            False
            sage: g = Q.geodesic(w)
            sage: g
            [3, 14]
            sage: Q.is_geodesic(g)
            True

        The boundary of an interior face is contractible and the boundary of
        a bigon is not::

            sage: q, bigons = OrientedMap(vp=[[0, 2, 1, 3]]).quad_system(punctured_faces=True, relabel=True)
            sage: Q = HyperbolicQuadrangulation(q, bigons)
            sage: q.faces(), bigons
            ([[0, 7, 2, 5], [1, 8], [3, 4, 9, 6]], [1])
            sage: Q.geodesic([0, 7, 2, 5])
            []
            sage: Q.geodesic([1, 8])
            [1, 8]

        The walk around a vertex goes, for each corner, along the sides of its
        face that do not touch the vertex. It is contractible when the vertex
        has no punctured corner. Here the vertex has one punctured corner,
        the one of the bigon: the walk goes around the vertex and the puncture
        and is homotopic to the boundary of the bigon::

            sage: vp, fp = q.vertex_permutation(), q.face_permutation()
            sage: walk = []
            sage: for x in q.vertices()[0]:
            ....:     y = fp[x]
            ....:     while y != vp[x] ^^ 1:
            ....:         walk.append(y)
            ....:         y = fp[y]
            sage: walk
            [7, 2, 9, 6, 5, 0, 3, 4]
            sage: Q.geodesic(walk)
            [8, 1]

        TESTS:

        Walks of length one and two, on the genus two surface with one vertex
        and its face punctured::

            sage: m = OrientedMap(vp=[[0, 2, 1, 3, 4, 6, 5, 7]])
            sage: P = HyperbolicQuadrangulation(m, punctured_faces=True)
            sage: P.geodesic([0]), P.geodesic([0, 0]), P.geodesic([0, 1])
            ([0], [0, 0], [])
            sage: P.geodesic([])
            []
            sage: Q.geodesic([0])
            Traceback (most recent call last):
            ...
            ValueError: non-closed walk [0]: 0 at position 0 ends at vertex 1 but 0 at position 0 starts at vertex 0
        """
        if check:
            w = self._m._check_walk(w, closed=True)
        if not len(w):
            return []
        return _WeightedTurnRuns(self, w).geodesic()


class _WeightedTurnRuns:
    r"""
    The weighted turns of a closed walk, run-length encoded in a circular
    doubly linked list, and their reduction to a geodesic.

    A run is a node with a value, the weighted turn, a length, the number of
    consecutive weighted turns equal to the value, and a first half-edge, the
    half-edge ``h`` such that the first weighted turn of the run is the one
    from ``h`` to the next half-edge of the walk. From the first half-edge of
    any run and the runs that follow, the walk is rebuilt in linear time.

    Only the runs of `\pm 2` are merged and can be longer than one: a bracket
    `s, 2s, \ldots, 2s, s` is replaced in constant time only when its turns
    `2s` form a single run, whereas a spur or a bracket without turn `\pm 2`
    costs a constant time anyway and no reduction looks inside a run of
    another value.
    """
    def __init__(self, Q, w):
        self._hv = Q._hv
        self._wpos = Q._wpos
        self._wdeg = Q._wdeg
        self._offset = Q._offset
        self._at = Q._at

        n = len(w)
        val = []
        length = []
        first = []
        for i in range(n):
            h = w[i]
            t = Q._weighted_turn(h, w[i + 1 if i + 1 < n else 0])
            if val and val[-1] == t and (t == 2 or t == -2):
                length[-1] += 1
            else:
                val.append(t)
                length.append(1)
                first.append(h)
        if len(val) > 1 and val[0] == val[-1] and (val[0] == 2 or val[0] == -2):
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
        v = self._hv[h]
        return self._at[self._offset[v] + (self._wpos[h] + t) % self._wdeg[v]]

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
        Merge the equal runs of `\pm 2` among the ``count`` adjacencies that
        follow the node ``start`` and return a live node among them.
        """
        val = self._val
        u = start
        for _ in range(count):
            v = self._next[u]
            if u != v and val[u] == val[v] and (val[u] == 2 or val[u] == -2):
                # u comes first and keeps its first half-edge
                self._len[u] += self._len[v]
                self._remove(v)
            else:
                u = v
        return start if self._alive[start] else u

    def _splice(self, P, middle, N, first_N, runs):
        r"""
        Replace the last weighted turn of ``P``, the runs ``middle`` and the
        first weighted turn of ``N`` by ``runs``, and return the node just
        before the new runs once merged with their neighbours.

        ``P`` and ``N`` may be the same run; ``first_N`` is the first
        half-edge of ``N`` once its first weighted turn is removed. The other
        runs keep their first half-edges, the half-edges before their first
        weighted turns being unchanged.
        """
        for u in middle:
            self._remove(u)
        new = self._insert_after(P, runs)
        self._len[P] -= 1
        self._len[N] -= 1
        self._first[N] = first_N
        if self._len[P] == 0:
            self._remove(P)
        if N != P and self._len[N] == 0:
            self._remove(N)
        self._num_turns -= 2
        return self._merge(self._prev[new], len(runs) + 1)

    def _replace_all(self, nodes, runs, num_turns):
        r"""
        Replace the whole cycle, made of ``nodes``, by ``runs`` and return a
        node of the new cycle.
        """
        for u in nodes:
            self._remove(u)
        a = self._new(*runs[0])
        for value, length, h in runs[1:]:
            self._insert_after(self._prev[a], [(value, length, h)])
        self._num_turns = num_turns
        return self._merge(a, len(runs))

    def _spur(self, z):
        r"""
        Remove the spur at the run ``z``, of value `0`.

        Weighted turns ``x, 0, y`` from ``h_{i-1}``, ``h_i = first[z]``,
        ``h_{i+1} = h_i ^ 1`` and ``h_{i+2}`` become the single weighted turn
        ``x + y`` from ``h_{i-1}`` to ``h_{i+2}``, at the vertex where ``h_i``
        starts. The new run gets the first half-edge ``h_{i-1}``, and the run
        after it loses its first weighted turn and gets the first half-edge
        ``h_{i+2}``; both are read in the position table of that vertex.

        Return ``None`` when the walk is ``(h_i, h_i ^ 1)``, which vanishes.
        """
        if self._num_turns == 2:
            return None
        P = self._prev[z]
        N = self._next[z]
        x = self._val[P]
        y = self._val[N]
        h = self._first[z]
        h_before = self._half_edge_at(h, -x) ^ 1
        h_after = self._half_edge_at(h, y)
        t = _normalize(x + y, self._wdeg[self._hv[h]])
        return self._splice(P, [z], N, h_after, [(t, 1, h_before)])

    def _bracket(self, a, s):
        r"""
        Remove the bracket that starts at the run ``a``, of value ``s = ±1``,
        if there is one, and return the node to resume from; return ``-1``
        when there is no bracket at ``a``.

        The weighted turns ``s, 2s, ..., 2s, s``, with `r` turns `2s`, go
        from ``h_i = first[a]`` along ``h_{i+1}, ..., h_{i+r+1}`` to ``L =
        h_{i+r+2}``, the first half-edge of the run after the bracket, all
        along one side of a strip of `r + 1` interior faces. The other side of
        the strip, ``f_0, ..., f_r``, starts with ``f_0``, at weighted
        position `-s` from ``h_i``, and ends with ``f_r``, whose reverse is at
        weighted position `s` from ``L ^ 1``. In ``x, s, 2s, ..., 2s, s, y``,
        the walk ``h_i, ..., L`` is replaced by ``f_0, ..., f_r`` and the
        weighted turns become ``x - s, -2s, ..., -2s, y - s``: the run of
        ``x - s`` gets the first half-edge ``h_{i-1}``, read at weighted
        position ``-x`` from ``h_i``, the run of `-2s` gets ``f_0``, the run
        of ``y - s`` gets ``f_r``, and the run after loses its first weighted
        turn and gets the first half-edge read at weighted position ``y`` from
        ``L ^ 1``.
        """
        val = self._val
        b = self._next[a]
        if val[b] == 2 * s:
            r = self._len[b]
            c = self._next[b]
            middle = [a, b, c]
        else:
            r = 0
            c = b
            middle = [a, c]
        if val[c] != s:
            return -1
        if c == a:
            if r == 0:
                return -1
            return self._cyclic_bracket(a, b, s, r)

        num_turns = self._num_turns
        if num_turns == r + 2:
            # the whole walk is s, 2s, ..., 2s, s: the two turns s are
            # consecutive and form a bracket that starts at c
            return -1
        h = self._first[a]
        f0 = self._half_edge_at(h, -s)
        if num_turns == r + 3:
            # near-cyclic: x, s, 2s, ..., 2s, s with x the weighted turn from
            # L = h_{i-1} to h_i; the whole walk becomes f_0, ..., f_r with
            # the weighted turns -2s, ..., -2s and x - 2s from f_r to f_0
            X = self._prev[a]
            L = self._first[X]
            fr = self._half_edge_at(L ^ 1, s) ^ 1
            t = _normalize(val[X] - 2 * s, self._wdeg[self._hv[h]])
            runs = [(t, 1, fr)]
            if r:
                runs.append((-2 * s, r, f0))
            return self._replace_all(middle + [X], runs, r + 1)

        P = self._prev[a]
        N = self._next[c]
        x = val[P]
        y = val[N]
        L = self._first[N]
        fr = self._half_edge_at(L ^ 1, s) ^ 1
        h_before = self._half_edge_at(h, -x) ^ 1
        h_after = self._half_edge_at(L ^ 1, y)
        runs = [(_normalize(x - s, self._wdeg[self._hv[h]]), 1, h_before)]
        if r:
            runs.append((-2 * s, r, f0))
        runs.append((_normalize(y - s, self._wdeg[self._hv[L ^ 1]]), 1, fr))
        return self._splice(P, middle, N, h_after, runs)

    def _cyclic_bracket(self, a, b, s, r):
        r"""
        Reduce the walk whose weighted turns are ``s, 2s, ..., 2s`` with `r`
        turns `2s`.

        The walk ``h_0 = first[a], h_1, ..., h_r`` goes around a ring of `r`
        interior faces ``Q_0, ..., Q_{r-1}``, with ``Q_0`` along ``h_0`` and
        ``h_1``, through the weighted turn ``s``, and ``Q_k`` along ``h_{k+1}``
        for `k \geq 1`. The other side of the ring is made of one side
        ``f_k`` of each ``Q_k`` for `k \geq 1` and gives the weighted turns
        ``-2s, ..., -2s, -3s``, the last one at the corner ``w`` of ``Q_0``
        opposite to the walk. The side of ``Q_0`` at weighted position `-s`
        from ``h_0`` goes to ``w``; from its reverse ``e``, ``f_1`` is at
        weighted position `-2s` and the reverse of ``f_{r-1}`` at weighted
        position `s`. The run of `-3s` gets the first half-edge ``f_{r-1}``
        and the run of `-2s` gets ``f_1``.
        """
        if r < 2:
            # with r = 1, the two sides of Q_0 away from the walk would be the
            # two half-edges of one edge, ending at a vertex of degree one
            raise RuntimeError("cyclic bracket with less than two weighted turns ±2; the quadrangulation is not hyperbolic")
        h = self._first[a]
        e = self._half_edge_at(h, -s) ^ 1
        f1 = self._half_edge_at(e, -2 * s)
        f_last = self._half_edge_at(e, s) ^ 1
        runs = [(_normalize(-3 * s, self._wdeg[self._hv[e]]), 1, f_last)]
        if r > 2:
            runs.append((-2 * s, r - 2, f1))
        return self._replace_all([a, b], runs, r - 1)

    def _mark_dirty(self, start):
        r"""
        Mark dirty the nodes around ``start`` that a reduction may have turned
        into the start of a spur or a bracket, and return the first of them.

        A spur or a bracket is looked for from its first run, of value `0` or
        `\pm 1`, and spans at most three runs `s, 2s, s`. A reduction changes
        the runs from ``start``, the node before the new runs, to the node
        after them, at most four nodes further, so that the spurs and
        brackets it creates start between two nodes before ``start`` and
        five nodes after it.
        """
        u = self._prev[self._prev[start]]
        x = u
        for _ in range(min(8, self._size)):
            if not self._dirty[x]:
                self._dirty[x] = True
                self._num_dirty += 1
            x = self._next[x]
        return u

    def geodesic(self):
        r"""
        Reduce the walk and return it.

        This is the scan of :ref:`erickson-whittlesey-2013`, Lemma 4.3. All
        the runs are dirty at first. At a dirty run where no spur nor bracket
        starts, the run is marked clean and the cursor moves on; a reduction
        marks dirty the runs around the new ones and moves the cursor back to
        the first of them. So the dirty runs always follow the cursor, each
        run is marked clean once plus once per nearby reduction, and each
        reduction shortens the walk by two: the total time is linear.
        """
        val = self._val
        dirty = self._dirty
        nxt = self._next
        cur = 0
        while self._num_dirty:
            if not dirty[cur]:
                cur = nxt[cur]
                continue
            t = val[cur]
            if t == 0:
                start = self._spur(cur)
                if start is None:
                    return []
            elif t == 1 or t == -1:
                start = self._bracket(cur, t)
            else:
                start = -1
            if start == -1:
                dirty[cur] = False
                self._num_dirty -= 1
                cur = nxt[cur]
            else:
                cur = self._mark_dirty(start)
        return self._walk(cur)

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
