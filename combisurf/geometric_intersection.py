r"""
Geometric intersection numbers of closed curves on punctured surfaces

The surface is given by an :class:`~combisurf.oriented_map.OrientedMap` ``m``
without folded edges, together with a non-empty list of its faces that are
punctures (or boundary components); the other faces are filled in. A closed
curve is a closed walk on ``m`` (see :mod:`~combisurf.oriented_map`).

The computation happens on a reduced map (see
:meth:`~combisurf.oriented_map.OrientedMap.reduced_map`) that has a single
vertex and whose faces are exactly the punctures. It is obtained from ``m`` by
contracting a spanning forest and deleting the edges of a coforest in which
the unpunctured faces get merged into the punctured ones. The fundamental group of
the surface is then free on the edges of the reduced map, and a closed curve
becomes a cyclically reduced word in its half-edges.

The algorithm uses the description from :ref:`cohen-lustig-1987` that
reduces the problem of geometric intersection to word combinatorics. To
make such procedure efficient, some specific data structures and algorithms
should be used:

- a little extension of the standard *suffix tree* that is called *conjugate
  tree* (see :mod:`~combisurf.conjugate_tree`). A conjugate tree is used to
  encode the ordering of endpoints of the lifts of the curve to the universal
  cover. A standard reference for suffix trees is :ref:`ukkonen-1995`.

- the Fenwick sweep (:ref:`fenwick-1994`) that is used to compute prefix sums
  and allows to compute the geometric intersections from the conjugate ordering
  obtained from the conjugate tree.

Different approaches with precise time complexities for computing intersections
were designed for closed surfaces by V. Despre and F. Lazarus
:ref:`despre-lazarus-2019` (quadratic in the length) and more recently
:ref:`dubois-2024` (linear in the length and linear in the map size). The functions
in this module make the geometric intersection computation linear in the length
and logarithmic in the map size.
"""

from array import array
from collections.abc import Sequence

from combisurf.word import word_init, word_cyclically_reduce, word_cancel_ends, word_apply_morphism
from combisurf.permutation import perm_dense_cycle_positions
from combisurf.oriented_map import OrientedMap
from combisurf.conjugate_tree import ConjugateTree
from combisurf.crossing_arcs import (crossing_arcs_sweep_sorted, _cyclically_sorted_leaf_arcs, _leaf_weights,
                                     startpoint_sweep_sorted, startpoint_sweep_weighted, _tree_add_with_inverse,
                                     _word_arcs)


class GeometricIntersectionPairing:
    r"""
    Geometric intersection pairing on a surface.

    The geometric intersection pairing associates to each pair of closed curves
    on a surface a non-negative integer. This function is actually defined on
    free homotopy classes of curves. Curves here are encoded as walks on an
    :class:`~combisurf.oriented_map.OrientedMap`, that is a sequence of
    half-edges.

    INPUT:

    - ``m`` -- an :class:`~combisurf.oriented_map.OrientedMap`

    - ``punctured_faces`` -- either a boolean (default ``False``) or a list of
      face indices. If set to ``False`` or ``True`` then respectively no face or
      all faces are considered as punctured. If a list of integers is provided,
      then the faces corresponding to them are considered as punctured.

    EXAMPLES:

    The pair of pants is the map with one vertex and three faces. Its
    half-edges ``0`` and ``2`` are the generators `a` and `b` of its
    fundamental group, and ``1 = ~0`` and ``3 = ~2`` their inverses::

        sage: from combisurf import OrientedMap
        sage: from combisurf.geometric_intersection import GeometricIntersectionPairing
        sage: from combisurf.word import word_init, word_cyclically_reduce
        sage: P = OrientedMap(fp="(0)(1)(~0,~1)")
        sage: P
        OrientedMap("(0,~0,1,~1)", "(0)(~0,~1)(1)")
        sage: P.num_vertices(), P.num_faces()
        (1, 3)
        sage: gi = GeometricIntersectionPairing(P, punctured_faces=True)

    The three boundaries, read along the faces, are disjoint simple curves::

        sage: boundaries = [[0], [2], [1, 3]]
        sage: gi.matrix(boundaries).matrix()
        [0 0 0]
        [0 0 0]
        [0 0 0]

    The figure eight `a b^{-1}` has one self-intersection::

        sage: gi([[0, 3]])
        1
        sage: gi([[0, 2, 0, 0, 2]])
        2

    The exchange of `a` and `b`, and the rotation of order three that sends
    `a` to `b` and `b` to `a^{-1} b^{-1}`, are mapping classes. A letter ``h``
    is sent to the word ``mor[h]``. The rotation permutes the boundaries, and
    both preserve the intersection numbers::

        sage: mor0 = [[2], [3], [0], [1]]
        sage: rot = [[2], [3], [1, 3], [2, 0]]
        sage: def image(mor, w):
        ....:     return list(word_cyclically_reduce(word_init(sum((mor[h] for h in w), []))))
        sage: [image(rot, b) for b in boundaries]
        [[2], [1, 3], [0]]
        sage: curves = [[0, 3], [0, 0, 3], [0, 2, 0, 0, 2]]
        sage: M = gi.matrix(curves).matrix()
        sage: M
        [2 2 2]
        [2 4 2]
        [2 2 4]
        sage: gi.matrix([image(mor0, w) for w in curves]).matrix() == M
        True
        sage: gi.matrix([image(rot, w) for w in curves]).matrix() == M
        True

    Works with multiple vertices::

        sage: m = OrientedMap(vp="(0,1)(~0,~1)")
        sage: m.num_vertices(), m.has_folded_edge()
        (2, False)
        sage: GeometricIntersectionPairing(m, punctured_faces=True)
        GeometricIntersectionPairing(OrientedMap("(0,1)(~0,~1)", "(0,~1)(~0,1)"), [0, 1])

    But at least one face has to be punctured::

        sage: m = OrientedMap(vp="(0,1)(~0,~1)")
        sage: m.num_vertices(), m.has_folded_edge()
        (2, False)
        sage: GeometricIntersectionPairing(m, punctured_faces=False)
        Traceback (most recent call last):
        ...
        NotImplementedError: geometric intersection is not implemented on closed surfaces

    On a sphere with four faces, a figure eight around boundaries 0 and 1 is a
    boundary component if either 0 or 1 is not punctured and has
    self-intersection 1 otherwise::

        sage: from combisurf.word import word_init
        sage: m = OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~4,~5)(~2,5,~3)")
        sage: gi = GeometricIntersectionPairing(m, punctured_faces=[0,1,2,3])
        sage: gi0 = GeometricIntersectionPairing(m, punctured_faces=[1,2,3])
        sage: gi1 = GeometricIntersectionPairing(m, punctured_faces=[0,2,3])
        sage: gi2 = GeometricIntersectionPairing(m, punctured_faces=[0,1,3])
        sage: gi3 = GeometricIntersectionPairing(m, punctured_faces=[0,1,2])
        sage: eight01 = word_init("0,1,2,0,~4,~3")
        sage: gi([eight01])
        1
        sage: gi0([eight01])
        0
        sage: gi1([eight01])
        0
        sage: gi2([eight01])
        1
        sage: gi3([eight01])
        1

    TESTS:

    Repeated face indices count once::

        sage: P = OrientedMap(fp="(0)(1)(~0,~1)")
        sage: GeometricIntersectionPairing(P, punctured_faces=[0, 0, 1])
        GeometricIntersectionPairing(OrientedMap("(0,~0,1,~1)", "(0)(~0,~1)(1)"), [0, 1])
        sage: GeometricIntersectionPairing(P, punctured_faces=[0, 0, 1])([[0, 0, 2, 2]])
        1
        sage: GeometricIntersectionPairing(P, punctured_faces=True)([[0, 0, 2, 2]])
        2

    A map with a folded edge is rejected::

        sage: m = OrientedMap(fp="(0,1,~0,~1,2)")
        sage: GeometricIntersectionPairing(m)
        Traceback (most recent call last):
        ...
        NotImplementedError: geometric intersection is not implemented for maps with folded edges
    """
    def __init__(self, m, punctured_faces=False):
        if not isinstance(m, OrientedMap):
            raise TypeError("m must be an oriented map")
        if not isinstance(punctured_faces, bool) and not isinstance(punctured_faces, Sequence):
            raise TypeError(f"punctured_faces must either be a boolean or a list of face indices (got {punctured_faces})")

        if m.has_folded_edge():
            raise NotImplementedError("geometric intersection is not implemented for maps with folded edges")

        self._cm = m.copy(mutable=False)

        if isinstance(punctured_faces, Sequence):
            punctured_faces = sorted(set(self._cm._check_face_index(f) for f in punctured_faces))
        elif punctured_faces:
            punctured_faces = list(range(self._cm.num_faces()))
        else:
            punctured_faces = []

        self._punctured_faces = punctured_faces

        if not self._punctured_faces:
            raise NotImplementedError("geometric intersection is not implemented on closed surfaces")

        if self._cm.num_vertices() != 1 or len(self._punctured_faces) != self._cm.num_faces():
            f, c, _ = self._cm.forest_coforest_decomposition(root_faces=self._punctured_faces)
            r, mor = self._cm.reduced_map(forest=f, coforest=c, relabel=True, mapping=True, mutable=False)
            self._r = r
            self._mor = mor
        else:
            self._r = self._cm
            self._mor = None

        # NOTE: angles are the absolute angles made by half-edges around the unique
        # vertex of self._r
        self._angles = perm_dense_cycle_positions(self._r._vp)

    def _reduce(self, w, check=False):
        if check:
            w = self._cm._check_walk(w, closed=True)

        if self._mor is None:
            return word_cyclically_reduce(w)

        # NOTE: the images of the letters under self._mor are freely reduced so
        # that the output of word_apply_morphism is freely reduced as well
        return word_cancel_ends(word_apply_morphism(w, self._mor, reduce=True))

    def __repr__(self):
        return f"GeometricIntersectionPairing({self._cm}, {self._punctured_faces})"

    def __call__(self, u, v=None):
        r"""
        Return the geometric intersection between the multiwalks ``u`` and ``v``.

        If ``v`` is not provided, return the self-intersection of ``u``.

        INPUT:

        - ``u`` -- a multiwalk of closed walks on the underlying map (see
          :mod:`~combisurf.oriented_map`)

        - ``v`` -- an optional multiwalk of closed walks on the underlying map

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.word import word_init
            sage: from combisurf.geometric_intersection import GeometricIntersectionPairing

            sage: torus = OrientedMap(vp="(0,1,~0,~1)")
            sage: gi = GeometricIntersectionPairing(torus, punctured_faces=True)

            sage: gi([[0]], [[2]])
            1
            sage: gi([[0, 0, 2, 2]])
            1

            sage: gi([[0]], [[0, 2]])
            1
            sage: gi([[0, 2]], [[0, 2, 0]])
            1
            sage: gi([[0, 2, 0]], [[0, 2, 0, 0, 2]])
            1

        Intersection is multilinear::

            sage: ulist = [[0], [0, 2], [0, 0, 2]]
            sage: vlist = [[0, 2, 2, 0, 2], [2]]
            sage: gi(ulist, vlist)
            12
            sage: gi(ulist * 2, vlist)
            24
            sage: gi(ulist, vlist * 3)
            36
            sage: gi(ulist * 5, vlist * 3)
            180

        For intersection of two multicurves, non-primitivity plays the same role as multiplicity::

            sage: u0 = [0, 0, 2, 0, 3]
            sage: u1 = [0, 0, 2, 2, 1, 1, 3, 3]
            sage: gi([u0, u0, u1], [u1, u1, u1])
            108
            sage: gi([u0 * 2, u1], [u1, u1 * 2])
            108
            sage: gi([u0 * 2, u1], [u1 * 3])
            108

        For self-intersection, the ``k``-th power of a primitive curve ``u`` has
        self-intersection `k^2 i(u, u) + k - 1`::

            sage: u = [0, 0, 2, 2]
            sage: gi([u])
            1
            sage: gi([u * 2])
            5
            sage: gi([u * 3])
            11

        Two examples in genus 2 following :ref:`birman-series1984`, pages
        336-337. Their words list the sides of a fundamental octagon crossed
        by the curve, so they are closed walks on the one-vertex map whose
        vertex sees these sides in cyclic order, hence the ``vp``::

            sage: octagon = OrientedMap(vp="(0,1,2,3,~0,~1,~2,~3)")
            sage: gi = GeometricIntersectionPairing(octagon, punctured_faces=True)
            sage: w = word_init("0,~1,3")
            sage: gi([w])
            0
            sage: w = word_init("0,1,1,~2,1,1,~2")
            sage: gi([w])
            4

        Testing the simplicity criterion of :ref:`lapointe2019` on positive
        words::

            sage: W = Words([0, 2, 4, 6])
            sage: for l in range(2, 7):
            ....:     for w in W.iterate_by_length(l):
            ....:         if not w.is_primitive():
            ....:             continue
            ....:         bwt = w.BWT()
            ....:         ans1 = all(bwt[i + 1] <= bwt[i] for i in range(l - 1))
            ....:         ans2 = gi([list(w)]) == 0
            ....:         assert ans1 == ans2

            sage: for u in [[0], [0, 2], [0, 2, 0], [0, 2, 0, 0, 2],  [0, 2, 0, 0, 2, 0, 2, 0]]:
            ....:     assert gi([u]) == 0
            sage: for u in [[0, 0, 2, 2], [0, 2, 0, 2, 0, 0], [0, 2, 0, 0, 2, 0, 0, 2, 0, 2],
            ....:           [0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0],
            ....:           [0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0, 2]]:
            ....:     assert gi([u]) == 1
        """
        uwords, umult = self._cm._check_multiwalk(u, closed=True)
        if v is None:
            vwords = vmult = None
        else:
            vwords, vmult = self._cm._check_multiwalk(v, closed=True)

        return self._call(uwords, umult, vwords, vmult)

    def _call(self, uwords, umult, vwords=None, vmult=None):
        r"""
        Return the geometric intersection of the multiwalks ``(uwords, umult)``
        and ``(vwords, vmult)``, without checking the input.

        This is the fast path of ``__call__``. The input must be as returned by
        :meth:`~combisurf.oriented_map.OrientedMap._check_multiwalk` with
        ``closed=True`` on the original map ``self._cm``: ``uwords`` and
        ``vwords`` are lists of closed walks given as ``array('i')`` and
        ``umult`` and ``vmult`` are lists of positive integers of the same
        lengths. Nothing is checked.

        TESTS::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersectionPairing
            sage: from combisurf.word import word_init
            sage: octagon = OrientedMap(vp="(0,1,2,3,~0,~1,~2,~3)")
            sage: gi = GeometricIntersectionPairing(octagon, punctured_faces=True)
            sage: uwords = [word_init([0, 2, 5, 1, 2]), word_init([3])]
            sage: umult = [2, 1]
            sage: vwords = [word_init([2]), word_init([3, 7]), word_init([0])]
            sage: vmult = [1, 1, 1]
            sage: gi._call(uwords, umult)
            6
            sage: gi._call(uwords, umult, vwords, vmult)
            14
        """
        # For general multicurves where u and v might have common components, each primitive
        # word (and hence each arc) has an associated u-multiplicity and v-multiplicity.
        intersections = 0  # result
        n = len(self._r._vp)
        T = ConjugateTree(n)
        u_multiplicities = []
        v_multiplicities = []
        # NOTE: _tree_add_with_inverse adds a new word together with its inverse,
        # the pair getting the indices (2 * slot, 2 * slot + 1) where slot is
        # the next free one
        for u, m in zip(uwords, umult):
            u = self._reduce(u)
            if not u:
                continue
            i, exponent = _tree_add_with_inverse(T, u)
            if (i >> 1) == len(u_multiplicities):
                # u added to T
                u_multiplicities.append(0)
                v_multiplicities.append(0)
            u_multiplicities[i >> 1] += m * exponent
            if vwords is None:
                # NOTE: non-primitive contribution to self-intersection
                intersections += m * 2 * (exponent - 1)

        if vwords is not None:
            self_intersection = False
            for v, m in zip(vwords, vmult):
                v = self._reduce(v)
                if not v:
                    continue
                i, exponent = _tree_add_with_inverse(T, v)
                if (i >> 1) == len(u_multiplicities):
                    # v added to T
                    u_multiplicities.append(0)
                    v_multiplicities.append(0)
                v_multiplicities[i >> 1] += m * exponent
        else:
            self_intersection = True
            v_multiplicities = u_multiplicities

        # NOTE: the words are read many times below, so they are taken out
        # of the tree once rather than through an accessor at every letter
        words = T.words()

        # Essential intersection coming from pairs of arcs with four distinct
        # endpoints. Arcs sharing both endpoints are merged, their weights
        # added up.
        angles = self._angles
        ukeys, ukey_weights = _word_arcs(n, angles, words, u_multiplicities)
        # NOTE: the sweep is O(n + (len(u) + len(v)) log(n)). The O(n^2)
        # double sum of the same number (test/test_geometric_intersection.py)
        # is slower at every n, so there is no threshold: on the one-vertex
        # 4g-gon with two random curves of length 8, the sweep takes 0.3 us
        # against 5.2 us at n = 4 and 1.0 us against 8000 us at n = 256; when
        # the arcs fill the n^2 / 2 possible pairs (n = 64, curves of length
        # 4000) it takes 260 us against 620 us.
        if self_intersection:
            intersections += crossing_arcs_sweep_sorted(n, ukeys, ukey_weights, check=False)
        else:
            vkeys, vkey_weights = _word_arcs(n, angles, words, v_multiplicities)
            intersections += crossing_arcs_sweep_sorted(n, ukeys, ukey_weights, vkeys, vkey_weights, check=False)
            intersections *= 2

        # Essential intersections coming from pairs of conjugates with identical
        # start, each leaf being described by its startpoint, the angle from
        # its startpoint to its endpoint and its two multiplicities. Total cost
        # is O(n + (len(u) + len(v)) log(n)), where the n comes from the tables
        # indexed by the angles and the log(n) factor from partial sums.
        word_index, starts, arc_angles = _cyclically_sorted_leaf_arcs(T, angles)
        uweights = _leaf_weights(word_index, u_multiplicities)
        if self_intersection:
            # with the v-weights equal to the u-weights, the sweep counts each
            # pair of leaves in both orders
            intersections += startpoint_sweep_weighted(n, starts, arc_angles, uweights) // 2
        else:
            vweights = _leaf_weights(word_index, v_multiplicities)
            intersections += startpoint_sweep_weighted(n, starts, arc_angles, uweights, vweights)

        # we got twice the geometric intersection because we register all arcs and their inverses
        assert intersections % 2 == 0
        return intersections // 2

    def matrix(self, curves, check=True):
        r"""
        Return a :class:`GeometricIntersectionMatrix` for the primitive curves
        ``curves``.

        The :class:`GeometricIntersectionMatrix` is built on the reduced map
        of this object. It is much faster than calling this object on each
        pair, at the cost of fixing the list of curves once and for all.

        INPUT:

        - ``curves`` -- a list of closed walks on the underlying map

        - ``check`` -- boolean (default: ``True``); whether to check that the
          curves are closed walks on the underlying map. With ``check=False``
          they must already be given as ``array('i')``

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersectionPairing
            sage: torus = OrientedMap(vp="(0,1,~0,~1)")
            sage: gi = GeometricIntersectionPairing(torus, punctured_faces=True)
            sage: gi.matrix([[0], [2], [0, 2]]).matrix()
            [0 1 1]
            [1 0 1]
            [1 1 0]
        """
        curves = [self._reduce(w, check=check) for w in curves]
        return GeometricIntersectionMatrix(self._r, curves, check=False)


class GeometricIntersectionMatrix:
    r"""
    The geometric intersection numbers of a fixed list of curves.

    The list of curves is given once and for all at construction time. All the
    conjugates of all the curves (and of their inverses) are stored in a single
    :class:`~combisurf.conjugate_tree.ConjugateTree` and ordered on the boundary
    at infinity once. An intersection number is then a merge of two rank-sorted
    lists rather than a fresh tree, which is what makes this much faster than
    calling a :class:`GeometricIntersectionPairing` on each pair.
    Each curve keeps data of size proportional to its length, and an entry
    costs `O((|u| + |v|) \log(n))` where `n` is the number of half-edges.

    Asking for a curve outside of ``curves`` means building another object.

    INPUT:

    - ``m`` -- an :class:`~combisurf.oriented_map.OrientedMap` or a
      :class:`GeometricIntersectionPairing` built on it

    - ``curves`` -- a list of closed walks on ``m``; each of them must be
      primitive

    - ``check`` -- boolean (default: ``True``); whether to cyclically reduce
      the curves in input. A curve that is not cyclically reduced raises a
      ``ValueError`` even with ``check=False``

    EXAMPLES::

        sage: from combisurf import OrientedMap
        sage: from combisurf.geometric_intersection import GeometricIntersectionPairing, GeometricIntersectionMatrix

        sage: torus = OrientedMap(vp="(0,1,~0,~1)")
        sage: I = GeometricIntersectionMatrix(torus, [[0], [2], [0, 2], [0, 2, 2]])
        sage: I
        GeometricIntersectionMatrix of 4 curves on OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)")
        sage: I.matrix()
        [0 1 1 2]
        [1 0 1 1]
        [1 1 0 1]
        [2 1 1 0]

    The same object is reachable from an existing
    :class:`GeometricIntersectionPairing`::

        sage: gi = GeometricIntersectionPairing(torus, punctured_faces=True)
        sage: gi.matrix([[0], [2]]).matrix()
        [0 1]
        [1 0]

    An entry is the intersection of the two corresponding one-element
    multicurves, so the diagonal is ``i(c, c) = 2 i(c)`` rather than the
    self-intersection ``i(c)``::

        sage: c = [0, 0, 2, 0, 3]
        sage: I = GeometricIntersectionMatrix(torus, [c])
        sage: I.entry(0, 0)
        4
        sage: gi([c], [c])
        4
        sage: gi([c])
        2

    Curves that are conjugate or inverse to one another share the same internal
    data, which is correct since they have the same intersection numbers with
    everything::

        sage: from combisurf.word import word_init, word_free_group_inverse
        sage: w = word_init([0, 0, 2, 0, 3])
        sage: I = GeometricIntersectionMatrix(torus, [w, w[2:] + w[:2], word_free_group_inverse(w), [0, 2]])
        sage: I.matrix()
        [4 4 4 3]
        [4 4 4 3]
        [4 4 4 3]
        [3 3 3 0]

    Non-primitive curves are not supported. They are detected both when the
    primitive root is already known and when it is not::

        sage: GeometricIntersectionMatrix(torus, [[0, 2], [0, 2, 0, 2]])
        Traceback (most recent call last):
        ...
        NotImplementedError: non-primitive curve at index 1
        sage: GeometricIntersectionMatrix(torus, [[0, 2, 0, 2]])
        Traceback (most recent call last):
        ...
        NotImplementedError: non-primitive curve at index 0

    TESTS:

    A map with a folded edge is rejected already by the underlying
    :class:`GeometricIntersectionPairing`::

        sage: m = OrientedMap(fp="(0,1,~0,~1,2)")
        sage: GeometricIntersectionMatrix(m, [[0]])
        Traceback (most recent call last):
        ...
        NotImplementedError: geometric intersection is not implemented for maps with folded edges
    """
    def __init__(self, m, curves, check=True):
        if check and not isinstance(m, OrientedMap) or m.num_vertices() != 1:
            raise ValueError("m must be an OrientedMap with a single vertex")

        if m.has_folded_edge():
            raise NotImplementedError("geometric intersection is not implemented for maps with folded edges")

        if m.num_vertices() != 1:
            raise ValueError("the map must have a single vertex")

        self._cm = m.copy(mutable=False)
        angles = perm_dense_cycle_positions(self._cm._vp)
        n = len(angles)
        self._n = n

        # 1. all the curves and their inverses in a single tree, a curve and
        # its inverse sitting at the consecutive indices (2 * slot, 2 * slot + 1)
        self._curves = []
        for j, c in enumerate(curves):
            w = word_init(c)
            if check:
                w = word_cyclically_reduce(w)
            if not w:
                raise ValueError(f"trivial curve at index {j}")
            self._curves.append(w)

        # NOTE: the tree holds every curve and its inverse, so it stores twice
        # as many letters as the curves have, and a conjugate tree over T
        # letters has at most 2 T + 1 nodes. Reserving that makes the
        # construction below allocation-free.
        total = sum(len(w) for w in self._curves)
        T = self._tree = ConjugateTree(n, reserve=4 * total + 1)
        self._slot = []
        for j, w in enumerate(self._curves):
            # 2. the slot of the curve, rejecting the non-primitive ones
            i, exponent = _tree_add_with_inverse(T, w)
            if exponent != 1:
                raise NotImplementedError(f"non-primitive curve at index {j}")
            self._slot.append(i >> 1)

        num_slots = T.num_words() // 2
        # NOTE: read once rather than through an accessor at every letter
        words = T.words()

        # 3. the cyclic order at infinity of all the leaves, once. For each
        # slot we keep its own leaves in increasing order of rank, each of them
        # described by its rank, its startpoint and the angle from its
        # startpoint to its endpoint.
        ranks = [array('q') for _ in range(num_slots)]
        starts = [array('q') for _ in range(num_slots)]
        arc_angles = [array('q') for _ in range(num_slots)]
        leaf_word, leaf_start, leaf_angle = _cyclically_sorted_leaf_arcs(T, angles)
        for rank in range(len(leaf_word)):
            slot = leaf_word[rank] >> 1
            ranks[slot].append(rank)
            starts[slot].append(leaf_start[rank])
            arc_angles[slot].append(leaf_angle[rank])
        self._ranks = ranks
        self._starts = starts
        self._arc_angles = arc_angles

        # 4. the arcs of each slot for the crossing arcs term of an entry: the
        # consecutive pairs of letters of the curve, each an arc between two
        # angles first < last, keyed by last * n + first. We keep the distinct
        # keys in increasing order and their multiplicities, in the layout
        # that crossing_arcs_sweep_sorted reads (see _double_sum).
        self._arc_keys = []
        self._arc_weights = []
        for slot in range(num_slots):
            keys, weights = _word_arcs(n, angles, [words[2 * slot]], [1])
            self._arc_keys.append(keys)
            self._arc_weights.append(weights)

        # Scratch space for the two sweeps of an entry, allocated once.
        # NOTE: the sweeps leave it filled with zeros, which saves an
        # allocation per entry: 0.32 us against 0.42 us per call of
        # crossing_arcs_sweep_sorted at n = 1000 with curves of length 8.
        self._arc_scratch = array('q', [0]) * (2 * (n + 1))

    def __repr__(self):
        return f"GeometricIntersectionMatrix of {len(self._curves)} curves on {self._cm}"

    def __len__(self):
        return len(self._curves)

    def curves(self):
        r"""
        Return the list of curves of this matrix, cyclically reduced.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersectionMatrix
            sage: torus = OrientedMap(vp="(0,1,~0,~1)")
            sage: GeometricIntersectionMatrix(torus, [[0, 0, 1, 2], [2]]).curves()
            [array('i', [0, 2]), array('i', [2])]
        """
        return [w[:] for w in self._curves]

    def _double_sum(self, sx, sy):
        r"""
        Return the contribution of the pairs of arcs whose four endpoints are
        pairwise distinct, for the slots ``sx`` and ``sy``.

        This is :func:`~combisurf.crossing_arcs.crossing_arcs_sweep_sorted` on
        the arcs of the two slots built at construction time, with the
        `u`-weights from ``sx`` and the `v`-weights from ``sy``. It costs
        `O((|u| + |v|) \log(n))` where `n` is the number of half-edges.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersectionMatrix
            sage: octagon = OrientedMap(vp="(0,1,2,3,~0,~1,~2,~3)")
            sage: I = GeometricIntersectionMatrix(octagon, [[0, 3, 6], [0, 2, 2, 5, 2, 2, 5]])
            sage: I._double_sum(0, 1), I._double_sum(1, 0), I._double_sum(1, 1)
            (1, 1, 6)
        """
        # NOTE: the O(n^2) double sum (test/test_geometric_intersection.py)
        # computes the same number, from two flat vectors of length
        # (n - 3)(n - 2) / 2 per curve. It is slower at every n, so there is
        # no threshold: with 100 random curves of length 8 on the one-vertex
        # 4g-gon, the sweep takes 0.25 us per pair against 0.46 us at n = 4
        # and 0.87 us against 350 us at n = 128. Evaluating all the double
        # sums at once as a float64 matrix product was faster per pair, but
        # not on the whole matrix(), and it needed O(n^2) memory per curve.
        keys = self._arc_keys
        weights = self._arc_weights
        if sx == sy:
            return crossing_arcs_sweep_sorted(self._n, keys[sx], weights[sx], None, None,
                                              self._arc_scratch, False)
        return crossing_arcs_sweep_sorted(self._n, keys[sx], weights[sx], keys[sy], weights[sy],
                                          self._arc_scratch, False)

    def _entry_from(self, sx, sy, double_sum):
        r"""
        Return the intersection number of the slots ``sx`` and ``sy``, given
        the value ``double_sum`` of ``self._double_sum(sx, sy)``.

        The other term, coming from the pairs of arcs with identical
        startpoint, is :func:`~combisurf.crossing_arcs.startpoint_sweep_sorted`
        on the leaves of the two slots built at construction time.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersectionMatrix
            sage: octagon = OrientedMap(vp="(0,1,2,3,~0,~1,~2,~3)")
            sage: I = GeometricIntersectionMatrix(octagon, [[0, 3, 6], [0, 2, 2, 5, 2, 2, 5]])
            sage: I._entry_from(0, 1, I._double_sum(0, 1)), I._entry_from(1, 1, I._double_sum(1, 1))
            (2, 8)
        """
        if sx == sy:
            # a slot against itself is swept once, merging it with a copy of
            # itself is wrong
            sweep = startpoint_sweep_sorted(self._n, self._ranks[sx], self._starts[sx], self._arc_angles[sx],
                                            None, None, None, self._arc_scratch, False)
        else:
            sweep = startpoint_sweep_sorted(self._n, self._ranks[sx], self._starts[sx], self._arc_angles[sx],
                                            self._ranks[sy], self._starts[sy], self._arc_angles[sy],
                                            self._arc_scratch, False)
        # the two arcs of a crossing are counted once in each direction
        ans = 2 * double_sum + sweep
        assert ans % 2 == 0
        return ans // 2

    def entry(self, x, y):
        r"""
        Return the geometric intersection number of the curves of indices ``x``
        and ``y``.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersectionPairing, GeometricIntersectionMatrix
            sage: from combisurf.word import word_init
            sage: octagon = OrientedMap(vp="(0,1,2,3,~0,~1,~2,~3)")
            sage: gi = GeometricIntersectionPairing(octagon, punctured_faces=True)
            sage: curves = [word_init("0"), word_init("~1"), word_init("0,~1,3"),
            ....:           word_init("0,1,1,~2,1,1,~2")]
            sage: I = GeometricIntersectionMatrix(octagon, curves)
            sage: [I.entry(3, y) for y in range(4)]
            [2, 3, 2, 8]
            sage: [gi([curves[3]], [v]) for v in curves]
            [2, 3, 2, 8]
        """
        sx = self._slot[x]
        sy = self._slot[y]
        return self._entry_from(sx, sy, self._double_sum(sx, sy))

    def row(self, x):
        r"""
        Return the list of the geometric intersection numbers of the curve of
        index ``x`` with all the curves.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersectionMatrix
            sage: torus = OrientedMap(vp="(0,1,~0,~1)")
            sage: I = GeometricIntersectionMatrix(torus, [[0], [2], [0, 2], [0, 2, 2]])
            sage: I.row(1)
            [1, 0, 1, 1]

        It agrees with the entries taken one at a time::

            sage: all(I.row(x) == [I.entry(x, y) for y in range(len(I))] for x in range(len(I)))
            True
        """
        sx = self._slot[x]
        entry_from = self._entry_from
        double_sum = self._double_sum
        return [entry_from(sx, sy, double_sum(sx, sy)) for sy in self._slot]

    def matrix(self):
        r"""
        Return the full symmetric matrix of geometric intersection numbers over
        the integers.

        Only the entries with ``x <= y`` are computed, the others being
        obtained by symmetry.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersectionMatrix
            sage: torus = OrientedMap(vp="(0,1,~0,~1)")
            sage: mat = GeometricIntersectionMatrix(torus, [[0], [2], [0, 2], [0, 2, 2]]).matrix()
            sage: mat
            [0 1 1 2]
            [1 0 1 1]
            [1 1 0 1]
            [2 1 1 0]
            sage: mat.is_symmetric()
            True
            sage: mat.base_ring()
            Integer Ring
        """
        from sage.matrix.constructor import matrix as sage_matrix
        from sage.rings.integer_ring import ZZ

        N = len(self._curves)
        slot = self._slot
        entry_from = self._entry_from
        double_sum = self._double_sum
        rows = [[0] * N for _ in range(N)]
        for x in range(N):
            sx = slot[x]
            for y in range(x, N):
                sy = slot[y]
                e = entry_from(sx, sy, double_sum(sx, sy))
                rows[x][y] = e
                rows[y][x] = e
        return sage_matrix(ZZ, rows)
