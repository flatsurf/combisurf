r"""
Geometric intersection of arcs and geodesics on punctured and closed surfaces

See :ref:`despre-lazarus2019` for the computation of geometric intersection
numbers of curves given as words.
"""

from array import array

from combisurf.word import word_init, word_is_cyclically_reduced, word_cyclically_reduce, word_free_group_inverse
from combisurf.oriented_map import OrientedMap
from combisurf.conjugate_tree import ConjugateTree
from combisurf.crossing_arcs import (crossing_arcs_sweep_sorted, cyclically_sorted_leaf_arcs, _leaf_weights,
                                     startpoint_sweep_sorted, startpoint_sweep_weighted, tree_add_with_inverse,
                                     word_arcs)

class GeometricIntersection:
    def __init__(self, m):
        if not isinstance(m, OrientedMap):
            raise ValueError("m must be an oriented map")

        # keeps an immutable copy
        self._cm = m.copy(mutable=False)

        # NOTE: for now we assume that we have a single vertex and that
        # all faces are punctured
        # TODO: implement reduction and buffering through reducing triangulations for closed surfaces
        if self._cm.num_vertices() != 1:
            raise NotImplementedError
        if self._cm.has_folded_edge():
            raise NotImplementedError("geometric intersection is not implemented for maps with folded edges")

        n = len(self._cm._vp)
        # NOTE: an array rather than a list, since the Cython functions it is
        # handed to read it as a C array: converting a list of n = 128 angles
        # costs 1.7 us, against 0.08 us for an array, and a call of
        # geometric_intersection makes three such conversions
        self._angles = array('i', [-1]) * n
        self._angles[0] = 0
        i = 0
        for _ in range(n - 1):
            j = self._cm._fp[i]
            self._angles[j] = self._angles[i] + 1
            i = j

    def __repr__(self):
        return f"GeometricIntersection({self._cm})"

    def __call__(self, w1, w2=None):
        r"""
        Return the self-intersection of ``w1`` or the intersection between
        ``w1`` and ``w2``.
        """
        raise NotImplementedError

    # could return
    # [not a permutation]
    # [periods]
    def conjugate_sort(self, words, check=True):
        r"""
        Given a list of distinct words return their cyclic ordering on the
        boundary at infinity.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersection
            sage: octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")
            sage: gi = GeometricIntersection(octagon)

        A simple example (that turns out to be equivalent to lexicographic sort
        of conjugates)::

            sage: w = [0, 2, 0, 0, 2]
            sage: words, shifts = gi.conjugate_sort([w])
            sage: words
            [0, 0, 0, 0, 0]
            sage: shifts
            [2, 0, 3, 1, 4]
            sage: for k in shifts:
            ....:     print(w[k:] + w[:k])
            [0, 0, 2, 0, 2]
            [0, 2, 0, 0, 2]
            [0, 2, 0, 2, 0]
            [2, 0, 0, 2, 0]
            [2, 0, 2, 0, 0]
        """
        T = ConjugateTree(len(self._angles))
        for i, w in enumerate(words):
            if not isinstance(w, array):
                w = array('i', w)
            if check and not word_is_cyclically_reduced(w):
                raise ValueError
            ans = T.process(list(w))
            if ans <= 0:
                # was already given in T
                raise ValueError(f"conjugate words at position {-ans} and {i}")

        # NOTE: below a node, the angles are measured from the half-edge
        # through which the curve came in, the inverse of the last letter read
        n = len(self._angles)
        pivot = array('i', [self._angles[b ^ 1] for b in range(n)])
        word_indices, word_shifts = T.sorted_leaves_as_conjugates(self._angles, pivot)

        return list(word_indices), list(word_shifts)

    def conjugate_plot(self, words):
        from sage.rings.complex_double import CDF
        from sage.plot.colors import rainbow
        from sage.plot.text import text
        from sage.plot.point import point2d
        from sage.plot.circle import circle
        from sage.plot.line import line2d

        n = len(words)
        words = [array('i', w) for w in words]
        words_with_inverse = list(words) + [word_free_group_inverse(w) for w in words]
        l = sum(len(w) for w in words_with_inverse)
        # print(f"words_with_inverse={words_with_inverse}")
        colors = rainbow(n, 'rgbtuple')
        word_indices, word_shifts = self.conjugate_sort(words_with_inverse)
        assert len(word_indices) == len(word_shifts) == l
        # print(f"word_indices={word_indices} word_shifts={word_shifts}")
        conj_to_pos = [[None] * len(w) for w in words_with_inverse]
        for pos, (i, k) in enumerate(zip(word_indices, word_shifts)):
            conj_to_pos[i][k] = pos

        G = circle((0,0),1,color='black')
        z = CDF.zeta(l)
        for i, w, positions, color in zip(range(2 * n), words_with_inverse, conj_to_pos, colors * 2):
            G += point2d([z ** pos for pos in positions], color=color, pointsize=50)
            for k, pos in enumerate(positions):
                zz = z**pos
                G += text(''.join(map(str, w[k:] + w[:k])), (1.2*zz.real(), 1.2*zz.imag()), rotation=360. * pos / l, color=color)

        for i, w in enumerate(words):
            for k in range(len(w)):
                endpoint = conj_to_pos[i][k]
                startpoint = conj_to_pos[n+i][-k]
                G += line2d([z**startpoint, z**endpoint], color=colors[i])
        G.set_aspect_ratio(1)
        G.axes(False)
        return G

    def geometric_intersection(self, ulist, vlist=None, check=True):
        r"""
        Return the geometric intersection between the multicurves ``ulist``
        and ``vlist`` given as a list of walks on the half-edges of the
        underlying map.

        If ``vlist`` is not provided, return the self-intersection of
        ``ulist``.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.word import word_init
            sage: from combisurf.geometric_intersection import GeometricIntersection

            sage: torus = OrientedMap(fp="(0,1,~0,~1)")
            sage: gi = GeometricIntersection(torus)

            sage: gi.geometric_intersection([[0]], [[2]])
            1
            sage: gi.geometric_intersection([[0, 0, 2, 2]])
            1

            sage: gi.geometric_intersection([[0]], [[0, 2]])
            1
            sage: gi.geometric_intersection([[0, 2]], [[0, 2, 0]])
            1
            sage: gi.geometric_intersection([[0, 2, 0]], [[0, 2, 0, 0, 2]])
            1

        Two examples in genus 2 following :ref:`birman-series1984`, pages
        336-337::

            sage: octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")
            sage: gi = GeometricIntersection(octagon)
            sage: w = word_init("0,~1,3")
            sage: gi.geometric_intersection([w])
            0
            sage: w = word_init("0,1,1,~2,1,1,~2")
            sage: gi.geometric_intersection([w])
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
            ....:         ans2 = gi.geometric_intersection([list(w)]) == 0
            ....:         assert ans1 == ans2

            sage: ulist = [[0], [0, 2], [0, 0, 2]]
            sage: vlist = [[0, 2, 2, 0, 2], [2]]
            sage: gi.geometric_intersection(ulist, vlist)
            12
            sage: gi.geometric_intersection([[0, 0, 2, 2]])
            1

            sage: gi.geometric_intersection([[0]], [[0, 2]])
            1
            sage: gi.geometric_intersection([[0, 2]], [[0, 2, 0]])
            1
            sage: gi.geometric_intersection([[0, 2, 0]], [[0, 2, 0, 0, 2]])
            1

            sage: for u in [[0], [0, 2], [0, 2, 0], [0, 2, 0, 0, 2],  [0, 2, 0, 0, 2, 0, 2, 0]]:
            ....:     assert gi.geometric_intersection([u]) == 0
            sage: for u in [[0, 0, 2, 2], [0, 2, 0, 2, 0, 0], [0, 2, 0, 0, 2, 0, 0, 2, 0, 2],
            ....:           [0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0],
            ....:           [0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0, 2]]:
            ....:     assert gi.geometric_intersection([u]) == 1

        Two examples in genus 2 following :ref:`birman-series1984`, pages
        336-337::

            sage: octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")
            sage: gi = GeometricIntersection(octagon)
            sage: w = word_init("0,~1,3")
            sage: gi.geometric_intersection([w])
            0
            sage: w = word_init("0,1,1,~2,1,1,~2")
            sage: gi.geometric_intersection([w])
            4

        Testing the simplicity criterion of :ref:`lapointe2019` on positive
        words::

            sage: W = Words([0, 2, 4, 6])
            sage: for l in range(3, 6):
            ....:     for w in W.iterate_by_length(l):
            ....:         if not w.is_primitive():
            ....:             continue
            ....:         bwt = w.BWT()
            ....:         ans1 = all(bwt[i + 1] <= bwt[i] for i in range(l - 1))
            ....:         ans2 = gi.geometric_intersection([list(w)]) == 0
            ....:         assert ans1 == ans2

        Intersection is multilinear::

            sage: torus = OrientedMap(fp="(0,1,~0,~1)")
            sage: gi = GeometricIntersection(torus)
            sage: ulist = [[0], [0, 2], [0, 0, 2]]
            sage: vlist = [[0, 2, 2, 0, 2], [2]]
            sage: gi.geometric_intersection(ulist, vlist)
            12
            sage: gi.geometric_intersection(ulist * 2, vlist)
            24
            sage: gi.geometric_intersection(ulist, vlist * 3)
            36
            sage: gi.geometric_intersection(ulist * 5, vlist * 3)
            180

        For intersection of two multicurves, non-primitivity plays the same role as multiplicity::

            sage: u0 = [0, 0, 2, 0, 3]
            sage: u1 = [0, 0, 2, 2, 1, 1, 3, 3]
            sage: gi.geometric_intersection([u0, u0, u1], [u1, u1, u1])
            108
            sage: gi.geometric_intersection([u0 * 2, u1], [u1, u1 * 2])
            108
            sage: gi.geometric_intersection([u0 * 2, u1], [u1 * 3])
            108

        For self-intersection, non-primitivty adds a factor equal to the exponent minus one::

            sage: u = [0, 0, 2, 2]
            sage: gi.geometric_intersection([u])
            1
            sage: gi.geometric_intersection([u * 2])
            5
            sage: gi.geometric_intersection([u * 3])
            11

        TESTS:

        A map with a folded edge is rejected at construction, before any
        curve is even looked at::

            sage: m = OrientedMap(fp="(0,1,~0,~1,2)")
            sage: GeometricIntersection(m)
            Traceback (most recent call last):
            ...
            NotImplementedError: geometric intersection is not implemented for maps with folded edges
        """
        # For general multicurves where u and v might have common components, each primitive
        # word (and hence each arc) has an associated u-multiplicity and v-multiplicity.
        intersections = 0  # result
        n = len(self._cm._vp)
        T = ConjugateTree(n)
        u_multiplicities = []
        v_multiplicities = []
        # NOTE: tree_add_with_inverse adds a new word together with its inverse,
        # the pair getting the indices (2 * slot, 2 * slot + 1) where slot is
        # the next free one
        for u in ulist:
            if check:
                u = word_cyclically_reduce(word_init(u))
            if not u:
                continue
            i, exponent = tree_add_with_inverse(T, u)
            if (i >> 1) == len(u_multiplicities):
                # u added to T
                u_multiplicities.append(0)
                v_multiplicities.append(0)
            u_multiplicities[i >> 1] += exponent
            if vlist is None:
                # NOTE: non-primitive contribution to self-intersection
                intersections += 2 * (exponent - 1)

        if vlist is not None:
            self_intersection = False
            for v in vlist:
                if check:
                    v = word_cyclically_reduce(word_init(v))
                if not v:
                    continue
                i, exponent = tree_add_with_inverse(T, v)
                if (i >> 1) == len(u_multiplicities):
                    # v added to T
                    u_multiplicities.append(0)
                    v_multiplicities.append(0)
                v_multiplicities[i >> 1] += exponent
        else:
            self_intersection = True
            v_multiplicities = u_multiplicities

        # print(f"u_multiplicities={u_multiplicities} v_multiplicities={v_multiplicities}")
        # NOTE: the words are read many times below, so they are taken out
        # of the tree once rather than through an accessor at every letter
        words = T.words()

        # Essential intersection coming from pairs of arcs with four distinct
        # endpoints. Arcs sharing both endpoints are merged, their weights
        # added up.
        angles = self._angles
        ukeys, ukey_weights = word_arcs(n, angles, words, u_multiplicities)
        # NOTE: the sweep is O((len(u) + len(v)) log(n)). The O(n^2) double
        # sum of the same number (test/test_geometric_intersection.py) is
        # slower at every n, so there is no threshold: on the one-vertex
        # 4g-gon with two random curves of length 8, the sweep takes 0.3 us
        # against 5.2 us at n = 4 and 1.0 us against 8000 us at n = 256; when the arcs fill the n^2 / 2 possible pairs (n = 64,
        # curves of length 4000) it takes 260 us against 620 us.
        if self_intersection:
            intersections += crossing_arcs_sweep_sorted(n, ukeys, ukey_weights, check=False)
        else:
            vkeys, vkey_weights = word_arcs(n, angles, words, v_multiplicities)
            intersections += crossing_arcs_sweep_sorted(n, ukeys, ukey_weights, vkeys, vkey_weights, check=False)
            intersections *= 2

        # Essential intersections coming from pairs of conjugates with identical
        # start, each leaf being described by its startpoint, the angle from
        # its startpoint to its endpoint and its two multiplicities. Total cost
        # is (len(u) + len(v)) * log(n) where the log(n) factor comes from
        # partial sums.
        word_index, starts, arc_angles = cyclically_sorted_leaf_arcs(T, angles)
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

    def intersection_matrix(self, curves, check=True):
        r"""
        Return the matrix of geometric intersection numbers of the primitive
        curves ``curves``.

        This is a :class:`GeometricIntersectionMatrix` sharing the angle table
        of this object. It is much faster than calling
        :meth:`geometric_intersection` on each pair, at the cost of fixing the
        list of curves once and for all.

        INPUT:

        - ``curves`` -- a list of walks on the half-edges of the underlying map

        - ``check`` -- boolean (default: ``True``); whether to cyclically
          reduce the curves in input

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersection
            sage: torus = OrientedMap(fp="(0,1,~0,~1)")
            sage: gi = GeometricIntersection(torus)
            sage: gi.intersection_matrix([[0], [2], [0, 2]]).matrix()
            [0 1 1]
            [1 0 1]
            [1 1 0]
        """
        return GeometricIntersectionMatrix(self, curves, check=check)


class GeometricIntersectionMatrix:
    r"""
    The geometric intersection numbers of a fixed list of curves.

    The list of curves is given once and for all at construction time. All the
    conjugates of all the curves (and of their inverses) are stored in a single
    :class:`~combisurf.conjugate_tree.ConjugateTree` and ordered on the boundary
    at infinity once. An intersection number is then a merge of two rank-sorted
    lists rather than a fresh tree, which is what makes this much faster than
    calling :meth:`GeometricIntersection.geometric_intersection` on each pair.
    Each curve keeps data of size proportional to its length, and an entry
    costs `O((|u| + |v|) \log(n))` where `n` is the number of half-edges.

    Asking for a curve outside of ``curves`` means building another object.

    INPUT:

    - ``m`` -- an :class:`~combisurf.oriented_map.OrientedMap` or a
      :class:`GeometricIntersection` built on it

    - ``curves`` -- a list of walks on the half-edges of ``m``; each of them
      must be primitive

    - ``check`` -- boolean (default: ``True``); whether to cyclically reduce
      the curves in input

    EXAMPLES::

        sage: from combisurf import OrientedMap
        sage: from combisurf.geometric_intersection import GeometricIntersection, GeometricIntersectionMatrix

        sage: torus = OrientedMap(fp="(0,1,~0,~1)")
        sage: I = GeometricIntersectionMatrix(torus, [[0], [2], [0, 2], [0, 2, 2]])
        sage: I
        GeometricIntersectionMatrix of 4 curves on OrientedMap("(0,1,~0,~1)", "(0,1,~0,~1)")
        sage: I.matrix()
        [0 1 1 2]
        [1 0 1 1]
        [1 1 0 1]
        [2 1 1 0]

    The same object is reachable from an existing
    :class:`GeometricIntersection`, in which case the angle table is shared::

        sage: gi = GeometricIntersection(torus)
        sage: gi.intersection_matrix([[0], [2]]).matrix()
        [0 1]
        [1 0]

    An entry is the intersection of the two corresponding one-element
    multicurves, so the diagonal is ``i(c, c) = 2 i(c)`` rather than the
    self-intersection ``i(c)``::

        sage: c = [0, 0, 2, 0, 3]
        sage: I = GeometricIntersectionMatrix(torus, [c])
        sage: I.entry(0, 0)
        4
        sage: gi.geometric_intersection([c], [c])
        4
        sage: gi.geometric_intersection([c])
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
    :class:`GeometricIntersection`::

        sage: m = OrientedMap(fp="(0,1,~0,~1,2)")
        sage: GeometricIntersectionMatrix(m, [[0]])
        Traceback (most recent call last):
        ...
        NotImplementedError: geometric intersection is not implemented for maps with folded edges
    """
    def __init__(self, m, curves, check=True):
        if isinstance(m, GeometricIntersection):
            gi = m
        elif isinstance(m, OrientedMap):
            gi = GeometricIntersection(m)
        else:
            raise ValueError("m must be an oriented map or a GeometricIntersection")

        self._gi = gi
        angles = gi._angles
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
            i, exponent = tree_add_with_inverse(T, w)
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
        leaf_word, leaf_start, leaf_angle = cyclically_sorted_leaf_arcs(T, angles)
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
            keys, weights = word_arcs(n, angles, [words[2 * slot]], [1])
            self._arc_keys.append(keys)
            self._arc_weights.append(weights)

        # Scratch space for the two sweeps of an entry, allocated once.
        # NOTE: the sweeps leave it filled with zeros, which saves an
        # allocation per entry: 0.32 us against 0.42 us per call of
        # crossing_arcs_sweep_sorted at n = 1000 with curves of length 8.
        self._arc_scratch = array('q', [0]) * (2 * (n + 1))

    def __repr__(self):
        return f"GeometricIntersectionMatrix of {len(self._curves)} curves on {self._gi._cm}"

    def __len__(self):
        return len(self._curves)

    def curves(self):
        r"""
        Return the list of curves of this matrix, cyclically reduced.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.geometric_intersection import GeometricIntersectionMatrix
            sage: torus = OrientedMap(fp="(0,1,~0,~1)")
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
            sage: octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")
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
            sage: octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")
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
            sage: from combisurf.geometric_intersection import GeometricIntersection, GeometricIntersectionMatrix
            sage: from combisurf.word import word_init
            sage: octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")
            sage: gi = GeometricIntersection(octagon)
            sage: curves = [word_init("0"), word_init("~1"), word_init("0,~1,3"),
            ....:           word_init("0,1,1,~2,1,1,~2")]
            sage: I = gi.intersection_matrix(curves)
            sage: [I.entry(3, y) for y in range(4)]
            [2, 3, 2, 8]
            sage: [gi.geometric_intersection([curves[3]], [v]) for v in curves]
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
            sage: torus = OrientedMap(fp="(0,1,~0,~1)")
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
            sage: torus = OrientedMap(fp="(0,1,~0,~1)")
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
