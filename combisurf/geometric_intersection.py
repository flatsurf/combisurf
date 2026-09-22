r"""
Geometric intersection of arcs and geodesics on punctured and closed surfaces
"""

from array import array

from combisurf.word import word_init, word_is_cyclically_reduced, word_cyclically_reduce, word_free_group_inverse
from combisurf.oriented_map import OrientedMap
from combisurf.conjugate_tree import ConjugateTree
from combisurf.partial_sums import PartialSums, PartialSumsNaive

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

        n = len(self._cm._vp)
        self._angles = [-1] * n
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
        T = ConjugateTree()
        for i, w in enumerate(words):
            if not isinstance(w, array):
                w = array('i', w)
            if check and not word_is_cyclically_reduced(w):
                raise ValueError
            ans = T.process(list(w))
            if ans <= 0:
                # was already given in T
                raise ValueError(f"conjugate words at position {-ans} and {i}")

        word_indices = []
        word_shifts = []
        for s in T.cyclically_sorted_leaves(self._angles):
            i, k = T.leaf_as_conjugate(s)
            word_indices.append(i)
            word_shifts.append(k)

        return word_indices, word_shifts

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

        Two examples in genus 2 following Birman-Series p336-337::

            sage: octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")
            sage: gi = GeometricIntersection(octagon)
            sage: w = word_init("0,~1,3")
            sage: gi.geometric_intersection([w])
            0
            sage: w = word_init("0,1,1,~2,1,1,~2")
            sage: gi.geometric_intersection([w])
            4

        Testing the simplicity criterion of Lapointe on positive words::

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

        Two examples in genus 2 following Birman-Series p336-337::

            sage: octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")
            sage: gi = GeometricIntersection(octagon)
            sage: w = word_init("0,~1,3")
            sage: gi.geometric_intersection([w])
            0
            sage: w = word_init("0,1,1,~2,1,1,~2")
            sage: gi.geometric_intersection([w])
            4

        Testing the simplicity criterion of Lapointe on positive words::

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
        """
        # For general multicurves where u and v might have common components, each primitive
        # word (and hence each arc) has an associated u-multiplicity and v-multiplicity.
        intersections = 0  # result
        T = ConjugateTree()
        u_multiplicities = []
        v_multiplicities = []
        for u in ulist:
            if check:
                u = word_init(u)
                u = word_cyclically_reduce(u)
            if not u:
                continue
            status = T.process(u)
            if status <= 0:
                # u (or conjugate) already present
                i = -status
                assert len(u) % len(T._words[i]) == 0
                exponent = len(u) // len(T._words[i])
            else:
                # u added to T
                u_multiplicities.append(0)
                if vlist is not None:
                    v_multiplicities.append(0)
                i = len(T._words) - 1
                exponent = status
                ans = T.process(word_free_group_inverse(T._words[i]))
                assert ans == 1
            u_multiplicities[i >> 1] += exponent
            if vlist is None:
                # NOTE: non-primitive contribution to self-intersection
                intersections += 2 * (exponent - 1)

        if vlist is not None:
            self_intersection = False
            for v in vlist:
                if check:
                    v = word_init(v)
                    v = word_cyclically_reduce(v)
                if not v:
                    continue
                status = T.process(v)
                if status <= 0:
                    # v (or conjugate) already present
                    i = -status
                    assert len(v) % len(T._words[i])== 0
                    exponent = len(v) // len(T._words[i])
                else:
                    # v added to T
                    u_multiplicities.append(0)
                    v_multiplicities.append(0)
                    i = len(T._words) - 1
                    exponent = status
                    ans = T.process(word_free_group_inverse(T._words[i]))
                    assert ans == 1
                v_multiplicities[i >> 1] += exponent
        else:
            self_intersection = True
            v_multiplicities = u_multiplicities

        # print(f"u_multiplicities={u_multiplicities} v_multiplicities={v_multiplicities}")
        n = len(self._cm._vp)

        # Essential intersection coming from pairs of conjugates with four
        # distinct 1-order intervals associated to their startpoints and endpoints
        # NOTE: O(n^2 + len(u) + len(v)) cost
        Nu = [[0] * n for _ in range(n)]
        Nv = [[0] * n for _ in range(n)]
        for i in range(0, len(T._words), 2):
            w = T._words[i]
            for p in range(len(w)):
                first = self._angles[w[p]]
                last = self._angles[w[(p - 1) % len(w)] ^ 1]
                assert first != last
                if last < first:
                    first, last = last, first
                Nu[first][last] += u_multiplicities[i >> 1]
                Nv[first][last] += v_multiplicities[i >> 1]

        # NOTE: below is a O(n^2) time version of the two following O(n^4) time sums
        #     sum(Nu[i0][j0] * Nv[i1][j1]
        #         for i0 in range(n)
        #         for j0 in range(i0 + 1, n)
        #         for i1 in range(i0 + 1, j0)
        #         for j1 in range(j0 + 1, n))
        #
        #     sum(Nu[i1][j1] * Nv[i0][j0]
        #         for i0 in range(n)
        #         for j0 in range(i0 + 1, n)
        #         for i1 in range(i0 + 1, j0)
        #         for j1 in range(j0 + 1, n))
        #
        # We optimize the computation of the first sum by transforming Nu and
        # Nv to contain partial sums in respectively i0 and j1 respectively
        # (O(n^2) time).  Then we do a double sum in i1, j0 (O(n^2) time). We
        # reverse the role of Nu and Nv to handle the second sum.
        Nu1 = [l[:] for l in Nu]
        for j in range(n):
            for i in range(j - 1):
                Nu1[i + 1][j] += Nu1[i][j]
        Nv1 = [l[:] for l in Nv]
        for i in range(n):
            for j in range(n - 1, i + 1, -1):
                Nv1[i][j - 1] += Nv1[i][j]

        Nv2 = [l[:] for l in Nv]
        for j in range(n):
            for i in range(j - 1):
                Nv2[i + 1][j] += Nv2[i][j]
        Nu2 = [l[:] for l in Nu]
        for i in range(n):
            for j in range(n - 1, i + 1, -1):
                Nu2[i][j - 1] += Nu2[i][j]

        intersections += sum(Nu1[i1 - 1][j0] * Nv1[i1][j0 + 1] + Nv2[i1 - 1][j0] * Nu2[i1][j0 + 1]
                             for i1 in range(1, n - 2) for j0 in range(i1 + 1, n - 1))
        if not self_intersection:
            intersections *= 2

        # Essential intersections coming from pairs of conjugates with identical
        # start. Total cost is (len(u) + len(v)) * log(n)
        # where the log(n) factor comes from partial sums
        word_indices = []
        word_shifts = []
        for s in T.cyclically_sorted_leaves(self._angles):
            i, k = T.leaf_as_conjugate(s)
            word_indices.append(i)
            word_shifts.append(k)

        letter = 0  # current letter that is looked at
        pos = 0     # pointer in the list cs
        Nu = PartialSums(n - 1)
        Nv = PartialSums(n - 1)
        while pos < len(word_indices):
            Nu.reset()
            Nv.reset()
            startpoint = T._words[word_indices[pos]][word_shifts[pos]]
            startangle = self._angles[startpoint]
            while pos < len(word_indices) and T._words[word_indices[pos]][word_shifts[pos]] == startpoint:
                # print(f"pos={pos} intersections={intersections} Nu={Nu} Nv={Nv}")
                i = word_indices[pos]
                w = T._words[i]
                k = word_shifts[pos]
                endpoint = w[(k - 1) % len(w)] ^ 1
                assert startpoint != endpoint
                endangle = self._angles[endpoint]
                angle = (endangle - startangle) % n
                assert angle >= 1
                angle -= 1
                intersections += u_multiplicities[i >> 1] * Nv.partial_sum(0, angle)
                if not self_intersection:
                    intersections += v_multiplicities[i >> 1] * Nu.partial_sum(0, angle)
                Nu.update(angle, u_multiplicities[i >> 1])
                Nv.update(angle, v_multiplicities[i >> 1])
                pos += 1

        # we got twice the geometric intersection because we register all arcs and their inverses
        assert pos == len(word_indices), (pos, len(word_indices))
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
        T = self._tree = ConjugateTree()
        self._curves = []
        self._slot = []
        for j, c in enumerate(curves):
            w = word_init(c)
            if check:
                w = word_cyclically_reduce(w)
            if not w:
                raise ValueError(f"trivial curve at index {j}")
            self._curves.append(w)

            # 2. the slot of the curve, rejecting the non-primitive ones
            status = T.process(w[:], check=False)
            if status > 0:
                if status != 1:
                    # a power of a word that was not in the tree
                    raise NotImplementedError(f"non-primitive curve at index {j}")
                i = len(T._words) - 1
                ans = T.process(word_free_group_inverse(T._words[i]), check=False)
                assert ans == 1
            else:
                i = -status
                if len(w) != len(T._words[i]):
                    # a power of a word that was already in the tree
                    raise NotImplementedError(f"non-primitive curve at index {j}")
            self._slot.append(i >> 1)

        num_slots = len(T._words) // 2

        # 3. the cyclic order at infinity of all the leaves, once. For each
        # slot we keep its own leaves in increasing order of rank, each of them
        # described by its rank, its startpoint and the angle from its
        # startpoint to its endpoint.
        ranks = [[] for _ in range(num_slots)]
        starts = [[] for _ in range(num_slots)]
        arc_angles = [[] for _ in range(num_slots)]
        for rank, s in enumerate(T.cyclically_sorted_leaves(angles)):
            i, k = T.leaf_as_conjugate(s)
            w = T._words[i]
            startpoint = w[k]
            endpoint = w[k - 1] ^ 1
            slot = i >> 1
            ranks[slot].append(rank)
            starts[slot].append(startpoint)
            arc_angles[slot].append((angles[endpoint] - angles[startpoint]) % n - 1)
        self._ranks = ranks
        self._starts = starts
        self._arc_angles = arc_angles

        # 4. the two per-curve halves of the O(n^2) term of an entry. The arc
        # matrix M of a curve has M[first][last] counting the arcs going from
        # the angle first to the angle last (first < last); P is M prefix
        # summed down each column and S is M suffix summed along each row. The
        # term only reads them at the pairs (i, j) with 1 <= i <= n - 3 and
        # i + 1 <= j <= n - 2, so we flatten P[i - 1][j] and S[i][j + 1] over
        # those pairs and the term becomes a dot product (see _double_sum).
        pairs = [(i, j) for i in range(1, n - 2) for j in range(i + 1, n - 1)]
        self._K = len(pairs)
        self._A = []
        self._B = []
        for slot in range(num_slots):
            w = T._words[2 * slot]
            M = [[0] * n for _ in range(n)]
            for p in range(len(w)):
                first = angles[w[p]]
                last = angles[w[p - 1] ^ 1]
                if last < first:
                    first, last = last, first
                M[first][last] += 1
            P = [row[:] for row in M]
            for j in range(n):
                for i in range(j - 1):
                    P[i + 1][j] += P[i][j]
            S = M  # M itself is not needed anymore
            for i in range(n):
                for j in range(n - 1, i + 1, -1):
                    S[i][j - 1] += S[i][j]
            self._A.append([P[i - 1][j] for i, j in pairs])
            self._B.append([S[i][j + 1] for i, j in pairs])

        # Scratch space for the sweeps, allocated once. The sweep does one
        # partial sum and two updates per arc, and on that mix the naive
        # structure (O(1) updates, partial sums done by a C level sum over a
        # slice) measures faster than the binary splitting one up to n in the
        # hundreds, which covers every map this class is used on in practice.
        cls = PartialSumsNaive if n <= 256 else PartialSums
        self._Nu = cls(n - 1)
        self._Nv = cls(n - 1)

        # The float64 copies of _A and _B that row() and matrix() multiply,
        # with the element-products asked of them so far: None while they have
        # not been built, () once it is known that they cannot be. See
        # _dot_arrays.
        self._dot = None
        self._dot_work = 0

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

        This is a dot product of the flat vectors built at construction time;
        it costs `O(n^2)` where `n` is the number of half-edges.
        """
        Au = self._A[sx]
        Bu = self._B[sx]
        if sx == sy:
            return 2 * sum(a * b for a, b in zip(Au, Bu))
        Av = self._A[sy]
        Bv = self._B[sy]
        return sum(a * b + c * d for a, b, c, d in zip(Au, Bv, Av, Bu))

    # Setting up the many-pairs paths means importing numpy and building two
    # num_slots x K arrays, which measures at about 30 ms on a cold process;
    # the Python dot products they replace run at about 0.05 µs per element of
    # A. So the setup pays for itself after some 600000 element-products, and
    # _dot_arrays counts them across calls before doing it, which is what
    # makes a first call that is just under the line correct itself on the
    # second one rather than pay for a numpy it will not use.
    _DOT_SETUP_WORK = 600000

    def _dot_arrays(self, work):
        r"""
        Return the two ``num_slots x K`` ``float64`` arrays holding the flat
        vectors of every slot, or the empty tuple when they are not to be used.

        ``work`` is the number of element-products the caller would otherwise
        run through :meth:`_double_sum`. The arrays are built once the calls
        seen so far add up to more than the setup costs, and reused from then
        on; until then the empty tuple is returned and nothing is cached, so
        that a later and larger call can still build them.

        The many-pairs paths of :meth:`row` and :meth:`matrix` evaluate the
        `O(n^2)` term of all the pairs at once as a matrix product. NumPy has
        no BLAS path for integer matrix products — an ``int64`` product falls
        back to a naive loop and is an order of magnitude slower than
        ``float64`` — so the product goes through ``float64``.

        That is exact here rather than approximate. Every entry of ``A`` and
        ``B`` is at most the length of its curve, since both are partial sums
        of a matrix whose entries sum to that length, so a dot product is at
        most ``K * L^2``; below `2^{53}` every integer is a ``float64``. The
        bound is checked rather than assumed, and the empty tuple is returned
        when it does not hold, which sends the callers back to the exact
        Python dot product of :meth:`_double_sum`. It cannot in fact be
        reached while this is worth doing — the `O(n^2)` term is already
        negligible by the time the curves are long enough to overflow, the two
        regimes being disjoint — but a single comparison is a cheap way not to
        rely on that.
        """
        if self._dot is None:
            self._dot_work += work
            if self._dot_work < self._DOT_SETUP_WORK or not self._A:
                return ()
            import numpy

            lmax = max(len(w) for w in self._curves)
            if self._K * lmax * lmax >= 2 ** 53:
                self._dot = ()
            else:
                self._dot = (numpy.array(self._A, dtype=numpy.float64),
                             numpy.array(self._B, dtype=numpy.float64))
        return self._dot

    @staticmethod
    def _exact_ints(values):
        r"""
        Return the ``float64`` array ``values`` as integers, checking that
        nothing was lost on the way there.
        """
        import numpy

        rounded = numpy.rint(values)
        assert numpy.array_equal(rounded, values), "float64 lost the O(n^2) term"
        return rounded.astype(numpy.int64)

    def _double_sum_row(self, sx):
        r"""
        Return ``[self._double_sum(sx, sy) for sy in range(num_slots)]``, as
        two matrix-vector products, or ``None`` when that path is unavailable.
        """
        arrays = self._dot_arrays(len(self._curves) * self._K)
        if not arrays:
            return None
        A, B = arrays
        return self._exact_ints(A[sx] @ B.T + B[sx] @ A.T).tolist()

    def _double_sum_table(self):
        r"""
        Return the ``num_slots x num_slots`` array ``D`` with
        ``D[x][y] = A[x] . B[y]``, as one matrix product, or ``None`` when that
        path is unavailable.

        ``self._double_sum(x, y)`` is ``D[x][y] + D[y][x]``; the two halves are
        kept apart so that the symmetrization costs nothing here and is done
        one row at a time by :meth:`matrix`.
        """
        N = len(self._curves)
        arrays = self._dot_arrays(N * (N + 1) // 2 * self._K)
        if not arrays:
            return None
        A, B = arrays
        return self._exact_ints(A @ B.T)

    def _entry_from(self, sx, sy, double_sum):
        r"""
        Return the intersection number of the slots ``sx`` and ``sy``, given
        the value ``double_sum`` of ``self._double_sum(sx, sy)``.
        """
        # the two arcs of a crossing are counted once in each direction
        ans = 2 * double_sum + self._sweep(sx, sy)
        assert ans % 2 == 0
        return ans // 2

    def _sweep(self, sx, sy):
        r"""
        Return the contribution of the pairs of arcs with identical startpoint,
        for the slots ``sx`` and ``sy``.

        The leaves of the two slots are merged by rank and the resulting list
        is swept by groups of equal startpoints. This costs
        `O((|u| + |v|) \log(n))`.

        The two :class:`~combisurf.partial_sums.PartialSums` are zero on entry
        and are restored to zero at the end of each group, by undoing the
        updates of the group rather than by clearing the whole vector.
        """
        Nu = self._Nu
        Nv = self._Nv
        ans = 0

        if sx == sy:
            # a slot against itself is swept once with both multiplicities
            # equal to one; merging it with a copy of itself is wrong
            starts = self._starts[sx]
            arc_angles = self._arc_angles[sx]
            update = Nu.update
            partial_sum = Nu.partial_sum
            l = len(starts)
            pos = 0
            while pos < l:
                startpoint = starts[pos]
                first = pos
                pos += 1
                while pos < l and starts[pos] == startpoint:
                    pos += 1
                if pos - first == 1:
                    # a single arc through this startpoint crosses nothing
                    continue
                update(arc_angles[first], 1)
                for t in range(first + 1, pos):
                    angle = arc_angles[t]
                    ans += 2 * partial_sum(0, angle)
                    update(angle, 1)
                for t in range(first, pos):
                    update(arc_angles[t], -1)
            return ans

        ru = self._ranks[sx]
        su = self._starts[sx]
        au = self._arc_angles[sx]
        rv = self._ranks[sy]
        sv = self._starts[sy]
        av = self._arc_angles[sy]
        lu = len(ru)
        lv = len(rv)
        iu = iv = 0
        while iu < lu or iv < lv:
            # the leaves sharing a startpoint are consecutive in the cyclic
            # order, so the head of smaller rank opens the group
            if iv == lv or (iu < lu and ru[iu] < rv[iv]):
                startpoint = su[iu]
            else:
                startpoint = sv[iv]
            iu0 = iu
            while iu < lu and su[iu] == startpoint:
                iu += 1
            iv0 = iv
            while iv < lv and sv[iv] == startpoint:
                iv += 1
            if iu0 == iu or iv0 == iv:
                # only one of the two curves goes through this startpoint
                continue

            ju = iu0
            jv = iv0
            while ju < iu or jv < iv:
                if jv == iv or (ju < iu and ru[ju] < rv[jv]):
                    angle = au[ju]
                    ans += Nv.partial_sum(0, angle)
                    Nu.update(angle, 1)
                    ju += 1
                else:
                    angle = av[jv]
                    ans += Nu.partial_sum(0, angle)
                    Nv.update(angle, 1)
                    jv += 1
            for t in range(iu0, iu):
                Nu.update(au[t], -1)
            for t in range(iv0, iv):
                Nv.update(av[t], -1)
        return ans

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

        The `O(n^2)` term of all the entries of the row is obtained at once,
        as two matrix-vector products; see :meth:`_dot_arrays`.

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
        slot = self._slot
        sx = slot[x]
        terms = self._double_sum_row(sx)
        if terms is None:
            return [self._entry_from(sx, sy, self._double_sum(sx, sy)) for sy in slot]
        return [self._entry_from(sx, sy, terms[sy]) for sy in slot]

    def matrix(self):
        r"""
        Return the full symmetric matrix of geometric intersection numbers over
        the integers.

        Only the entries with ``x <= y`` are computed, the others being
        obtained by symmetry, and the `O(n^2)` term of all of them is obtained
        at once as a single matrix product; see :meth:`_dot_arrays`.

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
        D = self._double_sum_table()
        rows = [[0] * N for _ in range(N)]
        for x in range(N):
            sx = slot[x]
            terms = None if D is None else (D[sx] + D[:, sx]).tolist()
            for y in range(x, N):
                sy = slot[y]
                if terms is None:
                    e = self._entry_from(sx, sy, self._double_sum(sx, sy))
                else:
                    e = self._entry_from(sx, sy, terms[sy])
                rows[x][y] = e
                rows[y][x] = e
        return sage_matrix(ZZ, rows)
