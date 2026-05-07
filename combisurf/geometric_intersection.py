r"""
Geometric intersection of arcs and geodesics on punctured and closed surfaces
"""

from array import array

from combisurf.word import word_init, word_is_cyclically_reduced, word_cyclically_reduce, word_free_group_inverse
from combisurf.oriented_map import OrientedMap
from combisurf.conjugate_tree import ConjugateTree
from combisurf.partial_sums import PartialSums

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

    def _check_closed_walk(self, w):
        if not isinstance(w, array):
            w = array("i", w)

        if w.typecode != "i":
            raise ValueError

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

        for j in range(n):
            for i in range(j - 1):
                Nu[i + 1][j] += Nu[i][j]
        for i in range(n):
            for j in range(n - 1, i + 1, -1):
                Nv[i][j - 1] += Nv[i][j]

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
