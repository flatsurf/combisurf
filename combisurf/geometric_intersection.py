r"""
Geometric intersection of arcs and geodesics on punctured and closed surfaces
"""

from array import array

from combisurf.word import word_init, fg_word_is_cyclically_reduced, fg_word_inverse
from combisurf.oriented_map import OrientedMap
from combisurf.conjugate_tree import ConjugateTree
from combisurf.partial_sums import PartialSums

class GeometricIntersection:
    def __init__(self, m):
        r"""

        """
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

    def __call__(self, w1, w2=None):
        r"""
        Return the self-intersection of ``w1`` or the intersection between
        ``w1`` and ``w2``.
        """
        raise NotImplementedError

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
            if check and not fg_word_is_cyclically_reduced(w):
                raise ValueError
            ans = T.process(list(w))
            if ans <= 0:
                # was already given in T
                raise ValueError(f"conjugate words at position {-ans} and {i}")

        n = len(self._cm._vp)
        word_indices = []
        word_shifts = []

        queue = [T._transitions[0][letter] for letter in sorted(T._transitions[0], key = lambda letter: self._angles[letter], reverse=True)]
        while queue:
            s = queue.pop()
            if T._transition_end[s] == -1:
                # leaf
                i, k = T.leaf_as_conjugate(s)
                word_indices.append(i)
                word_shifts.append(k)
            else:
                i = T._transition_word[s]
                p = T._transition_end[s]
                last_letter = T._letter(i, p - 1) ^ 1
                transitions = sorted(T._transitions[s], key = lambda letter: (self._angles[letter] - self._angles[last_letter]) %  n, reverse=True)
                queue.extend(T._transitions[s][letter] for letter in transitions)

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
        words_with_inverse = list(words) + [fg_word_inverse(w) for w in words]
        l = sum(len(w) for w in words_with_inverse)
        print(f"words_with_inverse={words_with_inverse}")
        colors = rainbow(n, 'rgbtuple')
        word_indices, word_shifts = self.conjugate_sort(words_with_inverse)
        assert len(word_indices) == len(word_shifts) == l
        print(f"word_indices={word_indices} word_shifts={word_shifts}")
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

    # TODO: generalize to make the input a list of (possibly non-primitive) walks
    # TODO: some more simple checks
    # - pick unicellular maps with non-trivial automorphisms and check automorphism invariance
    def primitive_curve_geometric_intersection(self, u, v=None, check=True):
        r"""
        Return the geometric intersection between the primitive curves ``u``
        and ``v`` given as cyclically reduced words on the half-edges of the
        underlying map. If ``v`` is not provided, return the self-intersection
        of ``u``.

        If both ``u`` and ``v`` are provided, they must be different.

        EXAMPLES::

            sage: from combisurf import OrientedMap
            sage: from combisurf.word import word_init
            sage: from combisurf.geometric_intersection import GeometricIntersection

            sage: torus = OrientedMap(fp="(0,1,~0,~1)")
            sage: gi = GeometricIntersection(torus)
            sage: gi.primitive_curve_geometric_intersection([0], [2])
            1
            sage: gi.primitive_curve_geometric_intersection([0, 0, 2, 2])
            1

            sage: gi.primitive_curve_geometric_intersection([0], [0, 2])
            1
            sage: gi.primitive_curve_geometric_intersection([0, 2], [0, 2, 0])
            1
            sage: gi.primitive_curve_geometric_intersection([0, 2, 0], [0, 2, 0, 0, 2])
            1

            sage: for u in [[0], [0, 2], [0, 2, 0], [0, 2, 0, 0, 2],  [0, 2, 0, 0, 2, 0, 2, 0]]:
            ....:     assert gi.primitive_curve_geometric_intersection(u) == 0
            sage: for u in [[0, 0, 2, 2], [0, 2, 0, 2, 0, 0], [0, 2, 0, 0, 2, 0, 0, 2, 0, 2],
            ....:           [0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0],
            ....:           [0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0, 2, 0, 0, 2, 0, 2, 0, 0, 2]]:
            ....:     assert gi.primitive_curve_geometric_intersection(u) == 1

        Two examples in genus 2 following Birman-Series p336-337::

            sage: octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")
            sage: gi = GeometricIntersection(octagon)
            sage: w = word_init("0,~1,3")
            sage: gi.primitive_curve_geometric_intersection(w)
            0
            sage: w = word_init("0,1,1,~2,1,1,~2")
            sage: gi.primitive_curve_geometric_intersection(w)
            4

        Testing the simplicity criterion of Lapointe on positive words::

            sage: W = Words([0, 2, 4, 6])
            sage: for l in range(2, 7):
            ....:     for w in W.iterate_by_length(l):
            ....:         if not w.is_primitive():
            ....:             continue
            ....:         bwt = w.BWT()
            ....:         ans1 = all(bwt[i + 1] <= bwt[i] for i in range(l - 1))
            ....:         ans2 = gi.primitive_curve_geometric_intersection(list(w)) == 0
            ....:         assert ans1 == ans2
        """
        if check:
            u = word_init(u)
        U = fg_word_inverse(u)
        words = [u, U]
        if v is None:
            self_intersection = True
            v = u
            V = U
        else:
            self_intersection = False
            if check:
                v = word_init(v)
            V = fg_word_inverse(v)
            words.append(v)
            words.append(V)
        word_indices, word_shifts = self.conjugate_sort(words)
        if len(word_indices) != sum(len(w) for w in words):
            raise ValueError(f"non primitive or identical curves u={u} v={v}")

        n = len(self._cm._vp)

        # Essential intersection coming from pairs of conjugates with four
        # distinct 1-order intervals associated to their startpoints and endpoints
        # NOTE: O(n^2 + len(u) + len(v)) cost
        Nu = [[0] * n for _ in range(n)]
        for i in range(len(u)):
            first = self._angles[u[i]]
            last = self._angles[u[(i - 1) % len(u)] ^ 1]
            assert first != last
            if last < first:
                first, last = last, first
            Nu[first][last] += 1
        Nv = [[0] * n for _ in range(n)]
        for i in range(len(v)):
            first = self._angles[v[i]]
            last = self._angles[v[(i - 1) % len(v)] ^ 1]
            assert first != last
            if last < first:
                first, last = last, first
            Nv[first][last] += 1

        # NOTE: below is a O(n^2) time version of the following O(n^4) time sum
        #     sum(Nu[i0][j0] * Nv[i1][j1] + Nu[i1][j1] * Nv[i0][j0]
        #         for i0 in range(n)
        #         for j0 in range(i0 + 1, n)
        #         for i1 in range(i0 + 1, j0)
        #         for j1 in range(j0 + 1, n))
        # We optimize the computation by transforming Nu and Nv to contain
        # partial sums in respectively i0 and j1 respectively (O(n^2) time).
        # Then we do a double sum in i1, j0 (O(n^2) time).
        for j in range(n):
            for i in range(j - 1):
                Nu[i + 1][j] += Nu[i][j]
        for i in range(n):
            for j in range(n - 1, i + 1, -1):
                Nv[i][j - 1] += Nv[i][j]

        intersections = sum(Nu[i1 - 1][j0] * Nv[i1][j0 + 1] for i1 in range(n - 2) for j0 in range(i1 + 1, n - 1))

        if v is not None:
            intersections *= 2

        # Essential intersections coming from pairs of conjugates with identical
        # start. Total cost is (len(u) + len(v)) * log(n)
        # where the log(n) factor comes from partial sums
        letter = 0  # current letter that is looked at
        pos = 0     # pointer in the list cs
        Nu = PartialSums(n - 1)
        if self_intersection:
            Nv = Nu
        else:
            Nv = PartialSums(n - 1)
        while pos < len(word_indices):
            Nu.reset()
            Nv.reset()
            startpoint = words[word_indices[pos]][word_shifts[pos]]
            startangle = self._angles[startpoint]
            while pos < len(word_indices) and words[word_indices[pos]][word_shifts[pos]] == startpoint:
                i = word_indices[pos]
                w = words[i]
                k = word_shifts[pos]
                endpoint = w[(k - 1) % len(w)] ^ 1
                assert startpoint != endpoint
                endangle = self._angles[endpoint]
                angle = (endangle - startangle) % n
                assert angle >= 1
                angle -= 1
                # print(f"pos={pos} i={i} arc ({startpoint},{endpoint}) angle={angle}")
                if i == 0 or i == 1:
                    # print(f"add intersection from Nv={Nv} and update Nu={Nu}")
                    intersections += Nv.partial_sum(0, angle)
                    Nu.update(angle, 1)
                    # print(f"after update Nu={Nu}")
                elif i == 2 or i == 3:
                    # print(f"add intersection from Nu={Nu} and update Nv={Nv}")
                    intersections += Nu.partial_sum(0, angle)
                    Nv.update(angle, 1)
                    # print(f"after update Nv={Nv}")
                pos += 1

        assert pos == len(word_indices), (pos, len(word_indices))
        assert intersections % 2 == 0
        return intersections // 2
