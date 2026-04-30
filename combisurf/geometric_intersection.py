r"""
Geometric intersection of arcs and geodesics on punctured and closed surfaces
"""

from array import array

from combisurf.word import fg_word_is_cyclically_reduced, fg_word_inverse
from combisurf.oriented_map import OrientedMap
from combisurf.conjugate_tree import ConjugateTree

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

    def primitive_curve_geometric_intersection(self, u, v=None):
        r"""
        Return the geometric intersection between the primitive curves ``u``
        and ``v`` given as cyclically reduced words on the half-edges of the
        underlying map. If ``v`` is not provided, return the self-intersection
        of ``u``.

        Note that ``u`` must be different from ``v``.

        EXAMPLES::

            sage: from combisurf import OrientedMap
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
        """
        if not isinstance(u, array):
            u = array('i', u)
        U = fg_word_inverse(u)
        words = [u, U]
        if v is not None:
            if not isinstance(v, array):
                v = array('i', v)
            V = fg_word_inverse(v)
            words.append(v)
            words.append(V)
        word_indices, word_shifts = self.conjugate_sort(words)
        if len(word_indices) != sum(len(w) for w in words):
            raise ValueError(f"non primitive or identical curves u={u} v={v} cs={cs}")

        n = len(self._cm._vp)
        Nu = [[0] * n for _ in range(n)]
        for i in range(len(u)):
            first = self._angles[u[i]]
            last = self._angles[u[(i - 1) % len(u)] ^ 1]
            assert first != last
            if last < first:
                first, last = last, first
            Nu[first][last] += 1
        if v is None:
            Nv = Nu
        else:
            Nv = [[0] * n for _ in range(n)]
            for i in range(len(v)):
                first = self._angles[v[i]]
                last = self._angles[v[(i - 1) % len(v)] ^ 1]
                assert first != last
                if last < first:
                    first, last = last, first
                Nv[first][last] += 1

        intersections = sum(Nu[i0][j0] * Nv[i1][j1] + Nu[i1][j1] * Nv[i0][j0]
                                                        for i0 in range(n)
                                                        for j0 in range(i0 + 1, n)
                                                        for i1 in range(i0 + 1, j0)
                                                        for j1 in range(j0 + 1, n))
        if v is not None:
            intersections *= 2

        # TODO: we should implement the log(g) vs g sweeps using binary encoding
        # pair of arcs with conjugacy
        # we run by batch of arcs with endpoint in a fixed letter
        letter = 0  # current letter that is looked at
        pos = 0     # pointer in the list cs
        while pos < len(word_indices):
            Nu = [0] * (n - 1)
            if v is None:
                Nv = Nu
            else:
                Nv = [0] * (n - 1)
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
                    for a in range(angle):
                        intersections += Nv[a]
                    Nu[angle] += 1
                    # print(f"after update Nu={Nu}")
                elif i == 2 or i == 3:
                    # print(f"add intersection from Nu={Nu} and update Nv={Nv}")
                    for a in range(angle):
                        intersections += Nu[a]
                    Nv[angle] += 1
                    # print(f"after update Nv={Nv}")
                pos += 1

        assert pos == len(word_indices), (pos, len(word_indices))
        assert intersections % 2 == 0
        return intersections // 2
