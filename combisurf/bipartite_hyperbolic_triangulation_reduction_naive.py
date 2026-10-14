r"""
Strongly reduced closed walks of bipartite hyperbolic triangulations,
reference implementation

This module holds :class:`TurnRuns`, the pure Python reduction of a closed
walk of a
:class:`~combisurf.bipartite_hyperbolic_triangulation.BipartiteHyperbolicTriangulation`
to its strongly reduced closed walk. It is the reference against which the
Cython function
:func:`~combisurf.bipartite_hyperbolic_triangulation_reduction.strongly_reduced`
of :mod:`combisurf.bipartite_hyperbolic_triangulation_reduction` is tested;
:meth:`~combisurf.bipartite_hyperbolic_triangulation.BipartiteHyperbolicTriangulation.strongly_reduced`
uses the fast one. Both apply the same moves in the same order to the same
run numbering, so that they return the same list of half-edges, rotation
included.

EXAMPLES::

    sage: from combisurf import OrientedMap
    sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
    sage: from combisurf.bipartite_hyperbolic_triangulation_reduction_naive import TurnRuns
    sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
    sage: T = BipartiteHyperbolicTriangulation(*m.triangle_system(relabel=True))
    sage: w = [4, 3, 0, 3, 2, 5]
    sage: TurnRuns(T, w).strongly_reduced() == T.strongly_reduced(w)
    True
"""

from combisurf.hyperbolic_quadrangulation import _mod_normalize

# the weighted turns that the end of a flip must avoid once multiplied by the
# side s of the flip
_NOT_FLIP_END = (-1, 0, 1, 2, 3)


class TurnRuns:
    r"""
    The weighted turns of a closed walk of a
    :class:`~combisurf.bipartite_hyperbolic_triangulation.BipartiteHyperbolicTriangulation`,
    run-length encoded in a circular doubly linked list, and their reduction
    to the strongly reduced closed walk.

    A run is a node with a value, the weighted turn, a length, the number of
    consecutive weighted turns equal to the value, and a first half-edge, the
    half-edge ``h`` such that the first weighted turn of the run is the one
    from ``h`` to the next half-edge of the walk. Only the runs of `\pm 3` are
    merged and can be longer than one. A weighted turn `\pm 3` keeps the
    parity of the half-edges, so that all the turns of a run have the colour
    of its first half-edge.

    INPUT:

    - ``T`` -- a
      :class:`~combisurf.bipartite_hyperbolic_triangulation.BipartiteHyperbolicTriangulation`

    - ``w`` -- a list of half-edges of ``T``, a closed walk of positive length

    EXAMPLES::

        sage: from combisurf import OrientedMap
        sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
        sage: from combisurf.bipartite_hyperbolic_triangulation_reduction_naive import TurnRuns
        sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
        sage: T = BipartiteHyperbolicTriangulation(*m.triangle_system(relabel=True))
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
        # This guards an invariant the moves keep: a punctured corner weighs
        # 4, so the three positions after its first one are empty, and every
        # value the moves read (0, +-1, +-2, +-3) has a side with no
        # punctured corner. The error is reached only by a direct call.
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
