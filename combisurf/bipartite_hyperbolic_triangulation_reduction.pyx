r"""
Strongly reduced closed walks of bipartite hyperbolic triangulations

This module holds :func:`strongly_reduced`, the reduction of a closed walk
of a
:class:`~combisurf.bipartite_hyperbolic_triangulation.BipartiteHyperbolicTriangulation`
to its strongly reduced closed walk, called by
:meth:`~combisurf.bipartite_hyperbolic_triangulation.BipartiteHyperbolicTriangulation.strongly_reduced`,
where the moves are described. It is the algorithm of
:class:`~combisurf.bipartite_hyperbolic_triangulation_reduction_naive.TurnRuns`,
its pure Python reference, on C arrays: the same moves are applied in the
same order to the same run numbering, so that both return the same list of
half-edges, rotation included.

The runs of weighted turns are the nodes of a circular doubly linked list
held in parallel C arrays (value, length, first half-edge, next, previous,
alive, dirty). A node is never reused: a move removes nodes and appends new
ones, at most seven, so that the arrays are grown by doubling their capacity
whenever fewer than eight free nodes are left before a move. The number of
moves is linear in the length of the walk, and so is the memory.

EXAMPLES::

    sage: from array import array
    sage: from combisurf import OrientedMap
    sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
    sage: from combisurf.bipartite_hyperbolic_triangulation_reduction import strongly_reduced
    sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
    sage: T = BipartiteHyperbolicTriangulation(*m.triangle_system(relabel=True))
    sage: strongly_reduced(T._h2v, T._wpos, T._wdeg, T._offset, T._at, [4, 3, 0, 3, 2, 5])
    [3, 0]
"""
# ****************************************************************************
#  This file is part of combisurf
#
#       Copyright (C) 2026 Vincent Delecroix
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# ****************************************************************************

from cpython cimport array
from libc.limits cimport INT_MAX
from libc.stdlib cimport free, malloc, realloc
from libc.string cimport memmove

from cysignals.signals cimport sig_check

cdef enum:
    # return codes of the moves besides a node: no move applies, the walk
    # vanishes, an exception is set
    NO_MOVE = -1
    VANISHED = -2
    ERROR = -3
    # a move appends at most seven nodes, and a rewrite fills at most seven
    # slots of a Pending
    MAX_NEW = 8


cdef struct Runs:
    int* val
    int* len
    int* first
    int* nxt
    int* prv
    char* alive
    char* dirty
    Py_ssize_t cap
    Py_ssize_t count
    Py_ssize_t size
    Py_ssize_t num_dirty
    Py_ssize_t num_turns
    const int* h2v
    const int* wpos
    const int* wdeg
    const int* offset
    const int* at


# runs waiting to be inserted: slot 0 and the slot after the inner runs are
# filled by the rewrite that surrounds them
cdef struct Pending:
    int val[MAX_NEW]
    int len[MAX_NEW]
    int first[MAX_NEW]
    int count


cdef inline int pos_mod(int t, int d) noexcept nogil:
    r"""
    Return ``t`` modulo ``d > 0`` in `[0, d)`, also for a negative ``t``,
    which the C ``%`` of ``cdivision`` truncates toward zero.
    """
    t %= d
    return t + d if t < 0 else t


cdef inline int mod_normalize(int t, int d) noexcept nogil:
    r"""
    Return the representative of ``t`` modulo ``d`` in `(-d/2, d/2]`.
    """
    t = pos_mod(t, d)
    return t - d if 2 * t > d else t


cdef inline int norm(Runs* R, int t, int h) noexcept nogil:
    return mod_normalize(t, R.wdeg[R.h2v[h]])


cdef int half_edge_at(Runs* R, int h, int t) except -1:
    r"""
    Return the half-edge at weighted position ``t`` counterclockwise from
    ``h`` around its vertex.
    """
    cdef int v = R.h2v[h]
    cdef int x = R.at[R.offset[v] + pos_mod(R.wpos[h] + t, R.wdeg[v])]
    if x == -1:
        raise RuntimeError(f"the weighted position {t} from the half-edge {h} is inside a punctured corner")
    return x


cdef int grow(Runs* R) except -1:
    r"""
    Make room for ``MAX_NEW`` more nodes.
    """
    if R.count + MAX_NEW <= R.cap:
        return 0
    cdef Py_ssize_t cap = 2 * R.cap + MAX_NEW
    cdef void* p
    p = realloc(R.val, cap * sizeof(int))
    if p == NULL:
        raise MemoryError
    R.val = <int*> p
    p = realloc(R.len, cap * sizeof(int))
    if p == NULL:
        raise MemoryError
    R.len = <int*> p
    p = realloc(R.first, cap * sizeof(int))
    if p == NULL:
        raise MemoryError
    R.first = <int*> p
    p = realloc(R.nxt, cap * sizeof(int))
    if p == NULL:
        raise MemoryError
    R.nxt = <int*> p
    p = realloc(R.prv, cap * sizeof(int))
    if p == NULL:
        raise MemoryError
    R.prv = <int*> p
    p = realloc(R.alive, cap * sizeof(char))
    if p == NULL:
        raise MemoryError
    R.alive = <char*> p
    p = realloc(R.dirty, cap * sizeof(char))
    if p == NULL:
        raise MemoryError
    R.dirty = <char*> p
    R.cap = cap
    return 0


cdef inline int new_node(Runs* R, int value, int length, int first) noexcept nogil:
    cdef int i = <int> R.count
    R.count += 1
    R.val[i] = value
    R.len[i] = length
    R.first[i] = first
    R.nxt[i] = i
    R.prv[i] = i
    R.alive[i] = 1
    R.dirty[i] = 1
    R.num_dirty += 1
    R.size += 1
    return i


cdef inline void remove(Runs* R, int u) noexcept nogil:
    cdef int p = R.prv[u]
    cdef int q = R.nxt[u]
    R.nxt[p] = q
    R.prv[q] = p
    R.alive[u] = 0
    if R.dirty[u]:
        R.dirty[u] = 0
        R.num_dirty -= 1
    R.size -= 1


cdef int insert_after(Runs* R, int a, Pending* runs, int start, int stop) noexcept nogil:
    r"""
    Insert the runs ``start`` to ``stop - 1`` of ``runs`` after the node
    ``a`` and return the node of the first.
    """
    cdef int b = R.nxt[a]
    cdef int first = -1
    cdef int i, u
    for i in range(start, stop):
        u = new_node(R, runs.val[i], runs.len[i], runs.first[i])
        R.nxt[a] = u
        R.prv[u] = a
        a = u
        if first == -1:
            first = u
    R.nxt[a] = b
    R.prv[b] = a
    return first


cdef int merge(Runs* R, int start, int count) noexcept nogil:
    r"""
    Merge the equal runs of `\pm 3` among the ``count`` adjacencies that
    follow the node ``start`` and return a live node among them.
    """
    cdef int u = start
    cdef int v, _i
    for _i in range(count):
        v = R.nxt[u]
        if u != v and R.val[u] == R.val[v] and (R.val[u] == 3 or R.val[u] == -3):
            # u comes first and keeps its first half-edge
            R.len[u] += R.len[v]
            remove(R, v)
        else:
            u = v
    return start if R.alive[start] else u


cdef int mark_dirty(Runs* R, int start, Py_ssize_t count) noexcept nogil:
    r"""
    Mark dirty ``count`` nodes from two nodes before ``start`` and return the
    first of them.
    """
    cdef int u = R.prv[R.prv[start]]
    cdef int x = u
    cdef Py_ssize_t _i
    if count > R.size:
        count = R.size
    for _i in range(count):
        if not R.dirty[x]:
            R.dirty[x] = 1
            R.num_dirty += 1
        x = R.nxt[x]
    return u


cdef int splice(Runs* R, int P, const int* middle, int nm, int N, int first_N, Pending* runs) noexcept nogil:
    r"""
    Replace the last weighted turn of ``P``, the ``nm`` runs ``middle`` and
    the first weighted turn of ``N`` by ``runs``, mark the runs around them
    dirty and return the node to resume the scan from.
    """
    cdef Py_ssize_t removed = 2
    cdef Py_ssize_t added = 0
    cdef int i, u, new
    for i in range(nm):
        u = middle[i]
        removed += R.len[u]
        remove(R, u)
    new = insert_after(R, P, runs, 0, runs.count)
    R.len[P] -= 1
    R.len[N] -= 1
    R.first[N] = first_N
    if R.len[P] == 0:
        remove(R, P)
    if N != P and R.len[N] == 0:
        remove(R, N)
    for i in range(runs.count):
        added += runs.len[i]
    R.num_turns += added - removed
    u = merge(R, R.prv[new], runs.count + 1)
    return mark_dirty(R, u, runs.count + 6)


cdef int replace_all(Runs* R, const int* nodes, int nn, Pending* runs) noexcept nogil:
    r"""
    Replace the whole cycle, made of the ``nn`` nodes ``nodes``, by ``runs``,
    mark it dirty and return a node of it.
    """
    cdef int i, a
    cdef Py_ssize_t total = runs.len[0]
    for i in range(nn):
        remove(R, nodes[i])
    a = new_node(R, runs.val[0], runs.len[0], runs.first[0])
    for i in range(1, runs.count):
        insert_after(R, R.prv[a], runs, i, i + 1)
        total += runs.len[i]
    R.num_turns = total
    a = merge(R, a, runs.count)
    return mark_dirty(R, a, R.size)


cdef inline Py_ssize_t frame(Runs* R, const int* middle, int nm, int* P, int* N) noexcept nogil:
    r"""
    Set ``P`` and ``N`` to the runs before and after the ``nm`` runs
    ``middle`` and return the number of weighted turns of the walk out of
    them.
    """
    cdef Py_ssize_t pattern = 0
    cdef int i
    for i in range(nm):
        pattern += R.len[middle[i]]
    P[0] = R.prv[middle[0]]
    N[0] = R.nxt[middle[nm - 1]]
    return R.num_turns - pattern


cdef int rewrite(Runs* R, int* middle, int nm, int s, Py_ssize_t rest, int P, int N,
                 Pending* runs, int g_last) except ERROR:
    r"""
    Replace the part of the walk made of the runs ``middle`` by a walk with
    the weighted turns ``runs.val[1:runs.count]`` inside, and return the
    node to resume the scan from; see ``TurnRuns._rewrite`` of the reference
    module.

    The slot ``0`` of ``runs`` is left free for the run before; ``middle``
    has room for one more node.
    """
    cdef int h0 = R.first[middle[0]]
    cdef int x = R.val[P]
    cdef int hm, y, first_N, i
    if rest == 1:
        runs.val[0] = norm(R, x - 2 * s, h0)
        runs.len[0] = 1
        runs.first[0] = g_last
        middle[nm] = P
        return replace_all(R, middle, nm + 1, runs)
    hm = R.first[N]
    y = R.val[N]
    runs.val[0] = norm(R, x - s, h0)
    runs.len[0] = 1
    runs.first[0] = half_edge_at(R, h0, -x) ^ 1
    i = runs.count
    runs.val[i] = norm(R, y - s, hm ^ 1)
    runs.len[i] = 1
    runs.first[i] = g_last
    runs.count = i + 1
    first_N = half_edge_at(R, hm ^ 1, y)
    return splice(R, P, middle, nm, N, first_N, runs)


cdef int spur(Runs* R, int u) except ERROR:
    r"""
    Remove the spur at the run ``u``, of value `0`, and return the node to
    resume from, ``VANISHED`` when the walk vanishes.
    """
    cdef int middle[1]
    cdef int P, N, h0, x, y, first_N
    cdef Pending runs
    middle[0] = u
    cdef Py_ssize_t rest = frame(R, middle, 1, &P, &N)
    if rest <= 1:
        return VANISHED
    h0 = R.first[u]
    x = R.val[P]
    y = R.val[N]
    runs.val[0] = norm(R, x + y, h0)
    runs.len[0] = 1
    runs.first[0] = half_edge_at(R, h0, -x) ^ 1
    runs.count = 1
    first_N = half_edge_at(R, h0, y)
    return splice(R, P, middle, 1, N, first_N, &runs)


cdef int spike(Runs* R, int u, int s) except ERROR:
    r"""
    Remove the spike at the run ``u``, of value ``s = ±1``, and return the
    node to resume from, ``NO_MOVE`` if there is none.
    """
    cdef int middle[2]
    cdef int P, N, f
    cdef Pending runs
    middle[0] = u
    cdef Py_ssize_t rest = frame(R, middle, 1, &P, &N)
    if rest < 1:
        return NO_MOVE
    f = half_edge_at(R, R.first[u], -s)
    runs.count = 1
    return rewrite(R, middle, 1, s, rest, P, N, &runs, f)


cdef int bracket(Runs* R, int u, int s) except ERROR:
    r"""
    Remove the bracket ``2s, (3s)^k, 2s`` that starts at the run ``u``, if
    there is one, and return the node to resume from, ``NO_MOVE`` if there is
    none.
    """
    cdef int middle[4]
    cdef int nm, k, b, c, P, N, g_first, g_last
    cdef Pending runs
    b = R.nxt[u]
    middle[0] = u
    if R.val[b] == 3 * s:
        k = R.len[b]
        c = R.nxt[b]
        middle[1] = b
        middle[2] = c
        nm = 3
    else:
        k = 0
        c = b
        middle[1] = c
        nm = 2
    if R.val[c] != 2 * s:
        return NO_MOVE
    cdef Py_ssize_t rest = frame(R, middle, nm, &P, &N)
    if rest < 1:
        return NO_MOVE
    g_first = half_edge_at(R, R.first[u], -s)
    g_last = half_edge_at(R, R.first[N] ^ 1, s) ^ 1
    runs.val[1] = -3 * s
    runs.len[1] = k + 1
    runs.first[1] = g_first
    runs.count = 2
    return rewrite(R, middle, nm, s, rest, P, N, &runs, g_last)


cdef inline bint not_flip_end(int t) noexcept nogil:
    r"""
    Return whether ``t``, a weighted turn at an end of a flip multiplied by
    the side of the flip, is one of `-1, 0, 1, 2, 3`.
    """
    return -1 <= t <= 3


cdef int flip(Runs* R, int u, int s) except ERROR:
    r"""
    Flip the run ``u``, of value ``2s``, if it is white and is the middle of
    a flip ``x, (3s)^k, 2s, (3s)^l, y``, and return the node to resume from,
    ``NO_MOVE`` otherwise.
    """
    cdef int middle[4]
    cdef int nm, k, l, p, q, P, N, g_first, g_last, c1, e_in, apex, i
    cdef Pending runs
    if R.first[u] & 1:
        return NO_MOVE
    p = R.prv[u]
    q = R.nxt[u]
    if R.val[p] == 3 * s:
        k = R.len[p]
        middle[0] = p
        middle[1] = u
        nm = 2
    else:
        k = 0
        middle[0] = u
        nm = 1
    if R.val[q] == 3 * s:
        l = R.len[q]
        middle[nm] = q
        nm += 1
    else:
        l = 0
    cdef Py_ssize_t rest = frame(R, middle, nm, &P, &N)
    if rest < 1 or not_flip_end(s * R.val[P]) or not_flip_end(s * R.val[N]):
        return NO_MOVE
    g_first = half_edge_at(R, R.first[middle[0]], -s)
    g_last = half_edge_at(R, R.first[N] ^ 1, s) ^ 1
    c1 = half_edge_at(R, R.first[u] ^ 1, s) ^ 1
    i = 1
    if k:
        runs.val[i] = -2 * s
        runs.len[i] = 1
        runs.first[i] = g_first
        i += 1
        if k > 1:
            runs.val[i] = -3 * s
            runs.len[i] = k - 1
            runs.first[i] = half_edge_at(R, g_first ^ 1, -2 * s)
            i += 1
        e_in = half_edge_at(R, c1, 2 * s) ^ 1
    else:
        e_in = g_first
    apex = -s * (2 + (k > 0) + (l > 0))
    runs.val[i] = norm(R, apex, c1)
    runs.len[i] = 1
    runs.first[i] = e_in
    i += 1
    if l:
        if l > 1:
            runs.val[i] = -3 * s
            runs.len[i] = l - 1
            runs.first[i] = half_edge_at(R, c1, -2 * s)
            i += 1
        runs.val[i] = -2 * s
        runs.len[i] = 1
        runs.first[i] = half_edge_at(R, g_last, 2 * s) ^ 1
        i += 1
    runs.count = i
    return rewrite(R, middle, nm, s, rest, P, N, &runs, g_last)


cdef int cylinder(Runs* R, int u) except ERROR:
    r"""
    Replace the walk, whose weighted turns are all ``3s`` and form the single
    run ``u``, by the other side of the cylinder of triangles on its side
    ``s``, and return its node.
    """
    # the value is 3 or -3, so that the division is exact
    cdef int s = R.val[u] // 3
    cdef int h = R.first[u]
    cdef int f = half_edge_at(R, half_edge_at(R, h ^ 1, s) ^ 1, -s)
    cdef int nodes[1]
    cdef Pending runs
    nodes[0] = u
    runs.val[0] = -3 * s
    runs.len[0] = <int> R.num_turns
    runs.first[0] = f
    runs.count = 1
    return replace_all(R, nodes, 1, &runs)


cdef int reduce(Runs* R) except ERROR:
    r"""
    Apply the moves until there is none and return a node of the walk, or
    ``VANISHED`` when the walk vanishes.
    """
    cdef Py_ssize_t limit = 8 * R.num_turns + 16
    cdef Py_ssize_t moves = 0
    cdef int cur = 0
    cdef int t, r, s
    while R.num_dirty:
        sig_check()
        if not R.dirty[cur]:
            cur = R.nxt[cur]
            continue
        grow(R)
        t = R.val[cur]
        if t == 0:
            r = spur(R, cur)
            if r == VANISHED:
                return VANISHED
        elif t == 1 or t == -1:
            r = spike(R, cur, t)
        elif t == 2 or t == -2:
            # t is 2 or -2, so that the division is exact
            s = t // 2
            r = bracket(R, cur, s)
            if r == NO_MOVE:
                r = flip(R, cur, s)
        else:
            r = NO_MOVE
        if r == NO_MOVE:
            R.dirty[cur] = 0
            R.num_dirty -= 1
            cur = R.nxt[cur]
        else:
            cur = r
            moves += 1
            if moves > limit:
                raise RuntimeError("the reduction does not terminate; the triangulation is not a bipartite hyperbolic triangulation")
    return cur


cdef int check_table(array.array a, str name) except -1:
    if a.ob_descr.typecode != b'i':
        raise TypeError(f"{name} must be an array of typecode 'i'")
    return 0


def strongly_reduced(array.array h2v not None, array.array wpos not None, array.array wdeg not None,
                     array.array offset not None, array.array at not None, w):
    r"""
    Return the strongly reduced closed walk freely homotopic to the closed
    walk ``w`` of a bipartite hyperbolic triangulation given by its tables.

    This is the reduction of
    :meth:`~combisurf.bipartite_hyperbolic_triangulation.BipartiteHyperbolicTriangulation.strongly_reduced`,
    which checks its input and calls this function. Here ``w`` is assumed to
    be a closed walk: only its half-edges are checked to be active
    half-edges of the tables.

    INPUT:

    - ``h2v``, ``wpos``, ``wdeg``, ``offset``, ``at`` -- arrays of typecode
      ``'i'``, the attributes ``_h2v``, ``_wpos``, ``_wdeg``, ``_offset`` and
      ``_at`` of a
      :class:`~combisurf.bipartite_hyperbolic_triangulation.BipartiteHyperbolicTriangulation`

    - ``w`` -- a closed walk, a list or an array of half-edges

    OUTPUT: a list of half-edges, empty when ``w`` is contractible; the same
    list, rotation included, as
    :meth:`~combisurf.bipartite_hyperbolic_triangulation_reduction_naive.TurnRuns.strongly_reduced`

    EXAMPLES::

        sage: from array import array
        sage: from combisurf import OrientedMap
        sage: from combisurf.bipartite_hyperbolic_triangulation import BipartiteHyperbolicTriangulation
        sage: from combisurf.bipartite_hyperbolic_triangulation_reduction import strongly_reduced
        sage: from combisurf.bipartite_hyperbolic_triangulation_reduction_naive import TurnRuns
        sage: m = OrientedMap(fp="(0,1,~0,~1,2,3,~2,~3)")
        sage: t, mor, _ = m.triangle_system(relabel=True, mapping=True)
        sage: T = BipartiteHyperbolicTriangulation(t)
        sage: tables = (T._h2v, T._wpos, T._wdeg, T._offset, T._at)

    The boundary of a triangle is contractible, a conjugate of the image of
    the edge ``0`` of ``m`` is not::

        sage: strongly_reduced(*tables, t.faces()[0])
        []
        sage: w = array('i', [x for h in [2, 0, 3] for x in mor[h]])
        sage: strongly_reduced(*tables, w)
        [3, 0]
        sage: strongly_reduced(*tables, w) == TurnRuns(T, w).strongly_reduced()
        True

    A walk whose weighted turns are all white `-3` is replaced by the other
    side of its cylinder of triangles::

        sage: strongly_reduced(*tables, [0, 18, 24])
        [24, 4, 12]

    TESTS:

    The empty walk, half-edges out of range and tables of a wrong type::

        sage: strongly_reduced(*tables, [])
        []
        sage: strongly_reduced(*tables, [0, 1000])
        Traceback (most recent call last):
        ...
        ValueError: 1000 is not an active half-edge
        sage: strongly_reduced(array('l', T._h2v), *tables[1:], [0, 1])
        Traceback (most recent call last):
        ...
        TypeError: h2v must be an array of typecode 'i'

    ALGORITHM:

    The one of
    :class:`~combisurf.bipartite_hyperbolic_triangulation_reduction_naive.TurnRuns`,
    on runs held in C arrays freed on every exit, an interruption by
    ``sig_check`` included.
    """
    check_table(h2v, "h2v")
    check_table(wpos, "wpos")
    check_table(wdeg, "wdeg")
    check_table(offset, "offset")
    check_table(at, "at")
    cdef array.array word
    if isinstance(w, array.array) and (<array.array> w).ob_descr.typecode == b'i':
        word = <array.array> w
    else:
        word = array.array('i', w)
    cdef Py_ssize_t n = len(word)
    if n == 0:
        return []
    if n >= INT_MAX // 16:
        raise OverflowError("the walk is too long")
    cdef const int* ww = word.data.as_ints
    cdef const int* hv = h2v.data.as_ints
    cdef Py_ssize_t num_half_edges = len(h2v)
    cdef Py_ssize_t i
    cdef int h, g, t, k
    for i in range(n):
        h = ww[i]
        if h < 0 or h >= num_half_edges or (h ^ 1) >= num_half_edges or hv[h] < 0 or hv[h ^ 1] < 0:
            raise ValueError(f"{h} is not an active half-edge")

    cdef Runs R
    R.val = R.len = R.first = R.nxt = R.prv = NULL
    R.alive = R.dirty = NULL
    R.h2v = hv
    R.wpos = wpos.data.as_ints
    R.wdeg = wdeg.data.as_ints
    R.offset = offset.data.as_ints
    R.at = at.data.as_ints
    cdef int* out = NULL
    cdef int u, j, _run, _turn
    try:
        R.cap = n + MAX_NEW
        R.val = <int*> malloc(R.cap * sizeof(int))
        R.len = <int*> malloc(R.cap * sizeof(int))
        R.first = <int*> malloc(R.cap * sizeof(int))
        R.nxt = <int*> malloc(R.cap * sizeof(int))
        R.prv = <int*> malloc(R.cap * sizeof(int))
        R.alive = <char*> malloc(R.cap * sizeof(char))
        R.dirty = <char*> malloc(R.cap * sizeof(char))
        if (R.val == NULL or R.len == NULL or R.first == NULL or R.nxt == NULL or R.prv == NULL
                or R.alive == NULL or R.dirty == NULL):
            raise MemoryError

        k = 0
        for i in range(n):
            h = ww[i]
            g = h ^ 1
            t = mod_normalize(R.wpos[ww[i + 1 if i + 1 < n else 0]] - R.wpos[g], R.wdeg[hv[g]])
            if k and R.val[k - 1] == t and (t == 3 or t == -3):
                R.len[k - 1] += 1
            else:
                R.val[k] = t
                R.len[k] = 1
                R.first[k] = h
                k += 1
        if k > 1 and R.val[0] == R.val[k - 1] and (R.val[0] == 3 or R.val[0] == -3):
            # the run at the end goes on at the start: it absorbs the first
            # one and keeps its own first half-edge, which comes before
            R.len[k - 1] += R.len[0]
            k -= 1
            memmove(R.val, R.val + 1, k * sizeof(int))
            memmove(R.len, R.len + 1, k * sizeof(int))
            memmove(R.first, R.first + 1, k * sizeof(int))
        for j in range(k):
            R.nxt[j] = j + 1
            R.prv[j] = j - 1
            R.alive[j] = 1
            R.dirty[j] = 1
        R.nxt[k - 1] = 0
        R.prv[0] = k - 1
        R.count = k
        R.size = k
        R.num_dirty = k
        R.num_turns = n

        u = reduce(&R)
        if u == VANISHED:
            return []
        if R.size == 1:
            t = R.val[u]
            # all white -3 or all black 3
            if (t == -3 or t == 3) and (R.first[u] & 1) == (t == 3):
                grow(&R)
                u = cylinder(&R, u)

        out = <int*> malloc(R.num_turns * sizeof(int))
        if out == NULL and R.num_turns:
            raise MemoryError
        j = 0
        for _run in range(R.size):
            h = R.first[u]
            t = R.val[u]
            for _turn in range(R.len[u]):
                out[j] = h
                j += 1
                h = half_edge_at(&R, h ^ 1, t)
            u = R.nxt[u]
        return [out[i] for i in range(j)]
    finally:
        free(out)
        free(R.val)
        free(R.len)
        free(R.first)
        free(R.nxt)
        free(R.prv)
        free(R.alive)
        free(R.dirty)
