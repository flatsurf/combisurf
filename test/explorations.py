r"""
Helpers to explore collections of maps up to isomorphism by local mutations.

The mutations are those of the automata of the sibling projects veerer
(flips of triangulations, flips of the reduced core automaton of veering
triangulations) and surface-dynamics (the generators ``L`` and ``I`` of
``GL(2, Z)`` acting on origamis), written as functions on lists, together
with :func:`explore` that counts the classes of the states reached from a
seed. They are not methods of :class:`OrientedMap`: they belong to those
automata.

A face permutation ``fp`` is a list over the half-edges, the edge ``e``
having the half-edges ``2e`` and ``2e + 1``; ``fp[2e + 1] == -1`` when ``e``
is folded. Every mutation returns new lists and leaves its arguments
unchanged.

The module is not collected by pytest (no ``test_`` prefix) and imports
neither veerer nor surface-dynamics.
"""
import collections
from array import array

# veerer's edge colours (veerer/constants.py)
RED = 1
BLUE = 2
PURPLE = 4
GREEN = 8


def parse_colours(text):
    r"""
    Return the list of veerer colours of a colour string such as ``"RRBB"``.
    """
    return [{'R': RED, 'B': BLUE, 'P': PURPLE, 'G': GREEN}[c] for c in text]


def origami_face_permutation(r, u):
    r"""
    Return the face permutation of the oriented map of the origami
    ``(r, u)`` (``r``, ``u`` the 0-based images of the ``n`` squares):
    square ``i`` has the half-edges ``4*i, ..., 4*i + 3``, and its face is
    ``(4i, 4i+2, 4*r[i]+1, 4*u[i]+3)``.
    """
    n = len(r)
    fp = [-1] * (4 * n)
    for i in range(n):
        fp[4 * i] = 4 * i + 2
        fp[4 * i + 2] = 4 * r[i] + 1
        fp[4 * r[i] + 1] = 4 * u[i] + 3
        fp[4 * u[i] + 3] = 4 * i
    return fp


def origami_decoration(n):
    r"""
    Return the half-edge decoration ``h % 4`` of an origami map of ``n``
    squares (the position of ``h`` in the face of its square): it restores
    the origami's orientation of the squares, which the map alone forgets.
    """
    return [h % 4 for h in range(4 * n)]


def _twin(fp, h):
    r"""
    Return the other half-edge of the edge of ``h``, or ``h`` itself on a
    folded edge (veerer ``Constellation._ep``).
    """
    return h if fp[h ^ 1] == -1 else h ^ 1


def is_flippable(fp, e):
    r"""
    Return whether the edge ``e`` of the triangulation ``fp`` is flippable:
    the two triangles on its sides make a quadrilateral (veerer
    ``edge_triangulation_is_flippable``, no boundary). A folded edge is
    flippable.
    """
    h = 2 * e
    H = _twin(fp, h)
    a = fp[h]
    b = fp[a]
    return a != H and b != H


def flippable_edges(fp):
    r"""
    Return the list of the flippable edges of the triangulation ``fp``
    (veerer ``Triangulation.flippable_edges``; the edges whose half-edges
    are both inactive are skipped, veerer has none).
    """
    return [e for e in range(len(fp) // 2) if fp[2 * e] != -1 and is_flippable(fp, e)]


def flip(fp, e):
    r"""
    Return the face permutation of the triangulation ``fp`` with the edge
    ``e`` flipped (veerer ``Triangulation.flip``).

    On a folded edge (``H == h``) the assignments are the same, in the
    same order: they overwrite one another and reverse the triangle of
    ``h``, as in veerer.
    """
    fp = list(fp)
    h = 2 * e
    H = _twin(fp, h)
    a = fp[h]
    b = fp[a]
    if a == H or b == H:
        raise ValueError('edge %d is not flippable' % e)
    c = fp[H]
    d = fp[c]
    fp[h] = b
    fp[b] = c
    fp[c] = h
    fp[a] = H
    fp[H] = d
    fp[d] = a
    return fp


def _square(fp, h):
    r"""
    Return the half-edges ``a, b, c, d`` of the quadrilateral whose diagonal
    is ``h`` (veerer ``Triangulation.square_about_half_edge``).
    """
    a = fp[h]
    b = fp[a]
    c = fp[_twin(fp, h)]
    d = fp[c]
    return a, b, c, d


def is_forward_flippable(fp, colours, e):
    r"""
    Return whether the edge ``e`` of the veering triangulation
    ``(fp, colours)`` is forward flippable (veerer
    ``edge_veering_triangulation_is_forward_flippable``, no boundary).
    """
    if colours[e] == GREEN or not is_flippable(fp, e):
        return False
    if colours[e] == PURPLE:
        return True
    a, b, c, d = _square(fp, 2 * e)
    return (bool(colours[a // 2] & (BLUE | GREEN)) and bool(colours[b // 2] & (RED | GREEN))
            and bool(colours[c // 2] & (BLUE | GREEN)) and bool(colours[d // 2] & (RED | GREEN)))


def edge_has_curve(fp, colours, e):
    r"""
    Return whether a curve of the dual train track crosses the edge ``e``
    of the veering triangulation ``(fp, colours)``: a breadth-first search
    in the train track (veerer ``VeeringTriangulation.edge_has_curve``).
    """
    a, b, c, d = _square(fp, 2 * e)
    if colours[a // 2] == BLUE or colours[b // 2] == RED:
        POS, NEG = BLUE, RED
    else:
        POS, NEG = RED, BLUE
    if colours[e] == NEG:
        start = b
        end = _twin(fp, d)
    else:
        start = a
        end = _twin(fp, c)
    if start == end:
        return True
    seen = [False] * len(fp)
    seen[start] = True
    q = collections.deque([start])
    while q:
        f = q.popleft()
        F = _twin(fp, f)
        r = fp[F]
        s = fp[r]
        if not seen[r] and not (colours[r // 2] == POS and colours[F // 2] == NEG):
            if r == end:
                return True
            seen[r] = True
            q.append(r)
        if not seen[s] and not (colours[s // 2] == NEG and colours[F // 2] == POS):
            if s == end:
                return True
            seen[s] = True
            q.append(s)
    return False


def reduced_core_seed(fp, colours):
    r"""
    Return the state ``(fp, colours)`` of the reduced core automaton of
    the veering triangulation ``(fp, colours)``: its forward flippable
    edges turned ``PURPLE`` (veerer
    ``VeeringTriangulation.forgot_forward_flippable_colour``, the
    ``_seed_setup`` of ``ReducedCoreAutomaton`` but its canonical labels).
    """
    new = list(colours)
    for e in range(len(colours)):
        if fp[2 * e] != -1 and is_forward_flippable(fp, colours, e):
            new[e] = PURPLE
    return list(fp), new


def reduced_core_neighbours(fp, colours):
    r"""
    Iterate over the out-neighbours ``(fp, colours)`` of a state of the
    reduced core automaton (veerer ``ReducedCoreAutomaton._out_neighbors``
    and ``_flip``): each purple edge flipped to ``BLUE`` and to ``RED``,
    the edges of its quadrilateral that became forward flippable turned
    ``PURPLE``, kept when a curve of the train track crosses the edge.
    """
    for e in range(len(colours)):
        if colours[e] != PURPLE:
            continue
        for col in (BLUE, RED):
            nfp = flip(fp, e)
            ncol = list(colours)
            ncol[e] = col
            if edge_has_curve(nfp, ncol, e):
                a, b, c, d = (x // 2 for x in _square(nfp, 2 * e))
                if col == BLUE:
                    if is_forward_flippable(nfp, ncol, b):
                        ncol[b] = PURPLE
                    if b != d and is_forward_flippable(nfp, ncol, d):
                        ncol[d] = PURPLE
                else:
                    if is_forward_flippable(nfp, ncol, a):
                        ncol[a] = PURPLE
                    if a != c and is_forward_flippable(nfp, ncol, c):
                        ncol[c] = PURPLE
                # veerer tests the curve again, after the recolouring
                if edge_has_curve(nfp, ncol, e):
                    yield nfp, ncol


def flip_neighbours(fp):
    r"""
    Return the list of the triangulations obtained from ``fp`` by one flip
    (veerer ``FlipGraph._out_neighbors``).
    """
    return [flip(fp, e) for e in flippable_edges(fp)]


def origami_l(r, u):
    r"""
    Return the image ``(r, u)`` of the origami ``(r, u)`` (0-based lists)
    under the generator ``L`` of ``GL(2, Z)``: ``(r, r o u)``, as in
    surface-dynamics ``Origami_dense_pyx._compute_gl2z_edges``.
    """
    return list(r), [r[x] for x in u]


def origami_i(r, u):
    r"""
    Return the image ``(u, r)`` of the origami ``(r, u)`` under the
    generator ``I`` of ``GL(2, Z)`` (the exchange of the two directions,
    surface-dynamics ``Origami_dense_pyx._compute_gl2z_edges``).
    """
    return list(u), list(r)


def gl2z_neighbours(state):
    r"""
    Return the images of the origami ``state = (r, u)`` under ``L`` and
    ``I``.
    """
    r, u = state
    return [origami_l(r, u), origami_i(r, u)]


def combisurf_key(fp, half_edge_data=None, edge_data=None):
    r"""
    Return a canonical key of the map ``fp`` with its decoration, as
    bytes: equal keys if and only if the decorated maps are isomorphic.

    The key is the face permutation relabelled by
    ``OrientedMap.canonical_labelling`` (``-1`` at the images of the
    inactive half-edges), followed by the relabelled ``edge_data`` then
    the relabelled ``half_edge_data``. The ``OrientedMap`` is built at
    each call: its construction is part of the cost of a query.
    """
    from combisurf import OrientedMap
    m = OrientedMap(fp=fp)
    lab = m.canonical_labelling(half_edge_data=half_edge_data, edge_data=edge_data)
    n = len(fp)
    key = [-1] * n
    for h in range(n):
        if fp[h] != -1:
            key[lab[h]] = lab[fp[h]]
    if edge_data is not None:
        ce = [0] * len(edge_data)
        for e in range(len(edge_data)):
            ce[lab[2 * e] // 2] = edge_data[e]
        key.extend(ce)
    if half_edge_data is not None:
        ch = [0] * n
        for h in range(n):
            ch[lab[h]] = half_edge_data[h]
        key.extend(ch)
    return array('i', key).tobytes()


def explore(seed, neighbours, key):
    r"""
    Return the number of classes of the states reached from ``seed`` by
    ``neighbours`` (a function from a state to an iterable of states),
    two states being in the same class when their ``key`` agree: a
    breadth-first search with a dict from key to state. The states kept
    for expansion are the mutated states themselves, not relabelled
    copies.
    """
    seen = {key(seed): seed}
    queue = collections.deque([seed])
    while queue:
        state = queue.popleft()
        for new in neighbours(state):
            k = key(new)
            if k not in seen:
                seen[k] = new
                queue.append(new)
    return len(seen)
