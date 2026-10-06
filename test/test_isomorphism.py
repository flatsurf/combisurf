r"""
Reference tests of isomorphism, automorphisms and canonical labels on
examples with known answers from veerer and surface-dynamics.

The maps are built with combisurf's own constructors from the encodings of
the two projects (face permutations in edge notation for veerer,
dart-indexed permutations and unicellular words for the fat graphs of
surface-dynamics, cylinder diagrams). Neither project is imported here: the
expected values are copied from their doctests, whose location is given
next to each test.

The functionality under test is reached only through the adapter functions
below. Subspaces are supported through the keyword ``key``; they raise
``NotImplementedError`` for what combisurf does not offer yet (tuples of
permutations), and the tests that need it are expected to fail with that
exception. Disconnected maps are supported,
among them the cylinder diagrams whose components are joined only by the
pairing of their cylinders.

Regular maps of ``S_5`` and ``PSL(2, 7)`` (built with Sage groups) and
origamis as maps (decorated by ``h % 4``) give automorphism counts that
follow from their construction. Every automorphism count is also compared
with a brute-force count that uses neither the canonical key of the
reference implementation nor :class:`OrientedMap`.

The sizes of the collections explored by the automata of veerer (flip graph
of triangulations, reduced core automaton) and surface-dynamics (the
GL(2, Z)-orbit of an origami, with the keys of
:func:`combisurf.isomorphism.canonical_key`) are checked with the mutations
of :mod:`explorations`.

Coloured maps raise ``NotImplementedError`` on the isomorphism methods: the
four methods of :class:`ColoredOrientedMap` and plain maps compared with
coloured maps in :meth:`OrientedMap.is_isomorphic`.

The last tests check the reference implementation
:mod:`combisurf.isomorphism_naive` against the defining properties of
isomorphisms and canonical forms on random maps, connected and
disconnected, with decorations matched exactly and up to relabelling of
their values (on a disconnected map, the classes shared by the components
tie them together). The fast implementation :mod:`combisurf.isomorphism`
is then compared with it by equality of the results, on random maps of up
to 400 half-edges.

The properties on random maps of up to 48 half-edges, the comparison of
the two implementations on them, the automorphism counts and the
isomorphism verdicts of cylinder diagrams are checked under each rule of
the roots (the keyword ``root_rule``) as well, with and without the flag
``ROOTS_DATA``; under every rule the automorphisms form the same group.
"""
import math
import random
import re
from array import array

import pytest

from combisurf import OrientedMap
from combisurf import isomorphism as isomorphism_fast
from combisurf import isomorphism_naive

from explorations import (BLUE, RED, combisurf_key, explore, flip_neighbours,
                          gl2z_neighbours, origami_decoration,
                          origami_face_permutation, origami_i, origami_l,
                          parse_colours, reduced_core_neighbours,
                          reduced_core_seed)

NOT_IMPLEMENTED = pytest.mark.xfail(raises=NotImplementedError, strict=True)

# the rules of the roots, by their names in OrientedMap and their integers
# in the modules isomorphism and isomorphism_naive
ROOT_RULE_NAMES = ["all", "half-edge-max", "half-edge-fewest", "edge-max", "edge-fewest",
                   "half-edge-max+data", "half-edge-fewest+data", "edge-max+data", "edge-fewest+data"]
BASE_ROOT_RULES = [isomorphism_naive.ROOTS_ALL, isomorphism_naive.ROOTS_HALF_EDGE_MAX,
                   isomorphism_naive.ROOTS_HALF_EDGE_FEWEST, isomorphism_naive.ROOTS_EDGE_MAX,
                   isomorphism_naive.ROOTS_EDGE_FEWEST]
ROOT_RULES = BASE_ROOT_RULES + [rule | isomorphism_naive.ROOTS_DATA for rule in BASE_ROOT_RULES[1:]]
WITH_ROOT_RULES = pytest.mark.parametrize("root_rule", ROOT_RULES, ids=ROOT_RULE_NAMES)
WITH_ROOT_RULE_NAMES = pytest.mark.parametrize("root_rule", ROOT_RULE_NAMES)


# Adapters

def decoration_keywords(half_edge_data=None, edge_data=None, face_data=None,
                        face_classes=None, edge_subspace=None):
    r"""
    Return the keywords of the decorations, as the methods of
    :class:`OrientedMap` take them: ``half_edge_data``, ``edge_data`` and
    ``face_data`` matched exactly, ``face_classes`` up to relabelling of its
    values; ``edge_subspace`` becomes the keyword ``key`` (see
    :func:`subspace_key`).
    """
    return {'half_edge_data': half_edge_data, 'edge_data': edge_data,
            'face_data': face_data, 'face_classes': face_classes,
            'key': None if edge_subspace is None else subspace_key(edge_subspace)}


def subspace_key(rows):
    r"""
    Return the ``key`` of a linear subspace of ``QQ^E`` spanned by ``rows``
    (lists indexed by edges): the function of a relabelling ``r`` that moves
    the column ``e`` to the edge ``r[2 * e] // 2`` (the sign of the half-edge
    is dropped) and returns the reduced echelon form over ``QQ`` as a tuple
    of tuples of its nonzero rows.
    """
    from sage.all import QQ, matrix

    def key(r):
        columns = len(rows[0])
        new = [[0] * columns for _ in rows]
        for i, row in enumerate(rows):
            for e, x in enumerate(row):
                new[i][r[2 * e] // 2] = x
        echelon = matrix(QQ, new).echelon_form()
        return tuple(tuple(row) for row in echelon.rows() if any(row))
    return key


def automorphism_count(m, half_edge_data=None, edge_data=None,
                       face_data=None, face_classes=None, edge_subspace=None):
    r"""
    Return the number of automorphisms of the decorated map ``m``.

    The decorations are ``None`` (absent) or:

    - ``half_edge_data``, ``face_data``: lists of length
      ``len(m.vertex_permutation())`` indexed by half-edges (any value on
      the inactive ones), preserved exactly by an automorphism;
    - ``edge_data``: a list indexed by edges, preserved exactly;
    - ``face_classes``: a list indexed by half-edges, preserved up to a
      relabelling of its values;
    - ``edge_subspace``: a list of rows indexed by edges spanning a
      subspace; an automorphism permutes the columns by the permutation it
      induces on edges (signs ignored) and must preserve the row space
      (matched through ``key``).
    """
    keywords = decoration_keywords(half_edge_data, edge_data, face_data,
                                   face_classes, edge_subspace)
    return len(m.automorphisms(**keywords))


def group_closure(gens, n):
    r"""
    Return the set of the elements (as tuples) of the group generated by the
    permutations ``gens`` of ``{0, ..., n-1}``.
    """
    identity = tuple(range(n))
    elements = {identity}
    todo = [identity]
    while todo:
        x = todo.pop()
        for g in gens:
            y = tuple(g[i] for i in x)
            if y not in elements:
                elements.add(y)
                todo.append(y)
    return elements


def automorphism_group_generated(m, gens, **decorations):
    r"""
    Return whether the automorphism group of ``m`` with ``decorations`` (as
    in :func:`automorphism_count`) is the group generated by ``gens``, a
    list of half-edge permutations given as lists.
    """
    keywords = decoration_keywords(**decorations)
    group = {tuple(g) for g in m.automorphisms(**keywords)}
    return group_closure(gens, len(m.vertex_permutation())) == group


def isomorphism(a, b, root_rule="all", **decorations):
    r"""
    Return a half-edge relabelling ``r`` (a list) such that ``a``
    relabelled by ``r`` is ``b`` with the decorations of ``a`` carried to
    those of ``b``, or ``None`` if there is none, searched under the rule
    of the roots ``root_rule``.

    Each decoration keyword (as in :func:`automorphism_count`) takes a pair
    ``(decoration of a, decoration of b)``, the subspaces included.
    """
    subspaces = decorations.pop('edge_subspace', None)
    keywords = decoration_keywords(**decorations)
    keywords = {key: value for key, value in keywords.items() if value is not None}
    if subspaces is not None:
        keywords['key'] = (subspace_key(subspaces[0]), subspace_key(subspaces[1]))
    ok, r = a.is_isomorphic(b, certificate=True, root_rule=root_rule, **keywords)
    return list(r) if ok else None


def canonical_form(m, root_rule="all", **decorations):
    r"""
    Return an object comparable with ``==`` that is equal for two decorated
    maps if and only if they are isomorphic, under the rule of the roots
    ``root_rule``.

    The decorations are as in :func:`automorphism_count`.
    """
    return m.canonical_form(root_rule=root_rule, **decoration_keywords(**decorations))


def tuple_automorphism_count(perms):
    r"""
    Return the number of permutations that commute with every permutation
    of ``perms`` (a list of permutations of ``{0, ..., n-1}`` given as
    lists), for tuples of permutations without an involution.
    """
    raise NotImplementedError


# Encodings

def map_from_fp(s):
    r"""
    Return the map with face permutation ``s``, a cycle string in edge
    notation (``i`` and ``~i`` are the two halves of edge ``i``), as veerer
    writes them.
    """
    return OrientedMap(fp=s)


def darts_perm(s, n=None):
    r"""
    Return the list of images of the dart-indexed cycle string ``s`` (plain
    integers, no ``~``), on ``{0, ..., n-1}`` with ``n`` one more than the
    largest dart by default.
    """
    cycles = []
    for c in s.replace(' ', '').split(')'):
        if not c:
            continue
        if c[0] != '(':
            raise ValueError('invalid cycle string %r' % s)
        cycles.append([int(x) for x in c[1:].split(',')])
    if n is None:
        n = 1 + max((max(c) for c in cycles), default=-1)
    p = list(range(n))
    for c in cycles:
        for i, x in enumerate(c):
            p[x] = c[(i + 1) % len(c)]
    return p


def map_from_darts(vp, fp):
    r"""
    Return the map of the fat graph with dart-indexed cycle strings ``vp``
    and ``fp``, as surface-dynamics writes them.

    Dart ``2e`` and ``2e+1`` form edge ``e``. The strings are converted to
    lists first: as strings, combisurf reads them in edge notation.
    """
    n = max(len(darts_perm(vp)), len(darts_perm(fp)))
    return OrientedMap(vp=darts_perm(vp, n), fp=darts_perm(fp, n))


def map_from_unicellular_word(w):
    r"""
    Return the one-face map of the word ``w`` on ``{0, ..., k-1}``, each
    letter appearing twice, with the dart labels of surface-dynamics'
    ``FatGraph.from_unicellular_word``: the first occurrence of letter
    ``e`` is dart ``2e``, the second ``2e+1``, and the face visits the
    darts in the order of the word.
    """
    n = len(w)
    fp = [-1] * n
    seen = [0] * (n // 2)
    previous = 2 * w[-1] + 1
    for k in w:
        if seen[k] == 0:
            h = 2 * k
        elif seen[k] == 1:
            h = 2 * k + 1
        else:
            raise ValueError('invalid unicellular word')
        seen[k] += 1
        fp[previous] = h
        previous = h
    if previous != 2 * w[-1] + 1:
        raise ValueError('invalid unicellular word')
    return OrientedMap(fp=fp)


def cylinder_diagram(s):
    r"""
    Return ``(m, colour, cylinder)`` for the cylinder diagram ``s`` written
    ``'(b...)-(t...) (b...)-(t...) ...'``.

    Separatrix ``i`` is edge ``i``. A bot cycle ``(b0, ..., bk)`` is the
    face ``2*b0 -> 2*b1 -> ... -> 2*b0`` and a top cycle ``(t0, ...)`` the
    face ``2*t0+1 -> 2*t1+1 -> ...``. ``colour[h]`` is ``0`` on even darts
    (bot faces) and ``1`` on odd darts (top faces); ``cylinder[h]`` is the
    index of the cylinder whose bot or top cycle contains the face of
    ``h``.
    """
    cylinders = []
    for c in s.split():
        bot, top = c.split('-')
        if bot[0] != '(' or bot[-1] != ')' or top[0] != '(' or top[-1] != ')':
            raise ValueError('invalid cylinder diagram %r' % s)
        bot = [int(x) for x in bot[1:-1].split(',')]
        top = [int(x) for x in top[1:-1].split(',')]
        cylinders.append((bot, top))
    nseps = 1 + max(max(b + t) for b, t in cylinders)
    n = 2 * nseps
    fp = [-1] * n
    cylinder = [-1] * n
    for i, (bot, top) in enumerate(cylinders):
        for cycle, shift in ((bot, 0), (top, 1)):
            for j, x in enumerate(cycle):
                h = 2 * x + shift
                if fp[h] != -1:
                    raise ValueError('invalid cylinder diagram %r' % s)
                fp[h] = 2 * cycle[(j + 1) % len(cycle)] + shift
                cylinder[h] = i
    if -1 in fp:
        raise ValueError('invalid cylinder diagram %r' % s)
    colour = [h % 2 for h in range(n)]
    return OrientedMap(fp=fp), colour, cylinder


def separatrix_to_half_edge(s):
    r"""
    Return the half-edge permutation ``h -> 2*s[h // 2] + (h % 2)`` of the
    permutation ``s`` of separatrices (edges).
    """
    return [2 * s[h // 2] + (h % 2) for h in range(2 * len(s))]


def relabel_map(m, r):
    r"""
    Return the map ``m`` relabelled by the half-edge relabelling ``r``: the
    image of ``r[h]`` is ``r[fp[h]]``.
    """
    fp = m.face_permutation()
    new_fp = [-1] * len(fp)
    for h, x in enumerate(fp):
        if x != -1:
            new_fp[r[h]] = r[x]
    return OrientedMap(fp=new_fp)


def random_relabelling(m, rng):
    r"""
    Return ``(r, sigma)`` where ``sigma`` is a random permutation of the
    edges of ``m`` and ``r`` the half-edge relabelling that sends edge
    ``e`` to ``sigma[e]``, swapping its two halves at random when it is not
    folded.
    """
    vp = m.vertex_permutation()
    ne = len(vp) // 2
    sigma = list(range(ne))
    rng.shuffle(sigma)
    r = [-1] * len(vp)
    for e in range(ne):
        swap = rng.randrange(2) if (vp[2 * e] != -1 and vp[2 * e + 1] != -1) else 0
        r[2 * e] = 2 * sigma[e] + swap
        r[2 * e + 1] = 2 * sigma[e] + 1 - swap
    return r, sigma


def carry_decorations(r, sigma, decorations):
    r"""
    Return the decorations carried by the half-edge relabelling ``r`` and
    the edge relabelling ``sigma`` (see :func:`automorphism_count` for the
    decorations).
    """
    carried = {}
    for key, value in decorations.items():
        if value is None:
            carried[key] = None
        elif key in ('half_edge_data', 'face_data', 'face_classes'):
            new = [None] * len(value)
            for h, x in enumerate(value):
                new[r[h]] = x
            carried[key] = new
        elif key == 'edge_data':
            new = [None] * len(value)
            for e, x in enumerate(value):
                new[sigma[e]] = x
            carried[key] = new
        elif key == 'edge_subspace':
            rows = []
            for row in value:
                new = [None] * len(row)
                for e, x in enumerate(row):
                    new[sigma[e]] = x
                rows.append(tuple(new))
            carried[key] = rows
        else:
            raise ValueError('unknown decoration %r' % key)
    return carried


# Data of the oracles

O1 = "(0,~1,2)(~0,1,~3)(4,~5,3)(~4,6,~2)(7,~6,8)(~7,5,~9)(10,~11,9)(~10,11,~8)"
O2 = "(0,8,~7)(1,9,~0)(2,10,~1)(3,11,~2)(4,12,~3)(5,13,~4)(6,14,~5)(7,15,~6)"
O3 = "(0,1,3)(2,4,~4)(~2,5,~5)(6,7,8)"
O4 = "(0,3,7,~6,~2,1)(2,5,~4,~3,~1,~0)(4,8,~5)(6,~8,~7)"
O5_EDGE_DATA = [2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 2, 1]  # "BRBBBRRBBBBR", R=1, B=2
O6 = "(0,1,2)(~0)(~1)(~2)"
O7 = "(0,1,2)(~0,~1,~2)"
O9_SUBSPACE = [(2, 0, 0, 0, 1, 0, 1, 0, 2),
               (0, 2, 0, 0, 0, 1, 0, 1, 2),
               (0, 0, 1, 1, 0, 0, 0, 0, 2)]

UNICELLULAR = {
    'O10': [0, 1, 0, 1, 2, 3, 2, 3],
    'O11': [0, 1, 0, 2, 1, 3, 2, 3],
    'O12': [0, 1, 0, 2, 3, 1, 2, 3],
    'O13': [0, 1, 2, 3, 0, 1, 2, 3],
    'O14': [0, 1, 2, 3, 0, 4, 1, 2, 3, 4],
}

DARTS = {
    'O15': ('(0,9,6,5,7,4,8,2,1,3)', '(0,2,1,3,8)(4,6,5,7,9)'),
    'O16': ('(0,5,4)(1,2,3)', '(0,3,1,4)(2)(5)'),
    'O17': ('(0,4,3)(5,2,1)', '(0,2,4,1,3,5)'),
    'O18': ('(0,2,1,3)', '(0,2,1,3)'),
    'O19': ('(0,8,6,4,3,7)(1,9,11,5,2,10)', '(0,10,9)(2,4,11)(1,7,8)(3,5,6)'),
}

# the triangulation with boundary of veerer/linear_family.py:745-747 as
# veerer stores it (_fp, _half_edges_data, _edges_data, _subspace)
O20_FP = [2, 3, 4, 6, 0, 13, 1, 12, 24, 14, 5, 7, 18, 16, 9, 20, 10, 26, 11,
          27, 22, 23, 15, 25, 8, 21, 28, 29, 17, 19]
O20_DATA = {
    'half_edge_data': [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0,
                       1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    'edge_data': [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2],
    'edge_subspace': [(1, 0, -1, -1, 0, 0, 0, -2, -1, -1, 2, 0, -2, 1, 0),
                      (0, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 2, 2, 0, 1),
                      (0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)],
}

O21 = "(0,2,3,1,5,6)(~0,~5,7,~8,~2,~6)(~1,4,~3)(~4,~7,8)"
O21_SUBSPACE = [(2, 1, 0, 0, 2, 0, 0, 0, 1),
                (0, 0, 1, 0, 2, 1, 0, 0, 0),
                (0, 0, 0, 1, 2, 0, 2, 1, 0)]

CD_I3_A = '(0,2,1)-(3,4,5) (3)-(1) (4)-(2) (5)-(0)'
CD_I3_B = '(0,2,1)-(3,5,4) (3)-(1) (4)-(2) (5)-(0)'
CD_I4_B = '(1,3,5)-(2,4,0) (2)-(5) (4)-(3) (0)-(1)'
CD_I5_B = '(4,2,3)-(1,5,0) (1)-(3) (0)-(2) (5)-(4)'
CD_I6_A = '(0,5)-(3,4) (1,4)-(2,5) (2)-(0) (3)-(1)'
CD_I6_B = '(0,5)-(3,4) (1,4)-(2,5) (2)-(1) (3)-(0)'
CD_I7 = ['(0,1)-(4,5) (2,4)-(0,3) (3,5)-(1,2)',
         '(0,1)-(2,3) (2,4)-(1,4) (3,5)-(0,5)',
         '(0,5)-(2,5) (1,4)-(0,4) (2,3)-(1,3)']
CD_I8 = '(0,1)-(0,2) (2,3)-(1,3)'


# Tests of the encodings

@pytest.mark.parametrize("fp, half_edges, edges, folded, connected", [
    (O1, 24, 12, 0, True),
    (O2, 24, 16, 8, True),
    (O3, 12, 9, 6, False),
    (O4, 18, 9, 0, True),
    (O6, 6, 3, 0, True),
    (O7, 6, 3, 0, True),
    (O21, 18, 9, 0, True),
], ids=["O1", "O2", "O3", "O4", "O6", "O7", "O21"])
def test_map_from_fp(fp, half_edges, edges, folded, connected):
    m = map_from_fp(fp)
    assert m.num_half_edges() == half_edges
    assert m.num_edges() == edges
    assert m.num_folded_edges() == folded
    assert m.is_connected() == connected


def test_map_from_fp_list_O20():
    m = OrientedMap(fp=O20_FP)
    assert m.num_half_edges() == 30
    assert m.num_edges() == 15
    assert not m.is_connected()
    assert len(O20_DATA['half_edge_data']) == 30
    assert len(O20_DATA['edge_data']) == 15
    assert all(len(row) == 15 for row in O20_DATA['edge_subspace'])


@pytest.mark.parametrize("name, half_edges", [
    ("O10", 8), ("O11", 8), ("O12", 8), ("O13", 8), ("O14", 10)])
def test_map_from_unicellular_word(name, half_edges):
    w = UNICELLULAR[name]
    m = map_from_unicellular_word(w)
    assert m.num_half_edges() == half_edges
    assert m.num_faces() == 1
    assert m.num_folded_edges() == 0
    # the face visits the darts in the order of the word
    fp = m.face_permutation()
    seen = set()
    darts = []
    for k in w:
        darts.append(2 * k + (k in seen))
        seen.add(k)
    for i, h in enumerate(darts):
        assert fp[h] == darts[(i + 1) % len(darts)]


def test_map_from_unicellular_word_invalid():
    with pytest.raises(ValueError):
        map_from_unicellular_word([0, 0, 0, 1])


@pytest.mark.parametrize("name, half_edges", [
    ("O15", 10), ("O16", 6), ("O17", 6), ("O18", 4), ("O19", 12)])
def test_map_from_darts(name, half_edges):
    vp, fp = DARTS[name]
    m = map_from_darts(vp, fp)
    assert m.num_half_edges() == half_edges
    assert m.num_folded_edges() == 0
    assert m.is_connected()
    assert list(m.vertex_permutation()) == darts_perm(vp)
    assert list(m.face_permutation()) == darts_perm(fp)


def test_darts_perm():
    assert darts_perm('(0,2,1)') == [2, 0, 1]
    assert darts_perm('(0,3)(1)', 5) == [3, 1, 2, 0, 4]
    assert darts_perm('(0,4)(1,5)(2,6)(3,7)(8,9)') == [4, 5, 6, 7, 0, 1, 2, 3, 9, 8]


@pytest.mark.parametrize("s, vertices, faces, connected", [
    (CD_I3_A, 2, 8, False),
    (CD_I3_B, 2, 8, False),
    (CD_I4_B, 2, 8, False),
    (CD_I5_B, 2, 8, False),
    (CD_I6_A, 2, 8, False),
    (CD_I6_B, 2, 8, False),
    (CD_I7[0], 2, 6, True),
    (CD_I7[1], 2, 6, True),
    (CD_I7[2], 2, 6, True),
    (CD_I8, 2, 4, True),
], ids=["I3a", "I3b", "I4b", "I5b", "I6a", "I6b", "I7a", "I7b", "I7c", "I8"])
def test_cylinder_diagram(s, vertices, faces, connected):
    m, colour, cylinder = cylinder_diagram(s)
    assert m.num_half_edges() == (8 if s == CD_I8 else 12)
    assert m.num_folded_edges() == 0
    assert m.num_vertices() == vertices
    assert m.num_faces() == faces
    assert m.is_connected() == connected
    assert len(colour) == len(cylinder) == m.num_half_edges()
    # every face has one colour and one cylinder, and every cylinder one bot
    # face and one top face
    faces_of_cylinder = {}
    for face in m.faces():
        assert len({colour[h] for h in face}) == 1
        assert len({cylinder[h] for h in face}) == 1
        faces_of_cylinder.setdefault(cylinder[face[0]], []).append(colour[face[0]])
    assert all(sorted(c) == [0, 1] for c in faces_of_cylinder.values())


def test_separatrix_to_half_edge():
    assert separatrix_to_half_edge([3, 2, 1, 0]) == [6, 7, 4, 5, 2, 3, 0, 1]


def test_random_relabelling():
    rng = random.Random(0)
    for fp in (O1, O3):
        m = map_from_fp(fp)
        r, sigma = random_relabelling(m, rng)
        mm = relabel_map(m, r)
        assert mm.num_half_edges() == m.num_half_edges()
        assert mm.num_folded_edges() == m.num_folded_edges()
        assert sorted(mm.face_profile()) == sorted(m.face_profile())
        assert sorted(mm.vertex_profile()) == sorted(m.vertex_profile())


# Automorphism counts, undecorated maps

@pytest.mark.parametrize("build, expected", [
    (lambda: map_from_fp(O1), 4),   # veerer/constellation.py:1769-1773
    (lambda: map_from_fp(O2), 8),   # veerer/constellation.py:1777-1779
    # veerer/constellation.py:1813-1815, disconnected
    (lambda: map_from_fp(O3), 36),
    (lambda: map_from_fp(O4), 2),   # veerer/constellation.py:1805-1807
    # surface-dynamics/fat_graph.py:2968-2982
    (lambda: map_from_unicellular_word(UNICELLULAR['O10']), 2),
    (lambda: map_from_unicellular_word(UNICELLULAR['O11']), 1),
    (lambda: map_from_unicellular_word(UNICELLULAR['O12']), 1),
    (lambda: map_from_unicellular_word(UNICELLULAR['O13']), 8),
    # surface-dynamics/fat_graph.py:2987-2989
    (lambda: map_from_unicellular_word(UNICELLULAR['O14']), 2),
    (lambda: map_from_darts(*DARTS['O16']), 2),  # surface-dynamics/fat_graph.py:3001-3005
    (lambda: map_from_darts(*DARTS['O17']), 6),  # surface-dynamics/fat_graph.py:3006-3011
    (lambda: map_from_darts(*DARTS['O18']), 4),  # surface-dynamics/fat_graph.py:3013-3015
    (lambda: map_from_darts(*DARTS['O19']), 4),  # surface-dynamics/fat_graph.py:2868-2880
], ids=["O1", "O2", "O3", "O4", "O10", "O11", "O12", "O13", "O14", "O16",
        "O17", "O18", "O19"])
def test_automorphism_count(build, expected):
    assert automorphism_count(build()) == expected


def regular_face_permutation(G, vg, eg, fg):
    r"""
    Return the face permutation of the regular map of the group ``G`` with
    vertex automorphism ``vg``, edge automorphism ``eg`` (an involution)
    and face automorphism ``fg = (vg * eg)^-1``.

    Half-edges are the elements of ``G`` in the order of ``list(G)``,
    acted on by right multiplication; the result is relabelled so that the
    ``c``-th cycle of the edge permutation (right multiplication by
    ``eg``) becomes the edge ``(2c, 2c+1)``. The left multiplication by
    ``G`` commutes with these right multiplications: it acts on the map by
    automorphisms, transitively on the half-edges.
    """
    from combisurf.permutation import perm_compose, perm_conjugate, perm_cycles, perm_is_one

    elements = list(G)
    n = len(elements)
    to_int = {g: i for i, g in enumerate(elements)}
    vgp = array('i', [to_int[elements[i] * vg] for i in range(n)])
    egp = array('i', [to_int[elements[i] * eg] for i in range(n)])
    fgp = array('i', [to_int[elements[i] * fg] for i in range(n)])
    if not perm_is_one(perm_compose(perm_compose(vgp, egp), fgp)):
        raise RuntimeError('vg * eg * fg is not the identity')
    r = array('i', [-1] * n)
    for c, (i, j) in enumerate(perm_cycles(egp)):
        r[i] = 2 * c
        r[j] = 2 * c + 1
    return perm_conjugate(fgp, r)


def regular_map_S5():
    r"""
    Return the regular map of ``S_5`` with 120 half-edges (a quadrangulation
    of a surface of genus 6).
    """
    from sage.all import SymmetricGroup

    S = SymmetricGroup(5)
    vg = S("(1,2,4)(3,5)")
    eg = S("(1,2)(3,4)")
    fg = (vg * eg).inverse()
    return OrientedMap(fp=list(regular_face_permutation(S, vg, eg, fg)))


def regular_map_PSL27():
    r"""
    Return the regular triangulation of genus 3 of ``PSL(2, 7)`` with 168
    half-edges (Klein's map: vertices of degree 7).
    """
    from sage.all import PSL

    G = PSL(2, 7)
    L = list(G)
    rng = random.Random(0)
    eg = fg = None
    while (eg is None or eg.order() != 2 or fg.order() != 3
           or (eg * fg).order() != 7 or G.subgroup([eg, fg]).order() != len(L)):
        eg = rng.choice(L)
        fg = rng.choice(L)
    vg = (eg * fg).inverse()
    return OrientedMap(fp=list(regular_face_permutation(G, vg, eg, fg)))


@pytest.mark.parametrize("build, expected", [
    (regular_map_S5, 120),      # the order of S_5
    (regular_map_PSL27, 168),   # the order of PSL(2, 7)
], ids=["O22", "O23"])
def test_automorphism_count_regular_map(build, expected):
    # the left multiplication by the group is a family of automorphisms
    # that is regular on the half-edges
    m = build()
    assert m.num_half_edges() == expected
    automorphisms = m.automorphisms()
    assert len(automorphisms) == expected
    h0 = next(m.half_edges())
    assert sorted(g[h0] for g in automorphisms) == sorted(m.half_edges())


@pytest.mark.parametrize("build", [regular_map_S5, regular_map_PSL27], ids=["O22", "O23"])
def test_fast_generators_regular_map(build):
    fp = list(build().face_permutation())
    fast = isomorphism_fast.automorphism_generators(fp)
    assert fast == isomorphism_naive.automorphism_generators(fp)
    assert 0 < len(fast) <= 7


def test_fast_generators_connected_without_closure(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("_group_closure called on a connected map")

    monkeypatch.setattr(isomorphism_naive, "_group_closure", forbidden)
    triangle = [2, 5, 4, 1, 0, 3]
    assert isomorphism_fast.automorphism_generators(triangle) == [
        array('i', [1, 0, 5, 4, 3, 2]), array('i', [2, 3, 4, 5, 0, 1])]
    assert isomorphism_fast.automorphism_generators([2, -1, 0, 3]) == []
    assert isomorphism_fast.automorphism_generators([]) == []
    assert len(isomorphism_fast.automorphism_generators(list(range(2)))) == 1
    for build in (regular_map_S5, regular_map_PSL27):
        fp = list(build().face_permutation())
        assert len(isomorphism_fast.automorphism_generators(fp)) > 0


def face_indices(m):
    r"""
    Return the list that gives to each half-edge of ``m`` the index of its
    face (distinct per face).
    """
    data = [-1] * len(m.vertex_permutation())
    for i, face in enumerate(m.faces()):
        for h in face:
            data[h] = i
    return data


@pytest.mark.parametrize("name, expected", [
    ("O16", 1),  # surface-dynamics/fat_graph.py:3001-3005
    ("O17", 6),  # surface-dynamics/fat_graph.py:3006-3011
], ids=["O16", "O17"])
def test_automorphism_count_faces_fixed(name, expected):
    m = map_from_darts(*DARTS[name])
    assert automorphism_count(m, face_data=face_indices(m)) == expected


# Automorphism groups

def test_automorphism_group_O15():
    # surface-dynamics/fat_graph.py:2993-2997
    m = map_from_darts(*DARTS['O15'])
    g = darts_perm('(0,4)(1,5)(2,6)(3,7)(8,9)', 10)
    assert automorphism_group_generated(m, [g])


def test_automorphism_group_I8():
    # surface-dynamics/separatrix_diagram.py:3058-3062
    m, colour, cylinder = cylinder_diagram(CD_I8)
    g = separatrix_to_half_edge([3, 2, 1, 0])
    assert automorphism_group_generated(m, [g], face_data=colour,
                                        face_classes=cylinder)


# Automorphism counts, decorated maps

@pytest.mark.parametrize("build, decorations, expected", [
    # veerer/constellation.py:1782-1788
    (lambda: map_from_fp(O1), {'edge_data': O5_EDGE_DATA}, 4),
    # veerer/constellation.py:1792-1794
    (lambda: map_from_fp(O6), {'half_edge_data': [0, 1, 0, 1, 0, 1]}, 3),
    # veerer/constellation.py:1795-1797
    (lambda: map_from_fp(O7), {'half_edge_data': [0, 1, 0, 1, 0, 1]}, 3),
    # veerer/constellation.py:1798-1800
    (lambda: map_from_fp(O7), {'half_edge_data': [0, 1, 0, 1, 0, 2]}, 1),
    # veerer/constellation.py:1805-1809
    (lambda: map_from_fp(O4), {'edge_subspace': O9_SUBSPACE}, 2),
    # veerer/linear_family.py:745-747, over the disconnected map
    (lambda: OrientedMap(fp=O20_FP), O20_DATA, 8),
    # veerer/linear_family.py:749-751
    (lambda: map_from_fp(O21), {'edge_subspace': O21_SUBSPACE}, 2),
], ids=["O5", "O6", "O7", "O8", "O9", "O20", "O21"])
def test_automorphism_count_decorated(build, decorations, expected):
    assert automorphism_count(build(), **decorations) == expected


def origami_map(r, u):
    r"""
    Return the map of the origami ``(r, u)`` and its half-edge decoration
    ``h % 4``.
    """
    n = len(r)
    return (OrientedMap(fp=origami_face_permutation(r, u)),
            origami_decoration(n))


# The decorated count is the number of translations of the origami:
# surface-dynamics/origami_dense.pyx:2838-2850 for the first two (in origami
# order), and the n translations of a torus of n squares for the last two.
# The undecorated counts (map automorphisms may rotate the squares) were
# measured with the naive implementation and are not independent.
@pytest.mark.parametrize("r, u, undecorated, decorated", [
    ([1, 0, 2], [2, 1, 0], 4, 1),
    ([1, 0, 3, 2], [0, 2, 1, 3], 4, 2),
    ([1, 2, 3, 4, 5, 0], [0, 1, 2, 3, 4, 5], None, 6),
    ([1, 2, 0, 4, 5, 3], [3, 4, 5, 1, 2, 0], None, 6),
], ids=["O24", "O25", "O26-cylinder", "O26-sheared"])
def test_automorphism_count_origami(r, u, undecorated, decorated):
    m, d = origami_map(r, u)
    assert automorphism_count(m, half_edge_data=d) == decorated
    if undecorated is not None:
        assert automorphism_count(m) == undecorated


def half_edge_array_of_edges(m, values):
    r"""
    Return the half-edge array of the edge array ``values``.
    """
    return [values[h // 2] for h in range(m.num_half_edges())]


BRUTE_FORCE_CASES = [
    (lambda: map_from_fp(O1), {}),
    (lambda: map_from_fp(O2), {}),
    (lambda: map_from_fp(O3), {}),
    (lambda: map_from_fp(O4), {}),
    (lambda: map_from_unicellular_word(UNICELLULAR['O10']), {}),
    (lambda: map_from_unicellular_word(UNICELLULAR['O11']), {}),
    (lambda: map_from_unicellular_word(UNICELLULAR['O12']), {}),
    (lambda: map_from_unicellular_word(UNICELLULAR['O13']), {}),
    (lambda: map_from_unicellular_word(UNICELLULAR['O14']), {}),
    (lambda: map_from_darts(*DARTS['O16']), {}),
    (lambda: map_from_darts(*DARTS['O17']), {}),
    (lambda: map_from_darts(*DARTS['O18']), {}),
    (lambda: map_from_darts(*DARTS['O19']), {}),
    (lambda: map_from_fp(O1), {'edge_data': O5_EDGE_DATA}),
    (lambda: map_from_fp(O6), {'half_edge_data': [0, 1, 0, 1, 0, 1]}),
    (lambda: map_from_fp(O7), {'half_edge_data': [0, 1, 0, 1, 0, 1]}),
    (lambda: map_from_fp(O7), {'half_edge_data': [0, 1, 0, 1, 0, 2]}),
    (regular_map_S5, {}),
    (regular_map_PSL27, {}),
    (lambda: origami_map([1, 0, 2], [2, 1, 0])[0], {}),
    (lambda: origami_map([1, 0, 3, 2], [0, 2, 1, 3])[0], {}),
    (lambda: origami_map([1, 2, 3, 4, 5, 0], [0, 1, 2, 3, 4, 5])[0], {}),
    (lambda: origami_map([1, 2, 0, 4, 5, 3], [3, 4, 5, 1, 2, 0])[0], {}),
    (lambda: origami_map([1, 0, 2], [2, 1, 0])[0], 'origami'),
    (lambda: origami_map([1, 0, 3, 2], [0, 2, 1, 3])[0], 'origami'),
    (lambda: origami_map([1, 2, 3, 4, 5, 0], [0, 1, 2, 3, 4, 5])[0], 'origami'),
    (lambda: origami_map([1, 2, 0, 4, 5, 3], [3, 4, 5, 1, 2, 0])[0], 'origami'),
]


@pytest.mark.parametrize("build, decorations", BRUTE_FORCE_CASES,
                         ids=list(range(len(BRUTE_FORCE_CASES))))
def test_brute_force_count_agrees(build, decorations):
    m = build()
    if decorations == 'origami':
        decorations = {'half_edge_data': origami_decoration(m.num_half_edges() // 4)}
    fp = list(m.face_permutation())
    data = []
    if 'half_edge_data' in decorations:
        data.append(decorations['half_edge_data'])
    if 'edge_data' in decorations:
        data.append(half_edge_array_of_edges(m, decorations['edge_data']))
    expected = automorphism_count(m, **decorations)
    assert brute_force_automorphism_count(fp, tuple(data)) == expected
    # the same group under every rule of the roots
    group = {tuple(g) for g in m.automorphisms(**decoration_keywords(**decorations))}
    for root_rule in ROOT_RULE_NAMES[1:]:
        automorphisms = m.automorphisms(root_rule=root_rule, **decoration_keywords(**decorations))
        assert len(automorphisms) == expected
        assert {tuple(g) for g in automorphisms} == group


# Isomorphism verdicts

def check_relabelling_between(a, b, r):
    r"""
    Check that the half-edge relabelling ``r`` restricted to the active
    half-edges of ``a`` is a bijection onto those of ``b`` that respects
    twins and sends folded edges to folded edges.

    :func:`relabel_map` conjugates the face permutation only, and relies on
    these properties.
    """
    vpa = a.vertex_permutation()
    vpb = b.vertex_permutation()
    active_a = [h for h in range(len(vpa)) if vpa[h] != -1]
    active_b = [h for h in range(len(vpb)) if vpb[h] != -1]
    image = [r[h] for h in active_a]
    assert len(set(image)) == len(image)
    assert sorted(image) == active_b
    for h in active_a:
        folded_a = vpa[h ^ 1] == -1
        folded_b = vpb[r[h] ^ 1] == -1
        assert folded_a == folded_b
        if not folded_a:
            assert r[h ^ 1] == r[h] ^ 1


def check_isomorphism(a, b, expected, colours=None, cylinders=None, root_rule="all"):
    r"""
    Check the verdict of :func:`isomorphism` on ``a`` and ``b`` under the
    rule of the roots ``root_rule``; when they are isomorphic, check that
    the relabelling sends ``a`` to ``b`` and carries the colours exactly and
    the cylinders up to relabelling.
    """
    decorations = {}
    if colours is not None:
        decorations['face_data'] = colours
    if cylinders is not None:
        decorations['face_classes'] = cylinders
    r = isomorphism(a, b, root_rule=root_rule, **decorations)
    if not expected:
        assert r is None
        return
    assert r is not None
    check_relabelling_between(a, b, r)
    assert relabel_map(a, r) == b
    if colours is not None:
        ca, cb = colours
        assert all(cb[r[h]] == ca[h] for h in a.half_edges())
    if cylinders is not None:
        ca, cb = cylinders
        matching = {}
        for h in a.half_edges():
            assert matching.setdefault(ca[h], cb[r[h]]) == cb[r[h]]
        assert len(set(matching.values())) == len(matching)


@pytest.mark.parametrize("fa, fb, expected", [
    ("(0,1,2)(~0,~2,~1)", "(0,2,1)(~0,~1,~2)", True),   # veerer/constellation.py:2165-2169
    ("(0,1,2)(~0,~2,~1)", "(0,1,2)(~0,~1,~2)", False),  # veerer/constellation.py:2167-2171
], ids=["I1", "I2"])
def test_isomorphism_fp(fa, fb, expected):
    check_isomorphism(map_from_fp(fa), map_from_fp(fb), expected)


@WITH_ROOT_RULE_NAMES
@pytest.mark.parametrize("sa, sb, expected", [
    # the maps of I3-I6 are disconnected, their components tied together
    # by the cylinders (face classes)
    # surface-dynamics/separatrix_diagram.py:2376-2379
    (CD_I3_A, CD_I3_B, False),
    # surface-dynamics/separatrix_diagram.py:2381-2385
    (CD_I3_A, CD_I4_B, True),
    # surface-dynamics/separatrix_diagram.py:2387-2391
    (CD_I3_B, CD_I5_B, True),
    # surface-dynamics/separatrix_diagram.py:2393-2396
    (CD_I6_A, CD_I6_B, False),
    (CD_I7[0], CD_I7[1], False),  # surface-dynamics/separatrix_diagram.py:2398-2406
    (CD_I7[0], CD_I7[2], False),
    (CD_I7[1], CD_I7[2], False),
], ids=["I3", "I4", "I5", "I6", "I7ab", "I7ac", "I7bc"])
def test_isomorphism_cylinder_diagram(sa, sb, expected, root_rule):
    a, colour_a, cylinder_a = cylinder_diagram(sa)
    b, colour_b, cylinder_b = cylinder_diagram(sb)
    check_isomorphism(a, b, expected, colours=(colour_a, colour_b),
                      cylinders=(cylinder_a, cylinder_b), root_rule=root_rule)


@pytest.mark.parametrize("s", CD_I7 + [CD_I8], ids=["I7a", "I7b", "I7c", "I8"])
def test_isomorphism_cylinder_diagram_relabelled(s):
    # a cylinder diagram is isomorphic to its relabellings, with the
    # cylinders renumbered
    a, colour, cylinder = cylinder_diagram(s)
    rng = random.Random(0)
    for _ in range(3):
        r, sigma = random_relabelling(a, rng)
        b = relabel_map(a, r)
        carried = carry_decorations(r, sigma, {'face_data': colour, 'face_classes': cylinder})
        names = list(range(len(set(cylinder))))
        rng.shuffle(names)
        renamed = [names[c] for c in carried['face_classes']]
        check_isomorphism(a, b, True, colours=(colour, carried['face_data']),
                          cylinders=(cylinder, renamed))


# Canonical forms: invariance under relabelling, the property encoded by the
# iso_sig doctests of veerer/constellation.py:2075-2133

@pytest.mark.parametrize("build, decorations", [
    (lambda: map_from_fp(O1), {}),
    (lambda: map_from_fp(O2), {}),
    (lambda: map_from_fp(O3), {}),
    (lambda: map_from_fp(O4), {}),
    (lambda: map_from_fp(O1), {'edge_data': O5_EDGE_DATA}),
    (lambda: OrientedMap(fp=O20_FP), O20_DATA),
], ids=["O1", "O2", "O3", "O4", "O5", "O20"])
def test_canonical_form_invariance(build, decorations):
    m = build()
    c = canonical_form(m, **decorations)
    rng = random.Random(0)
    for _ in range(3):
        r, sigma = random_relabelling(m, rng)
        mm = relabel_map(m, r)
        assert canonical_form(mm, **carry_decorations(r, sigma, decorations)) == c


@pytest.mark.parametrize("root_rule", ROOT_RULE_NAMES[1:])
def test_canonical_form_invariance_root_rule(root_rule):
    # O1-O5 under the other rules; O3 has folded edges and two non-folded
    # edges whose invariants differ, so that the classes differ by rule
    for build, decorations in [(lambda: map_from_fp(O1), {}), (lambda: map_from_fp(O2), {}),
                               (lambda: map_from_fp(O3), {}), (lambda: map_from_fp(O4), {}),
                               (lambda: map_from_fp(O1), {'edge_data': O5_EDGE_DATA})]:
        m = build()
        c = canonical_form(m, root_rule=root_rule, **decorations)
        rng = random.Random(0)
        for _ in range(3):
            r, sigma = random_relabelling(m, rng)
            mm = relabel_map(m, r)
            assert canonical_form(mm, root_rule=root_rule, **carry_decorations(r, sigma, decorations)) == c
    m = map_from_fp(O3)
    assert m.is_isomorphic(relabel_map(m, random_relabelling(m, random.Random(1))[0]), root_rule=root_rule)


# Tuples of permutations: surface-dynamics origami_dense.pyx:2838-2850,
# Origami('(1,2)','(1,3)') and Origami('(1,2)(3,4)','(2,3)'), 0-based here

@NOT_IMPLEMENTED
@pytest.mark.parametrize("perms, expected", [
    ([[1, 0, 2], [2, 1, 0]], 1),
    ([[1, 0, 3, 2], [0, 2, 1, 3]], 2),
], ids=["origami-3", "origami-4"])
def test_tuple_automorphism_count(perms, expected):
    assert tuple_automorphism_count(perms) == expected


# Sizes of collections up to isomorphism, explored by local mutations

def test_explore_flip_graph():
    # veerer/automata.py, triangulation_flip_graph_small: 236 states
    seed = "(0,1,2)(~0,~3,~8)(3,5,4)(~4,~1,~5)(6,7,8)(~6,9,~2)"
    start = list(OrientedMap(fp=seed)._fp)
    assert explore(start, flip_neighbours, combisurf_key) == 236


def test_explore_reduced_core_automaton():
    # veerer/automata.py, core_automaton_small: 234 states
    seed = "(0,8,~7)(~0,~6,7)(1,11,~10)(~1,~11,4)(2,10,~9)(~2,~4,5)(3,9,~8)(~3,~5,6)"
    colours = parse_colours("RRRRBBBBBBBB")
    assert colours[:4] == [RED] * 4 and colours[4:] == [BLUE] * 8
    start = reduced_core_seed(list(OrientedMap(fp=seed)._fp), colours)
    count = explore(start, lambda s: reduced_core_neighbours(*s),
                    lambda s: combisurf_key(s[0], edge_data=s[1]))
    assert count == 234


def test_explore_gl2z_orbit_O31():
    # surface-dynamics automata.py, gl2z_orbit_small: the orbit of
    # Origami("(1,2,...,49)", "(1,50)") has 32400 states; each state has two
    # images, hence 2 * 32400 + 1 keys computed
    n = 50
    r = list(range(1, n - 1)) + [0, n - 1]
    u = list(range(n))
    u[0], u[n - 1] = n - 1, 0
    decoration = (origami_decoration(n),)
    calls = 0

    def key(state):
        nonlocal calls
        calls += 1
        return isomorphism_fast.canonical_key(origami_face_permutation(*state), decoration)

    assert explore((r, u), gl2z_neighbours, key) == 32400
    assert calls == 64801


def test_origami_mutations():
    r, u = [1, 0, 2], [2, 1, 0]
    for image in (origami_l(r, u), origami_i(r, u)):
        assert len(image) == 2
        for p in image:
            assert sorted(p) == list(range(len(r)))
    assert len(gl2z_neighbours((r, u))) == 2


# Coloured maps


@pytest.mark.parametrize("call", [
    lambda m: m.canonical_labelling(),
    lambda m: m.canonical_form(),
    lambda m: m.automorphisms(),
    lambda m: m.is_isomorphic(m),
    lambda m: m.canonical_roots(),
    lambda m: m.num_automorphisms(),
    lambda m: m.automorphism_vertex_action([0, 1, 2, 3, 4, 5]),
    lambda m: m.automorphism_face_action([0, 1, 2, 3, 4, 5]),
    lambda m: m.automorphism_generators(),
    lambda m: m.canonical_labelling(root_rule="all"),
    lambda m: m.canonical_labelling(half_edge_classes=[0] * 6),
    lambda m: m.canonical_roots(key=len),
    lambda m: m.automorphism_generators(root_rule="all"),
], ids=["canonical_labelling", "canonical_form", "automorphisms",
        "is_isomorphic", "canonical_roots", "num_automorphisms",
        "automorphism_vertex_action", "automorphism_face_action",
        "automorphism_generators", "root_rule", "half_edge_classes",
        "key", "generators_root_rule"])
def test_coloured_map_methods(call):
    r"""
    Coloured maps raise NotImplementedError on the isomorphism methods.
    """
    from combisurf import ColoredOrientedMap
    m = ColoredOrientedMap(fp="(0,1,2)(~0,~2,~1)", ecolors=[0, 1, 1])
    with pytest.raises(NotImplementedError):
        call(m)


def test_oriented_map_isomorphic_with_coloured():
    r"""
    OrientedMap.is_isomorphic raises NotImplementedError when other is a coloured map.
    """
    from combisurf import ColoredOrientedMap
    m = OrientedMap(fp="(0,1,2)(~0,~2,~1)")
    m_coloured = ColoredOrientedMap(fp="(0,1,2)(~0,~2,~1)", ecolors=[0, 1, 1])
    with pytest.raises(NotImplementedError):
        m.is_isomorphic(m_coloured)


# Properties of the reference implementation on random maps, connected and
# disconnected

MAX_ATTEMPTS = 1000


def is_connected_fp(fp):
    r"""
    Return whether the face permutation ``fp`` (``-1`` on inactive
    half-edges) is connected.
    """
    active = [h for h in range(len(fp)) if fp[h] != -1]
    if not active:
        return True
    ne = sum(1 for e in range(len(fp) // 2) if fp[2 * e] != -1 or fp[2 * e + 1] != -1)
    _, k = isomorphism_naive.relabelling_from(fp, active[0])
    return k == 2 * ne


def random_fp(rng, ne, folded_probability):
    r"""
    Return a random face permutation with ``ne`` edges, each folded with
    probability ``folded_probability``.
    """
    active = []
    for e in range(ne):
        active.append(2 * e)
        if rng.random() >= folded_probability:
            active.append(2 * e + 1)
    images = active[:]
    rng.shuffle(images)
    fp = [-1] * (2 * ne)
    for h, x in zip(active, images):
        fp[h] = x
    return fp


def cyclic_cover(fp, k, shifts):
    r"""
    Return the face permutation of the ``k``-fold cyclic cover of ``fp``
    in which the half-edge ``h`` of the sheet ``i`` goes to ``fp[h]`` in
    the sheet ``i + shifts[h]``; the edge ``e`` of the sheet ``i`` is the
    edge ``k * e + i``. The rotation of the sheets is an automorphism.
    """
    n = len(fp)
    cover = [-1] * (k * n)
    for h in range(n):
        if fp[h] == -1:
            continue
        for i in range(k):
            j = (i + shifts[h]) % k
            cover[2 * (k * (h // 2) + i) + h % 2] = 2 * (k * (fp[h] // 2) + j) + fp[h] % 2
    return cover


def insert_gap(fp, g):
    r"""
    Return ``fp`` with an inactive edge inserted before the edge ``g``.
    """
    new = [-1] * (len(fp) + 2)
    for h, x in enumerate(fp):
        if x != -1:
            new[h if h < 2 * g else h + 2] = x if x < 2 * g else x + 2
    return new


def random_connected_fp(rng, max_half_edges=40):
    r"""
    Return a random connected face permutation of at most
    ``max_half_edges`` half-edge labels, sometimes with folded edges,
    sometimes a cyclic cover (with automorphisms), sometimes with a label
    gap that is not the last edge.
    """
    for _ in range(MAX_ATTEMPTS):
        k = rng.choice([1, 1, 2, 3, 4])
        ne = rng.randint(1, max(1, (max_half_edges - 2) // (2 * k)))
        fp = random_fp(rng, ne, rng.choice([0, 0.2, 0.5]))
        if k > 1:
            fp = cyclic_cover(fp, k, [rng.randrange(k) for _ in fp])
        if not is_connected_fp(fp):
            continue
        if len(fp) >= 4 and rng.random() < 0.5:
            fp = insert_gap(fp, rng.randrange(len(fp) // 2))
        assert len(fp) <= max_half_edges
        return fp
    raise RuntimeError("no connected map found in {} attempts".format(MAX_ATTEMPTS))


def random_disconnected_fp(rng):
    r"""
    Return the face permutation of the disjoint union of one to three
    random connected maps of at most 16 half-edge labels each (at most 48
    in all), a component being sometimes a relabelled copy of an earlier
    one, with the labels of the union shuffled by edges.
    """
    components = []
    for _ in range(rng.randint(1, 3)):
        if components and rng.random() < 0.5:
            fp = rng.choice(components)
            fp = relabel_arrays(fp, (), random_edge_relabelling(rng, fp))[0]
        else:
            fp = random_connected_fp(rng, max_half_edges=16)
        components.append(fp)
    union = []
    for fp in components:
        shift = len(union)
        union.extend(-1 if x == -1 else x + shift for x in fp)
    assert len(union) <= 48
    return relabel_arrays(union, (), random_edge_relabelling(rng, union))[0]


def random_data(rng, fp):
    r"""
    Return a random tuple of zero, one or two half-edge arrays for ``fp``,
    with few values so that automorphisms may survive.
    """
    return tuple([rng.randrange(2) for _ in fp] for _ in range(rng.randrange(3)))


def random_classes(rng, fp):
    r"""
    Return a random tuple of zero, one or two half-edge arrays for ``fp``,
    with values in ``range(3)``, to be matched up to relabelling.
    """
    return tuple([rng.randrange(3) for _ in fp] for _ in range(rng.randrange(3)))


def rename_values(rng, classes):
    r"""
    Return the arrays of ``classes`` (values in ``range(m)``, ``m`` at least
    3), each with its values renamed by its own random permutation of
    ``range(m)``.
    """
    renamed = []
    for c in classes:
        names = list(range(max(3, max(c, default=0) + 1)))
        rng.shuffle(names)
        renamed.append([names[x] for x in c])
    return tuple(renamed)


def check_classes_carried(fp, c1, c2, g):
    r"""
    Check that ``c1[h] -> c2[g[h]]`` on the active half-edges ``h`` of
    ``fp`` is a well-defined bijection of values.
    """
    forward = {}
    for h in range(len(fp)):
        if fp[h] != -1:
            assert forward.setdefault(c1[h], c2[g[h]]) == c2[g[h]]
    assert len(set(forward.values())) == len(forward)


def random_edge_relabelling(rng, fp):
    r"""
    Return a random twin-compatible relabelling of ``fp`` (a permutation of
    ``range(len(fp))``) that swaps the halves of non-folded edges at random
    and keeps the last edge active.
    """
    ne = len(fp) // 2
    for _ in range(MAX_ATTEMPTS):
        sigma = list(range(ne))
        rng.shuffle(sigma)
        if fp[2 * sigma.index(ne - 1)] != -1:
            break
    else:
        raise RuntimeError("no relabelling found in {} attempts".format(MAX_ATTEMPTS))
    r = [-1] * len(fp)
    for e in range(ne):
        swap = rng.randrange(2) if (fp[2 * e] != -1 and fp[2 * e + 1] != -1) else 0
        r[2 * e] = 2 * sigma[e] + swap
        r[2 * e + 1] = 2 * sigma[e] + 1 - swap
    return r


def relabel_arrays(fp, data, r):
    r"""
    Return ``fp`` and ``data`` relabelled by ``r``.
    """
    new_fp = [-1] * len(fp)
    for h, x in enumerate(fp):
        if x != -1:
            new_fp[r[h]] = r[x]
    new_data = []
    for d in data:
        new = [None] * len(d)
        for h, x in enumerate(d):
            new[r[h]] = x
        new_data.append(new)
    return new_fp, tuple(new_data)


def check_automorphism(fp, data, g):
    r"""
    Check that ``g`` is an automorphism of ``fp`` decorated by ``data``.
    """
    n = len(fp)
    assert sorted(g) == list(range(n))
    for h in range(n):
        if fp[h] == -1 and fp[h ^ 1] == -1:
            assert g[h] == h
        elif fp[h] != -1:
            assert g[fp[h]] == fp[g[h]]
            if fp[h ^ 1] != -1:
                assert g[h ^ 1] == g[h] ^ 1
            else:
                assert fp[g[h] ^ 1] == -1
            for d in data:
                assert d[g[h]] == d[h]


def propagate(fp, data, classes, state, h0, x0):
    r"""
    Return the partial map ``state`` extended by ``h0 -> x0`` and every image
    it forces, or ``None`` on a contradiction; ``state`` is not modified.

    A state is a triple ``(g, inverse, value_maps)``: the partial map of
    half-edges as a dictionary, its inverse, and for each array of
    ``classes`` the value map ``c[h] -> c[g[h]]`` and its inverse, shared by
    all the components. The images are forced along ``fp`` and the twins
    of non-folded edges (a folded edge must go to a folded edge), so the
    whole component of ``h0`` gets mapped, injectively, onto the component
    of ``x0``. Each image keeps every array of ``data``, and the value maps
    must stay well defined: they are then bijections.
    """
    g = dict(state[0])
    inverse = dict(state[1])
    value_maps = [(dict(forward), dict(backward)) for forward, backward in state[2]]
    todo = [(h0, x0)]
    while todo:
        h, x = todo.pop()
        if h in g:
            if g[h] != x:
                return None
            continue
        if x in inverse or fp[x] == -1:
            return None
        if any(d[h] != d[x] for d in data):
            return None
        for c, (forward, backward) in zip(classes, value_maps):
            if forward.setdefault(c[h], c[x]) != c[x]:
                return None
            if backward.setdefault(c[x], c[h]) != c[h]:
                return None
        g[h] = x
        inverse[x] = h
        todo.append((fp[h], fp[x]))
        if fp[h ^ 1] != -1:
            if fp[x ^ 1] == -1:
                return None
            todo.append((h ^ 1, x ^ 1))
        elif fp[x ^ 1] != -1:
            return None
    return g, inverse, value_maps


def count_extensions(fp, data, classes, state):
    r"""
    Return the number of automorphisms of ``fp`` decorated by ``data`` and
    ``classes`` that extend the partial map ``state`` (as in
    :func:`propagate`), a union of components mapped onto components.

    The least active half-edge not mapped yet is sent to every active
    half-edge not yet an image in turn, the images its component forces
    are propagated, and the complete bijections are counted; each is
    checked to commute with ``fp``.
    """
    g, inverse, _ = state
    active = [h for h in range(len(fp)) if fp[h] != -1]
    unmapped = [h for h in active if h not in g]
    if not unmapped:
        assert sorted(g.values()) == active
        return int(all(g[fp[h]] == fp[g[h]] for h in active))
    h0 = unmapped[0]
    total = 0
    for x in active:
        if x not in inverse:
            extended = propagate(fp, data, classes, state, h0, x)
            if extended is not None:
                total += count_extensions(fp, data, classes, extended)
    return total


def empty_state(classes):
    r"""
    Return the state of :func:`propagate` of the empty partial map.
    """
    return {}, {}, [({}, {}) for _ in classes]


def brute_force_automorphism_count(fp, data, classes=()):
    r"""
    Return the number of automorphisms of the face permutation ``fp``,
    connected or not, decorated by the tuple of half-edge arrays ``data``
    (kept exactly) and ``classes`` (kept up to relabelling of their values,
    by one bijection for all components).

    The automorphisms are enumerated by a backtracking over the components
    (see :func:`count_extensions`): on a connected map, every active
    half-edge is tried as the image of the least one. This uses neither
    the canonical key of :mod:`combisurf.isomorphism_naive` nor
    :class:`OrientedMap`.
    """
    return count_extensions(fp, data, classes, empty_state(classes))


def check_labelling_gives_form(fp, data, classes, r, form):
    r"""
    Check that the canonical labelling ``r`` relabels ``fp``, ``data`` and
    ``classes`` into the canonical form ``form``: the face permutation on
    the labels of the active edges, then the arrays read in the order of
    the labels of the active half-edges, the classes renumbered by first
    appearance.
    """
    assert sorted(r) == list(range(len(fp)))
    assert all(r[h ^ 1] == r[h] ^ 1 for h in range(len(fp)))
    k = len(form[0])
    fp_c, arrays_c = relabel_arrays(fp, data + classes, r)
    assert tuple(fp_c[:k]) == form[0]
    read = [tuple(a[x] for x in range(k) if fp_c[x] != -1) for a in arrays_c]
    renumbered = []
    for values in read[len(data):]:
        names = {}
        renumbered.append(tuple(names.setdefault(v, len(names)) for v in values))
    assert tuple(read[:len(data)]) + tuple(renumbered) == form[1:]


def test_random_connected_fp():
    rng = random.Random(0)
    folded = gaps = symmetric = 0
    for _ in range(100):
        fp = random_connected_fp(rng)
        assert is_connected_fp(fp)
        OrientedMap(fp=fp)
        folded += any(fp[2 * e] != -1 and fp[2 * e + 1] == -1 for e in range(len(fp) // 2))
        gaps += any(fp[2 * e] == -1 for e in range(len(fp) // 2))
        symmetric += len(isomorphism_naive.automorphisms(fp)) > 1
    assert folded and gaps and symmetric


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_naive_canonical_form_invariance(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_connected_fp(rng)
        data = random_data(rng, fp)
        c = isomorphism_naive.canonical_form(fp, data, root_rule=root_rule)
        r = isomorphism_naive.canonical_labelling(fp, data, root_rule=root_rule)
        check_labelling_gives_form(fp, data, (), r, c)
        fp_c, data_c = relabel_arrays(fp, data, r)
        for _ in range(3):
            s = random_edge_relabelling(rng, fp)
            fp2, data2 = relabel_arrays(fp, data, s)
            assert isomorphism_naive.canonical_form(fp2, data2, root_rule=root_rule) == c
            r2 = isomorphism_naive.canonical_labelling(fp2, data2, root_rule=root_rule)
            assert relabel_arrays(fp2, data2, r2)[0][:len(c[0])] == fp_c[:len(c[0])]


def check_same_group(fp, data, classes, automorphisms, root_rule):
    r"""
    Check that the list ``automorphisms`` under the rule ``root_rule`` is,
    as a set, the list under ``ROOTS_ALL``, and that the generators under
    the rule generate it.
    """
    listed = {tuple(g) for g in automorphisms}
    if root_rule != isomorphism_naive.ROOTS_ALL:
        assert listed == {tuple(g) for g in isomorphism_naive.automorphisms(fp, data, classes)}
    gens = isomorphism_naive.automorphism_generators(fp, data, classes, root_rule=root_rule)
    assert group_closure(gens, len(fp)) == listed


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_naive_automorphisms(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_connected_fp(rng)
        data = random_data(rng, fp)
        automorphisms = isomorphism_naive.automorphisms(fp, data, root_rule=root_rule)
        assert list(automorphisms[0]) == list(range(len(fp)))
        assert len({tuple(g) for g in automorphisms}) == len(automorphisms)
        for g in automorphisms:
            check_automorphism(fp, data, g)
        assert len(automorphisms) == brute_force_automorphism_count(fp, data)
        check_same_group(fp, data, (), automorphisms, root_rule)


@pytest.mark.parametrize("seed", range(5))
def test_is_isomorphic_certificate(seed):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_connected_fp(rng)
        data = [rng.randrange(3) for _ in fp]
        a = OrientedMap(fp=fp)
        s = random_edge_relabelling(rng, fp)
        fp2, (data2,) = relabel_arrays(fp, (data,), s)
        b = OrientedMap(fp=fp2)
        ok, r = a.is_isomorphic(b, certificate=True, half_edge_data=(data, data2))
        assert ok
        check_relabelling_between(a, b, r)
        assert relabel_map(a, r) == b
        assert all(data2[r[h]] == data[h] for h in a.half_edges())
        aa = a.copy(mutable=True)
        aa.relabel(r)
        assert aa == b
        data3 = data2[:]
        h = next(b.half_edges())
        data3[h] = 3
        assert not a.is_isomorphic(b, half_edge_data=(data, data3))


@pytest.mark.parametrize("seed", range(5))
def test_edge_data_as_half_edge_data(seed):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_connected_fp(rng)
        m = OrientedMap(fp=fp)
        colours = [rng.randrange(3) for _ in range(len(fp) // 2)]
        half_edge_colours = [colours[h // 2] for h in range(len(fp))]
        r = m.canonical_labelling(edge_data=colours)
        assert r == m.canonical_labelling(half_edge_data=half_edge_colours)
        assert m.canonical_form(edge_data=colours) == m.canonical_form(half_edge_data=half_edge_colours)
        automorphisms = m.automorphisms(edge_data=colours)
        assert automorphisms == m.automorphisms(half_edge_data=half_edge_colours)
        assert isinstance(r, array)


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_naive_automorphisms_classes(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_connected_fp(rng)
        data = random_data(rng, fp)
        classes = random_classes(rng, fp)
        automorphisms = isomorphism_naive.automorphisms(fp, data, classes, root_rule=root_rule)
        assert list(automorphisms[0]) == list(range(len(fp)))
        assert len({tuple(g) for g in automorphisms}) == len(automorphisms)
        for g in automorphisms:
            check_automorphism(fp, data, g)
            for c in classes:
                check_classes_carried(fp, c, c, g)
        assert len(automorphisms) == brute_force_automorphism_count(fp, data, classes)
        check_same_group(fp, data, classes, automorphisms, root_rule)


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_naive_canonical_form_classes(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_connected_fp(rng)
        data = random_data(rng, fp)
        classes = random_classes(rng, fp)
        c = isomorphism_naive.canonical_form(fp, data, classes, root_rule=root_rule)
        assert len(c) == 1 + len(data) + len(classes)
        check_labelling_gives_form(fp, data, classes,
                                   isomorphism_naive.canonical_labelling(fp, data, classes, root_rule=root_rule), c)
        for _ in range(3):
            s = random_edge_relabelling(rng, fp)
            fp2, data2 = relabel_arrays(fp, data, s)
            _, classes2 = relabel_arrays(fp, classes, s)
            classes2 = rename_values(rng, classes2)
            assert isomorphism_naive.canonical_form(fp2, data2, classes2, root_rule=root_rule) == c
            r = isomorphism_naive.isomorphism(fp, data, fp2, data2, classes, classes2, root_rule=root_rule)
            assert r is not None
            for c1, c2 in zip(classes, classes2):
                check_classes_carried(fp, c1, c2, r)
        # an array kept exactly is also kept up to relabelling
        exact = {tuple(g) for g in isomorphism_naive.automorphisms(fp, data + classes)}
        relabelled = {tuple(g) for g in isomorphism_naive.automorphisms(fp, data, classes)}
        assert exact <= relabelled


@pytest.mark.parametrize("seed", range(5))
def test_is_isomorphic_classes(seed):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_connected_fp(rng)
        a = OrientedMap(fp=fp)
        half_edge_classes = [rng.randrange(3) for _ in fp]
        face_classes = [0] * len(fp)
        for face in a.faces():
            value = rng.randrange(3)
            for h in face:
                face_classes[h] = value
        s = random_edge_relabelling(rng, fp)
        fp2, arrays2 = relabel_arrays(fp, (half_edge_classes, face_classes), s)
        half_edge_classes2, face_classes2 = rename_values(rng, arrays2)
        b = OrientedMap(fp=fp2)
        ok, r = a.is_isomorphic(b, certificate=True,
                                half_edge_classes=(half_edge_classes, half_edge_classes2),
                                face_classes=(face_classes, face_classes2))
        assert ok
        check_relabelling_between(a, b, r)
        assert relabel_map(a, r) == b
        check_classes_carried(fp, half_edge_classes, half_edge_classes2, r)
        check_classes_carried(fp, face_classes, face_classes2, r)
        # a new value on one half-edge splits its class, unless the class
        # had that half-edge only
        h = next(b.half_edges())
        size = sum(1 for x in b.half_edges() if half_edge_classes2[x] == half_edge_classes2[h])
        half_edge_classes3 = half_edge_classes2[:]
        half_edge_classes3[h] = 3
        verdict = a.is_isomorphic(b, half_edge_classes=(half_edge_classes, half_edge_classes3),
                                  face_classes=(face_classes, face_classes2))
        assert verdict == (size == 1)


def test_random_disconnected_fp():
    rng = random.Random(0)
    disconnected = symmetric = 0
    for _ in range(50):
        fp = random_disconnected_fp(rng)
        assert len(fp) <= 48
        OrientedMap(fp=fp)
        disconnected += not is_connected_fp(fp)
        symmetric += brute_force_automorphism_count(fp, ()) > 1
    assert disconnected and symmetric


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_naive_automorphisms_disconnected(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_disconnected_fp(rng)
        data = random_data(rng, fp)
        classes = random_classes(rng, fp)
        automorphisms = isomorphism_naive.automorphisms(fp, data, classes, root_rule=root_rule)
        assert list(automorphisms[0]) == list(range(len(fp)))
        assert len({tuple(g) for g in automorphisms}) == len(automorphisms)
        for g in automorphisms:
            check_automorphism(fp, data, g)
            for c in classes:
                check_classes_carried(fp, c, c, g)
        assert len(automorphisms) == brute_force_automorphism_count(fp, data, classes)
        if root_rule != isomorphism_naive.ROOTS_ALL:
            check_same_group(fp, data, classes, automorphisms, root_rule)


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_naive_canonical_form_disconnected(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_disconnected_fp(rng)
        data = random_data(rng, fp)
        classes = random_classes(rng, fp)
        c = isomorphism_naive.canonical_form(fp, data, classes, root_rule=root_rule)
        assert len(c) == 1 + len(data) + len(classes)
        check_labelling_gives_form(fp, data, classes,
                                   isomorphism_naive.canonical_labelling(fp, data, classes, root_rule=root_rule), c)
        for _ in range(3):
            s = random_edge_relabelling(rng, fp)
            fp2, data2 = relabel_arrays(fp, data, s)
            _, classes2 = relabel_arrays(fp, classes, s)
            classes2 = rename_values(rng, classes2)
            assert isomorphism_naive.canonical_form(fp2, data2, classes2, root_rule=root_rule) == c


@pytest.mark.parametrize("seed", range(5))
def test_is_isomorphic_disconnected(seed):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_disconnected_fp(rng)
        a = OrientedMap(fp=fp)
        data = [rng.randrange(2) for _ in fp]
        classes = [rng.randrange(3) for _ in fp]
        s = random_edge_relabelling(rng, fp)
        fp2, (data2, classes2) = relabel_arrays(fp, (data, classes), s)
        (classes2,) = rename_values(rng, (classes2,))
        b = OrientedMap(fp=fp2)
        ok, r = a.is_isomorphic(b, certificate=True, half_edge_data=(data, data2),
                                half_edge_classes=(classes, classes2))
        assert ok
        check_relabelling_between(a, b, r)
        assert relabel_map(a, r) == b
        assert all(data2[r[h]] == data[h] for h in a.half_edges())
        check_classes_carried(fp, classes, classes2, r)


# The fast implementation against the reference one

def random_large_connected_fp(rng):
    r"""
    Return a random connected face permutation of 200 to 400 half-edge
    labels, as :func:`random_connected_fp`.
    """
    for _ in range(MAX_ATTEMPTS):
        fp = random_connected_fp(rng, max_half_edges=400)
        if len(fp) >= 200:
            return fp
    raise RuntimeError("no large map found in {} attempts".format(MAX_ATTEMPTS))


def flatten(form):
    r"""
    Return the canonical form ``form`` flattened as the bytes of C ``int``,
    the layout of :func:`combisurf.isomorphism.canonical_key`.
    """
    return array('i', [x for t in form for x in t]).tobytes()


def check_generators(gens, expected, automorphisms, fp, classes):
    r"""
    Check the generators ``gens`` of the map ``fp``: equal to the generators
    ``expected`` of the reference on a connected map; on a disconnected
    map, each one in the list ``automorphisms`` and generating its group.
    """
    if is_connected_fp(fp):
        assert gens == expected
        return
    listed = {tuple(g) for g in automorphisms}
    assert all(tuple(g) in listed for g in gens)
    assert group_closure(gens, len(fp)) == listed


def check_fast_against_naive(rng, fp, data, classes, root_rule=isomorphism_naive.ROOTS_ALL):
    r"""
    Check that :mod:`combisurf.isomorphism` and
    :mod:`combisurf.isomorphism_naive` agree on the decorated map ``fp``,
    on a relabelled copy of it and on a map that is not isomorphic to it,
    under the rule of the roots ``root_rule``.
    The arrays are passed as lists or as ``array('i')`` at random.
    """
    def arrays(fp, data):
        if rng.random() < 0.5:
            return array('i', fp), tuple(array('i', d) for d in data)
        return fp, data

    fast, naive = isomorphism_fast, isomorphism_naive
    rule = {'root_rule': root_rule}
    fp1, data1 = arrays(fp, data)
    r = naive.canonical_labelling(fp, data, classes, **rule)
    assert fast.canonical_labelling(fp1, data1, classes, **rule) == r
    form = naive.canonical_form(fp, data, classes, **rule)
    assert fast.canonical_form(fp1, data1, classes, **rule) == form
    automorphisms = naive.automorphisms(fp, data, classes, **rule)
    assert fast.automorphisms(fp1, data1, classes, **rule) == automorphisms
    check_generators(fast.automorphism_generators(fp1, data1, classes, **rule),
                     naive.automorphism_generators(fp, data, classes, **rule), automorphisms, fp, classes)
    key, r1, count = fast.canonical_key(fp1, data1, classes, relabelling=True, automorphism_count=True, **rule)
    assert key == flatten(form)
    assert r1 == r
    assert count == len(automorphisms)
    assert fast.canonical_key(fp1, data1, classes, **rule) == key
    assert fast.canonical_key(fp1, data1, classes, relabelling=True, **rule) == (key, r)
    assert fast.canonical_key(fp1, data1, classes, automorphism_count=True, **rule) == (key, count)

    s = random_edge_relabelling(rng, fp)
    fp2, data2 = relabel_arrays(fp, data, s)
    _, classes2 = relabel_arrays(fp, classes, s)
    classes2 = rename_values(rng, classes2)
    fp2a, data2a = arrays(fp2, data2)
    assert fast.canonical_key(fp2a, data2a, classes2, **rule) == key
    iso = fast.isomorphism(fp1, data1, fp2a, data2a, classes, classes2, **rule)
    assert iso is not None
    assert iso == naive.isomorphism(fp, data, fp2, data2, classes, classes2, **rule)

    # one more array, marking one half-edge in the copy only
    marked = [0] * len(fp2)
    marked[next(h for h in range(len(fp2)) if fp2[h] != -1)] = 1
    zeros = [0] * len(fp)
    assert fast.isomorphism(fp1, data1 + (zeros,), fp2a, data2a + (marked,), classes, classes2, **rule) is None
    assert naive.isomorphism(fp, data + (zeros,), fp2, data2 + (marked,), classes, classes2, **rule) is None


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_fast_against_naive(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(20):
        fp = random_connected_fp(rng)
        check_fast_against_naive(rng, fp, random_data(rng, fp), random_classes(rng, fp), root_rule)


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_fast_against_naive_large(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(2):
        fp = random_large_connected_fp(rng)
        check_fast_against_naive(rng, fp, random_data(rng, fp), random_classes(rng, fp), root_rule)


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_fast_against_naive_disconnected(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(10):
        fp = random_disconnected_fp(rng)
        check_fast_against_naive(rng, fp, random_data(rng, fp), random_classes(rng, fp), root_rule)


def check_same_verdict(fp1, data1, fp2, data2, classes1=(), classes2=(), expected=None):
    r"""
    Check that :func:`combisurf.isomorphism.isomorphism` returns the result
    of :func:`combisurf.isomorphism_naive.isomorphism` on the two decorated
    maps, in both orders; with ``expected="none"`` that both results are
    ``None``.
    """
    for a, b in [((fp1, data1, classes1), (fp2, data2, classes2)),
                 ((fp2, data2, classes2), (fp1, data1, classes1))]:
        naive = isomorphism_naive.isomorphism(a[0], a[1], b[0], b[1], a[2], b[2])
        assert isomorphism_fast.isomorphism(a[0], a[1], b[0], b[1], a[2], b[2]) == naive
        if expected == "none":
            assert naive is None


def test_isomorphism_rejects_active_counts():
    # the same number of half-edges, edges active or not: two edges against
    # one, and one folded edge against one edge with two monogons; then the
    # same counts of active edges and half-edges, placed differently
    check_same_verdict([1, 0, 2, 3], (), [1, 0, -1, -1], (), expected="none")
    check_same_verdict([0, -1], (), [0, 1], (), expected="none")
    check_same_verdict([0, -1, 2, 3], (), [0, 1, 2, -1], ())
    assert isomorphism_fast.isomorphism([0, -1, 2, 3], (), [0, 1, 2, -1], ()) is not None


def test_isomorphism_rejects_degrees():
    # three edges and two triangles each: the sphere (three vertices of
    # degree 2) against the torus (one vertex of degree 6), then their duals
    # (the same vertex degrees, faces of degrees 2, 2, 2 against 6)
    sphere = OrientedMap(fp=[2, 5, 4, 1, 0, 3])
    torus = OrientedMap(fp=[2, 3, 4, 5, 0, 1])
    check_same_verdict(list(sphere._fp), (), list(torus._fp), (), expected="none")
    check_same_verdict(list(sphere.dual()._fp), (), list(torus.dual()._fp), (), expected="none")


@pytest.mark.parametrize("seed", range(5))
def test_isomorphism_rejects_data_multiset(seed):
    # a relabelled copy whose data differs by one value on one half-edge, at
    # the extremes of a C int as well, against the copy whose values are
    # moved between half-edges (the same multiset)
    rng = random.Random(seed)
    for _ in range(10):
        fp = random_connected_fp(rng)
        values = [0, 1, -2**31, 2**31 - 1]
        data = ([rng.choice(values) for _ in fp], [rng.randrange(2) for _ in fp])
        fp2, data2 = relabel_arrays(fp, data, random_edge_relabelling(rng, fp))
        check_same_verdict(fp, data, fp2, data2)
        assert isomorphism_fast.isomorphism(fp, data, fp2, data2) is not None
        active = [x for x in range(len(fp2)) if fp2[x] != -1]
        h = rng.choice(active)
        for j in range(2):
            changed = list(data2)
            changed[j] = changed[j][:]
            changed[j][h] = values[(values.index(changed[j][h]) + 1) % 4]
            check_same_verdict(fp, data, fp2, tuple(changed), expected="none")
        if len(active) < 2:
            continue
        a, b = rng.sample(active, 2)
        moved = data2[0][:]
        moved[a], moved[b] = moved[b], moved[a]
        check_same_verdict(fp, data, fp2, (moved, data2[1]))


@pytest.mark.parametrize("seed", range(5))
def test_isomorphism_rejects_class_sizes(seed):
    # a relabelled copy whose classes are renamed, against one where a
    # half-edge leaves its class (other sizes) and one where two half-edges
    # of different classes swap them (the same sizes)
    rng = random.Random(seed)
    for _ in range(10):
        fp = random_connected_fp(rng)
        classes = ([rng.randrange(3) for _ in fp], [0] * len(fp))
        fp2, classes2 = relabel_arrays(fp, classes, random_edge_relabelling(rng, fp))
        classes2 = rename_values(rng, classes2)
        check_same_verdict(fp, (), fp2, (), classes, classes2)
        assert isomorphism_fast.isomorphism(fp, (), fp2, (), classes, classes2) is not None
        active = [x for x in range(len(fp2)) if fp2[x] != -1]
        if len(active) < 2:
            continue
        h = rng.choice(active)
        split = classes2[1][:]
        split[h] = max(split) + 1
        check_same_verdict(fp, (), fp2, (), classes, (classes2[0], split), expected="none")
        a, b = rng.sample(active, 2)
        swapped = classes2[0][:]
        swapped[a], swapped[b] = swapped[b], swapped[a]
        check_same_verdict(fp, (), fp2, (), classes, (swapped, classes2[1]))


def random_large_disconnected_fp(rng):
    r"""
    Return the face permutation of the disjoint union of two or three
    random connected maps of 60 to 150 half-edge labels each, a component
    after the first being a relabelled copy of an earlier one with
    probability one half, with the labels of the union shuffled by edges.
    """
    components = []
    for _ in range(rng.choice([2, 3])):
        if components and rng.random() < 0.5:
            fp = rng.choice(components)
            fp = relabel_arrays(fp, (), random_edge_relabelling(rng, fp))[0]
        else:
            for _ in range(MAX_ATTEMPTS):
                fp = random_connected_fp(rng, max_half_edges=150)
                if len(fp) >= 60:
                    break
            else:
                raise RuntimeError("no component found in {} attempts".format(MAX_ATTEMPTS))
        components.append(fp)
    union = []
    for fp in components:
        shift = len(union)
        union.extend(-1 if x == -1 else x + shift for x in fp)
    return relabel_arrays(union, (), random_edge_relabelling(rng, union))[0]


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_fast_against_naive_large_disconnected(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(3):
        fp = random_large_disconnected_fp(rng)
        check_fast_against_naive(rng, fp, random_data(rng, fp), random_classes(rng, fp), root_rule)


def test_default_root_rule():
    # omitting root_rule is the rule half-edge-fewest+data in each layer
    # (ROOTS_HALF_EDGE_FEWEST | ROOTS_DATA, the integer 10); the decoration
    # of the OrientedMap part makes the flag matter
    rng = random.Random(0)
    fewest = isomorphism_naive.ROOTS_HALF_EDGE_FEWEST | isomorphism_naive.ROOTS_DATA
    assert fewest == 10
    for fp in [random_connected_fp(rng) for _ in range(10)] + [random_disconnected_fp(rng) for _ in range(10)]:
        data, classes = random_data(rng, fp), random_classes(rng, fp)
        for module in (isomorphism_naive, isomorphism_fast):
            assert module.canonical_labelling(fp, data, classes) == module.canonical_labelling(fp, data, classes, root_rule=fewest)
            assert module.canonical_form(fp, data, classes) == module.canonical_form(fp, data, classes, root_rule=fewest)
            assert module.automorphisms(fp, data, classes) == module.automorphisms(fp, data, classes, root_rule=fewest)
            assert module.isomorphism(fp, data, fp, data, classes, classes) == \
                module.isomorphism(fp, data, fp, data, classes, classes, root_rule=fewest)
        assert isomorphism_fast.canonical_key(fp, data, classes) == \
            isomorphism_fast.canonical_key(fp, data, classes, root_rule=fewest)
        m = OrientedMap(fp=fp)
        rule = "half-edge-fewest+data"
        options = {"half_edge_data": data[0]} if data else {}
        assert m.canonical_labelling(**options) == m.canonical_labelling(root_rule=rule, **options)
        assert m.canonical_form(**options) == m.canonical_form(root_rule=rule, **options)
        assert m.automorphisms(**options) == m.automorphisms(root_rule=rule, **options)
        assert m.automorphism_generators(**options) == m.automorphism_generators(root_rule=rule, **options)
        pair = {"half_edge_data": (data[0], data[0])} if data else {}
        assert m.is_isomorphic(m, certificate=True, **pair) == m.is_isomorphic(m, certificate=True, root_rule=rule, **pair)


def random_origami(rng, n):
    r"""
    Return a random connected origami ``(r, u)`` of ``n`` squares.
    """
    for _ in range(MAX_ATTEMPTS):
        r = list(range(n))
        u = list(range(n))
        rng.shuffle(r)
        rng.shuffle(u)
        if is_connected_fp(origami_face_permutation(r, u)):
            return r, u
    raise RuntimeError("no connected origami found in {} attempts".format(MAX_ATTEMPTS))


@pytest.mark.parametrize("seed", range(3))
def test_fast_against_naive_origamis(seed):
    # the decoration h % 4 makes many roots tie on the face permutation, so
    # that the comparison of roots is decided by the data or not at all
    rng = random.Random(seed)
    origamis = [([1, 2, 3, 4, 5, 0], [0, 1, 2, 3, 4, 5]),
                ([1, 2, 0, 4, 5, 3], [3, 4, 5, 1, 2, 0]),
                ([(i + 1) % 12 for i in range(12)], [(i + 5) % 12 for i in range(12)])]
    origamis.extend(random_origami(rng, rng.randrange(2, 20)) for _ in range(5))
    for r, u in origamis:
        fp = origami_face_permutation(r, u)
        d = origami_decoration(len(r))
        check_fast_against_naive(rng, fp, (d,), ())
        check_fast_against_naive(rng, fp, (), (d,))
        check_fast_against_naive(rng, fp, (), ())


@pytest.mark.parametrize("build", [regular_map_S5, regular_map_PSL27], ids=["O22", "O23"])
def test_fast_against_naive_regular_map(build):
    # every root ties: no root is rejected early
    rng = random.Random(0)
    fp = list(build().face_permutation())
    check_fast_against_naive(rng, fp, (), ())


def test_fast_four_edges():
    fp = list(range(8))
    c = ([0, 1, 0, 1, 2, 3, 2, 3],)
    assert len(isomorphism_fast.automorphisms(fp)) == 384
    assert len(isomorphism_fast.automorphisms(fp, classes=c)) == 32
    assert isomorphism_fast.canonical_key(fp, automorphism_count=True)[1] == 384
    assert isomorphism_fast.canonical_key(fp, classes=c, automorphism_count=True)[1] == 32
    gens = isomorphism_fast.automorphism_generators(fp)
    assert len(gens) == 3
    check_generators(gens, isomorphism_naive.automorphism_generators(fp),
                     isomorphism_naive.automorphisms(fp), fp, ())
    gens = isomorphism_fast.automorphism_generators(fp, classes=c)
    check_generators(gens, isomorphism_naive.automorphism_generators(fp, classes=c),
                     isomorphism_naive.automorphisms(fp, classes=c), fp, c)


def wreath_fp(k, component):
    r"""
    Return the face permutation of ``k`` disjoint copies of ``component``.
    """
    return [-1 if x == -1 else x + len(component) * i for i in range(k) for x in component]


def is_automorphism(g, fp, data=()):
    r"""
    Whether ``g`` commutes with ``fp`` on the active half-edges, preserves
    the pairs ``{2e, 2e+1}``, the arrays ``data`` and is the identity on the
    inactive edges.
    """
    for h in range(len(fp)):
        if g[h ^ 1] != g[h] ^ 1:
            return False
        if fp[h] == -1:
            if fp[h ^ 1] == -1 and g[h] != h:
                return False
            continue
        if g[fp[h]] != fp[g[h]]:
            return False
        if any(d[g[h]] != d[h] for d in data):
            return False
    return sorted(g) == list(range(len(fp)))


def test_fast_wreath_generators():
    fast = isomorphism_fast
    # 20 isolated edges: the group is too large to be listed
    fp = list(range(40))
    gens = fast.automorphism_generators(fp)
    assert len(gens) == 3
    assert all(is_automorphism(g, fp) for g in gens)
    assert fast.canonical_key(fp, automorphism_count=True)[1] == 2**20 * math.factorial(20)

    triangle = [2, 5, 4, 1, 0, 3]
    fp = wreath_fp(3, triangle)
    gens = fast.automorphism_generators(fp)
    assert len(gens) == 4
    assert all(is_automorphism(g, fp) for g in gens)
    assert len(group_closure(gens, len(fp))) == 1296 == fast.canonical_key(fp, automorphism_count=True)[1]

    # data separating two of three isomorphic components
    data = ([0] * 6 + [1] * 6 + [0] * 6,)
    automorphisms = isomorphism_naive.automorphisms(fp, data)
    gens = fast.automorphism_generators(fp, data)
    assert all(is_automorphism(g, fp, data) for g in gens)
    assert group_closure(gens, len(fp)) == {tuple(g) for g in automorphisms}
    assert fast.canonical_key(fp, data, automorphism_count=True)[1] == len(automorphisms)

    # a folded edge and an inactive edge in the components, components interleaved
    fp = [2, -1, 0, 3, -1, -1, 8, -1, 6, 9, 10, 11, -1, -1]
    automorphisms = isomorphism_naive.automorphisms(fp)
    gens = fast.automorphism_generators(fp)
    assert all(is_automorphism(g, fp) for g in gens)
    assert group_closure(gens, len(fp)) == {tuple(g) for g in automorphisms}
    assert fast.canonical_key(fp, automorphism_count=True)[1] == len(automorphisms)


@WITH_ROOT_RULES
def test_fast_wreath_without_group_closure(monkeypatch, root_rule):
    def unusable(*args, **kwds):
        raise AssertionError("the group closure was called")

    monkeypatch.setattr(isomorphism_naive, "_group_closure", unusable)
    monkeypatch.setattr(isomorphism_naive, "automorphisms", unusable)
    monkeypatch.setattr(isomorphism_naive, "automorphism_generators", unusable)
    fast = isomorphism_fast
    triangle = [2, 5, 4, 1, 0, 3]
    for fp, data in [(list(range(14)), ()), (wreath_fp(3, triangle), ()),
                     (wreath_fp(3, triangle), ([0] * 6 + [1] * 6 + [0] * 6,)),
                     ([0, 1, 2, 3], ())]:
        assert fast.automorphism_generators(fp, data, root_rule=root_rule)
        assert fast.canonical_key(fp, data, automorphism_count=True, root_rule=root_rule)[1] > 1


@WITH_ROOT_RULES
def test_fast_disconnected_without_naive(monkeypatch, root_rule):
    rng = random.Random(0)
    cases = [([0, 1, 2, 3], (), ()),
             ([0, 1, 2, 3], (), ([0, 1, 1, 2],)),
             (list(range(8)), (), ([0, 1, 0, 1, 2, 3, 2, 3],)),
             ([2, 5, 4, 1, 0, 3, 6, -1], ([0, 1, 0, 0, 1, 1, 0, 0],), ())]
    for _ in range(3):
        fp = random_disconnected_fp(rng)
        while is_connected_fp(fp):
            fp = random_disconnected_fp(rng)
        cases.append((fp, random_data(rng, fp), random_classes(rng, fp)))

    # the results of the reference, then the reference made unusable
    rule = {'root_rule': root_rule}
    expected = []
    for fp, data, classes in cases:
        s = random_edge_relabelling(rng, fp)
        fp2, data2 = relabel_arrays(fp, data, s)
        _, classes2 = relabel_arrays(fp, classes, s)
        naive = isomorphism_naive
        expected.append((naive.canonical_labelling(fp, data, classes, **rule),
                         naive.canonical_form(fp, data, classes, **rule),
                         naive.automorphisms(fp, data, classes, **rule),
                         naive.automorphism_generators(fp, data, classes, **rule),
                         fp2, data2, classes2,
                         naive.isomorphism(fp, data, fp2, data2, classes, classes2, **rule)))

    def unusable(*args, **kwds):
        raise AssertionError("the reference implementation was called")

    for name in ("_search", "_roots", "_root_invariants", "_group_closure", "canonical_labelling",
                 "canonical_form", "automorphisms", "automorphism_generators", "isomorphism"):
        monkeypatch.setattr(isomorphism_naive, name, unusable)

    fast = isomorphism_fast
    for (fp, data, classes), (r, form, automorphisms, gens, fp2, data2, classes2, iso) in zip(cases, expected):
        assert fast.canonical_labelling(fp, data, classes, **rule) == r
        assert fast.canonical_form(fp, data, classes, **rule) == form
        assert fast.automorphisms(fp, data, classes, **rule) == automorphisms
        check_generators(fast.automorphism_generators(fp, data, classes, **rule), gens, automorphisms, fp, classes)
        assert fast.isomorphism(fp, data, fp2, data2, classes, classes2, **rule) == iso
        key = flatten(form)
        n = len(automorphisms)
        assert fast.canonical_key(fp, data, classes, **rule) == key
        assert fast.canonical_key(fp, data, classes, relabelling=True, **rule) == (key, r)
        assert fast.canonical_key(fp, data, classes, automorphism_count=True, **rule) == (key, n)
        assert fast.canonical_key(fp, data, classes, relabelling=True, automorphism_count=True, **rule) == (key, r, n)


def check_generated_group(fp, data, classes):
    r"""
    Check that the generators of :func:`combisurf.isomorphism.automorphism_generators`
    are automorphisms of the reference list and generate it, and that the
    count of :func:`combisurf.isomorphism.canonical_key` is its length.
    """
    automorphisms = isomorphism_naive.automorphisms(fp, data, classes)
    listed = {tuple(g) for g in automorphisms}
    gens = isomorphism_fast.automorphism_generators(fp, data, classes)
    assert all(tuple(g) in listed for g in gens)
    assert group_closure(gens, len(fp)) == listed
    assert isomorphism_fast.canonical_key(fp, data, classes, automorphism_count=True)[1] == len(automorphisms)


def cylinder_diagram_case(s):
    r"""
    Return the face permutation, the data (the colours) and the classes
    (the cylinders) of the cylinder diagram ``s``.
    """
    m, colour, cylinder = cylinder_diagram(s)
    return list(m._fp), (colour,), (cylinder,)


CLASSES_CASES = [
    (list(range(8)), (), ([0, 1, 0, 1, 2, 3, 2, 3],)),
    (list(range(14)), (), ([(h // 4) * 2 + h % 2 for h in range(14)],)),
    (list(range(12)), (), ([0] * 12,)),
    (wreath_fp(5, [2, 5, 4, 1, 0, 3]), (), ([h % 3 for h in range(30)],)),
] + [cylinder_diagram_case(s) for s in (CD_I3_A, CD_I3_B, CD_I6_A, CD_I6_B)]


@pytest.mark.parametrize("fp, data, classes", CLASSES_CASES,
                         ids=["4-edges", "7-edges-pairs", "6-edges-one-value",
                              "5-triangles", "I3a", "I3b", "I6a", "I6b"])
def test_fast_generators_classes(fp, data, classes):
    check_generated_group(fp, data, classes)


def test_fast_generators_classes_twenty_edges():
    from sage.all import PermutationGroup

    # 2^20 * 20! automorphisms: too many for a closure; the order of the
    # group generated is computed by Sage
    fp = list(range(40))
    c = [0] * 40
    gens = isomorphism_fast.automorphism_generators(fp, classes=(c,))
    for g in gens:
        assert is_automorphism(g, fp)
        check_classes_carried(fp, c, c, g)
    G = PermutationGroup([[x + 1 for x in g] for g in gens], domain=list(range(1, 41)))
    assert G.order() == 2**20 * math.factorial(20)


def test_fast_canonical_key_classes_twenty_edges():
    # the search of the key follows every tied root, and returns from each
    # leaf equal to the best one to the first level where their roots
    # differ, instead of visiting the 2^20 * 20! leaves of the minimum
    rng = random.Random(0)
    fp = list(range(40))
    c = [0] * 40
    key, count = isomorphism_fast.canonical_key(fp, classes=(c,), automorphism_count=True)
    assert count == 2**20 * math.factorial(20)
    assert key == flatten((tuple(range(40)), tuple([0] * 40)))
    for _ in range(5):
        fp2, classes2 = relabel_arrays(fp, (c,), random_edge_relabelling(rng, fp))
        assert isomorphism_fast.canonical_key(fp2, (), rename_values(rng, classes2)) == key


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(5))
def test_fast_against_naive_classes_copies(seed, root_rule):
    # two or three copies of a component, its classes repeated on each copy
    # so that the copies tie: the search of the key backjumps from the leaves
    # equal to the best one, the reference visits them all
    rng = random.Random(seed)
    fast, naive = isomorphism_fast, isomorphism_naive
    rule = {'root_rule': root_rule}
    for _ in range(4):
        component = random_connected_fp(rng, max_half_edges=12)
        k = rng.randint(2, 3)
        values = rng.randint(1, 2)
        classes = tuple([rng.randrange(values) for _ in component] * k for _ in range(rng.randint(1, 2)))
        union = wreath_fp(k, component)
        fp, classes = relabel_arrays(union, classes, random_edge_relabelling(rng, union))
        form = naive.canonical_form(fp, (), classes, **rule)
        r = naive.canonical_labelling(fp, (), classes, **rule)
        assert fast.canonical_form(fp, (), classes, **rule) == form
        assert fast.canonical_labelling(fp, (), classes, **rule) == r
        assert fast.canonical_key(fp, (), classes, relabelling=True, **rule) == (flatten(form), r)
        s = random_edge_relabelling(rng, fp)
        fp2, classes2 = relabel_arrays(fp, classes, s)
        assert fast.canonical_key(fp2, (), rename_values(rng, classes2), **rule) == flatten(form)


@pytest.mark.parametrize("seed", range(5))
def test_fast_generators_classes_copies(seed):
    rng = random.Random(seed)
    for _ in range(3):
        fp = wreath_fp(rng.randint(2, 4), random_connected_fp(rng, max_half_edges=16))
        for values in (2, 3):
            classes = tuple([rng.randrange(values) for _ in fp] for _ in range(rng.randint(1, 2)))
            fp2, classes2 = relabel_arrays(fp, classes, random_edge_relabelling(rng, fp))
            check_generated_group(fp2, (), classes2)


@pytest.mark.parametrize("seed", range(5))
def test_fast_generators_classes_mixed_copies(seed):
    # copies of two components of different sizes, the data repeated on the
    # copies so that they stay isomorphic, the classes of each copy drawn
    # from two patterns so that they cut the wreath product
    rng = random.Random(seed)
    for _ in range(3):
        union, data, classes = [], [], []
        for _ in range(2):
            component = random_connected_fp(rng, max_half_edges=12)
            d = [rng.randrange(2) for _ in component]
            patterns = [[rng.randrange(2) for _ in component] for _ in range(2)]
            for _ in range(rng.randint(1, 3)):
                shift = len(union)
                union.extend(-1 if x == -1 else x + shift for x in component)
                data.extend(d)
                classes.extend(rng.choice(patterns))
        r = random_edge_relabelling(rng, union)
        fp2, data2 = relabel_arrays(union, (data,), r)
        _, classes2 = relabel_arrays(union, (classes,), r)
        check_generated_group(fp2, data2, classes2)


@WITH_ROOT_RULES
def test_fast_generators_classes_without_naive(monkeypatch, root_rule):
    assert not hasattr(isomorphism_fast, "isomorphism_naive")
    rng = random.Random(1)
    cases = list(CLASSES_CASES)
    for _ in range(3):
        fp = wreath_fp(3, random_connected_fp(rng, max_half_edges=16))
        cases.append((fp, (), ([rng.randrange(2) for _ in fp],)))
    expected = [{tuple(g) for g in isomorphism_naive.automorphisms(fp, data, classes)}
                for fp, data, classes in cases]

    def unusable(*args, **kwds):
        raise AssertionError("the reference implementation was called")

    for name in ("_search", "_roots", "_root_invariants", "_group_closure", "canonical_labelling",
                 "canonical_form", "automorphisms", "automorphism_generators", "isomorphism"):
        monkeypatch.setattr(isomorphism_naive, name, unusable)
    for (fp, data, classes), listed in zip(cases, expected):
        gens = isomorphism_fast.automorphism_generators(fp, data, classes, root_rule=root_rule)
        assert group_closure(gens, len(fp)) == listed
        count = isomorphism_fast.canonical_key(fp, data, classes, automorphism_count=True, root_rule=root_rule)[1]
        assert count == len(listed)


@WITH_ROOT_RULES
@pytest.mark.parametrize("seed", range(3))
def test_fast_generators_rules_disconnected(seed, root_rule):
    # the generators and the count on disconnected maps do not depend on the
    # rule: they generate the group listed by the reference under the rule,
    # and the count is its length
    rng = random.Random(seed)
    triangle = [2, 5, 4, 1, 0, 3]
    cases = [(wreath_fp(3, triangle) + [18, -1], (), ()),
             (wreath_fp(3, triangle) + [18, -1], (), ([h % 2 for h in range(20)],)),
             ([2, 5, 4, 1, 0, 3, 6, -1], ([0, 1, 0, 0, 1, 1, 0, 0],), ())]
    while len(cases) < 9:
        fp = random_disconnected_fp(rng)
        if is_connected_fp(fp):
            continue
        classes = random_classes(rng, fp) if len(cases) % 2 else ()
        cases.append((fp, random_data(rng, fp), classes))
    for fp, data, classes in cases:
        automorphisms = isomorphism_naive.automorphisms(fp, data, classes, root_rule=root_rule)
        listed = {tuple(g) for g in automorphisms}
        assert len(listed) == len(automorphisms)
        gens = isomorphism_fast.automorphism_generators(fp, data, classes, root_rule=root_rule)
        assert all(tuple(g) in listed for g in gens)
        assert group_closure(gens, len(fp)) == listed
        count = isomorphism_fast.canonical_key(fp, data, classes, automorphism_count=True, root_rule=root_rule)[1]
        assert count == len(automorphisms)


def test_fast_sizes():
    isomorphism_fast._check_sizes(2**20, 2, 2)
    for n, nd, nc in [(2**31, 0, 0), (2**28, 0, 0), (2**26, 5, 5)]:
        with pytest.raises(OverflowError, match="map too large"):
            isomorphism_fast._check_sizes(n, nd, nc)


def test_fast_malformed_input():
    for fp, data, classes in [([0, -1, 2], (), ()), ([0, -1], [[0, 0]], ()),
                              ([0, -1], ([0],), ()), ([0, -1], (), [[0, 0]]),
                              ([0, -1], (), ([0],)),
                              ([1, -1], (), ()), ([1, 1], (), ()), ([0, 0], (), ()),
                              ([-2, 0], (), ()), ([0, 5], (), ()), ([2, 0], (), ())]:
        for name in ("canonical_labelling", "canonical_form", "automorphisms",
                     "automorphism_generators", "canonical_key"):
            with pytest.raises((ValueError, TypeError)) as fast_error:
                getattr(isomorphism_fast, name)(fp, data, classes)
            if name == "canonical_key":
                continue
            with pytest.raises((ValueError, TypeError)) as naive_error:
                getattr(isomorphism_naive, name)(fp, data, classes)
            assert fast_error.type is naive_error.type
            assert str(fast_error.value) == str(naive_error.value)
    with pytest.raises(ValueError, match="permutation"):
        isomorphism_fast.canonical_labelling([1, 1])
    with pytest.raises(OverflowError):
        isomorphism_fast.canonical_labelling([1, 0], ([2**40, 0],))


def test_root_rule_constants():
    for name in ("ROOTS_ALL", "ROOTS_HALF_EDGE_MAX", "ROOTS_HALF_EDGE_FEWEST",
                 "ROOTS_EDGE_MAX", "ROOTS_EDGE_FEWEST", "ROOTS_DATA"):
        assert getattr(isomorphism_fast, name) == getattr(isomorphism_naive, name)
    assert isomorphism_naive.ROOTS_DATA == 8
    assert ROOT_RULES == list(range(5)) + list(range(9, 13))
    assert [OrientedMap._root_rule(name) for name in ROOT_RULE_NAMES] == ROOT_RULES


def test_root_rule_restricts():
    # the random maps of the tests above exercise the rules: on some of them
    # a rule tries fewer roots than the active half-edges, and the rules
    # do not all try the same roots
    rng = random.Random(0)
    restricted = {rule: 0 for rule in ROOT_RULES[1:]}
    distinct = 0
    for _ in range(50):
        fp = random_connected_fp(rng)
        unlabelled = [-1] * len(fp)
        active = sum(1 for x in fp if x != -1)
        roots = {}
        for rule in ROOT_RULES:
            invariants = isomorphism_naive._root_invariants(fp, rule)
            roots[rule] = isomorphism_naive._roots(fp, unlabelled, invariants, rule)
            assert 0 < len(roots[rule]) <= active
        assert len(roots[isomorphism_naive.ROOTS_ALL]) == active
        for rule in ROOT_RULES[1:]:
            restricted[rule] += len(roots[rule]) < active
        distinct += len({tuple(r) for r in roots.values()}) > 2
    assert all(restricted.values())
    assert distinct


@pytest.mark.parametrize("base", BASE_ROOT_RULES[1:], ids=ROOT_RULE_NAMES[1:5])
def test_root_rule_data_without_data(base):
    # without data the flag adds nothing to the invariants: the results are
    # those of the base rule, in both modules and in OrientedMap
    rng = random.Random(1)
    flagged = base | isomorphism_naive.ROOTS_DATA
    name = ROOT_RULE_NAMES[base]
    for fp in [random_connected_fp(rng) for _ in range(10)] + [random_disconnected_fp(rng) for _ in range(10)]:
        classes = random_classes(rng, fp)
        for module in (isomorphism_naive, isomorphism_fast):
            for f in ("canonical_labelling", "canonical_form", "automorphisms", "automorphism_generators"):
                assert getattr(module, f)(fp, (), classes, root_rule=flagged) == \
                    getattr(module, f)(fp, (), classes, root_rule=base)
            assert module.isomorphism(fp, (), fp, (), classes, classes, root_rule=flagged) == \
                module.isomorphism(fp, (), fp, (), classes, classes, root_rule=base)
        assert isomorphism_fast.canonical_key(fp, (), classes, relabelling=True, automorphism_count=True,
                                              root_rule=flagged) == \
            isomorphism_fast.canonical_key(fp, (), classes, relabelling=True, automorphism_count=True,
                                           root_rule=base)
        m = OrientedMap(fp=fp)
        assert m.canonical_labelling(root_rule=name + "+data") == m.canonical_labelling(root_rule=name)
        assert m.canonical_form(root_rule=name + "+data") == m.canonical_form(root_rule=name)


def test_root_rule_data_restricts():
    # on decorated maps the flag splits the classes of the base rule: under
    # *_MAX its roots are among those of the base rule, under *_FEWEST its
    # class is no larger than the smallest class of the base rule, and on
    # some maps it has fewer roots; the fast module follows the reference
    # with values negative and beyond the number of half-edges
    rng = random.Random(2)
    wide = [-2**31, 2**31 - 1, -3, 1000]
    fewer = {rule: 0 for rule in BASE_ROOT_RULES[1:]}
    for _ in range(50):
        fp = random_connected_fp(rng)
        data = tuple([rng.choice(wide) for _ in fp] for _ in range(1 + rng.randrange(2)))
        unlabelled = [-1] * len(fp)
        for base in BASE_ROOT_RULES[1:]:
            flagged = base | isomorphism_naive.ROOTS_DATA
            roots = isomorphism_naive._roots(fp, unlabelled, isomorphism_naive._root_invariants(fp, base, data), base)
            roots_data = isomorphism_naive._roots(fp, unlabelled,
                                                  isomorphism_naive._root_invariants(fp, flagged, data), flagged)
            assert 0 < len(roots_data) <= len(roots)
            if base in (isomorphism_naive.ROOTS_HALF_EDGE_MAX, isomorphism_naive.ROOTS_EDGE_MAX):
                assert set(roots_data) <= set(roots)
            fewer[base] += len(roots_data) < len(roots)
            check_fast_against_naive(rng, fp, data, (), flagged)
    assert all(fewer.values())


@pytest.mark.parametrize("root_rule", [5, -1, 8, 13, 2**70, "all", "edge-max", 1.0, True, None])
def test_root_rule_invalid_integer(root_rule):
    fp, fp2 = [2, 5, 4, 1, 0, 3], [4, 3, 0, 5, 2, 1]
    for module in (isomorphism_fast, isomorphism_naive):
        for name in ("canonical_labelling", "canonical_form", "automorphisms",
                     "automorphism_generators", "canonical_key"):
            if name == "canonical_key" and module is isomorphism_naive:
                continue
            with pytest.raises(ValueError, match="root_rule must be one of the integers"):
                getattr(module, name)(fp, (), (), root_rule=root_rule)
        with pytest.raises(ValueError, match="root_rule must be one of the integers"):
            module.isomorphism(fp, (), fp2, (), root_rule=root_rule)


@pytest.mark.parametrize("root_rule", ["edge", "EDGE-MAX", "edge_max", 0, 3, None, ["all"], "all+data", "data"])
def test_root_rule_invalid_string(root_rule):
    m = OrientedMap(fp="(0,1,2)(~0,~2,~1)")
    for check in (True, False):
        for call in (m.canonical_labelling, m.canonical_form, m.automorphisms,
                     m.automorphism_generators):
            with pytest.raises(ValueError, match="root_rule must be one of 'all'"):
                call(root_rule=root_rule, check=check)
        with pytest.raises(ValueError, match="root_rule must be one of 'all'"):
            m.is_isomorphic(m, root_rule=root_rule, check=check)


# faces (0,2,6,8)(1,7,3); the edge (4,5) is inactive and the half-edge 9 is
# the inactive twin of the folded half-edge 8
DECORATED_FP = [2, 7, 6, 1, -1, -1, 8, 3, 0, -1]
INACTIVE = (4, 5, 9)


class IndexOnly:
    r"""
    An integer that defines ``__index__`` but not ``__int__``.
    """
    def __init__(self, value):
        self.value = value

    def __index__(self):
        return self.value


def decoration_with(values, positions, base):
    r"""
    Return a copy of the list ``base`` with ``values`` written at
    ``positions``.
    """
    d = list(base)
    for p, v in zip(positions, values):
        d[p] = v
    return d


HALF_EDGE_BASE = [0, 1, 0, 1, 0, 0, 0, 1, 0, 0]
FACE_BASE = [0, 1, 0, 1, 0, 0, 0, 1, 0, 0]
VERTEX_BASE = [1, 0, 0, 0, 0, 0, 0, 1, 1, 0]
EDGE_BASE = [0, 1, 0, 3, 4]


def decorated_labelling(**decorations):
    r"""
    Return the canonical labelling of the map of ``DECORATED_FP`` with
    the given decorations of :meth:`OrientedMap.canonical_labelling`
    converted by hand to half-edge arrays.
    """
    data = []
    for keyword in ("half_edge_data", "edge_data", "vertex_data", "face_data"):
        d = decorations.get(keyword)
        if d is None:
            continue
        if keyword == "edge_data":
            d = [d[h // 2] for h in range(len(DECORATED_FP))]
        data.append(array('i', [0 if DECORATED_FP[h] == -1 else int(d[h]) for h in range(len(d))]))
    return isomorphism_fast.canonical_labelling(DECORATED_FP, tuple(data))


@pytest.mark.parametrize("check", [True, False])
@pytest.mark.parametrize("value", [None, "a", 0.5, 2**31, -2**31 - 1])
@pytest.mark.parametrize("container", [list, tuple])
def test_decoration_inactive_values_ignored(value, container, check):
    # values on inactive half-edges are neither checked nor converted
    m = OrientedMap(fp=DECORATED_FP)
    he = container(decoration_with([value] * 3, INACTIVE, HALF_EDGE_BASE))
    face = container(decoration_with([value] * 3, INACTIVE, FACE_BASE))
    vertex = container(decoration_with([value] * 3, INACTIVE, VERTEX_BASE))
    edge = container(decoration_with([value], [2], EDGE_BASE))
    expected = decorated_labelling(half_edge_data=HALF_EDGE_BASE, edge_data=EDGE_BASE,
                                   vertex_data=VERTEX_BASE, face_data=FACE_BASE)
    assert m.canonical_labelling(half_edge_data=he, edge_data=edge, vertex_data=vertex,
                                 face_data=face, check=check) == expected
    assert m.canonical_labelling(half_edge_classes=he, check=check) == \
        m.canonical_labelling(half_edge_classes=HALF_EDGE_BASE)
    assert m.is_isomorphic(m, half_edge_data=(he, HALF_EDGE_BASE), check=check)
    assert len(m.automorphisms(half_edge_data=he, check=check)) == \
        len(m.automorphisms(half_edge_data=HALF_EDGE_BASE))


@pytest.mark.parametrize("check", [True, False])
def test_decoration_containers(check):
    m = OrientedMap(fp=DECORATED_FP)
    expected = decorated_labelling(half_edge_data=HALF_EDGE_BASE)
    for d in (HALF_EDGE_BASE, tuple(HALF_EDGE_BASE), array('i', HALF_EDGE_BASE), array('l', HALF_EDGE_BASE),
              decoration_with([2**40] * 3, INACTIVE, HALF_EDGE_BASE)):
        assert m.canonical_labelling(half_edge_data=d, check=check) == expected
        assert m.canonical_labelling(half_edge_data=(d, HALF_EDGE_BASE), check=check) == expected
    expected = decorated_labelling(edge_data=EDGE_BASE)
    for d in (EDGE_BASE, tuple(EDGE_BASE), array('i', EDGE_BASE), array('q', EDGE_BASE)):
        assert m.canonical_labelling(edge_data=d, check=check) == expected


@pytest.mark.parametrize("value", [None, "a", 0.5])
@pytest.mark.parametrize("kind", ["data", "classes"])
def test_decoration_non_integer_active(value, kind):
    m = OrientedMap(fp=DECORATED_FP)
    for keyword, base, position in (("half_edge", HALF_EDGE_BASE, 3), ("edge", EDGE_BASE, 4),
                                    ("vertex", VERTEX_BASE, 0), ("face", FACE_BASE, 8)):
        name = keyword + "_" + kind
        for container in (list, tuple):
            d = container(decoration_with([value], [position], base))
            with pytest.raises(ValueError, match="^{} must contain integers, got {}$".format(
                    name, re.escape(repr(value)))):
                m.canonical_labelling(**{name: d})


def test_decoration_non_integer_active_unchecked():
    # without the check, the conversion to C int refuses a float and a
    # string as operator.index does
    m = OrientedMap(fp=DECORATED_FP)
    d = decoration_with([0.5], [0], HALF_EDGE_BASE)
    with pytest.raises(TypeError, match="^'float' object cannot be interpreted as an integer$"):
        m.canonical_labelling(half_edge_data=d, check=False)
    d = decoration_with(["a"], [0], HALF_EDGE_BASE)
    with pytest.raises(TypeError, match="^'str' object cannot be interpreted as an integer$"):
        m.canonical_labelling(half_edge_data=d, check=False)


def test_decoration_integer_types():
    from sage.rings.integer_ring import ZZ
    m = OrientedMap(fp=DECORATED_FP)
    kinds = [ZZ, bool]
    try:
        import numpy
    except ImportError:
        pass
    else:
        kinds.append(numpy.int64)
    for check in (True, False):
        for kind in kinds:
            d = [kind(x) for x in HALF_EDGE_BASE]
            e = [kind(x) for x in [0, 1, 0, 1, 1]]
            assert m.canonical_labelling(half_edge_data=d, check=check) == \
                decorated_labelling(half_edge_data=HALF_EDGE_BASE)
            assert m.canonical_labelling(edge_data=e, check=check) == \
                decorated_labelling(edge_data=[0, 1, 0, 1, 1])


@pytest.mark.parametrize("check", [True, False])
def test_decoration_index_only(check):
    # objects with __index__ are accepted by the check and by the conversion
    m = OrientedMap(fp=DECORATED_FP)
    d = decoration_with([IndexOnly(1)], [0], HALF_EDGE_BASE)
    expected = m.canonical_labelling(half_edge_data=decoration_with([1], [0], HALF_EDGE_BASE), check=check)
    assert m.canonical_labelling(half_edge_data=d, check=check) == expected
    assert m.canonical_labelling(half_edge_data=tuple(d), check=check) == expected


@pytest.mark.parametrize("check", [True, False])
@pytest.mark.parametrize("value", [2**31, -2**31 - 1])
def test_decoration_overflow(value, check):
    m = OrientedMap(fp=DECORATED_FP)
    for keyword, base, position in (("half_edge_data", HALF_EDGE_BASE, 3), ("edge_data", EDGE_BASE, 4),
                                    ("vertex_data", [value] * 10, 0), ("face_data", [value] * 10, 0)):
        d = decoration_with([value], [position], base)
        with pytest.raises(OverflowError, match="value too large to convert to int"):
            m.canonical_labelling(**{keyword: d}, check=check)
    # classes are renumbered before the conversion
    d = decoration_with([value], [3], HALF_EDGE_BASE)
    assert m.canonical_labelling(half_edge_classes=d, check=check) == \
        m.canonical_labelling(half_edge_classes=HALF_EDGE_BASE)


def test_decoration_order_of_errors():
    m = OrientedMap(fp=DECORATED_FP)
    big = decoration_with([2**31], [0], HALF_EDGE_BASE)
    with pytest.raises(ValueError, match="^face_data must have length 10$"):
        m.canonical_labelling(half_edge_data=big, face_data=[0])
    with pytest.raises(ValueError, match="^vertex_data must be constant on vertices$"):
        m.canonical_labelling(half_edge_data=big, vertex_data=HALF_EDGE_BASE)
    with pytest.raises(ValueError, match="^half_edge_data must contain integers, got 0.5$"):
        m.canonical_labelling(half_edge_data=[0.5] * 10, face_data=[0])
    with pytest.raises(ValueError, match="^vertex_data must contain integers, got 0.5$"):
        m.canonical_labelling(vertex_data=decoration_with([0.5], [0], HALF_EDGE_BASE))
    with pytest.raises(ValueError, match="^face_data must be constant on faces$"):
        m.canonical_labelling(face_data=VERTEX_BASE)
    with pytest.raises(ValueError, match="^edge_data must have length 5$"):
        m.canonical_labelling(edge_data=[0.5] * 4)
    assert m.canonical_labelling(face_data=VERTEX_BASE, check=False) == \
        decorated_labelling(face_data=VERTEX_BASE)


def test_canonical_key_containers():
    expected = isomorphism_fast.canonical_key(array('i', DECORATED_FP), (), relabelling=True)
    for fp in (DECORATED_FP, tuple(DECORATED_FP), array('i', DECORATED_FP), array('l', DECORATED_FP),
               [IndexOnly(x) for x in DECORATED_FP]):
        assert isomorphism_fast.canonical_key(fp, (), relabelling=True) == expected
    data = array('i', HALF_EDGE_BASE)
    expected = isomorphism_fast.canonical_key(DECORATED_FP, (data,), relabelling=True)
    for value in (None, "a", 0.5, 2**31, -2**31 - 1):
        d = decoration_with([value] * 3, INACTIVE, HALF_EDGE_BASE)
        for fp in (DECORATED_FP, array('i', DECORATED_FP)):
            assert isomorphism_fast.canonical_key(fp, (d,), relabelling=True) == expected
            assert isomorphism_fast.canonical_key(fp, (tuple(d),), relabelling=True) == expected
    assert isomorphism_fast.canonical_key(DECORATED_FP, ([True] + HALF_EDGE_BASE[1:],), relabelling=True) == \
        isomorphism_fast.canonical_key(DECORATED_FP, ([1] + HALF_EDGE_BASE[1:],), relabelling=True)


@pytest.mark.parametrize("value, error, message", [
    (0.5, TypeError, "'float' object cannot be interpreted as an integer"),
    ("a", TypeError, "'str' object cannot be interpreted as an integer"),
    (None, TypeError, "'NoneType' object cannot be interpreted as an integer"),
    (2**31, OverflowError, "signed integer is greater than maximum"),
    (-2**31 - 1, OverflowError, "signed integer is less than minimum"),
])
def test_canonical_key_fp_list_invalid(value, error, message):
    for position in (0, 4, 9):
        fp = decoration_with([value], [position], DECORATED_FP)
        with pytest.raises(error, match="^{}$".format(re.escape(message))):
            isomorphism_fast.canonical_key(fp, ())


@pytest.mark.parametrize("value, error, message", [
    ("a", TypeError, "'str' object cannot be interpreted as an integer"),
    (None, TypeError, "'NoneType' object cannot be interpreted as an integer"),
    (2**31, OverflowError, "value too large to convert to int"),
    (-2**31 - 1, OverflowError, "value too large to convert to int"),
])
def test_canonical_key_data_invalid_active(value, error, message):
    for container in (list, tuple):
        d = container(decoration_with([value], [8], HALF_EDGE_BASE))
        with pytest.raises(error, match="^{}$".format(re.escape(message))):
            isomorphism_fast.canonical_key(DECORATED_FP, (d,))


def test_canonical_key_data_index_active():
    d = decoration_with([IndexOnly(1)], [8], HALF_EDGE_BASE)
    expected = isomorphism_fast.canonical_key(DECORATED_FP, (decoration_with([1], [8], HALF_EDGE_BASE),))
    for container in (list, tuple):
        assert isomorphism_fast.canonical_key(DECORATED_FP, (container(d),)) == expected


def test_canonical_key_data_float_active():
    # a float on an active half-edge is refused, not truncated to C int
    for container in (list, tuple):
        d = container(decoration_with([1.5], [8], HALF_EDGE_BASE))
        with pytest.raises(TypeError, match="^'float' object cannot be interpreted as an integer$"):
            isomorphism_fast.canonical_key(DECORATED_FP, (d,))
    with pytest.raises(TypeError, match="^'float' object cannot be interpreted as an integer$"):
        isomorphism_fast.canonical_key([1, 0], ([0.5, 0.25],), ())


def test_canonical_key_namedtuple_containers():
    from collections import namedtuple
    Data = namedtuple("Data", ["a"])
    d = array('i', HALF_EDGE_BASE)
    expected = isomorphism_fast.canonical_key(DECORATED_FP, (d,), (d,))
    assert isomorphism_fast.canonical_key(DECORATED_FP, Data(d), Data(d)) == expected
    assert isomorphism_fast.canonical_labelling(DECORATED_FP, Data(d), Data(d)) == \
        isomorphism_fast.canonical_labelling(DECORATED_FP, (d,), (d,))
    assert (isomorphism_fast.isomorphism(DECORATED_FP, Data(d), DECORATED_FP, Data(d), Data(d), Data(d))
            == isomorphism_fast.isomorphism(DECORATED_FP, (d,), DECORATED_FP, (d,), (d,), (d,)))


# The hook key (a function of the relabelling, here a linear subspace) and
# the method canonical_roots

def preserves_subspace(rows, a):
    r"""
    Return whether the half-edge permutation ``a`` preserves the row space
    of ``rows``: the columns are permuted by the edges ``a[2 * e] // 2``, and
    the rank of the union of the two row spaces is that of ``rows``.
    """
    from sage.all import QQ, matrix
    columns = len(rows[0])
    moved = [[0] * columns for _ in rows]
    for i, row in enumerate(rows):
        for e, x in enumerate(row):
            moved[i][a[2 * e] // 2] = x
    rank = matrix(QQ, rows).rank()
    return matrix(QQ, rows + moved).rank() == rank == matrix(QQ, moved).rank()


def random_subspace(rng, fp):
    r"""
    Return one or two random rows indexed by the edges of ``fp``, small
    integers on the active edges (the labels of the inactive edges are not
    determined by the isomorphism class, so the rows vanish there), not all
    zero.
    """
    active = [fp[2 * e] != -1 or fp[2 * e + 1] != -1 for e in range(len(fp) // 2)]
    while True:
        rows = [[rng.randint(-1, 2) if a else 0 for a in active] for _ in range(rng.randint(1, 2))]
        if all(any(row) for row in rows):
            return rows


def carry_subspace(rows, r):
    r"""
    Return ``rows`` with the column of the edge ``e`` moved to the edge
    ``r[2 * e] // 2``.
    """
    carried = []
    for row in rows:
        new = [0] * len(row)
        for e, x in enumerate(row):
            new[r[2 * e] // 2] = x
        carried.append(new)
    return carried


@pytest.mark.parametrize("build, edge_count, expected", [
    (lambda: map_from_fp(O7), 3, 2),
    (lambda: map_from_fp(O1), 12, 2),
], ids=["O7", "O1"])
def test_subspace_removes_automorphisms(build, edge_count, expected):
    m = build()
    rows = [[1] + [0] * (edge_count - 1)]
    direct = [a for a in m.automorphisms() if preserves_subspace(rows, a)]
    assert len(direct) == expected < len(m.automorphisms())
    assert automorphism_count(m, edge_subspace=rows) == expected
    assert sorted(tuple(a) for a in direct) == sorted(tuple(a) for a in m.automorphisms(**decoration_keywords(edge_subspace=rows)))


def random_fps(rng):
    r"""
    Return a random face permutation, connected or disconnected.
    """
    return random_connected_fp(rng) if rng.random() < 0.5 else random_disconnected_fp(rng)


@pytest.mark.parametrize("seed", range(6))
def test_key_random_maps(seed):
    rng = random.Random(seed)
    for _ in range(10):
        fp = random_fps(rng)
        m = OrientedMap(fp=fp)
        rows = random_subspace(rng, fp)
        key = subspace_key(rows)
        stabiliser = m.automorphisms(key=key)
        everything = m.automorphisms()
        assert stabiliser[0] == everything[0]
        assert [a for a in everything if tuple(a) in {tuple(b) for b in stabiliser}] == stabiliser
        r = m.canonical_labelling(key=key)
        assert {tuple(a) for a in stabiliser} == {tuple(a) for a in everything if key(array('i', [r[x] for x in a])) == key(r)}
        # the labelling is a minimiser of the key among those of the canonical form
        r0 = m.canonical_labelling()
        assert key(r) == min(key(array('i', [r0[x] for x in a])) for a in everything)
        assert m.canonical_form(key=key) == m.canonical_form() + (key(r),)
        generators = m.automorphism_generators(key=key)
        assert group_closure(generators, len(fp)) == {tuple(a) for a in stabiliser}
        assert (generators == []) == (len(stabiliser) == 1)
        assert all(preserves_subspace(rows, a) for a in stabiliser)
        # carried by a relabelling
        s = random_edge_relabelling(rng, fp)
        fp2 = relabel_arrays(fp, (), s)[0]
        m2 = OrientedMap(fp=fp2)
        rows2 = carry_subspace(rows, s)
        key2 = subspace_key(rows2)
        assert m2.canonical_form(key=key2) == m.canonical_form(key=key)
        ok, cert = m.is_isomorphic(m2, certificate=True, key=(key, key2))
        assert ok
        assert relabel_map(m, cert) == m2
        assert key(cert) == key2(array('i', range(len(fp2))))
        assert m.is_isomorphic(m2, key=(key, key2))
        assert len(m2.automorphisms(key=key2)) == len(stabiliser)


def test_key_is_isomorphic_errors():
    m = map_from_fp(O7)
    key = subspace_key([[1, 0, 0]])
    with pytest.raises(ValueError, match="key must be given for both maps or for none"):
        m.is_isomorphic(m, key=(key, None))
    with pytest.raises(ValueError, match="key must be a pair"):
        m.is_isomorphic(m, key=key)
    with pytest.raises(TypeError, match="key must be a pair of callables"):
        m.is_isomorphic(m, key=(key, 3))
    assert not m.is_isomorphic(m, key=(key, subspace_key([[0, 1, 1]])))
    assert m.is_isomorphic(m, key=(None, None))


@pytest.mark.parametrize("call", [
    lambda m, key: m.canonical_labelling(key=key),
    lambda m, key: m.canonical_form(key=key),
    lambda m, key: m.automorphisms(key=key),
    lambda m, key: m.automorphism_generators(key=key),
    lambda m, key: m.canonical_roots(key=key),
], ids=["canonical_labelling", "canonical_form", "automorphisms",
        "automorphism_generators", "canonical_roots"])
def test_key_not_callable(call):
    m = map_from_fp(O7)
    with pytest.raises(TypeError, match="key must be callable"):
        call(m, 3)


def test_key_none_is_unchanged():
    m = map_from_fp(O1)
    assert m.canonical_labelling(key=None) == m.canonical_labelling()
    assert m.canonical_form(key=None) == m.canonical_form()
    assert m.automorphisms(key=None) == m.automorphisms()
    assert m.automorphism_generators(key=None) == m.automorphism_generators()
    assert m.canonical_roots(key=None) == m.canonical_roots()


@pytest.mark.parametrize("fp, decorations, expected", [
    (O7, {}, [0, 1, 2, 3, 4, 5]),
    (O7, {'half_edge_data': [0, 1, 0, 1, 0, 2]}, [5]),
    (O7, {'half_edge_data': [0, 1, 0, 1, 0, 1]}, [1, 3, 5]),
    ("(0,1,2)(~0,~2,~1)(3)", {}, [6]),
    ("(0)(~0)(1)(~1)", {}, [0, 1, 2, 3]),
    ("(0)(~0)(1)(~1)", {'face_classes': [0, 1, 1, 0]}, [0, 1, 2, 3]),
], ids=["O7", "O7-marked", "O7-half", "triangle-folded", "two-edges", "two-edges-classes"])
def test_canonical_roots_examples(fp, decorations, expected):
    assert map_from_fp(fp).canonical_roots(**decorations) == expected


def test_canonical_roots_key():
    m = map_from_fp(O7)
    assert m.canonical_roots(key=lambda r: r[0] // 2) == [0, 1]
    assert m.canonical_roots(key=lambda r: 0) == m.canonical_roots()
    assert OrientedMap(fp=[]).canonical_roots() == []


@pytest.mark.parametrize("root_rule", ROOT_RULE_NAMES)
@pytest.mark.parametrize("seed", range(4))
def test_canonical_roots_random(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(8):
        fp = random_fps(rng)
        m = OrientedMap(fp=fp)
        data = [rng.randrange(2) for _ in fp] if rng.random() < 0.5 else None
        keywords = {} if data is None else {'half_edge_data': data}
        roots = m.canonical_roots(root_rule=root_rule, **keywords)
        r = m.canonical_labelling(root_rule=root_rule, **keywords)
        group = m.automorphisms(root_rule=root_rule, **keywords)
        assert roots == [h for h in range(len(fp)) if any(r[a[h]] == 0 for a in group)]
        s = random_edge_relabelling(rng, fp)
        fp2, carried = relabel_arrays(fp, () if data is None else (data,), s)
        keywords2 = {} if data is None else {'half_edge_data': carried[0]}
        roots2 = OrientedMap(fp=fp2).canonical_roots(root_rule=root_rule, **keywords2)
        assert roots2 == sorted(s[h] for h in roots)


# Sibling API: what veerer and surface-dynamics call

def sibling_case(rng):
    r"""
    Return a random map, connected or disconnected, a half-edge decoration
    and an edge decoration with two values each.
    """
    fp = random_fps(rng)
    return (OrientedMap(fp=fp), [rng.randrange(2) for _ in fp],
            [rng.randrange(2) for _ in range(len(fp) // 2)])


def sibling_relabelling(rng, m, hd, ed):
    r"""
    Return ``(r, carried, m2)``: a random relabelling of ``m`` that keeps
    its last edge active (a trailing inactive edge is not a valid map), the
    decorations ``hd`` and ``ed`` carried by it as keywords, and the
    relabelled map.
    """
    r = random_edge_relabelling(rng, list(m.face_permutation()))
    sigma = [r[2 * e] // 2 for e in range(len(r) // 2)]
    carried = carry_decorations(r, sigma, {'half_edge_data': hd, 'edge_data': ed})
    return r, carried, relabel_map(m, r)


def relabelled_fp(fp, r):
    r"""
    Return the face permutation ``fp`` relabelled by ``r``, without its
    trailing inactive edges, as the first entry of a canonical form.
    """
    new = relabel_arrays(fp, (), r)[0]
    while new and new[-1] == -1 and new[-2] == -1:
        del new[-2:]
    return new


def isolated_edges(k):
    r"""
    Return the map of ``k`` isolated edges.
    """
    return OrientedMap(fp="".join("(%d)(~%d)" % (i, i) for i in range(k)))


@pytest.mark.parametrize("seed", range(4))
def test_sibling_veerer_1(seed):
    r"""
    veerer 1: the canonical labelling of a disconnected map with
    ``half_edge_data`` and ``edge_data``: the relabelled map is the map of
    the canonical form, and the form does not change under a random
    relabelling.
    """
    rng = random.Random(seed)
    for _ in range(6):
        m, hd, ed = sibling_case(rng)
        form = m.canonical_form(half_edge_data=hd, edge_data=ed)
        r = m.canonical_labelling(half_edge_data=hd, edge_data=ed)
        fp = list(m.face_permutation())
        assert relabelled_fp(fp, r) == list(form[0])
        r2, carried, m2 = sibling_relabelling(rng, m, hd, ed)
        assert m2.canonical_form(**carried) == form
        fp2 = list(m2.face_permutation())
        assert relabelled_fp(fp2, m2.canonical_labelling(**carried)) == list(form[0])


@pytest.mark.parametrize("seed", range(4))
def test_sibling_veerer_2(seed):
    r"""
    veerer 2: ``automorphisms`` and ``automorphism_generators`` with the same
    decorations: the group generated is the list, and ``num_automorphisms``
    is its length.
    """
    rng = random.Random(seed)
    for _ in range(6):
        m, hd, ed = sibling_case(rng)
        for keywords in ({}, {'half_edge_data': hd, 'edge_data': ed}):
            automorphisms = m.automorphisms(**keywords)
            generators = m.automorphism_generators(**keywords)
            assert group_closure(generators, len(hd)) == {tuple(a) for a in automorphisms}
            assert m.num_automorphisms(**keywords) == len(automorphisms)


@pytest.mark.parametrize("seed", range(4))
def test_sibling_veerer_3(seed):
    r"""
    veerer 3: ``is_isomorphic`` with a certificate and the same decorations,
    as pairs, between a map and a relabelling.
    """
    rng = random.Random(seed)
    for _ in range(6):
        m, hd, ed = sibling_case(rng)
        r, carried, m2 = sibling_relabelling(rng, m, hd, ed)
        pairs = {'half_edge_data': (hd, carried['half_edge_data']), 'edge_data': (ed, carried['edge_data'])}
        ok, cert = m.is_isomorphic(m2, certificate=True, **pairs)
        assert ok
        assert relabel_map(m, cert) == m2
        for h in m.half_edges():
            assert carried['half_edge_data'][cert[h]] == hd[h]
            assert carried['edge_data'][cert[h] // 2] == ed[h // 2]
        changed = list(carried['half_edge_data'])
        h0 = next(m.half_edges())
        changed[cert[h0]] = 1 - changed[cert[h0]]
        assert not m.is_isomorphic(m2, half_edge_data=(hd, changed))


def test_sibling_veerer_4():
    r"""
    veerer 4: the subspace key (see ``test_automorphism_count_decorated``
    O9, O20, O21, ``test_subspace_removes_automorphisms`` and, for the
    minimum over ``r0 o a`` and the stabiliser, ``test_key_random_maps``),
    here through the enumeration of all the minimisers ``r0 o a``.
    """
    m = map_from_fp(O7)
    key = subspace_key([[1, 0, 0]])
    r0 = m.canonical_labelling()
    minimisers = [array('i', [r0[x] for x in a]) for a in m.automorphisms()]
    assert len({tuple(r) for r in minimisers}) == 6
    best = min(key(r) for r in minimisers)
    kept = [r for r in minimisers if key(r) == best]
    assert m.canonical_labelling(key=key) in kept
    assert len(kept) == m.num_automorphisms(key=key) == 2


@pytest.mark.parametrize("seed", range(4))
def test_sibling_veerer_5(seed):
    r"""
    veerer 5: the entry point ``isomorphism.canonical_key`` on ``array('i')``
    face permutation and decoration arrays, with ``relabelling=True`` and
    ``automorphism_count=True``, agrees with the methods of ``OrientedMap``
    (speed is measured by the benchmark).
    """
    rng = random.Random(seed)
    for _ in range(6):
        fp = random_fps(rng)
        m = OrientedMap(fp=fp)
        hd = [rng.randrange(2) for _ in fp]
        cl = [rng.randrange(3) for _ in fp]
        key, r, count = isomorphism_fast.canonical_key(
            array('i', fp), (array('i', hd),), (array('i', cl),), relabelling=True, automorphism_count=True)
        keywords = {'half_edge_data': hd, 'half_edge_classes': cl}
        assert r == m.canonical_labelling(**keywords)
        assert count == m.num_automorphisms(**keywords) == len(m.automorphisms(**keywords))
        s = random_edge_relabelling(rng, fp)
        fp2, (hd2, cl2) = relabel_arrays(fp, (hd, cl), s)
        assert isomorphism_fast.canonical_key(array('i', fp2), (array('i', hd2),), (array('i', cl2),)) == key


def test_sibling_surface_dynamics_1():
    r"""
    surface-dynamics 1: the canonical labelling of a disconnected fat graph
    and of a separatrix diagram (a face 2-colouring as ``face_data``).
    """
    rng = random.Random(0)
    m = map_from_fp(O3)
    assert not is_connected_fp(list(m.face_permutation()))
    a, colour, _ = cylinder_diagram(CD_I8)
    for m, decorations in ((m, {}), (a, {'face_data': colour})):
        form = m.canonical_form(**decorations)
        assert relabelled_fp(list(m.face_permutation()), m.canonical_labelling(**decorations)) == list(form[0])
        for _ in range(3):
            r, sigma = random_relabelling(m, rng)
            carried = carry_decorations(r, sigma, decorations)
            assert relabel_map(m, r).canonical_form(**carried) == form


def test_sibling_surface_dynamics_2():
    r"""
    surface-dynamics 2: cylinder diagrams through ``face_classes`` (see
    ``test_isomorphism_cylinder_diagram``).
    """
    forms = {}
    for name, s in [("I3a", CD_I3_A), ("I3b", CD_I3_B), ("I4b", CD_I4_B), ("I5b", CD_I5_B)]:
        m, colour, cylinder = cylinder_diagram(s)
        forms[name] = m.canonical_form(face_data=colour, face_classes=cylinder)
    assert forms["I3a"] == forms["I4b"]
    assert forms["I3b"] == forms["I5b"]
    assert forms["I3a"] != forms["I3b"]


@pytest.mark.parametrize("seed", range(4))
def test_sibling_surface_dynamics_3(seed):
    r"""
    surface-dynamics 3: ``fix_vertices`` and ``fix_faces`` as ``vertex_data``
    and ``face_data`` by distinct values (see
    ``test_automorphism_count_faces_fixed``), ``num_automorphisms``, and the
    actions of the automorphisms on the vertices and the faces.
    """
    rng = random.Random(seed)
    for _ in range(6):
        fp = random_fps(rng)
        m = OrientedMap(fp=fp)
        n = len(fp)
        nv, nf = m.num_vertices(), m.num_faces()
        automorphisms = m.automorphisms()
        vertex_actions = [m.automorphism_vertex_action(a) for a in automorphisms]
        face_actions = [m.automorphism_face_action(a) for a in automorphisms]
        for va, fa in zip(vertex_actions, face_actions):
            assert sorted(va) == list(range(nv))
            assert sorted(fa) == list(range(nf))
        for _ in range(5):
            i = rng.randrange(len(automorphisms))
            j = rng.randrange(len(automorphisms))
            a, b = automorphisms[i], automorphisms[j]
            ab = array('i', [a[b[h]] for h in range(n)])
            assert m.automorphism_vertex_action(ab) == array('i', [vertex_actions[i][x] for x in vertex_actions[j]])
            assert m.automorphism_face_action(ab) == array('i', [face_actions[i][x] for x in face_actions[j]])
        to_vertex = list(m.half_edge_to_vertex())
        to_face = list(m.half_edge_to_face())
        fixed = m.automorphisms(vertex_data=to_vertex)
        assert m.num_automorphisms(vertex_data=to_vertex) == len(fixed)
        assert all(m.automorphism_vertex_action(a) == array('i', range(nv)) for a in fixed)
        fixed = m.automorphisms(face_data=to_face)
        assert m.num_automorphisms(face_data=to_face) == len(fixed)
        assert all(m.automorphism_face_action(a) == array('i', range(nf)) for a in fixed)


def test_sibling_surface_dynamics_4():
    r"""
    surface-dynamics 4: ``h in m.canonical_roots()`` as the test of McKay
    augmentation (checked against the automorphisms in
    ``test_canonical_roots_random``), here on 8 isolated edges (`2^8 8!`
    automorphisms, none listed).
    """
    m = isolated_edges(8)
    assert m.canonical_roots() == list(range(16))
    assert m.canonical_roots(face_classes=[0] * 16) == list(range(16))


def test_num_automorphisms_examples():
    assert map_from_fp(O7).num_automorphisms() == 6
    assert map_from_fp(O7).num_automorphisms(key=lambda r: r[0] // 2) == 2
    assert OrientedMap(fp=[]).num_automorphisms() == 1
    count = isolated_edges(8).num_automorphisms()
    assert count == 2**8 * math.factorial(8) == 10321920
    assert type(count) is int


@pytest.mark.parametrize("root_rule", ROOT_RULE_NAMES)
@pytest.mark.parametrize("seed", range(3))
def test_num_automorphisms_random(seed, root_rule):
    rng = random.Random(seed)
    for _ in range(6):
        fp = random_fps(rng)
        m = OrientedMap(fp=fp)
        hd = [rng.randrange(2) for _ in fp]
        cl = [rng.randrange(3) for _ in fp]
        for keywords in ({}, {'half_edge_data': hd}, {'half_edge_classes': cl},
                         {'half_edge_data': hd, 'half_edge_classes': cl}):
            assert m.num_automorphisms(root_rule=root_rule, **keywords) == len(m.automorphisms(root_rule=root_rule, **keywords))
        key = subspace_key(random_subspace(rng, fp))
        assert m.num_automorphisms(root_rule=root_rule, key=key) == len(m.automorphisms(root_rule=root_rule, key=key))


def test_automorphism_actions_theta():
    m = map_from_fp("(0,1,2)(~0,~2,~1)")
    assert m.vertices() == [[0, 5], [1, 2], [3, 4]]
    assert m.faces() == [[0, 2, 4], [1, 5, 3]]
    for a, vertices, faces in [([1, 0, 5, 4, 3, 2], [1, 0, 2], [1, 0]),
                               ([2, 3, 4, 5, 0, 1], [1, 2, 0], [0, 1]),
                               ([5, 4, 3, 2, 1, 0], [0, 2, 1], [1, 0])]:
        assert m.automorphism_vertex_action(a) == array('i', vertices)
        assert m.automorphism_face_action(array('i', a)) == array('i', faces)
    empty = OrientedMap(fp=[])
    assert empty.automorphism_vertex_action(array('i')) == array('i', [0])
    assert empty.automorphism_face_action(array('i')) == array('i', [0])
    for bad in ([0, 1, 4, 5, 2, 3], [0, 1, 2], [0, 0, 1, 2, 3, 4], [0, 1, 2, 3, 4, 6], [1.0, 0, 5, 4, 3, 2]):
        with pytest.raises(ValueError):
            m.automorphism_vertex_action(bad)
        with pytest.raises(ValueError):
            m.automorphism_face_action(bad)
    assert m.automorphism_vertex_action([0, 1, 4, 5, 2, 3], check=False) is not None
