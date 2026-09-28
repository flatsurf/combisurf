r"""
Cyclic ordering of conjugates on the boundary at infinity and its picture

These helpers used to be methods of
:class:`~combisurf.geometric_intersection.GeometricIntersectionPairing`. They
take such a pairing ``gi`` as first argument and work on its reduced map
``gi._r``, a map with a single vertex.
"""

from array import array

from combisurf.word import word_free_group_inverse
from combisurf.conjugate_tree import ConjugateTree


def _conjugate_sort(gi, words):
    r"""
    Return the cyclic ordering on the boundary at infinity of the conjugates of
    ``words``, given as cyclically reduced words on the reduced map ``gi._r``.

    Nothing is checked. See :func:`conjugate_sort`.
    """
    angles = gi._angles
    T = ConjugateTree(len(angles))
    for i, w in enumerate(words):
        ans = T.process(list(w))
        if ans <= 0:
            # was already given in T
            raise ValueError(f"conjugate words at position {-ans} and {i}")

    # NOTE: below a node, the angles are measured from the half-edge
    # through which the curve came in, the inverse of the last letter read
    n = len(angles)
    pivot = array('i', [angles[b ^ 1] for b in range(n)])
    word_indices, word_shifts = T.sorted_leaves_as_conjugates(angles, pivot)

    return list(word_indices), list(word_shifts)


def conjugate_sort(gi, words, check=True):
    r"""
    Return the cyclic ordering on the boundary at infinity of the conjugates of
    the closed walks ``words``.

    INPUT:

    - ``gi`` -- a :class:`~combisurf.geometric_intersection.GeometricIntersectionPairing`

    - ``words`` -- a list of closed walks on the map of ``gi``, pairwise non
      conjugate once reduced

    - ``check`` -- boolean (default: ``True``); whether to check that the
      words are closed walks on the map of ``gi``. With ``check=False`` they
      must already be given as ``array('i')``

    OUTPUT: a pair of lists ``(indices, shifts)`` of the same length, the total
    length of the reduced words. The conjugate at position ``p`` in the cyclic
    ordering is the reduced word ``indices[p]`` rotated by ``shifts[p]``.

    The words are first reduced onto the reduced map of ``gi`` (cyclically
    reduced and, if the map of ``gi`` has more than one vertex or unpunctured
    faces, mapped to the one-vertex reduced map). The shifts refer to these
    reduced words, not to the input.

    EXAMPLES::

        sage: from combisurf import OrientedMap
        sage: from combisurf.geometric_intersection import GeometricIntersectionPairing
        sage: octagon = OrientedMap(vp="(0,1,2,3,~0,~1,~2,~3)")
        sage: gi = GeometricIntersectionPairing(octagon, punctured_faces=True)

    A simple example (that turns out to be equivalent to lexicographic sort
    of conjugates)::

        sage: w = [0, 2, 0, 0, 2]
        sage: words, shifts = conjugate_sort(gi, [w])
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

    A more involved example::

        sage: W = [[0], [2], [0, 2, 0, 3]]
        sage: words, shifts = conjugate_sort(gi, W)
        sage: words
        [2, 0, 2, 2, 1, 2]
        sage: shifts
        [2, 0, 0, 1, 0, 3]
    """
    words = [gi._reduce(w, check=check) for w in words]
    return _conjugate_sort(gi, words)


def conjugate_plot(gi, words):
    r"""
    Return a picture of the conjugates of ``words`` and of their inverses on the
    boundary at infinity, each conjugate being joined to the one of the inverse
    at the other end of the same lift.

    INPUT:

    - ``gi`` -- a :class:`~combisurf.geometric_intersection.GeometricIntersectionPairing`

    - ``words`` -- a list of closed walks on the map of ``gi``

    The labels are the conjugates of the words reduced on the reduced map of
    ``gi`` (see :func:`conjugate_sort`).

    EXAMPLES::

        sage: from combisurf import OrientedMap
        sage: from combisurf.geometric_intersection import GeometricIntersectionPairing
        sage: m = OrientedMap(fp="(0,1,2)(~0,3,4)(~1,~4,~5)(~2,5,~3)")
        sage: gi = GeometricIntersectionPairing(m, punctured_faces=True)
        sage: conjugate_plot(gi, ["0,1,2,0,~4,~3"])
        Graphics object consisting of ... graphics primitives
    """
    from sage.rings.complex_double import CDF
    from sage.plot.colors import rainbow
    from sage.plot.text import text
    from sage.plot.point import point2d
    from sage.plot.circle import circle
    from sage.plot.line import line2d

    n = len(words)
    words = [gi._reduce(w, check=True) for w in words]
    words_with_inverse = list(words) + [word_free_group_inverse(w) for w in words]
    l = sum(len(w) for w in words_with_inverse)
    colors = rainbow(n, 'rgbtuple')
    word_indices, word_shifts = _conjugate_sort(gi, words_with_inverse)
    assert len(word_indices) == len(word_shifts) == l
    conj_to_pos = [[None] * len(w) for w in words_with_inverse]
    for pos, (i, k) in enumerate(zip(word_indices, word_shifts)):
        conj_to_pos[i][k] = pos

    G = circle((0, 0), 1, color='black')
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
