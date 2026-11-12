References
==========

.. _mosher-mcg-guide:

Mosher
------

Lee Mosher. *A user's guide to the mapping class group: once punctured
surfaces*.

.. _massuyeau-mcg-lectures:

Massuyeau
---------

Gwénaël Massuyeau. *Lectures on mapping class groups, braid groups and
formality*.

.. _tutte1963:

Tutte (1963)
------------

W. T. Tutte. *How to Draw a Graph*. Proceedings of the London
Mathematical Society, s3-13(1):743-767, 1963.

Establishes that a triconnected planar graph, drawn with its outer face
fixed as a convex polygon and every other vertex placed at a positive
weighted average of its neighbors, always yields a planar straight-line
drawing with every face convex. Used throughout
:mod:`combisurf.layout.tutte_barycentric` as the correctness argument
behind :func:`~combisurf.layout.tutte_barycentric.tutte_layout`, and
throughout :mod:`combisurf.layout.gluing` as the guarantee that keeps a
composed piece's entire drawing within the convex hull of its own fixed
boundary.

.. _hopcroft-kahn1992:

Hopcroft & Kahn (1992)
-----------------------

John E. Hopcroft and Peter C. Kahn. *A Paradigm for Robust Geometric
Algorithms*. Algorithmica, 7(1-6):339-380, 1992.

Cited in :mod:`combisurf.layout.tutte_barycentric` for the generalization
of Tutte's theorem to asymmetric edge weights (the two directions of an
edge may carry different weights, so long as both are positive).

.. _eades1992:

Eades (1992)
------------

Peter Eades. *Drawing Free Trees*. Bulletin of the Institute for
Combinatorics and its Applications, 5:10-36, 1992.

Source of the radial tree-drawing algorithm (``DrawSubTree1``)
implemented as :func:`~combisurf.layout.tree_pieces.eades_radial_layout`.

.. _eades-garvan1996:

Eades & Garvan (1996)
----------------------

Peter Eades and Patrick Garvan. *Drawing Stressed Planar Graphs in Three
Dimensions*. In Graph Drawing (GD 1995), Lecture Notes in Computer
Science 1027, pages 212-223. Springer, 1996.

Cited in :mod:`combisurf.layout.tutte_barycentric` for the observation
that plain uniform-weight Tutte embeddings can have an edge-length ratio
growing exponentially in the number of vertices, motivating the
alternative weighting schemes implemented there.

.. _floater2003:

Floater (2003)
---------------

Michael S. Floater. *Mean Value Coordinates*. Computer Aided Geometric
Design, 20(1):19-27, 2003.

Source of the mean value coordinate weights implemented as
:func:`~combisurf.layout.tutte_barycentric.mean_value_weights`.

.. _chiu-eppstein-goodrich2023:

Chiu, Eppstein & Goodrich (2023)
---------------------------------

Man-Kwun Chiu, David Eppstein, and Michael T. Goodrich. *Manipulating
Weights to Improve Stress-Graph Drawings of 3-Connected Planar Graphs*.
In Graph Drawing and Network Visualization (GD 2023). arXiv:2307.10527.

Source of the edge-length-ratio drawing-quality metric
(:func:`~combisurf.layout.tutte_barycentric.edge_length_ratio`) and the
BFS-depth weighting scheme
(:func:`~combisurf.layout.tutte_barycentric.bfs_depth_weights`,
:func:`~combisurf.layout.tutte_barycentric.bfs_depth_layout`).

.. _fenwick-1994:

Fenwick (1994)
--------------

P. M. Fenwick. *A new data structure for cumulative frequency tables*.
Software: Practice and Experience, 24(3):327-336, 1994.

The binary indexed tree behind the sweeps
:func:`~combisurf.crossing_arcs.crossing_arcs_sweep_sorted`,
:func:`~combisurf.crossing_arcs.startpoint_sweep_sorted` and
:func:`~combisurf.crossing_arcs.startpoint_sweep_weighted`, where it answers
each prefix sum and each update in logarithmic time.

.. _chazelle1986:

Chazelle (1986)
---------------

B. Chazelle. *Reporting and counting segment intersections*. Journal of
Computer and System Sciences, 32:156-182, 1986.

Cited in :mod:`combisurf.crossing_arcs` for where its method comes from:
counting the crossing chords of a circle is the easy case of counting
segment intersections, because the cyclic order of the endpoints already
gives the sweep order.

.. _despre-lazarus-2019:

Despré & Lazarus (2019)
-----------------------

V. Despré and F. Lazarus. *Computing the geometric intersection number of
curves*. Journal of the ACM, 66(6), Article 45, 2019. arXiv:1511.09327.

Directly on the problem solved by
:meth:`~combisurf.geometric_intersection.GeometricIntersection.geometric_intersection`,
the geometric intersection numbers of curves on a surface given as words.

.. _birman-series1984:

Birman & Series (1984)
----------------------

J. S. Birman and C. Series. *An algorithm for simple curves on surfaces*.
Journal of the London Mathematical Society (2), 29:331-342, 1984.

Source of the two genus 2 examples, pages 336-337, in the doctests of
:meth:`~combisurf.geometric_intersection.GeometricIntersection.geometric_intersection`.

.. _lapointe2019:

Lapointe (2019)
---------------

M. Lapointe. *Number of orbits of discrete interval exchanges*. Discrete
Mathematics & Theoretical Computer Science, 21(3), Paper No. 13, 16 pages,
2019.

Source of the simplicity criterion checked in the doctests of
:meth:`~combisurf.geometric_intersection.GeometricIntersection.geometric_intersection`:
a primitive positive word gives a simple curve exactly when its
Burrows-Wheeler transform is non-increasing.

.. _ukkonen-1995:

Ukkonen (1995)
--------------

E. Ukkonen. *On-line construction of suffix trees*. Algorithmica,
14(3):249-260, 1995.

Source of the construction of :class:`~combisurf.conjugate_tree.ConjugateTree`
(the active point, ``canonize`` and ``test_and_split``), here applied to
words read cyclically.

.. _cohen-lustig-1987:

Cohen Lustig (1987)
-------------------

M. Cohen, M. Lustig. *Paths of geodesics and geometric intersection numbers. I*.
Combinatorial group theory and topology, Sel. Pap. Conf., Alta/Utah 1984, Ann. Math. Stud. (111):479-500 1987.

.. _dubois-2024:

Dubois (2024)
-------------

L. Dubois. *Making Multicurves Cross Minimally on Surfaces*
ESA 2024: 50:1-50:15
