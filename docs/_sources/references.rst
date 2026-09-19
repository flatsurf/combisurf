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
