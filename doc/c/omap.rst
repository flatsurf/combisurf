``omap.h``: oriented maps in C
==============================

Files: ``combisurf/src/omap.{h,c}``

An oriented map on ``n`` half-edges (``n`` even) is given by its vertex and
face permutations ``vp`` and ``fp``, arrays of ``n`` ints with ``-1`` on the
inactive half-edges. The twin of the half-edge ``h`` is ``h ^ 1``. An edge is
active when at least one of its half-edges is; a folded edge has one active
half-edge. ``vp`` and ``fp`` must be permutations of the active half-edges,
related by ``vp[fp[t]] = h`` for every active half-edge ``h``, where ``t`` is
its twin, or ``h`` itself on a folded edge. The functions of this header do
not check them (the Cython layer does).

A map can be decorated by two kinds of integer arrays on its half-edges,
*data* and *classes*, which behave differently under isomorphism. The values
of a data array are kept by an isomorphism, while the values of a class array
are kept up to a relabelling: only the partition
:math:`\{h : c[h] = i\}_i` of the half-edges of a class array ``c`` matters,
whatever the order of its parts. A coloring of the half-edges is a typical
data array; a coloring up to a permutation of the colors is a typical class
array. A decorated map has ``nd`` data arrays, given one after the other in
an array of ``nd * n`` ints (the array ``i`` at ``data + i * n``), and ``nc``
class arrays, given likewise in an array of ``nc * n`` ints. The values on
the inactive half-edges are never read. The values of a class array must lie
in :math:`0, \ldots, n - 1` (renumber them before the call, for instance by
first appearance).

The functions that need a work buffer come in pairs, as in :doc:`perm`: the
raw form ``_X(..., work)`` takes its buffer from the caller, sized by its
``*_WORK_SIZE`` macro, and allocates nothing; the public form ``X(...)``
allocates its buffer through :c:func:`combisurf_malloc` and frees it before
returning. The output arrays are always handed in by the caller. The Cython
module calls only the raw forms, under ``sig_on()``, so that an interrupt,
which jumps out of the C call, leaves no buffer allocated by the C code.

The algorithms are described in :doc:`omap_algorithms`.

Types
-----

.. c:struct:: omap_struct

   An oriented map on ``n`` half-edges.

   .. c:member:: int *vp
   .. c:member:: int *fp

      The vertex and face permutations, ``n`` ints each, ``-1`` on the
      inactive half-edges.

   .. c:member:: int n

      The number of half-edges, even.

.. c:type:: omap_struct omap_t[1]

   An :c:struct:`omap_struct` as an array of one element, so that an
   ``omap_t`` declared on the stack is passed by pointer.

.. c:struct:: domap_struct

   An oriented map decorated by ``nd`` data arrays and ``nc`` class arrays.
   Its first member is ``map``, so that a pointer to a ``domap_struct`` is
   also a pointer to an :c:struct:`omap_struct` (C11 6.7.2.1) and the
   ``omap_*`` functions apply to ``(const omap_struct *) m`` for a
   :c:type:`domap_t` ``m``.

   .. c:member:: omap_struct map

      The map.

   .. c:member:: const int *data

      The ``nd`` data arrays of ``map.n`` ints, one after the other.

   .. c:member:: const int *classes

      The ``nc`` class arrays of ``map.n`` ints, one after the other.

   .. c:member:: int nd

      The number of data arrays.

   .. c:member:: int nc

      The number of class arrays.

.. c:type:: domap_struct domap_t[1]

   A :c:struct:`domap_struct` as an array of one element, so that a
   ``domap_t`` declared on the stack is passed by pointer.

Basic functions
---------------

.. c:function:: int _omap_num_edges(const int *p, int n)
.. c:function:: int omap_num_edges(const omap_t m)

   Return the number of active edges of the map ``m``: the edges with at least
   one active half-edge.

   The raw form takes a permutation ``p`` of ``n`` ints, and accepts either
   ``vp`` or ``fp``, which have the same active half-edges.

.. c:function:: void _omap_vp_from_fp(int *vp, const int *fp, int n)

   Write into ``vp`` (``n`` ints) the vertex permutation of the map of face
   permutation ``fp``: ``vp[fp[t]] = h`` for every active half-edge ``h``,
   ``t`` being its twin, or ``h`` itself on a folded edge, and ``-1`` on the
   inactive half-edges.

.. c:function:: int omap_next_at_vertex(const omap_t m, int h)
.. c:function:: int omap_previous_at_vertex(const omap_t m, int h)
.. c:function:: int omap_next_in_face(const omap_t m, int h)
.. c:function:: int omap_previous_in_face(const omap_t m, int h)

   Return the half-edge after ``h`` around its vertex (``vp[h]``), the one
   before ``h`` around its vertex (the preimage of ``h`` under ``vp``,
   ``fp[t]`` for the twin ``t`` of ``h``, or ``h`` itself on a folded edge),
   the one after ``h`` in its face (``fp[h]``), and the one before ``h`` in
   its face (the preimage of ``h`` under ``fp``, the twin of ``vp[h]``, or
   ``vp[h]`` itself on a folded edge).

   Return ``-1`` when ``h`` is invalid: out of :math:`0, \ldots, n - 1` or
   inactive.

Isomorphisms
------------

The Cython module :mod:`combisurf.isomorphism` wraps the functions below.

.. _omap-relabelling:

Relabelling from a root
^^^^^^^^^^^^^^^^^^^^^^^

The building block of the isomorphism code is the relabelling of a map from
an active half-edge, its root: a breadth-first walk along the faces. The
half-edges are taken from a queue that initially holds the root; the face of
the half-edge taken is walked from it, and each half-edge met on an edge
without labels gets the next label ``k`` and its twin the label ``k + 1``,
``k`` going up by 2 each time (starting at 0); the twin is pushed when it is
active. The walk reaches the edges of the component of the root.

.. c:function:: int _omap_relabelling_from(int *r, int *rinv, const int *p, int n, int root, int *queue)
.. c:function:: int omap_relabelling_from(int *r, int *rinv, const omap_t m, int root)

   Compute the relabelling of ``m`` from the active half-edge ``root``.

   Write into ``r`` (``n`` ints) the label of each half-edge, ``-1`` on the
   half-edges of the edges not reached. Return the number ``k`` of labels
   used: the labels are given in pairs to the two half-edges of an edge, so
   that ``k`` is even, twice the number of edges reached, folded edges
   included, and the labels used are :math:`0, \ldots, k - 1`. Write into
   ``rinv`` (``n`` ints) the inverse of ``r`` on these labels only:
   ``rinv[i]`` is the half-edge of label ``i``, or ``-1`` when that half-edge
   is inactive (the twin of a folded edge); the entries of ``rinv`` from
   ``k`` on are not written. The map is connected if and only if ``k`` is
   twice its number of active edges.

   The raw form walks the cycles of the permutation ``p`` of ``n`` ints it is
   given, and accepts either ``vp`` or ``fp``: with ``fp`` it walks the faces
   as the public form does, with ``vp`` the vertices. ``queue`` is a work
   buffer of ``OMAP_RELABELLING_WORK_SIZE(n)`` ints.

.. c:macro:: OMAP_RELABELLING_WORK_SIZE(n)

   The size of the work buffer of :c:func:`_omap_relabelling_from`, in ints,
   ``n``: the queue holds the root and one half-edge per edge reached, at
   most ``n / 2 + 1 <= n`` entries.

.. _omap-keys:

Keys
^^^^

The relabelling ``(r, rinv)`` of a decorated map from a root, using ``k``
labels, gives its key: the face permutation, the data and the classes read
in the new labels, packed into one array. Let ``na`` be the number of active
half-edges among the labels used, the labels ``i < k`` with
``rinv[i] != -1``, that is ``k`` minus the number of folded edges reached.
The key is made of ``k + (nd + nc) na`` ints:

- the face permutation in the new labels (``-1`` at the labels of inactive
  half-edges), ``k`` ints;
- then each data array read in the order of the new labels of the active
  half-edges, ``na`` ints each;
- then each class array read in the same order and renumbered
  ``0, 1, 2, ...`` by first appearance, ``na`` ints each.

Keys are compared as sequences of signed ints.

.. c:function:: int _domap_iso_key(int *key, const domap_t m, const int *r, const int *rinv, int k, int *names)
.. c:function:: int domap_iso_key(int *key, const domap_t m, const int *r, const int *rinv, int k)

   Write into ``key`` the key of ``m`` under the relabelling ``(r, rinv)``
   using ``k`` labels, as computed by :c:func:`omap_relabelling_from`.

   ``key`` has room for ``DOMAP_ISO_KEY_SIZE(n, nd, nc)`` ints.

   Return the length of the key, ``k + (nd + nc) na``.

   In the raw form, ``names`` is a work buffer of
   ``DOMAP_ISO_KEY_WORK_SIZE(n)`` ints, which must hold ``-1`` everywhere on
   entry and holds ``-1`` everywhere again on return, so that it serves
   several calls without being reset; the public form fills its own buffer
   with ``-1``.

.. c:macro:: DOMAP_ISO_KEY_SIZE(n, nd, nc)

   An upper bound on the length of a key, in ints, ``n (1 + nd + nc)``: the
   length ``k + (nd + nc) na`` of a key has ``k <= n`` (two labels per edge
   reached, at most ``n / 2`` edges) and ``na <= n``.

.. c:macro:: DOMAP_ISO_KEY_WORK_SIZE(n)

   The size of the work buffer of :c:func:`_domap_iso_key`, in ints, ``n``:
   the new number of each class value, which lies in
   :math:`0, \ldots, n - 1`.

.. _omap-root-rules:

Rules of the roots
^^^^^^^^^^^^^^^^^^

The roots tried by :c:func:`domap_canonical` and :c:func:`domap_iso_search`
are restricted by a rule (:c:enum:`domap_iso_root_rule`), the argument
``root_rule`` of the functions of :mod:`combisurf.isomorphism`.

For an active half-edge ``h``, ``vdeg(h)`` and ``fdeg(h)`` are the degrees of
its vertex and of its face, and its twin ``t`` is ``h ^ 1``, or ``h`` itself
on a folded edge. A rule gives each active half-edge an invariant:

- the pair ``(vdeg(h), fdeg(h))`` under the half-edge rules;
- the tuple ``(vdeg(h), fdeg(h), vdeg(t), fdeg(t), loop, monogon)`` under the
  edge rules, ``loop`` (``monogon``) being 1 when ``h`` and ``t`` are at the
  same vertex (in the same face) and 0 otherwise.

The half-edges of equal invariant form a class, and only the half-edges of
one class are tried: under ``*_MAX`` the class of the largest invariant
(tuples compared lexicographically), under ``*_FEWEST`` the class with the
fewest half-edges, ties broken by the largest invariant. The invariants are
computed once per call on the whole map. The base rules read ``vp`` and
``fp`` alone, never the decorations.

The flag :c:enumerator:`DOMAP_ISO_ROOTS_DATA`, combined with a rule other
than :c:enumerator:`DOMAP_ISO_ROOTS_ALL` by ``|``, appends to the invariant
the values of the data arrays (signed ints, never the class arrays, whose
values are renumbered and so not invariant): ``data_0[h], ...,
data_{nd-1}[h]`` under the half-edge rules, followed by ``data_0[t], ...,
data_{nd-1}[t]`` under the edge rules; with ``nd = 0`` it changes nothing.

Under :c:enumerator:`DOMAP_ISO_ROOTS_ALL` every active half-edge is a root,
and ``vp`` is not read.

.. c:enum:: domap_iso_root_rule

   The rules restricting the roots: a base rule (0 to 4), possibly combined
   with the flag :c:enumerator:`DOMAP_ISO_ROOTS_DATA` by ``|``. The values of
   the enumerators below are those of the constants ``ROOTS_ALL``, ...,
   ``ROOTS_EDGE_FEWEST`` and ``ROOTS_DATA`` of :mod:`combisurf.isomorphism`.
   The functions of this header read the base rule as
   ``rule & ~DOMAP_ISO_ROOTS_DATA``; another value of the base rule is not
   checked.

   .. c:enumerator:: DOMAP_ISO_ROOTS_ALL = 0

      Every active half-edge is a root.

   .. c:enumerator:: DOMAP_ISO_ROOTS_HALF_EDGE_MAX = 1

      The class of the largest half-edge invariant.

   .. c:enumerator:: DOMAP_ISO_ROOTS_HALF_EDGE_FEWEST = 2

      The class with the fewest half-edges under the half-edge invariant,
      ties broken by the largest invariant.

   .. c:enumerator:: DOMAP_ISO_ROOTS_EDGE_MAX = 3

      The class of the largest edge invariant.

   .. c:enumerator:: DOMAP_ISO_ROOTS_EDGE_FEWEST = 4

      The class with the fewest half-edges under the edge invariant, ties
      broken by the largest invariant.

   .. c:enumerator:: DOMAP_ISO_ROOTS_DATA = 8

      The flag appending the data values to the invariant of a base rule.

Connected maps
^^^^^^^^^^^^^^

.. c:function:: int _domap_canonical(int *key, int *relabelling, int *roots, int *nroots, const domap_t m, domap_iso_root_rule rule, int *work)
.. c:function:: int domap_canonical(int *key, int *relabelling, int *roots, int *nroots, const domap_t m, domap_iso_root_rule rule)

   Compute the canonical key of the connected decorated map ``m``: try every
   active half-edge of ``m`` of the class of the rule as a root (every active
   half-edge under :c:enumerator:`DOMAP_ISO_ROOTS_ALL`), in increasing order,
   and write into ``key`` (``DOMAP_ISO_KEY_SIZE(n, nd, nc)`` ints) the
   minimal key.

   When ``relabelling`` is not ``NULL`` (``n`` ints), write there the
   relabelling from the first root of minimal key (``-1`` on the half-edges
   of the inactive edges). When ``roots`` is not ``NULL`` (``n`` ints), write
   there the roots of minimal key in increasing order. When ``nroots`` is not
   ``NULL``, store there their number, which is the number of automorphisms
   of the decorated map under every rule (0 when the map has no active
   half-edge): an automorphism of a connected map is determined by the image
   of one half-edge and preserves the invariants, so the roots of minimal key
   form one orbit of the group.

   Return the length of the key (0 when the map has no active half-edge), or
   :c:enumerator:`DOMAP_ISO_DISCONNECTED` when the walk from the first root
   tried does not reach every active edge; nothing is written into the output
   buffers in that case. The first root tried is the smallest active
   half-edge under :c:enumerator:`DOMAP_ISO_ROOTS_ALL`, and the smallest
   half-edge of the class of the rule under the other rules; connectedness
   does not depend on the root, so that this walk alone decides it.

   In the raw form, ``work`` is a buffer of
   ``DOMAP_CANONICAL_WORK_SIZE(n, nd, nc)`` ints.

   The rejection of the roots, the skipped roots and the cost are in
   :ref:`omap-algorithms-canonical`.

.. c:macro:: DOMAP_CANONICAL_WORK_SIZE(n, nd, nc)

   The size of the work buffer of :c:func:`_domap_canonical`, in ints,
   ``9 n + 1 + DOMAP_ISO_KEY_SIZE(n, nd, nc)``:

   - the relabelling, its inverse, the queue and the numbering of the class
     values (``4 n``);
   - a candidate key (:c:macro:`DOMAP_ISO_KEY_SIZE`);
   - under a rule other than :c:enumerator:`DOMAP_ISO_ROOTS_ALL`, the active
     half-edges sorted by invariant, the class of each half-edge and the
     start of each class in that order (``3 n + 1``);
   - then the relabelling of the best root and the parents of the
     union-find of the orbits of the automorphisms found (``2 n``).

   Under :c:enumerator:`DOMAP_ISO_ROOTS_DATA` the order-preserving numbers of
   the data values (``nd n``) live in the candidate key, which is not used
   before the classes are computed, so that the flag needs no more room.

.. c:function:: int _domap_invariants_differ(const domap_t m1, const domap_t m2, int *work)
.. c:function:: int domap_invariants_differ(const domap_t m1, const domap_t m2)

   Return 1 when an invariant of the decorated maps ``m1`` and ``m2``
   differs, so that they are not isomorphic, and 0 otherwise.

   The two maps must have as many data arrays (``m1->nd == m2->nd``), as many
   class arrays (``m1->nc == m2->nc``) and as many active half-edges. An
   isomorphism is a bijection of the active half-edges that keeps the data
   values, and that relabels the values of each class array, sending a class
   to a class of the same size.

   In the raw form, ``work`` is a buffer of ``DOMAP_INVARIANTS_WORK_SIZE(n)``
   ints, with ``n = max(m1->map.n, m2->map.n)``, read only when
   ``m1->nc > 0`` (it may be ``NULL`` otherwise); the public form allocates
   it only then.

   The invariants compared, and why the degrees are not, are in
   :ref:`omap-algorithms-invariants`.

.. c:macro:: DOMAP_INVARIANTS_WORK_SIZE(n)

   The size of the work buffer of :c:func:`_domap_invariants_differ`, in
   ints, ``2 n + 1``: the size of each class and the number of classes of
   each size.

.. c:enumerator:: DOMAP_ISO_DISCONNECTED = -1

   The return value of :c:func:`domap_canonical` on a map that is not
   connected.

.. _omap-blocks:

Maps that may be disconnected
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A map that may be disconnected is labelled component after component: each
component is relabelled from a root as above, its labels following those of
the components labelled before, and gives a block. The block of a component
labelled from a root, with ``k`` labels, is made of:

- ``k``, its number of labels;
- its face permutation in its local labels :math:`0, \ldots, k - 1` (``-1``
  at the labels of inactive half-edges);
- each data array read in the order of the local labels of its active
  half-edges;
- then each class array read in the same order and renumbered by first
  appearance, continuing the numbering of the components labelled before.

A sequence of blocks is stored flat, the blocks one after the other; since
equal ``k`` and equal face permutations give blocks of equal lengths, the
lexicographic order of the flat sequences of ints is the order of the
sequences of blocks. A connected map gives a single block, which is its key
preceded by ``k``.

.. c:function:: int _domap_iso_search(int *sequence, int *length, int *relabelling, const domap_t m, domap_iso_root_rule rule, int all_leaves, domap_iso_leaf_callback callback, void *arg, int *work)
.. c:function:: int domap_iso_search(int *sequence, int *length, int *relabelling, const domap_t m, domap_iso_root_rule rule, int all_leaves, domap_iso_leaf_callback callback, void *arg)

   Search the lexicographically minimal sequence of blocks over the
   labellings of ``m`` component after component: from the components
   labelled so far, the roots of the next component are the half-edges not
   labelled yet of the class of the rule (every active half-edge not labelled
   yet under :c:enumerator:`DOMAP_ISO_ROOTS_ALL`), and the component gets the
   labels following those used. When several roots give the minimal block,
   every one of them is followed if ``all_leaves`` is nonzero; otherwise only
   the first one if ``nc == 0``, and the first one of each extension of the
   numbering of the class values if ``nc > 0``. Tied roots of different
   extensions are all followed, even when their subtrees are images of each
   other. The search and its cuts are in :ref:`omap-algorithms-iso-search`.

   Any map is accepted, connected or not. The empty map gives the empty
   sequence and one leaf.

   The minimal sequence is written into ``sequence``
   (``DOMAP_ISO_SEQUENCE_SIZE(n, nd, nc)`` ints) and its length into
   ``*length``. When ``relabelling`` is not ``NULL`` (``n`` ints), the
   relabelling of the first leaf reaching the minimal sequence is written
   there (``-1`` on the half-edges of the inactive edges). When ``callback``
   is not ``NULL``, it is called with ``arg`` on every leaf that is smaller
   than the best one so far (``first = 1``) and, when ``all_leaves`` is
   nonzero, on every leaf equal to the best one so far (``first = 0``). The
   leaves reported since the last call with ``first = 1`` are then the leaves
   reaching the minimal sequence, in the order of the search, as many as the
   automorphisms when ``all_leaves`` is nonzero.

   Return 0, or the nonzero value returned by the callback, which stops the
   search and leaves the outputs undefined.

   In the raw form, ``work`` is a buffer of
   ``DOMAP_ISO_SEARCH_WORK_SIZE(n, nd, nc)`` ints.

.. c:macro:: DOMAP_ISO_SEQUENCE_SIZE(n, nd, nc)

   An upper bound on the length of a sequence of blocks, in ints,
   ``DOMAP_ISO_KEY_SIZE(n, nd, nc) + n / 2``. Let ``ne`` be the number of
   active edges of the map and ``na`` its number of active half-edges. A
   component of ``k`` labels, ``na_c`` of which are on active half-edges,
   gives a block of ``1 + k + (nd + nc) na_c`` ints; the ``k`` of the
   components add up to ``2 ne``, their ``na_c`` to ``na``, and there are at
   most ``ne`` components, each having an active edge. A sequence thus has at
   most ``ne + 2 ne + (nd + nc) na`` ints, which ``ne <= n / 2`` and
   ``na <= n`` bound by ``n (1 + nd + nc) + n / 2``.

.. c:macro:: DOMAP_ISO_SEARCH_WORK_SIZE(n, nd, nc)

   The size of the work buffer of :c:func:`_domap_iso_search`, in ints:

   - the global and local relabellings with their inverses and the queue
     (``5 n``);
   - the numberings of the class values (``nc n``);
   - the state of each of the at most ``n / 2 + 1`` levels of the search
     (its root, the root of the best leaf at that level, its first label, the
     position and the length of its block, and the counts of class values
     before it: ``(5 + nc) (n / 2 + 1)``);
   - the counts of class values (``nc``);
   - the current sequence (:c:macro:`DOMAP_ISO_SEQUENCE_SIZE`);
   - one candidate block (``DOMAP_ISO_KEY_SIZE + 1``);
   - for the rules other than :c:enumerator:`DOMAP_ISO_ROOTS_ALL`, the active
     half-edges sorted by invariant, the class of each half-edge, the start
     of each class in that order and the number of half-edges not labelled
     yet in each class (``4 n + 1``).

   Under :c:enumerator:`DOMAP_ISO_ROOTS_DATA` the order-preserving numbers of
   the data values (``nd n``) live in the candidate block, which is not used
   before the classes are computed, so that the flag needs no more room.

   The sizes of this header are computed in ``int``. This one, the largest,
   overflows first, for about :math:`2^{31} / (14 + 2 nd + 3.5 nc)`
   half-edges; the caller checks the sizes first, and the Cython module
   raises an ``OverflowError`` on a map too large.

.. c:type:: domap_iso_leaf_callback

   The callback of :c:func:`domap_iso_search`, called on a leaf:

   .. code-block:: c

      typedef int (*domap_iso_leaf_callback)(const int *leaf, int first, void *arg);

   ``leaf`` is the relabelling of the leaf (``n`` ints, ``-1`` on the
   half-edges of the inactive edges), valid during the call only. ``first``
   is 1 when the leaf is smaller than every leaf reported before, which the
   caller discards, and 0 when it is equal to the last leaf reported with
   ``first = 1``. ``arg`` is the argument given to
   :c:func:`domap_iso_search`.

   A nonzero return value stops the search, which returns it; the Cython
   layer returns -1 with a Python exception set.

Maps with classes
^^^^^^^^^^^^^^^^^

Maps that may be disconnected, with classes: the automorphisms of ``fp`` and
``data`` that preserve the partitions of the half-edges by the values of each
class array.

.. c:function:: int _domap_stabiliser(int *orbits, const domap_t m, const domap_wreath_t w, domap_stabiliser_callback callback, void *arg, int *work)
.. c:function:: int domap_stabiliser(int *orbits, const domap_t m, const domap_wreath_t w, domap_stabiliser_callback callback, void *arg)

   Compute the group :math:`H` of the automorphisms of the map ``m``, laid
   out by ``w`` (a :c:struct:`domap_wreath_struct`), that preserve the
   partitions defined by its class arrays.

   The base is :math:`b_0, \ldots, b_{ncomp-1}`, the half-edge of canonical
   label 0 of each level (``lab[off[i]]``). Write into ``orbits``
   (``w->ncomp`` ints) the size of the orbit of :math:`b_i` right after its
   level ``i`` (1 at a skipped level), so that the order of :math:`H` is
   their product: it can exceed every C integer type (``k`` isolated edges of
   equal values have :math:`2^k k!` automorphisms), and the caller multiplies
   the sizes in its own integers.

   When ``callback`` is not ``NULL``, call it with ``arg`` on each generator,
   from the last level to the first and, at a level ``i``, in the order of
   the candidates :math:`x`, the points tried as images of :math:`b_i`
   (:ref:`omap-algorithms-stabiliser`): there are up to ``n - 1`` generators
   of ``n`` ints each, which no buffer of the caller holds. The generators
   are the identity on the half-edges of the inactive edges.

   Return 0, or the nonzero value returned by the callback, which stops the
   search and leaves ``orbits`` undefined.

   In the raw form, ``work`` is a buffer of
   ``DOMAP_STABILISER_WORK_SIZE(n, nc)`` ints.

   The subgroup search, the levels skipped and the placements cut are in
   :ref:`omap-algorithms-stabiliser`.

.. c:macro:: DOMAP_STABILISER_WORK_SIZE(n, nc)

   The size of the work buffer of :c:func:`_domap_stabiliser`, in ints, for
   at most ``n / 2`` components of at most ``n`` half-edges each,
   ``3 nc n + 5 (n / 2) + 8 n``:

   - the forward and backward maps of the values of each class array and the
     trail of the entries set (``3 nc n``);
   - the arrays of the levels used, ``target``, ``choice``, ``mark`` and
     ``dead`` (``5 (n / 2)``);
   - the buffers of the size of the largest component ``rb``, ``rbinv``,
     ``queue``, ``perm`` and ``ident`` (``5 n``);
   - the parents and the sizes of the union-find of the orbits (``2 n``);
   - and the generator being built (``n``).

   This size is computed in ``int``; the caller checks it first.

.. c:type:: domap_stabiliser_callback

   The callback of :c:func:`domap_stabiliser`, called on each generator:

   .. code-block:: c

      typedef int (*domap_stabiliser_callback)(const int *g, void *arg);

   ``g`` is the generator (``n`` ints, ``g[h]`` the image of ``h``), valid
   during the call only. ``arg`` is the argument given to
   :c:func:`domap_stabiliser`.

   A nonzero return value stops the search, which returns it; the Cython
   layer returns -1 with a Python exception set.

.. c:struct:: domap_wreath_struct

   The components of the map, grouped into kinds of isomorphic components
   (``fp`` and ``data`` alone), an input of :c:func:`domap_stabiliser`. The
   arrays it points into are owned by the caller. The levels are the
   ``ncomp`` components, kinds in order, copies in order. The isomorphism
   from a copy to another sends the half-edge of label ``i`` to the one of
   label ``i``.

   .. c:member:: const int *kind

      Per level ``j``: ``kind[j]`` its kind.

   .. c:member:: const int *off

      Per level ``j``: ``off[j]`` the offset of its half-edges in ``lab``.

   .. c:member:: const int *lab

      The half-edges of each level by canonical label: ``lab[off[j] + i]`` is
      the half-edge of label ``i`` of the level ``j``.

   .. c:member:: const int *start

      Per kind ``t``: its copies are the levels ``start[t], ...,
      start[t + 1] - 1``.

   .. c:member:: const int *coff

      Per kind ``t``: the offset ``coff[t]`` of its entries in ``cfp``,
      ``cra`` and ``crainv``.

   .. c:member:: const int *csize

      Per kind ``t``: the number ``csize[t]`` of half-edges of a copy.

   .. c:member:: const int *cfp

      Per kind ``t``: the face permutation of its first copy as a map of its
      own, ``csize[t]`` entries at the offset ``coff[t]``.

   .. c:member:: const int *cra

      Per kind ``t``: the canonical relabelling of its first copy (local
      index to label), ``csize[t]`` entries at the offset ``coff[t]``.

   .. c:member:: const int *crainv

      Per kind ``t``: the inverse of ``cra``, ``csize[t]`` entries at the
      offset ``coff[t]``.

   .. c:member:: const int *loff

      Per kind ``t``: the labels of its roots of minimal key are
      ``labels[loff[t]], ..., labels[loff[t + 1] - 1]``.

   .. c:member:: const int *labels

      The labels of the roots of minimal key of each kind, in increasing
      order.

   .. c:member:: int ncomp

      The number of components, which are the levels.

.. c:type:: domap_wreath_struct domap_wreath_t[1]

   A :c:struct:`domap_wreath_struct` as an array of one element.

.. c:macro:: OMAP_INLINE

   The storage class of the functions defined in the header,
   ``static inline``.
