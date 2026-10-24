``misc.h``: helpers in C
========================

Files: ``combisurf/src/misc.h``

Small helpers shared by the C code, all defined in the header.

Memory
------

.. c:function:: void *combisurf_malloc(size_t size)

   Allocate ``size`` bytes, at least one, since ``malloc(0)`` may return
   ``NULL`` and ``size = 0`` is a valid size (a buffer for an empty map).
   Call ``abort()`` when the allocation fails, as ``flint_malloc`` does, so
   that the result is never ``NULL``. The memory is released with ``free``.

Union-find
----------

A partition of the points :math:`0, \ldots, n - 1` is kept in a union-find
forest, given by the array ``parent`` of the points: the root of a tree is
its own parent, and it is the representative of its part. Merging
always puts the larger root under the smaller one, so that the representative
of a part is its smallest point. ``parent`` starts as the identity (for
instance by :c:func:`_perm_id`), every point in a part of its own.

The caller finds the roots with the lookup of its choice, then links them:
the two lookups shorten the paths differently, and the faster one depends on
the use.

.. c:function:: int uf_find_halving(int *parent, int h)

   Return the root of the tree of ``h``, with path halving: every point met
   on the way gets its grandparent as its parent.

.. c:function:: int uf_find_compression(int *parent, int h)

   Return the root of the tree of ``h``, with path compression: a first pass
   finds the root, a second pass sets the parent of every point of the path
   from ``h`` to the root to the root.

.. c:function:: void uf_link(int *parent, int a, int b)

   Merge the trees of the roots ``a`` and ``b``: the larger one gets the
   smaller one as its parent, so that the root of the merged tree is again
   its smallest point. Nothing changes when ``a == b``.

.. c:macro:: MISC_INLINE

   The storage class of the functions defined in the header,
   ``static inline``.
