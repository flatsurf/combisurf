``perm.h``: partial permutations in C
=====================================

Files: ``combisurf/src/perm.{h,c}``

*Partial permutations* are bijections on subsets of :math:`\{0, \ldots, n -
1\}`, given as arrays of ``n`` ints with ``-1`` on the inactive elements of
the domain.

.. c:struct:: perm_struct

   A partial permutation: the image ``p[h]`` of each point ``h`` of
   :math:`\{0, \ldots, n - 1\}`, ``-1`` when ``h`` is inactive.

   .. c:member:: int *p

      The images, an array of ``n`` ints.

   .. c:member:: int n

      The number of points.

.. c:type:: perm_struct perm_t[1]

   A :c:struct:`perm_struct` as an array of one element, so that a
   ``perm_t`` declared on the stack is passed by pointer.

The functions on a whole permutation come in pairs: the raw form
``_perm_func`` takes the array and its length, and the public form
``perm_func`` takes a :c:type:`perm_t`. A raw form that needs a work buffer
takes it from the caller, sized by its ``*_WORK_SIZE`` macro, and allocates
nothing; its public form allocates the buffer through
:c:func:`combisurf_malloc` and frees it before returning. The output arrays
are always handed in by the caller.

Constructors
------------

.. c:function:: void _perm_id(int *res, int n)
.. c:function:: void perm_id(perm_t res)

   Set ``res`` (``n`` ints) to the identity: ``res[i] = i`` for every ``i``.

Checks
------

.. c:function:: int _perm_check(const int *p, int n, int *seen)
.. c:function:: int perm_check(const perm_t p)

   Check that ``p`` is a permutation of its active points.

   Return ``1`` when ``p`` is a permutation of its active points, ``0`` when
   it is injective and in range but sends an active point to an inactive
   one, and ``-1`` when an image is out of :math:`\{-1, 0, \ldots, n - 1\}`
   or two active points have the same image.

   In the raw form, ``seen`` is a work buffer of ``PERM_CHECK_WORK_SIZE(n)``
   ints; its content on return is undefined.

.. c:macro:: PERM_CHECK_WORK_SIZE(n)

   The size of the work buffer of :c:func:`_perm_check`, in ints, ``n``: a
   mark per point, set when the point is met as an image.

Cycles
------

.. c:function:: int _perm_dense_cycles(int *id, int *size, const int *p, int n)
.. c:function:: int perm_dense_cycles(int *id, int *size, const perm_t p)

   Number the cycles of ``p`` ``0, 1, 2, ...`` in the increasing order of their
   smallest points and record their lengths: ``id[i] = k`` and ``size[i] = s``
   when the point ``i`` lies in the cycle of number ``k`` and this cycle has
   length ``s``.

   Write into ``id`` (``n`` ints) the number of the cycle of each active
   point, and into ``size`` (``n`` ints) the length of that cycle; both are
   ``-1`` on the inactive points. ``size`` may be ``NULL``, and is then not
   written.

   Return the number of cycles, or ``-1`` when the walk along a cycle leaves
   :math:`\{0, \ldots, n - 1\}` (an image ``-1`` or out of range on an active
   point); ``id`` and ``size`` are then partly written. A walk stops at a
   point already numbered, so that a non-injective ``p`` does not loop.

.. c:macro:: PERM_INLINE

   The storage class of the functions defined in the header,
   ``static inline``.
