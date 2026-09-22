r"""
Data structures for efficient consecutive partial sums, reference implementation

This module holds :class:`PartialSumsNaive` and
:class:`PartialSumsBinarySplitting`, the pure Python data structures for a
vector v of fixed size n on which we allow two operations
- updates of the form v[i] += x
- computation of the partial sum sum(v[i:j])

They are the reference against which the fast Cython classes of the same
names in :mod:`combisurf.partial_sums` are tested; everything else in the
package uses those. See :mod:`combisurf.partial_sums` for what the two
structures are and for the :func:`~combisurf.partial_sums.PartialSums`
factory that picks between them.

A plain array gives O(1) updates and O(n) partial sums; this direct approach
is :class:`PartialSumsNaive`. :class:`PartialSumsBinarySplitting` does both in
O(\log(n)) time instead, at the cost of a larger constant.

Both structures answer the same queries and are checked against each other,
and against their Cython counterparts, in ``test/test_partial_sums.py``::

    sage: from combisurf.partial_sums_naive import PartialSumsNaive, PartialSumsBinarySplitting
    sage: P0 = PartialSumsNaive(5)
    sage: P1 = PartialSumsBinarySplitting(5)
    sage: for i, x in enumerate([3, 1, 4, 1, 5]):
    ....:     P0.update(i, x)
    ....:     P1.update(i, x)
    sage: P0.partial_sum(1, 4) == P1.partial_sum(1, 4)
    True
"""


class PartialSumsNaive(object):
    r"""
    A vector of ``n`` integers, initially zero, updated in `O(1)` and summed
    over a range in `O(n)`.

    INPUT:

    - ``n`` -- non-negative integer, the size of the vector

    EXAMPLES::

        sage: from combisurf.partial_sums_naive import PartialSumsNaive
        sage: P = PartialSumsNaive(5)
        sage: P
        PartialSumsNaive([0, 0, 0, 0, 0])
        sage: P.update(2, 3)
        sage: P.update(4, 1)
        sage: P
        PartialSumsNaive([0, 0, 3, 0, 1])
        sage: P.partial_sum(0, 5)
        4
        sage: P.partial_sum(2, 4)
        3
        sage: P.reset()
        sage: P
        PartialSumsNaive([0, 0, 0, 0, 0])
    """
    def __init__(self, n):
        r"""
        TESTS::

            sage: from combisurf.partial_sums_naive import PartialSumsNaive
            sage: PartialSumsNaive(3)
            PartialSumsNaive([0, 0, 0])
            sage: PartialSumsNaive(0)
            PartialSumsNaive([])
        """
        self._values = [0] * n

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsNaive
            sage: PartialSumsNaive(2)
            PartialSumsNaive([0, 0])
        """
        return f"PartialSumsNaive({self._values})"

    def reset(self):
        r"""
        Set every entry back to zero.

        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsNaive
            sage: P = PartialSumsNaive(3)
            sage: P.update(1, 7)
            sage: P
            PartialSumsNaive([0, 7, 0])
            sage: P.reset()
            sage: P
            PartialSumsNaive([0, 0, 0])
        """
        for i in range(len(self._values)):
            self._values[i] = 0

    def update(self, i, x):
        r"""
        Add ``x`` to the entry at position ``i``.

        INPUT:

        - ``i`` -- integer with `0 \le i < n`, the position to update

        - ``x`` -- integer, the value added to that position

        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsNaive
            sage: P = PartialSumsNaive(3)
            sage: P.update(0, 2)
            sage: P.update(0, 5)
            sage: P
            PartialSumsNaive([7, 0, 0])
        """
        self._values[i] += x

    def partial_sum(self, start, end):
        r"""
        Return the partial sum from ``start`` (included) to ``end`` (excluded).

        INPUT:

        - ``start``, ``end`` -- integers with `0 \le \text{start} \le
          \text{end} \le n`, the half-open range to sum over

        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsNaive
            sage: P = PartialSumsNaive(5)
            sage: for i, x in enumerate([1, 2, 3, 4, 5]):
            ....:     P.update(i, x)
            sage: P.partial_sum(0, 5)
            15
            sage: P.partial_sum(1, 3)
            5

        The range is half-open, so an empty range sums to zero::

            sage: P.partial_sum(2, 2)
            0
        """
        return sum(self._values[start: end])


class PartialSumsBinarySplitting(object):
    r"""
    A vector of ``n`` integers, initially zero, updated and summed over a
    range in `O(\log(n))`.

    INPUT:

    - ``n`` -- positive integer, the size of the vector

    EXAMPLES::

        sage: from combisurf.partial_sums_naive import PartialSumsBinarySplitting
        sage: P = PartialSumsBinarySplitting(5)
        sage: P
        PartialSumsBinarySplitting([0, 0, 0, 0, 0, 0, 0, 0])
        sage: P.update(2, 3)
        sage: P.update(4, 1)
        sage: P
        PartialSumsBinarySplitting([0, 0, 3, 0, 1, 0, 0, 0])
        sage: P.partial_sum(0, 5)
        4
        sage: P.partial_sum(2, 4)
        3
        sage: P.reset()
        sage: P
        PartialSumsBinarySplitting([0, 0, 0, 0, 0, 0, 0, 0])

    ALGORITHM:

    Given `n`, let `b` be the number of bits of `n - 1`, so `m = 2^b \ge n`.
    The vector is padded to size `m` and stored as a complete binary tree with
    `m` leaves: the partial sum over every dyadic interval `[a 2^l, (a+1)
    2^l)` with `0 \le l \le b` is kept up to date. There are `2m - 1` such
    intervals in total, stored in one flat list of that length with the usual
    implicit heap layout — the root is at index `1`, the leaf for position `i`
    is at index `m + i`, and the parent of index `j` is at index `j // 2`. An
    ``update`` touches the `b + 1` dyadic intervals containing the position,
    and a ``partial_sum`` decomposes its range into at most `2b` of them, both
    root-to-leaf walks.

    TESTS:

    ``n`` must be positive, since a size-``0`` tree has no root to walk::

        sage: PartialSumsBinarySplitting(0)
        Traceback (most recent call last):
        ...
        ValueError: n must be a positive integer
        sage: PartialSumsBinarySplitting(-1)
        Traceback (most recent call last):
        ...
        ValueError: n must be a positive integer
    """
    def __init__(self, n):
        r"""
        TESTS::

            sage: from combisurf.partial_sums_naive import PartialSumsBinarySplitting
            sage: PartialSumsBinarySplitting(3)
            PartialSumsBinarySplitting([0, 0, 0, 0])
        """
        if n <= 0:
            raise ValueError("n must be a positive integer")
        self._b = (n - 1).bit_length()
        self._values = [0] * (2 ** (self._b + 1))

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsBinarySplitting
            sage: PartialSumsBinarySplitting(2)
            PartialSumsBinarySplitting([0, 0])
        """
        m = 1 << self._b
        return f"PartialSumsBinarySplitting({[self._values[m + i] for i in range(m)]})"

    def reset(self):
        r"""
        Set every entry back to zero.

        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsBinarySplitting
            sage: P = PartialSumsBinarySplitting(3)
            sage: P.update(1, 7)
            sage: P
            PartialSumsBinarySplitting([0, 7, 0, 0])
            sage: P.reset()
            sage: P
            PartialSumsBinarySplitting([0, 0, 0, 0])
        """
        for i in range(len(self._values)):
            self._values[i] = 0

    def update(self, i, x):
        r"""
        Add ``x`` to the entry at position ``i``.

        INPUT:

        - ``i`` -- integer with `0 \le i < 2^b`, the position to update, where
          `2^b` is the padded size reported by :meth:`index_to_interval`

        - ``x`` -- integer, the value added to that position

        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsBinarySplitting
            sage: P = PartialSumsBinarySplitting(3)
            sage: P.update(0, 2)
            sage: P.update(0, 5)
            sage: P
            PartialSumsBinarySplitting([7, 0, 0, 0])
        """
        ii = 1 << self._b
        assert 0 <= i < ii
        ii += i
        while ii:
            self._values[ii] += x
            ii >>= 1

    def index_to_interval(self, i):
        r"""
        Return the dyadic interval represented by node ``i`` of the tree.

        INPUT:

        - ``i`` -- integer with `1 \le i < 2^{b+1}`, an index into the
          underlying array, where `b` is as in :meth:`__init__`

        OUTPUT: a pair ``(a, l)`` such that node ``i`` represents the interval
        `[a \cdot 2^l, (a+1) \cdot 2^l)`. This is the inverse of
        :meth:`interval_to_index`.

        EXAMPLES:

        Node ``1`` is the root and represents the whole padded range::

            sage: from combisurf.partial_sums_naive import PartialSumsBinarySplitting
            sage: P = PartialSumsBinarySplitting(5)
            sage: P._b
            3
            sage: P.index_to_interval(1)
            (0, 3)

        The two children of the root split that range in half, and the
        leaves, at the far end of the array, each cover a single position::

            sage: P.index_to_interval(2)
            (0, 2)
            sage: P.index_to_interval(3)
            (1, 2)
            sage: P.index_to_interval(8)
            (0, 0)
            sage: P.index_to_interval(9)
            (1, 0)

        It is inverse to :meth:`interval_to_index`::

            sage: all(P.interval_to_index(*P.index_to_interval(i)) == i
            ....:      for i in range(1, 2 * 2**P._b))
            True
        """
        l = 0
        j = i
        while j:
            j >>= 1
            l += 1
        a = i - 2**(l - 1)  # remove the highest bit weight
        return (a, self._b - l + 1)

    def interval_to_index(self, a, l):
        r"""
        Return the index of the node representing the dyadic interval
        `[a \cdot 2^l, (a+1) \cdot 2^l)`.

        INPUT:

        - ``a``, ``l`` -- integers with `0 \le l \le b` and `0 \le a < 2^{b -
          l}`, describing the interval `[a \cdot 2^l, (a+1) \cdot 2^l)`, where
          `b` is as in :meth:`__init__`

        OUTPUT: the index into the underlying array of the node representing
        that interval. This is the inverse of :meth:`index_to_interval`.

        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsBinarySplitting
            sage: P = PartialSumsBinarySplitting(5)
            sage: P.interval_to_index(0, 3)
            1
            sage: P.interval_to_index(0, 0)
            8
            sage: P.interval_to_index(1, 0)
            9

        It is inverse to :meth:`index_to_interval`::

            sage: all(P.index_to_interval(P.interval_to_index(a, l)) == (a, l)
            ....:      for l in range(P._b + 1) for a in range(2**(P._b - l)))
            True
        """
        return a + 2**(self._b - l)

    def partial_sum(self, start, end):
        r"""
        Return the partial sum from ``start`` (included) to ``end`` (excluded).

        INPUT:

        - ``start``, ``end`` -- integers with `0 \le \text{start} \le
          \text{end} \le 2^b`, the half-open range to sum over, where `2^b`
          is the padded size reported by :meth:`index_to_interval`

        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsBinarySplitting
            sage: P = PartialSumsBinarySplitting(5)
            sage: for i, x in enumerate([1, 2, 3, 4, 5]):
            ....:     P.update(i, x)
            sage: P.partial_sum(0, 5)
            15
            sage: P.partial_sum(1, 3)
            5

        The range is half-open, so an empty range sums to zero::

            sage: P.partial_sum(2, 2)
            0
        """
        m = 1 << self._b
        assert 0 <= start < m
        assert 0 <= end <= m
        s = 0
        while start < end:
            if start % 2:
                s += self._values[m + start]
                start += 1
            if end % 2:
                end -= 1
                s += self._values[m + end]
            start >>= 1
            end >>= 1
            m >>= 1

        return s
