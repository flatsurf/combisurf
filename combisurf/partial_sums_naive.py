r"""
Data structures for efficient consecutive partial sums, reference implementation

This module holds :class:`PartialSumsNaive` and :class:`PartialSumsFenwick`,
the pure Python data structures for a vector v of fixed size n on which we
allow two operations
- updates of the form v[i] += x
- computation of the partial sum sum(v[i:j])

They are the reference against which the fast Cython classes of the same
names in :mod:`combisurf.partial_sums` are tested; everything else in the
package uses those. See :mod:`combisurf.partial_sums` for what the two
structures are and for the :func:`~combisurf.partial_sums.PartialSums`
factory that picks between them.

A plain array gives O(1) updates and O(n) partial sums; this direct approach
is :class:`PartialSumsNaive`. :class:`PartialSumsFenwick` does both in
O(\log(n)) time instead, at the cost of a larger constant.

Both structures answer the same queries and are checked against each other,
and against their Cython counterparts, in ``test/test_partial_sums.py``::

    sage: from combisurf.partial_sums_naive import PartialSumsNaive, PartialSumsFenwick
    sage: P0 = PartialSumsNaive(5)
    sage: P1 = PartialSumsFenwick(5)
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


class PartialSumsFenwick(object):
    r"""
    A vector of ``n`` integers, initially zero, updated and summed over a
    range in `O(\log(n))`.

    INPUT:

    - ``n`` -- positive integer, the size of the vector

    EXAMPLES::

        sage: from combisurf.partial_sums_naive import PartialSumsFenwick
        sage: P = PartialSumsFenwick(5)
        sage: P
        PartialSumsFenwick([0, 0, 0, 0, 0])
        sage: P.update(2, 3)
        sage: P.update(4, 1)
        sage: P
        PartialSumsFenwick([0, 0, 3, 0, 1])
        sage: P.partial_sum(0, 5)
        4
        sage: P.partial_sum(2, 4)
        3
        sage: P.reset()
        sage: P
        PartialSumsFenwick([0, 0, 0, 0, 0])

    ALGORITHM:

    A Fenwick tree, or binary indexed tree (:ref:`fenwick1994`): the vector is
    stored 1-indexed as ``tree[1..n]``, where ``tree[i]`` is the sum of the
    entries ending at position `i - 1` (0-based, excluded) whose count is the
    lowest set bit of `i`. An ``update`` at position `i` walks the nodes `i +
    1, i + 1 + \text{lowbit}, \ldots` up to `n` that cover it, and a
    ``partial_sum`` up to some position walks the nodes `i, i - \text{lowbit},
    \ldots` down to `0`; both walks have at most `\lfloor \log_2(n) \rfloor +
    1` steps.

    TESTS:

    ``n`` must be positive::

        sage: PartialSumsFenwick(0)
        Traceback (most recent call last):
        ...
        ValueError: n must be a positive integer
        sage: PartialSumsFenwick(-1)
        Traceback (most recent call last):
        ...
        ValueError: n must be a positive integer
    """
    def __init__(self, n):
        r"""
        TESTS::

            sage: from combisurf.partial_sums_naive import PartialSumsFenwick
            sage: PartialSumsFenwick(3)
            PartialSumsFenwick([0, 0, 0])
        """
        if n <= 0:
            raise ValueError("n must be a positive integer")
        self._n = n
        self._tree = [0] * (n + 1)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsFenwick
            sage: PartialSumsFenwick(2)
            PartialSumsFenwick([0, 0])
        """
        return f"PartialSumsFenwick({[self.partial_sum(i, i + 1) for i in range(self._n)]})"

    def reset(self):
        r"""
        Set every entry back to zero.

        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsFenwick
            sage: P = PartialSumsFenwick(3)
            sage: P.update(1, 7)
            sage: P
            PartialSumsFenwick([0, 7, 0])
            sage: P.reset()
            sage: P
            PartialSumsFenwick([0, 0, 0])
        """
        for i in range(len(self._tree)):
            self._tree[i] = 0

    def update(self, i, x):
        r"""
        Add ``x`` to the entry at position ``i``.

        INPUT:

        - ``i`` -- integer with `0 \le i < n`, the position to update

        - ``x`` -- integer, the value added to that position

        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsFenwick
            sage: P = PartialSumsFenwick(3)
            sage: P.update(0, 2)
            sage: P.update(0, 5)
            sage: P
            PartialSumsFenwick([7, 0, 0])

        Out-of-range positions raise::

            sage: P.update(3, 1)
            Traceback (most recent call last):
            ...
            IndexError: position out of range
            sage: P.update(-1, 1)
            Traceback (most recent call last):
            ...
            IndexError: position out of range
        """
        if not 0 <= i < self._n:
            raise IndexError("position out of range")
        i += 1
        while i <= self._n:
            self._tree[i] += x
            i += i & (-i)

    def _prefix(self, i):
        s = 0
        while i > 0:
            s += self._tree[i]
            i -= i & (-i)
        return s

    def partial_sum(self, start, end):
        r"""
        Return the partial sum from ``start`` (included) to ``end`` (excluded).

        INPUT:

        - ``start``, ``end`` -- integers with `0 \le \text{start} \le
          \text{end} \le n`, the half-open range to sum over

        EXAMPLES::

            sage: from combisurf.partial_sums_naive import PartialSumsFenwick
            sage: P = PartialSumsFenwick(5)
            sage: for i, x in enumerate([1, 2, 3, 4, 5]):
            ....:     P.update(i, x)
            sage: P.partial_sum(0, 5)
            15
            sage: P.partial_sum(1, 3)
            5

        The range is half-open, so an empty range sums to zero::

            sage: P.partial_sum(2, 2)
            0

        Out-of-range positions raise::

            sage: P.partial_sum(0, 6)
            Traceback (most recent call last):
            ...
            IndexError: range out of bounds
            sage: P.partial_sum(-1, 5)
            Traceback (most recent call last):
            ...
            IndexError: range out of bounds
        """
        if not 0 <= start <= end <= self._n:
            raise IndexError("range out of bounds")
        if start == 0:
            return self._prefix(end)
        return self._prefix(end) - self._prefix(start)
