r"""
Data structures for efficient consecutive partial sums

We consider a vector v of fixed size n on which we allow two operations
- updates of the form v[i] += x
- computation of the partial sum sum(v[i:j])

A plain array gives O(1) updates and O(n) partial sums; this direct approach
is :class:`PartialSumsNaive`. :class:`PartialSumsFenwick` does both in
O(\log(n)) time instead, at the cost of a larger constant. :func:`PartialSums`
picks whichever of the two is faster for a given size, so that a caller that
does not want to make that choice itself does not have to::

    sage: from combisurf.partial_sums import PartialSums, PartialSumsNaive, PartialSumsFenwick
    sage: type(PartialSums(3)) is PartialSumsNaive
    True
    sage: type(PartialSums(10000)) is PartialSumsFenwick
    True

Both structures answer the same queries and are checked against each other
and against a plain Python list, in ``test/test_partial_sums.py``::

    sage: P0 = PartialSumsNaive(5)
    sage: P1 = PartialSumsFenwick(5)
    sage: for i, x in enumerate([3, 1, 4, 1, 5]):
    ....:     P0.update(i, x)
    ....:     P1.update(i, x)
    sage: P0.partial_sum(1, 4) == P1.partial_sum(1, 4)
    True
"""
# ****************************************************************************
#  This file is part of combisurf
#
#       Copyright (C) 2026 Vincent Delecroix
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# ****************************************************************************

from cpython cimport array
from libc.string cimport memset


cdef class PartialSumsNaive:
    r"""
    A vector of ``n`` integers, initially zero, updated in `O(1)` and summed
    over a range in `O(n)`.

    INPUT:

    - ``n`` -- non-negative integer, the size of the vector

    EXAMPLES::

        sage: from combisurf.partial_sums import PartialSumsNaive
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
    def __cinit__(self, *args, **kwds):
        r"""
        Set up an empty buffer; ``__init__`` allocates the real one.
        """
        self.a_values = array.array('i', [])
        self.values = self.a_values.data.as_ints
        self.n = 0

    def __init__(self, n):
        r"""
        TESTS::

            sage: from combisurf.partial_sums import PartialSumsNaive
            sage: PartialSumsNaive(3)
            PartialSumsNaive([0, 0, 0])
            sage: PartialSumsNaive(0)
            PartialSumsNaive([])
        """
        self.a_values = array.array('i', [0] * n)
        self.values = self.a_values.data.as_ints
        self.n = len(self.a_values)

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.partial_sums import PartialSumsNaive
            sage: PartialSumsNaive(2)
            PartialSumsNaive([0, 0])
        """
        return f"PartialSumsNaive({list(self.a_values)})"

    def reset(self):
        r"""
        Set every entry back to zero.

        EXAMPLES::

            sage: from combisurf.partial_sums import PartialSumsNaive
            sage: P = PartialSumsNaive(3)
            sage: P.update(1, 7)
            sage: P
            PartialSumsNaive([0, 7, 0])
            sage: P.reset()
            sage: P
            PartialSumsNaive([0, 0, 0])
        """
        memset(self.values, 0, self.n * sizeof(int))

    def update(self, int i, int x):
        r"""
        Add ``x`` to the entry at position ``i``.

        INPUT:

        - ``i`` -- integer with `0 \le i < n`, the position to update

        - ``x`` -- integer, the value added to that position

        EXAMPLES::

            sage: from combisurf.partial_sums import PartialSumsNaive
            sage: P = PartialSumsNaive(3)
            sage: P.update(0, 2)
            sage: P.update(0, 5)
            sage: P
            PartialSumsNaive([7, 0, 0])
        """
        self.values[i] += x

    def partial_sum(self, int start, int end):
        r"""
        Return the partial sum from ``start`` (included) to ``end`` (excluded).

        INPUT:

        - ``start``, ``end`` -- integers with `0 \le \text{start} \le
          \text{end} \le n`, the half-open range to sum over

        EXAMPLES::

            sage: from combisurf.partial_sums import PartialSumsNaive
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
        cdef long long s = 0
        cdef int k
        for k in range(start, end):
            s += self.values[k]
        return s


cdef class PartialSumsFenwick:
    r"""
    A vector of ``n`` integers, initially zero, updated and summed over a
    range in `O(\log(n))`.

    INPUT:

    - ``n`` -- positive integer, the size of the vector

    EXAMPLES::

        sage: from combisurf.partial_sums import PartialSumsFenwick
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
    `\text{lowbit}(i)` entries at the 0-based positions `i -
    \text{lowbit}(i), \ldots, i - 1`, `\text{lowbit}(i)` being the lowest
    set bit of `i`. An ``update`` at position `i` walks the nodes `i +
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
    def __cinit__(self, *args, **kwds):
        r"""
        Set up an empty buffer; ``__init__`` allocates the real one.
        """
        self.a_values = array.array('q', [])
        self.values = self.a_values.data.as_longlongs
        self.n = 0

    def __init__(self, n):
        r"""
        TESTS::

            sage: from combisurf.partial_sums import PartialSumsFenwick
            sage: PartialSumsFenwick(3)
            PartialSumsFenwick([0, 0, 0])
        """
        if n <= 0:
            raise ValueError("n must be a positive integer")
        self.a_values = array.array('q', [0] * (n + 1))
        self.values = self.a_values.data.as_longlongs
        self.n = n

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from combisurf.partial_sums import PartialSumsFenwick
            sage: PartialSumsFenwick(2)
            PartialSumsFenwick([0, 0])
        """
        return f"PartialSumsFenwick({[self.partial_sum(i, i + 1) for i in range(self.n)]})"

    def reset(self):
        r"""
        Set every entry back to zero.

        EXAMPLES::

            sage: from combisurf.partial_sums import PartialSumsFenwick
            sage: P = PartialSumsFenwick(3)
            sage: P.update(1, 7)
            sage: P
            PartialSumsFenwick([0, 7, 0])
            sage: P.reset()
            sage: P
            PartialSumsFenwick([0, 0, 0])
        """
        memset(self.values, 0, (self.n + 1) * sizeof(long long))

    def update(self, int i, int x):
        r"""
        Add ``x`` to the entry at position ``i``.

        INPUT:

        - ``i`` -- integer with `0 \le i < n`, the position to update

        - ``x`` -- integer, the value added to that position

        EXAMPLES::

            sage: from combisurf.partial_sums import PartialSumsFenwick
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
        if not 0 <= i < self.n:
            raise IndexError("position out of range")
        fenwick_add(self.values, self.n, i, x)

    def partial_sum(self, int start, int end):
        r"""
        Return the partial sum from ``start`` (included) to ``end`` (excluded).

        INPUT:

        - ``start``, ``end`` -- integers with `0 \le \text{start} \le
          \text{end} \le n`, the half-open range to sum over

        EXAMPLES::

            sage: from combisurf.partial_sums import PartialSumsFenwick
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
        if not 0 <= start <= end <= self.n:
            raise IndexError("range out of bounds")
        if start == 0:
            return fenwick_prefix(self.values, end)
        return fenwick_prefix(self.values, end) - fenwick_prefix(self.values, start)


# Crossover of naive/Fenwick timings (interleaved) on two updates and one partial_sum per step.
cdef int _PARTIAL_SUMS_NAIVE_THRESHOLD = 256


def PartialSums(n):
    r"""
    Return a :class:`PartialSumsNaive` or :class:`PartialSumsFenwick`
    of size ``n``, whichever is faster at that size.

    INPUT:

    - ``n`` -- positive integer, the size of the vector

    EXAMPLES::

        sage: from combisurf.partial_sums import PartialSums, PartialSumsNaive, PartialSumsFenwick
        sage: type(PartialSums(256)) is PartialSumsNaive
        True
        sage: type(PartialSums(257)) is PartialSumsFenwick
        True

    Both answer the same queries, so a caller does not need to know which one
    it got::

        sage: P = PartialSums(5)
        sage: P.update(2, 3)
        sage: P.partial_sum(0, 5)
        3
    """
    if n <= _PARTIAL_SUMS_NAIVE_THRESHOLD:
        return PartialSumsNaive(n)
    return PartialSumsFenwick(n)
