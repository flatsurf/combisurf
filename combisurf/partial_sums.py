r"""
Data structure for efficient consecutive partial sums

We have a vector v of size 2^n on which we consider two kinds of operations
- updates of the form v[i] += x
- computation of the partial sums sum(v[:i])

Doing it the naive way, the first operation costs O(1) arithmetic step and the
second O(2^n) arithmetic steps. We implement a data structure that do both in
O(n) time.
"""

from sage.rings.integer_ring import ZZ


class PartialSumsNaive:
    r"""
    Update in O(1) and partial sum in O(n)
    """
    def __init__(self, n):
        self._values = [0] * n

    def reset(self):
        r"""
        Reset the values to zero.
        """
        for i in range(len(self._values)):
            self._values[i] = 0

    def update(self, i, x):
        r"""
        Add ``x`` to the entry at position ``i``.
        """
        self._values[i] += x

    def partial_sum(self, start, end):
        r"""
        Return the partial sum from ``start`` (included) to ``end`` (excluded).
        """
        return sum(self._values[start: end])


class PartialSumsBinarySplitting:
    r"""
    Update and partial sums in O(log(n))

    Given n=2^b, we store the 2^{b+1} - 1 partial sums [a2^l, (a+1)2^l) for l in {0,1,\ldots,b-1}.
    We store these partial sums on a plain list and the tree structure is implicit.
    The leaf [a,a+1) is stored at position 2^b a. The parent of i is i//2 (or i >> 1).
    """
    def __init__(self, n):
        if n <= 0:
            raise ValueError("n must be a positive integer")
        self._b = ZZ(n - 1).nbits()
        self._values = [0] * (2 ** (self._b + 1))

    def reset(self):
        r"""
        Reset the values to zero.
        """
        for i in range(len(self._values)):
            self._values[i] = 0

    def update(self, i, x):
        r"""
        Add ``x`` to the entry at position ``i``.
        """
        ii = 1 << self._b
        assert 0 <= i < ii
        ii += i
        while ii:
            self._values[ii] += x
            ii >>= 1

    def index_to_interval(self, i):
        l = 0
        j = ZZ(i)
        while j:
            j >>= 1
            l += 1
        a = i - 2**(l-1)  # remove the highest bit weight
        return (a, self._b - l + 2)

    def interval_to_index(self, a, l):
        return a + 2**(self._b - l)

    def partial_sum(self, start, end):
        r"""
        Return the partial sum from ``start`` (included) to ``end`` (excluded).
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


PartialSums = PartialSumsBinarySplitting
