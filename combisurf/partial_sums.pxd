from cpython cimport array


cdef inline void fenwick_add(long long *tree, Py_ssize_t size, Py_ssize_t i, long long x) noexcept:
    # add x at position i (0-based, 0 <= i < size; a negative i loops forever)
    # of the Fenwick tree tree[1..size]
    i += 1
    while i <= size:
        tree[i] += x
        i += i & (-i)


cdef inline long long fenwick_prefix(long long *tree, Py_ssize_t i) noexcept:
    # sum of the positions 0, ..., i - 1 of the Fenwick tree
    cdef long long s = 0
    while i > 0:
        s += tree[i]
        i -= i & (-i)
    return s


cdef inline void fenwick_clear(long long *tree, Py_ssize_t size, Py_ssize_t i) noexcept:
    # zero the cells of the Fenwick tree that fenwick_add(tree, size, i, x)
    # touches, so that a tree can be restored to zero in the time it took to
    # fill it rather than in O(size)
    i += 1
    while i <= size:
        tree[i] = 0
        i += i & (-i)


cdef class PartialSumsNaive:
    cdef array.array a_values  # the vector, handed out to Python as a list
    cdef int *values           # C view on it
    cdef int n                 # its length


cdef class PartialSumsFenwick:
    cdef array.array a_values  # the Fenwick tree, one flat array of length n + 1
    cdef long long *values     # C view on it
    cdef int n                 # the size of the vector
