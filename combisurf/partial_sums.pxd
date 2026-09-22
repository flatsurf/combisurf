from cpython cimport array


cdef class PartialSumsNaive:
    cdef array.array a_values  # the vector, handed out to Python as a list
    cdef int *values           # C view on it
    cdef int n                 # its length


cdef class PartialSumsBinarySplitting:
    cdef array.array a_values  # the segment tree, one flat array
    cdef int *values           # C view on it
    cdef int b                 # padded size is 2^b, at least n
