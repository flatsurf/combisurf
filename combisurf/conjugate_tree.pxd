from cpython cimport array
from libc.stdint cimport int64_t


cdef extern from "conjugate_tree.h":
    enum:
        CT_OK
        CT_ENOMEM
        CT_EEMPTY
        CT_ENEGATIVE
        CT_EALPHABET
        CT_ETOOLARGE
        CT_EINVALID
        CT_EINTERNAL
        CT_EBROKEN

    enum:
        CT_LAYOUT_DEFAULT
        CT_LAYOUT_SPARSE
        CT_LAYOUT_DENSE
        CT_LAYOUT_ROWS

    ctypedef struct ct_tree:
        int alphabet_size
        int layout
        int max_letter
        int broken
        int *wbuf
        int wbuf_size
        int wbuf_capacity
        int *wstart
        int *wlen
        int nwords
        int words_capacity
        int nstates
        int capacity
        int *dep
        int *sl
        int *parent
        int *tword
        int *tstart
        int *tend
        int *trans
        int *fchild
        int *nsib
        int *flet
        int promote
        int *row
        int *rows
        int nrows
        int rows_capacity

    int ct_init(ct_tree *T, int alphabet, int reserve, int layout) nogil
    void ct_free(ct_tree *T) nogil
    int ct_process(ct_tree *T, const int *w, int len, int *result) nogil
    int ct_reserve(ct_tree *T, int words, int letters) nogil
    int ct_sorted_leaves(const ct_tree *T, const int *order, const int *pivot, int n, int *out, int *num) nogil
    int ct_check(const ct_tree *T) nogil
    int ct_letter(const ct_tree *T, int i, int k) nogil
    int ct_child(const ct_tree *T, int s, int letter) nogil
    int ct_canonize(const ct_tree *T, int *s, int i, int *k, int p) nogil
    int ct_leaf_as_conjugate(const ct_tree *T, int s, int *i, int *k) nogil
    int64_t ct_size(const ct_tree *T) nogil
    const char *ct_strerror(int code) nogil


cdef class ConjugateTree:
    cdef ct_tree T

    cdef int _raise(self, int err, array.array w) except -1
    cdef int _reserve(self, int words, int letters) except -1
    cdef int _process(self, const int *w, int length, int *result) except -1
    cdef array.array _order_array(self, order)
    cdef array.array _pivot_array(self, pivot, int n)
    cdef int _sorted_leaves(self, array.array order, array.array pivot, int *out) except -1
