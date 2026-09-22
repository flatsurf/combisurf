from cpython cimport array


cdef class ConjugateTree:
    cdef int alphabet_size   # size of the alphabet, 0 when it is not known
    cdef bint dense          # whether transitions are a flat alphabet-indexed table
    cdef int max_letter      # largest letter seen so far, -1 when none

    # the words, laid end to end in one buffer addressed by offsets
    cdef array.array a_wbuf
    cdef int *wbuf
    cdef int wbuf_size       # letters in use
    cdef int wbuf_capacity   # letters allocated
    cdef array.array a_wstart
    cdef int *wstart         # word -> its offset in wbuf
    cdef array.array a_wlen
    cdef int *wlen           # word -> its length
    cdef int nwords
    cdef int words_capacity

    # one entry per node in each of these
    cdef int nstates
    cdef int capacity
    cdef array.array a_dep
    cdef int *dep            # node -> length of the word it reads
    cdef array.array a_sl
    cdef int *sl             # node -> suffix link
    cdef array.array a_parent
    cdef int *parent          # node -> parent
    cdef array.array a_tword
    cdef int *tword          # node -> word of the label of the edge into it
    cdef array.array a_tstart
    cdef int *tstart         # node -> start of that label
    cdef array.array a_tend
    cdef int *tend           # node -> end of that label, -1 for a leaf

    # transitions, dense or sibling lists; only one of the two is allocated
    cdef array.array a_trans
    cdef int *trans          # node * alphabet_size + letter -> child, -1 when absent
    cdef array.array a_fchild
    cdef int *fchild         # node -> first child, -1 when none
    cdef array.array a_nsib
    cdef int *nsib           # node -> next sibling, -1 when last

    cdef int _reserve_nodes(self, int capacity) except -1
    cdef int _reserve_words(self, int num_words, int num_letters) except -1
    cdef int _add_node(self) except -1
    cdef int _add_word(self, w) except -1
    cdef void _pop_word(self) noexcept nogil
    cdef void _truncate_word(self, int i, int size) noexcept nogil

    cdef inline int _letter(self, int i, int k) noexcept nogil
    cdef inline int _child(self, int s, int letter) noexcept nogil
    cdef inline void _add_child(self, int s, int letter, int t) noexcept nogil
    cdef void _replace_child(self, int s, int letter, int old, int new) noexcept nogil

    cdef int _test_and_split(self, int s, int i, int k, int p, int letter) except -2
    cdef void _canonize(self, int *s, int i, int *k, int p) noexcept nogil
    cdef int _update(self, int *s, int i, int *k, int p) except -1
