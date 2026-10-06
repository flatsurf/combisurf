cdef extern from "isomorphism.h":
    enum:
        ISO_DISCONNECTED

    ctypedef enum iso_root_rule:
        ISO_ROOTS_ALL
        ISO_ROOTS_HALF_EDGE_MAX
        ISO_ROOTS_HALF_EDGE_FEWEST
        ISO_ROOTS_EDGE_MAX
        ISO_ROOTS_EDGE_FEWEST

    int ISO_KEY_SIZE(int n, int nd, int nc)
    int ISO_WORK_SIZE(int n, int nd, int nc)
    int ISO_SEQUENCE_SIZE(int n, int nd, int nc)
    int ISO_SEARCH_WORK_SIZE(int n, int nd, int nc)

    # the callbacks written in Cython return -1 with a Python exception set,
    # which stops the search and is returned by iso_search
    ctypedef int (*iso_leaf_callback)(const int *leaf, int first, void *arg) except -1

    int iso_num_edges(int n, const int *fp)
    int iso_relabelling_from(int n, const int *fp, int root, int *r, int *rinv, int *queue)
    int iso_key(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                const int *r, const int *rinv, int k, int *key, int *names)
    int iso_canonical(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                      iso_root_rule rule, int *key, int *relabelling, int *roots, int *nroots,
                      int *work)
    int iso_search(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                   iso_root_rule rule, int all_leaves, int *sequence, int *length, int *relabelling,
                   iso_leaf_callback callback, void *arg, int *work) except -1
