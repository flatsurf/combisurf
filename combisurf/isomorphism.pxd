cdef extern from "omap.h":
    enum:
        DOMAP_DISCONNECTED

    ctypedef enum domap_root_rule:
        DOMAP_ROOTS_ALL
        DOMAP_ROOTS_HALF_EDGE_MAX
        DOMAP_ROOTS_HALF_EDGE_FEWEST
        DOMAP_ROOTS_EDGE_MAX
        DOMAP_ROOTS_EDGE_FEWEST

    int DOMAP_KEY_SIZE(int n, int nd, int nc)
    int DOMAP_WORK_SIZE(int n, int nd, int nc)
    int DOMAP_SEQUENCE_SIZE(int n, int nd, int nc)
    int DOMAP_SEARCH_WORK_SIZE(int n, int nd, int nc)

    # the callbacks written in Cython return -1 with a Python exception set,
    # which stops the search and is returned by domap_search
    ctypedef int (*domap_leaf_callback)(const int *leaf, int first, void *arg) except -1

    int omap_num_edges(int n, const int *fp)
    int omap_relabelling_from(int n, const int *fp, int root, int *r, int *rinv, int *queue)
    int domap_key(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                const int *r, const int *rinv, int k, int *key, int *names)
    int domap_canonical(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                      domap_root_rule rule, int *key, int *relabelling, int *roots, int *nroots,
                      int *work)
    int domap_search(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                   domap_root_rule rule, int all_leaves, int *sequence, int *length, int *relabelling,
                   domap_leaf_callback callback, void *arg, int *work) except -1
