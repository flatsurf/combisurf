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

    int DOMAP_INVARIANTS_WORK_SIZE(int n)

    ctypedef struct omap_struct:
        int *fp
        int n

    ctypedef omap_struct omap_t[1]

    ctypedef struct domap_struct:
        omap_struct map
        int nd
        const int *data
        int nc
        const int *classes

    ctypedef domap_struct domap_t[1]

    int _omap_num_edges(const int *fp, int n)
    int omap_num_edges(const omap_t m)
    int _omap_relabelling_from(int *r, int *rinv, const int *fp, int n, int root, int *queue)
    int omap_relabelling_from(int *r, int *rinv, const omap_t m, int root, int *queue)
    int domap_key(int *key, const domap_t m, const int *r, const int *rinv, int k, int *names)
    int domap_invariants_differ(const domap_t m1, const domap_t m2, int *work)
    int domap_canonical(int *key, int *relabelling, int *roots, int *nroots,
                      const domap_t m, domap_root_rule rule, int *work)
    int domap_search(int *sequence, int *length, int *relabelling,
                   const domap_t m, domap_root_rule rule, int all_leaves,
                   domap_leaf_callback callback, void *arg, int *work) except -1
