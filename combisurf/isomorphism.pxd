cdef extern from "omap.h":
    enum:
        DOMAP_ISO_DISCONNECTED

    ctypedef enum domap_iso_root_rule:
        DOMAP_ISO_ROOTS_ALL
        DOMAP_ISO_ROOTS_HALF_EDGE_MAX
        DOMAP_ISO_ROOTS_HALF_EDGE_FEWEST
        DOMAP_ISO_ROOTS_EDGE_MAX
        DOMAP_ISO_ROOTS_EDGE_FEWEST

    int DOMAP_ISO_KEY_SIZE(int n, int nd, int nc)
    int DOMAP_CANONICAL_WORK_SIZE(int n, int nd, int nc)
    int DOMAP_ISO_SEQUENCE_SIZE(int n, int nd, int nc)
    int DOMAP_ISO_SEARCH_WORK_SIZE(int n, int nd, int nc)

    # the callbacks written in Cython return -1 with a Python exception set,
    # which stops the search and is returned by domap_iso_search
    ctypedef int (*domap_iso_leaf_callback)(const int *leaf, int first, void *arg) except -1

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
    int domap_iso_key(int *key, const domap_t m, const int *r, const int *rinv, int k, int *names)
    int domap_invariants_differ(const domap_t m1, const domap_t m2, int *work)
    int domap_canonical(int *key, int *relabelling, int *roots, int *nroots,
                      const domap_t m, domap_iso_root_rule rule, int *work)
    int domap_iso_search(int *sequence, int *length, int *relabelling,
                       const domap_t m, domap_iso_root_rule rule, int all_leaves,
                       domap_iso_leaf_callback callback, void *arg, int *work) except -1

    ctypedef struct domap_wreath_struct:
        int ncomp
        const int *kind
        const int *off
        const int *lab
        const int *start
        const int *coff
        const int *csize
        const int *cfp
        const int *cra
        const int *crainv
        const int *loff
        const int *labels

    ctypedef domap_wreath_struct domap_wreath_t[1]

    # as domap_iso_leaf_callback: -1 with a Python exception set stops the
    # search and is returned by domap_stabiliser
    ctypedef int (*domap_stabiliser_callback)(const int *g, void *arg) except -1

    int DOMAP_STABILISER_WORK_SIZE(int n, int nc)
    int domap_stabiliser(int *orbits, const domap_t m, const domap_wreath_t w,
                         domap_stabiliser_callback callback, void *arg, int *work) except -1
