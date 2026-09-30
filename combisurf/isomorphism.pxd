cdef extern from "isomorphism.h":
    enum:
        ISO_DISCONNECTED

    int ISO_KEY_SIZE(int n, int nd, int nc)
    int ISO_WORK_SIZE(int n, int nd, int nc)

    int iso_num_edges(int n, const int *fp)
    int iso_relabelling_from(int n, const int *fp, int root, int *r, int *rinv, int *queue)
    int iso_key(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                const int *r, const int *rinv, int k, int *key, int *names)
    int iso_canonical(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                      int *key, int *relabelling, int *roots, int *nroots, int *work)
