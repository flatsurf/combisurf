cdef extern from "perm.h":
    int _perm_check(const int *p, int n, int *seen)
    void _perm_id(int *res, int n)
    int _perm_dense_cycles(int *id, int *size, const int *p, int n)
