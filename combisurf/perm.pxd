cdef extern from "perm.h":
    int _perm_check(const int *p, int n, int *seen)
    void _perm_id(int *res, int n)
    int _perm_dense_cycles(int *id, int *size, const int *p, int n)
    int _perm_orbits_find(int *parent, int h)
    void _perm_orbits_union(int *parent, int a, int b)
