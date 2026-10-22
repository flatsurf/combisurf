/*
 * Partial permutations of {0, ..., n - 1}, with -1 on the inactive points.
 * Documented in doc/c/perm.rst.
 *
 *       Copyright (C) 2026 Vincent Delecroix
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 */
#ifndef COMBISURF_PERM_H
#define COMBISURF_PERM_H

#define PERM_INLINE static inline

/* Types *********************************************************************/

/* A partial permutation: the image p[h] of each of the n points. */
typedef struct
{
    int *p;
    int n;
}
perm_struct;

typedef perm_struct perm_t[1];

/* Checks ********************************************************************/

/* 1 for a permutation of its active points, 0 if one is sent to an inactive
   point, -1 if not injective or out of range; seen has n ints. */
int _perm_check(const int *p, int n, int *seen);

PERM_INLINE int
perm_check(const perm_t p, int *seen)
{
    return _perm_check(p->p, p->n, seen);
}

/* Cycles ********************************************************************/

/* Set res to the identity. */
void _perm_id(int *res, int n);

PERM_INLINE void
perm_id(perm_t res)
{
    _perm_id(res->p, res->n);
}

/* Number the cycles; return their number, or -1 if a walk leaves [0, n). */
int _perm_dense_cycles(int *id, int *size, const int *p, int n);

PERM_INLINE int
perm_dense_cycles(int *id, int *size, const perm_t p)
{
    return _perm_dense_cycles(id, size, p->p, p->n);
}

/* Orbits ********************************************************************/

/* The representative of the orbit of h, with path halving. */
PERM_INLINE int
_perm_orbits_find(int *parent, int h)
{
    while (parent[h] != h)
    {
        parent[h] = parent[parent[h]];
        h = parent[h];
    }
    return h;
}

/* Merge the orbits of a and b, the larger representative under the smaller. */
PERM_INLINE void
_perm_orbits_union(int *parent, int a, int b)
{
    a = _perm_orbits_find(parent, a);
    b = _perm_orbits_find(parent, b);
    if (a < b)
        parent[b] = a;
    else if (b < a)
        parent[a] = b;
}

#endif
