/*
 * Small inline helpers of the C code: union-find.
 *
 *       Copyright (C) 2026 Vincent Delecroix
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 */
#ifndef COMBISURF_MISC_H
#define COMBISURF_MISC_H

#define MISC_INLINE static inline

/* Union-find ****************************************************************/

/* The root of the tree of h, with path halving. */
MISC_INLINE int
uf_find_halving(int *parent, int h)
{
    while (parent[h] != h)
    {
        parent[h] = parent[parent[h]];
        h = parent[h];
    }
    return h;
}

/* Merge the trees of the roots a and b, the larger under the smaller. */
MISC_INLINE void
uf_link(int *parent, int a, int b)
{
    if (a < b)
        parent[b] = a;
    else if (b < a)
        parent[a] = b;
}

#endif
