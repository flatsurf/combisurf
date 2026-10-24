/*
 * Small inline helpers of the C code: memory, union-find.
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

#include <stddef.h>
#include <stdlib.h>

#define MISC_INLINE static inline

/* Memory ********************************************************************/

/* Allocate size bytes, at least one, and abort() if the allocation fails. */
MISC_INLINE void *
combisurf_malloc(size_t size)
{
    void *ptr = malloc(size ? size : 1);

    if (ptr == NULL)
        abort();
    return ptr;
}

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

/* The root of the tree of h, every element of its path set to the root. */
MISC_INLINE int
uf_find_compression(int *parent, int h)
{
    int root = h;
    int next;

    while (parent[root] != root)
        root = parent[root];
    while (parent[h] != root && h != root)
    {
        next = parent[h];
        parent[h] = root;
        h = next;
    }
    return root;
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
