/*
 * Partial permutations of {0, ..., n - 1}. See doc/c/perm.rst for the
 * conventions and the sizes of the buffers.
 *
 * This file is part of combisurf.
 *
 *       Copyright (C) 2026 Vincent Delecroix
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 */
#include "perm.h"
#include "misc.h"

#include <stddef.h>

int
_perm_check(const int *p, int n, int *seen)
{
    int h, v, res = 1;

    for (h = 0; h < n; h++)
        seen[h] = 0;
    for (h = 0; h < n; h++)
    {
        v = p[h];
        if (v == -1)
            continue;
        if (v < 0 || v >= n || seen[v])
            return -1;
        seen[v] = 1;
    }
    /* the image is injective and in range, so p[h] is a valid index */
    for (h = 0; h < n; h++)
        if (p[h] != -1 && p[p[h]] == -1)
            res = 0;
    return res;
}

int
perm_check(const perm_t p)
{
    int *seen = combisurf_malloc((size_t) PERM_CHECK_WORK_SIZE(p->n)
                                 * sizeof(int));
    int res = _perm_check(p->p, p->n, seen);

    free(seen);
    return res;
}

void
_perm_id(int *res, int n)
{
    int i;

    for (i = 0; i < n; i++)
        res[i] = i;
}

int
_perm_dense_cycles(int *id, int *size, const int *p, int n)
{
    int h, i, x, len, c = 0;

    for (h = 0; h < n; h++)
        id[h] = -1;
    if (size != NULL)
        for (h = 0; h < n; h++)
            size[h] = -1;
    for (h = 0; h < n; h++)
    {
        if (p[h] == -1 || id[h] != -1)
            continue;
        len = 0;
        x = h;
        /* stop at a numbered point, so that a non-injective input ends too */
        while (id[x] == -1)
        {
            id[x] = c;
            len++;
            x = p[x];
            if (x < 0 || x >= n)
                return -1;
        }
        if (size != NULL)
        {
            /* len steps exactly: the walk need not return to h */
            x = h;
            for (i = 0; i < len; i++)
            {
                size[x] = len;
                x = p[x];
            }
        }
        c++;
    }
    return c;
}
