/*
 * Canonical keys, canonical labellings and automorphisms of connected maps.
 * See isomorphism.h for the conventions and the sizes of the buffers.
 *
 * Every active half-edge is tried as a root: the relabelling from it and
 * its key are computed in full and compared with the best key so far, which
 * costs O(n) per root and O(n^2) in all for n half-edges.
 *
 * With ISO_DEBUG defined, the relabellings and keys are checked as they are
 * computed, and a failed check aborts with a message.
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
#include "isomorphism.h"

#include <string.h>

#ifdef ISO_DEBUG
#include <stdio.h>
#include <stdlib.h>

#define ISO_CHECK(cond) do { \
        if (!(cond)) { \
            fprintf(stderr, "isomorphism.c:%d: check failed: %s\n", __LINE__, #cond); \
            abort(); \
        } \
    } while (0)

/* Check that (r, rinv) is a relabelling of the reached edges by k labels. */
static void debug_relabelling(int n, const int *fp, const int *r, const int *rinv, int k)
{
    int h, i;
    ISO_CHECK(k % 2 == 0 && 0 <= k && k <= n);
    for (h = 0; h < n; h++) {
        ISO_CHECK(-1 <= r[h] && r[h] < k);
        ISO_CHECK((r[h] == -1) == (r[h ^ 1] == -1));
        if (r[h] != -1)
            ISO_CHECK(r[h ^ 1] == (r[h] ^ 1));
        if (r[h] != -1 && fp[h] != -1)
            ISO_CHECK(r[fp[h]] != -1);
    }
    for (i = 0; i < k; i++) {
        if (rinv[i] == -1) {
            ISO_CHECK(i % 2 == 1 && rinv[i ^ 1] != -1 && fp[rinv[i ^ 1] ^ 1] == -1);
        } else {
            ISO_CHECK(0 <= rinv[i] && rinv[i] < n && fp[rinv[i]] != -1 && r[rinv[i]] == i);
        }
    }
}

/* Check that names is all -1. */
static void debug_names(int n, const int *names)
{
    int v;
    for (v = 0; v < n; v++)
        ISO_CHECK(names[v] == -1);
}
#endif

int iso_num_edges(int n, const int *fp)
{
    int e, ne = 0;
    for (e = 0; 2 * e < n; e++)
        if (fp[2 * e] != -1 || fp[2 * e + 1] != -1)
            ne++;
    return ne;
}

int iso_relabelling_from(int n, const int *fp, int root, int *r, int *rinv, int *queue)
{
    int head = 0, tail = 0, k = 0;
    int h, x;

    for (h = 0; h < n; h++)
        r[h] = -1;
    queue[tail++] = root;
    while (head < tail) {
        h = queue[head++];
        x = h;
        do {
            if (r[x] == -1) {
                r[x] = k;
                r[x ^ 1] = k + 1;
                rinv[k] = x;
                if (fp[x ^ 1] != -1) {
                    rinv[k + 1] = x ^ 1;
                    /* one push per edge, when it gets its labels, and the
                     * root: at most n / 2 + 1 <= n entries */
                    queue[tail++] = x ^ 1;
                } else {
                    rinv[k + 1] = -1;
                }
                k += 2;
            }
            x = fp[x];
        } while (x != h);
    }
#ifdef ISO_DEBUG
    debug_relabelling(n, fp, r, rinv, k);
#endif
    return k;
}

int iso_key(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
            const int *r, const int *rinv, int k, int *key, int *names)
{
    int i, j, h, pos, next;
    const int *a;

#ifdef ISO_DEBUG
    debug_names(n, names);
#endif
    for (i = 0; i < k; i++) {
        h = rinv[i];
        key[i] = h == -1 ? -1 : r[fp[h]];
    }
    pos = k;
    for (j = 0; j < nd; j++) {
        a = data + (long) j * n;
        for (i = 0; i < k; i++) {
            h = rinv[i];
            if (h != -1)
                key[pos++] = a[h];
        }
    }
    for (j = 0; j < nc; j++) {
        a = classes + (long) j * n;
        next = 0;
        for (i = 0; i < k; i++) {
            h = rinv[i];
            if (h != -1) {
#ifdef ISO_DEBUG
                ISO_CHECK(0 <= a[h] && a[h] < n);
#endif
                if (names[a[h]] == -1)
                    names[a[h]] = next++;
                key[pos++] = names[a[h]];
            }
        }
        for (i = 0; i < k; i++) {
            h = rinv[i];
            if (h != -1)
                names[a[h]] = -1;
        }
    }
#ifdef ISO_DEBUG
    debug_names(n, names);
#endif
    return pos;
}

int iso_canonical(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                  int *key, int *relabelling, int *roots, int *nroots, int *work)
{
    int *r = work;
    int *rinv = work + n;
    int *queue = work + 2 * n;
    int *names = work + 3 * n;
    int *best = key;
    int *cur = work + 4 * n;
    int *tmp;
    int ne, first, root, k, len = 0, count = 0, i, cmp;

    ne = iso_num_edges(n, fp);
    for (first = 0; first < n && fp[first] == -1; first++)
        ;
    if (first == n) {
        if (relabelling != NULL)
            for (i = 0; i < n; i++)
                relabelling[i] = -1;
        if (nroots != NULL)
            *nroots = 0;
        return 0;
    }

    for (i = 0; i < n; i++)
        names[i] = -1;
    for (root = first; root < n; root++) {
        if (fp[root] == -1)
            continue;
        k = iso_relabelling_from(n, fp, root, r, rinv, queue);
        if (k != 2 * ne) {
            /* connectedness does not depend on the root: only the first one
             * can fail */
            return ISO_DISCONNECTED;
        }
        if (count == 0) {
            len = iso_key(n, fp, nd, data, nc, classes, r, rinv, k, best, names);
            cmp = -1;
        } else {
            int l = iso_key(n, fp, nd, data, nc, classes, r, rinv, k, cur, names);
            (void) l;
#ifdef ISO_DEBUG
            ISO_CHECK(l == len);
#endif
            cmp = 0;
            for (i = 0; i < len; i++) {
                if (cur[i] != best[i]) {
                    cmp = cur[i] < best[i] ? -1 : 1;
                    break;
                }
            }
            if (cmp < 0) {
                tmp = best;
                best = cur;
                cur = tmp;
            }
        }
        if (cmp < 0) {
            count = 0;
            if (relabelling != NULL)
                memcpy(relabelling, r, (size_t) n * sizeof(int));
        }
        if (cmp <= 0) {
            if (roots != NULL)
                roots[count] = root;
            count++;
        }
    }
    if (best != key)
        memcpy(key, best, (size_t) len * sizeof(int));
    if (nroots != NULL)
        *nroots = count;
    return len;
}
