/*
 * Canonical keys, canonical labellings and automorphisms of maps. See
 * doc/c/omap.rst for the conventions and the sizes of the buffers, and
 * doc/c/omap_algorithms.rst for the algorithms and their costs.
 *
 * On a connected map (_domap_canonical), the roots of the rule are tried in
 * increasing order. The first root is walked and keyed in full; for the
 * others, the walk is compared with the best key so far as it goes and stops
 * at the first entry larger, so that a rejected root costs the labels it
 * gave (they alone are reset afterwards). A tie gives an automorphism, and
 * the orbits of the automorphisms found are kept in a union-find, so that a
 * root of an orbit already met costs one find.
 *
 * On any map (_domap_iso_search), the search goes component after component.
 * Its memory does not grow with the depth beyond a few ints per level: the
 * state of a level is its root, its first label, the position and length
 * of its block in the current sequence and the counts of class values
 * before it. A level is undone by walking the labels it gave, and its next
 * tied root is found again by scanning the roots after the current one for
 * a block equal to the one of the level, which is still in the sequence.
 * With classes and without all_leaves, a tied root found so is skipped when
 * an earlier tied root of the level extends the numbering of the class
 * values alike (see "Extensions of tied roots" in doc/c/omap_algorithms.rst);
 * the earlier roots are found again by the same scan, from the start of the
 * class, so that this needs no memory either.
 *
 * Under a rule other than DOMAP_ISO_ROOTS_ALL, _domap_canonical needs the
 * chosen class only, in increasing order. packed_class packs the entries of
 * the invariant of each active half-edge, shifted by their minimum, into
 * one int in mixed radix over their ranges, whose order is the
 * lexicographic order of the invariants; when the product of the ranges is
 * at most 2 n, it counts the half-edges per packed value, chooses the class
 * from these counts and collects its half-edges by one scan, in
 * O(n (1 + nd)) with one cycle walk for the vertices and one for the faces.
 * Otherwise, and always in _domap_iso_search, which needs every class, the
 * classes are computed by root_classes: a radix sort of the active
 * half-edges by invariant, O(n) per entry. Under DOMAP_ISO_ROOTS_DATA,
 * whose invariants end with data values, arbitrary ints, root_classes first
 * replaces each data array by order-preserving numbers of 0, ..., n - 1
 * (O(n) when its values span fewer than n values, O(n log n) otherwise), so
 * that every entry of the radix sort lies in 0, ..., n - 1. The flag is
 * masked wherever the base rule is tested. _domap_canonical then tries the
 * half-edges of the chosen class only. _domap_iso_search keeps the number of
 * half-edges not labelled yet in each class, updated as components are
 * labelled and undone, and chooses the class at each step from these
 * numbers (O(number of classes)); its scans run over the half-edges of that
 * class. Under DOMAP_ISO_ROOTS_ALL the loops over the roots are those
 * without classes.
 *
 * With OMAP_DEBUG defined, the relabellings and keys are checked as they are
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
#include "omap.h"
#include "misc.h"
#include "perm.h"

#include <string.h>

#ifdef OMAP_DEBUG
#include <stdio.h>
#include <stdlib.h>

#define OMAP_CHECK(cond)                                                    \
    do                                                                      \
    {                                                                       \
        if (!(cond))                                                        \
        {                                                                   \
            fprintf(stderr, "omap.c:%d: check failed: %s\n",                \
                    __LINE__, #cond);                                       \
            abort();                                                        \
        }                                                                   \
    } while (0)

/* Check that (r, rinv) is a relabelling of the reached edges by k labels. */
static void
debug_relabelling(const int *fp, int n, const int *r, const int *rinv, int k)
{
    int h, i;
    OMAP_CHECK(k % 2 == 0 && 0 <= k && k <= n);
    for (h = 0; h < n; h++)
    {
        OMAP_CHECK(-1 <= r[h] && r[h] < k);
        OMAP_CHECK((r[h] == -1) == (r[h ^ 1] == -1));
        if (r[h] != -1)
            OMAP_CHECK(r[h ^ 1] == (r[h] ^ 1));
        if (r[h] != -1 && fp[h] != -1)
            OMAP_CHECK(r[fp[h]] != -1);
    }
    for (i = 0; i < k; i++)
    {
        if (rinv[i] == -1)
        {
            OMAP_CHECK(i % 2 == 1 && rinv[i ^ 1] != -1
                       && fp[rinv[i ^ 1] ^ 1] == -1);
        }
        else
        {
            OMAP_CHECK(0 <= rinv[i] && rinv[i] < n && fp[rinv[i]] != -1
                       && r[rinv[i]] == i);
        }
    }
}

/* Check that names is all -1. */
static void
debug_names(const int *names, int n)
{
    int v;
    for (v = 0; v < n; v++)
        OMAP_CHECK(names[v] == -1);
}
#endif

/*
 * The inverse of the vertex permutation is h -> fp[t], t the twin of h (h
 * itself on a folded edge): vp[fp[t]] = h, one write per active h.
 */
void
_omap_vp_from_fp(int *vp, const int *fp, int n)
{
    int h, t;

    for (h = 0; h < n; h++)
        vp[h] = -1;
    for (h = 0; h < n; h++)
    {
        if (fp[h] == -1)
            continue;
        t = fp[h ^ 1] != -1 ? h ^ 1 : h;
        vp[fp[t]] = h;
    }
}

int
_omap_num_edges(const int *p, int n)
{
    int e, ne = 0;
    for (e = 0; 2 * e < n; e++)
        if (p[2 * e] != -1 || p[2 * e + 1] != -1)
            ne++;
    return ne;
}

int
_omap_relabelling_from(int *r, int *rinv, const int *p, int n, int root,
                       int *queue)
{
    int head = 0, tail = 0, k = 0;
    int h, x;

    for (h = 0; h < n; h++)
        r[h] = -1;
    queue[tail++] = root;
    while (head < tail)
    {
        h = queue[head++];
        x = h;
        do
        {
            if (r[x] == -1)
            {
                r[x] = k;
                r[x ^ 1] = k + 1;
                rinv[k] = x;
                if (p[x ^ 1] != -1)
                {
                    rinv[k + 1] = x ^ 1;
                    /* one push per edge, when it gets its labels, and the
                     * root: at most n / 2 + 1 <= n entries */
                    queue[tail++] = x ^ 1;
                }
                else
                {
                    rinv[k + 1] = -1;
                }
                k += 2;
            }
            x = p[x];
        } while (x != h);
    }
#ifdef OMAP_DEBUG
    debug_relabelling(p, n, r, rinv, k);
#endif
    return k;
}

int
omap_relabelling_from(int *r, int *rinv, const omap_t m, int root)
{
    /* allocates the queue of the raw form */
    int *queue = combisurf_malloc(
        (size_t) OMAP_RELABELLING_WORK_SIZE(m->n) * sizeof(int));
    int k = _omap_relabelling_from(r, rinv, m->fp, m->n, root, queue);

    free(queue);
    return k;
}

int
_domap_iso_key(int *key, const domap_t m, const int *r, const int *rinv, int k,
               int *names)
{
    const int *fp = m->map.fp;
    int n = m->map.n, nd = m->nd, nc = m->nc;
    const int *data = m->data, *classes = m->classes;
    int i, j, h, pos, next;
    const int *a;

#ifdef OMAP_DEBUG
    debug_names(names, n);
#endif
    for (i = 0; i < k; i++)
    {
        h = rinv[i];
        key[i] = h == -1 ? -1 : r[fp[h]];
    }
    pos = k;
    for (j = 0; j < nd; j++)
    {
        a = data + (long) j * n;
        for (i = 0; i < k; i++)
        {
            h = rinv[i];
            if (h != -1)
                key[pos++] = a[h];
        }
    }
    for (j = 0; j < nc; j++)
    {
        a = classes + (long) j * n;
        next = 0;
        for (i = 0; i < k; i++)
        {
            h = rinv[i];
            if (h != -1)
            {
#ifdef OMAP_DEBUG
                OMAP_CHECK(0 <= a[h] && a[h] < n);
#endif
                if (names[a[h]] == -1)
                    names[a[h]] = next++;
                key[pos++] = names[a[h]];
            }
        }
        for (i = 0; i < k; i++)
        {
            h = rinv[i];
            if (h != -1)
                names[a[h]] = -1;
        }
    }
#ifdef OMAP_DEBUG
    debug_names(names, n);
#endif
    return pos;
}

int
domap_iso_key(int *key, const domap_t m, const int *r, const int *rinv, int k)
{
    /* allocates names, which the raw form needs filled with -1 */
    int n = m->map.n;
    int *names = combisurf_malloc((size_t) DOMAP_ISO_KEY_WORK_SIZE(n)
                                  * sizeof(int));
    int i, ret;

    for (i = 0; i < n; i++)
        names[i] = -1;
    ret = _domap_iso_key(key, m, r, rinv, k, names);
    free(names);
    return ret;
}

/* A bijective mixing of the 32 bits of x into 64 bits (multiplications by
 * odd constants and xorshifts, each invertible). */
static inline unsigned long long
mix_value(int x)
{
    unsigned long long y =
        (unsigned long long) (unsigned int) x * 0x9E3779B97F4A7C15ULL;
    y ^= y >> 29;
    y *= 0xBF58476D1CE4E5B9ULL;
    return y ^ (y >> 32);
}

/*
 * The hash of the multiset of the values of a on the active half-edges of
 * fp: the sum modulo 2^64 of their mixed values, independent of their order.
 * Equal multisets have equal hashes, so different hashes prove different
 * multisets. One pass without a buffer or a branch: comparing the sorted
 * values instead, counted per value or heap sorted, needs a buffer and
 * branches.
 */
static unsigned long long
values_hash(const int *fp, int n, const int *a)
{
    unsigned long long s = 0;
    int h;
    for (h = 0; h < n; h++)
        s += fp[h] != -1 ? mix_value(a[h]) : 0;
    return s;
}

/*
 * Return whether the class arrays a1 on fp1 and a2 on fp2 (values in 0,
 * ..., n1 - 1 and 0, ..., n2 - 1) have different multisets of class sizes
 * on the active half-edges, as many on both sides. sz (n ints) and cnt
 * (n + 1 ints) are zero on entry and, when they are equal, on return.
 */
static int
class_sizes_differ(const int *fp1, int n1, const int *a1, const int *fp2,
                   int n2, const int *a2, int *sz, int *cnt)
{
    int h;

    for (h = 0; h < n1; h++)
        if (fp1[h] != -1)
            sz[a1[h]]++;
    /* each class counted at its first half-edge, its size then reset */
    for (h = 0; h < n1; h++)
    {
        if (fp1[h] != -1 && sz[a1[h]] != 0)
        {
            cnt[sz[a1[h]]]++;
            sz[a1[h]] = 0;
        }
    }
    for (h = 0; h < n2; h++)
        if (fp2[h] != -1)
            sz[a2[h]]++;
    /* the sizes sum to the same number of half-edges on both sides, so no
     * count going negative means equal counts */
    for (h = 0; h < n2; h++)
    {
        if (fp2[h] != -1 && sz[a2[h]] != 0)
        {
            if (--cnt[sz[a2[h]]] < 0)
                return 1;
            sz[a2[h]] = 0;
        }
    }
    return 0;
}

/*
 * Compare a hash of the multiset of the values of each data array, then the
 * multiset of the class sizes of each class array, O(n) each, stopping at
 * the first difference. The degrees of the cells are not compared: one more
 * pass over each map, for few pairs rejected.
 * See doc/c/omap_algorithms.rst, section "Invariants of two maps".
 */
int
_domap_invariants_differ(const domap_t m1, const domap_t m2, int *work)
{
    const int *fp1 = m1->map.fp, *fp2 = m2->map.fp;
    int n1 = m1->map.n, n2 = m2->map.n;
    int nd = m1->nd, nc = m1->nc;
    int n = n1 > n2 ? n1 : n2;
    int j;

    for (j = 0; j < nd; j++)
        if (values_hash(fp1, n1, m1->data + (long) j * n1)
            != values_hash(fp2, n2, m2->data + (long) j * n2))
            return 1;
    if (nc == 0)
        return 0;
    for (j = 0; j < 2 * n + 1; j++)
        work[j] = 0;
    /* each check leaves the buffers zero when it finds no difference */
    for (j = 0; j < nc; j++)
        if (class_sizes_differ(fp1, n1, m1->classes + (long) j * n1, fp2, n2,
                               m2->classes + (long) j * n2, work, work + n))
            return 1;
    return 0;
}

int
domap_invariants_differ(const domap_t m1, const domap_t m2)
{
    int n = m1->map.n > m2->map.n ? m1->map.n : m2->map.n;
    int *work = NULL;
    int ret;

    /* the raw form reads its work buffer only with classes */
    if (m1->nc > 0)
        work = combisurf_malloc((size_t) DOMAP_INVARIANTS_WORK_SIZE(n)
                                * sizeof(int));
    ret = _domap_invariants_differ(m1, m2, work);
    free(work);
    return ret;
}

/*
 * The entry j of the invariant of the active half-edge h (see the rules of
 * the roots in doc/c/omap.rst), shifted into 0, ..., n - 1: the degrees
 * minus 1, and loop and monogon as 0 or 1. vid and fid are the indices of
 * the vertices and of the faces, vdeg and fdeg their degrees.
 */
static int
invariant_entry(const int *fp, const int *vid, const int *vdeg,
                const int *fid, const int *fdeg, int j, int h)
{
    int t = fp[h ^ 1] != -1 ? h ^ 1 : h;
    switch (j)
    {
        case 0:
            return vdeg[h] - 1;
        case 1:
            return fdeg[h] - 1;
        case 2:
            return vdeg[t] - 1;
        case 3:
            return fdeg[t] - 1;
        case 4:
            return vid[h] == vid[t];
        default:
            return fid[h] == fid[t];
    }
}

/*
 * The entries of the invariants of a rule other than DOMAP_ISO_ROOTS_ALL,
 * each in 0, ..., n - 1: the base entries of invariant_entry, then, under
 * DOMAP_ISO_ROOTS_DATA, the order-preserving numbers of data_ranks for the
 * data values at h (nd entries), and under the edge rules those at the twin
 * t (nd more). rank holds the nd arrays of numbers, nd being 0 without the
 * flag.
 */
struct invariants
{
    const int *fp, *vid, *vdeg, *fid, *fdeg, *rank;
    int n, nd, base, edge;
};

/* The entry j of the invariant of the active half-edge h. */
static inline int
entry(const struct invariants *iv, int j, int h)
{
    if (j < iv->base)
        return invariant_entry(iv->fp, iv->vid, iv->vdeg, iv->fid, iv->fdeg,
                               j, h);
    j -= iv->base;
    if (j >= iv->nd)
    {
        j -= iv->nd;
        if (iv->fp[h ^ 1] != -1)
            h ^= 1;
    }
    return iv->rank[(long) j * iv->n + h];
}

/*
 * Sift down the entry i of the max-heap a[0:len] ordered by the values
 * v[a[.]].
 */
static void
sift_down(int *a, const int *v, int i, int len)
{
    int c, x = a[i];
    while ((c = 2 * i + 1) < len)
    {
        if (c + 1 < len && v[a[c + 1]] > v[a[c]])
            c++;
        if (v[a[c]] <= v[x])
            break;
        a[i] = a[c];
        i = c;
    }
    a[i] = x;
}

/*
 * Write into rank (n ints), at each active half-edge of fp, a value of
 * 0, ..., n - 1 that increases with the value of the array a there (not
 * written on the inactive half-edges), so that equal values of a, and only
 * those, get equal numbers. When the values of a on the active half-edges
 * span fewer than n values (max - min < n, compared in long) the number is
 * a - min, in O(n); the numbers are then not dense. Otherwise it is the rank
 * among the distinct values, in 0, ..., k - 1: tmp (n ints) receives the
 * active half-edges, sorted by value by a heap sort, which needs no buffer
 * beyond them and runs in O(n log n) whatever the values.
 */
static void
data_ranks(int *rank, const int *fp, int n, const int *a, int *tmp)
{
    int h, i, x, na = 0, k = 0, lo = 0, hi = 0;

    for (h = 0; h < n; h++)
    {
        if (fp[h] == -1)
            continue;
        if (na == 0 || a[h] < lo)
            lo = a[h];
        if (na == 0 || a[h] > hi)
            hi = a[h];
        na++;
    }
    /* hi - lo reaches 2^32 - 1, which a 32-bit long does not hold */
    if ((long long) hi - (long long) lo < (long long) n)
    {
        for (h = 0; h < n; h++)
            if (fp[h] != -1)
                rank[h] = (int) ((long long) a[h] - (long long) lo);
        return;
    }
    na = 0;
    for (h = 0; h < n; h++)
        if (fp[h] != -1)
            tmp[na++] = h;
    for (i = na / 2 - 1; i >= 0; i--)
        sift_down(tmp, a, i, na);
    for (i = na - 1; i > 0; i--)
    {
        x = tmp[0];
        tmp[0] = tmp[i];
        tmp[i] = x;
        sift_down(tmp, a, 0, i);
    }
    for (i = 0; i < na; i++)
    {
        if (i > 0 && a[tmp[i]] != a[tmp[i - 1]])
            k++;
        rank[tmp[i]] = k;
    }
}

#ifdef OMAP_DEBUG
/*
 * Compare the invariants of the active half-edges h1 and h2 computed
 * directly, from the degrees and the data values themselves rather than
 * their ranks: -1, 0 or 1.
 */
static int
debug_compare_invariants(const struct invariants *iv, const int *data, int h1,
                         int h2)
{
    int j, x, y, side, t1, t2;
    for (j = 0; j < iv->base; j++)
    {
        x = invariant_entry(iv->fp, iv->vid, iv->vdeg, iv->fid, iv->fdeg, j,
                            h1);
        y = invariant_entry(iv->fp, iv->vid, iv->vdeg, iv->fid, iv->fdeg, j,
                            h2);
        if (x != y)
            return x < y ? -1 : 1;
    }
    for (side = 0; side <= iv->edge; side++)
    {
        t1 = side && iv->fp[h1 ^ 1] != -1 ? h1 ^ 1 : h1;
        t2 = side && iv->fp[h2 ^ 1] != -1 ? h2 ^ 1 : h2;
        for (j = 0; j < iv->nd; j++)
        {
            x = data[(long) j * iv->n + t1];
            y = data[(long) j * iv->n + t2];
            if (x != y)
                return x < y ? -1 : 1;
        }
    }
    return 0;
}

/*
 * Check the class pos_cls[i] of each position i of order (na active
 * half-edges) against the invariants computed directly: the classes start
 * at 0 and grow by at most 1 from a position to the next, the half-edges of
 * one class have equal invariants and come in increasing order, and the
 * invariants increase from a class to the next.
 */
static void
debug_classes(const struct invariants *iv, const int *data,
              const int *order, const int *pos_cls, int na)
{
    int i, cmp;
    OMAP_CHECK(pos_cls[0] == 0);
    for (i = 1; i < na; i++)
    {
        cmp = debug_compare_invariants(iv, data, order[i - 1], order[i]);
        if (pos_cls[i] == pos_cls[i - 1])
        {
            OMAP_CHECK(cmp == 0 && order[i - 1] < order[i]);
        }
        else
        {
            OMAP_CHECK(pos_cls[i] == pos_cls[i - 1] + 1 && cmp < 0);
        }
    }
}
#endif

/*
 * The classes of the rule (not DOMAP_ISO_ROOTS_ALL) on the active
 * half-edges of fp, from the invariants of the whole map. Write into order
 * the active half-edges sorted by invariant, the half-edges of equal
 * invariant in increasing order, into cls (n ints) the class of each
 * half-edge, -1 on the inactive ones, the classes being numbered 0, 1, ...
 * in increasing order of invariant, and into start the start of each class
 * in order, followed by the number of active half-edges (the class c is
 * order[start[c]], ..., order[start[c + 1] - 1]). order holds n ints, start
 * n + 1 and tmp 4 n. Under DOMAP_ISO_ROOTS_DATA, the nd data arrays (nd n
 * ints) are read and ranks (nd n ints) receives their numbers of
 * data_ranks. Return the number of classes.
 *
 * The sort is a least significant digit radix sort: one stable counting
 * sort per entry of the invariant (2 or 6 of them, plus nd or 2 nd data
 * entries under DOMAP_ISO_ROOTS_DATA, each in 0, ..., n - 1), from the last
 * to the first, O(n) each.
 */
static int
root_classes(int *order, int *cls, int *start, int *ranks, const int *vp,
             const int *fp, int n, const int *data, int nd,
             domap_iso_root_rule rule, int *tmp)
{
    int *vdeg = tmp, *fdeg = tmp + n, *aux = tmp + 2 * n, *count = tmp + 3 * n;
    /* the indices of the vertices and of the faces, until the classes are
     * written over them */
    int *vid = cls, *fid = start;
    int base_rule = rule & ~DOMAP_ISO_ROOTS_DATA;
    int edge = base_rule == DOMAP_ISO_ROOTS_EDGE_MAX
        || base_rule == DOMAP_ISO_ROOTS_EDGE_FEWEST;
    struct invariants iv;
    int entries;
    int h, t, i, j, v, sum, na = 0, m;
    int *src = order, *dst = aux, *swap;

    if (!(rule & DOMAP_ISO_ROOTS_DATA))
        nd = 0;
    iv.fp = fp;
    iv.vid = vid;
    iv.vdeg = vdeg;
    iv.fid = fid;
    iv.fdeg = fdeg;
    iv.rank = ranks;
    iv.n = n;
    iv.nd = nd;
    iv.base = edge ? 6 : 2;
    iv.edge = edge;
    entries = iv.base + (edge ? 2 * nd : nd);

    /* the data numbers first, while aux is free */
    for (j = 0; j < nd; j++)
        data_ranks(ranks + (long) j * n, fp, n, data + (long) j * n, aux);

    _perm_dense_cycles(vid, vdeg, vp, n);
    _perm_dense_cycles(fid, fdeg, fp, n);

    for (h = 0; h < n; h++)
        if (fp[h] != -1)
            order[na++] = h;
    if (na == 0)
    {
        for (h = 0; h < n; h++)
            cls[h] = -1;
        start[0] = 0;
        return 0;
    }

    for (j = entries - 1; j >= 0; j--)
    {
        for (v = 0; v < n; v++)
            count[v] = 0;
        for (i = 0; i < na; i++)
            count[entry(&iv, j, src[i])]++;
        sum = 0;
        for (v = 0; v < n; v++)
        {
            t = count[v];
            count[v] = sum;
            sum += t;
        }
        for (i = 0; i < na; i++)
        {
            h = src[i];
            dst[count[entry(&iv, j, h)]++] = h;
        }
        swap = src;
        src = dst;
        dst = swap;
    }
    /* an even number of passes leaves the result in order; an odd one (2 +
     * nd passes with nd odd) leaves it in aux, copied back */
    if (src != order)
        memcpy(order, src, (size_t) na * sizeof(int));

    /* the class of each position, in aux */
    aux[0] = 0;
    for (i = 1; i < na; i++)
    {
        aux[i] = aux[i - 1];
        for (j = 0; j < entries; j++)
        {
            if (entry(&iv, j, order[i - 1]) != entry(&iv, j, order[i]))
            {
                aux[i]++;
                break;
            }
        }
    }
#ifdef OMAP_DEBUG
    /* before the classes are written over vid and fid */
    debug_classes(&iv, data, order, aux, na);
#endif
    m = aux[na - 1] + 1;
    for (h = 0; h < n; h++)
        cls[h] = -1;
    for (i = 0; i < na; i++)
    {
        cls[order[i]] = aux[i];
        if (i == 0 || aux[i] != aux[i - 1])
            start[aux[i]] = i;
    }
    start[m] = na;
    return m;
}

/*
 * Write into deg (n ints, 0 on entry) the length of the cycle of each active
 * half-edge in vp (vertex nonzero) or in fp, and when id is not NULL the
 * index of that cycle into id (neither written on the inactive half-edges);
 * store the smallest and the largest length into *lo and *hi (fp has an
 * active half-edge). The cycles are indexed in the order of their smallest
 * half-edge. Each cycle is walked once, its half-edges kept in stack (n
 * ints) for the writes.
 */
static inline void
cycle_degrees(int *id, int *deg, int *lo, int *hi, const int *vp,
              const int *fp, int n, int vertex, int *stack)
{
    int h, x, i, len, c = 0;

    for (h = 0; h < n; h++)
    {
        if (fp[h] == -1 || deg[h] != 0)
            continue;
        len = 0;
        x = h;
        do
        {
            stack[len++] = x;
            x = vertex ? vp[x] : fp[x];
        } while (x != h);
        for (i = 0; i < len; i++)
            deg[stack[i]] = len;
        if (id != NULL)
            for (i = 0; i < len; i++)
                id[stack[i]] = c;
        if (c == 0 || len < *lo)
            *lo = len;
        if (c == 0 || len > *hi)
            *hi = len;
        c++;
    }
}

/*
 * Multiply *prod by range and return whether the product is at most bound.
 * With *prod at most bound < 2^31 and range at most 2^32, the product fits
 * in long long. A test by division, range > bound / *prod, costs an integer
 * division where a multiplication in long long suffices.
 */
static inline int
radix_fits(long long *prod, long long range, long long bound)
{
    *prod *= range;
    return *prod <= bound;
}

static int choose_class(const int *cnt, int m, domap_iso_root_rule rule);

/*
 * The class of the rule (not DOMAP_ISO_ROOTS_ALL) chosen by choose_class among
 * the classes of root_classes, without sorting: every entry of the
 * invariant is shifted by its minimum over the active half-edges (the
 * degrees, the data values; loop and monogon are 0 or 1) and the entries
 * are packed in mixed radix over their ranges, the first entry most
 * significant, so that the order of the packed values is the lexicographic
 * order of the invariants. When the product of the ranges is at most
 * bound, the half-edges are counted per packed value in count (bound
 * ints), the class is chosen from these counts, and a scan writes its
 * half-edges into out in increasing order; their number is returned. When
 * the product exceeds bound, or when 2 nd > n + 1, 0 is returned and out,
 * which serves as a work buffer before, holds nothing. tmp holds 4 n ints
 * (the degrees, and the indices of the vertices and of the faces under the
 * edge rules), pv n ints (the packed values) and aux n + 1 ints (the
 * minimum and the range of each data array).
 *
 * The cost is O(n (1 + nd) + bound): one walk of each cycle of the vertices
 * and of the faces, one scan per data array for its range, one scan packing
 * and counting, one scan collecting the class; root_classes runs one
 * counting sort of n buckets per entry of the invariant (2 or 6, plus nd or
 * 2 nd data entries) and its data ranks first.
 */
static int
packed_class(int *out, const int *vp, const int *fp, int n, const int *data,
             int nd, domap_iso_root_rule rule, long long bound, int *tmp,
             int *pv, int *aux, int *count)
{
    int *vdeg = tmp, *fdeg = tmp + n, *vid = tmp + 2 * n, *fid = tmp + 3 * n;
    int *dmin = aux, *drange = aux + nd;
    int base_rule = rule & ~DOMAP_ISO_ROOTS_DATA;
    int edge = base_rule == DOMAP_ISO_ROOTS_EDGE_MAX
        || base_rule == DOMAP_ISO_ROOTS_EDGE_FEWEST;
    int h, t, j, c, v, k, vmin = 0, vmax = 0, fmin = 0,
        fmax = 0, lo, hi, rv, rf;
    long long prod = 1;
    const int *a;

    if (!(rule & DOMAP_ISO_ROOTS_DATA))
        nd = 0;
    if (2 * (long) nd > (long) n + 1)
        return 0;
    /* out serves as the stack of the cycle walks */
    for (h = 0; h < n; h++)
    {
        vdeg[h] = 0;
        fdeg[h] = 0;
    }
    cycle_degrees(edge ? vid : NULL, vdeg, &vmin, &vmax, vp, fp, n, 1, out);
    cycle_degrees(edge ? fid : NULL, fdeg, &fmin, &fmax, vp, fp, n, 0, out);
    /* the degrees are at most n, their ranges fit in int; under the edge
     * rules, the degrees at the twins span the same values, the twins of
     * the active half-edges being the active half-edges */
    rv = vmax - vmin + 1;
    rf = fmax - fmin + 1;
    if (!radix_fits(&prod, rv, bound) || !radix_fits(&prod, rf, bound))
        return 0;
    if (edge && (!radix_fits(&prod, rv, bound) || !radix_fits(&prod, rf, bound)
                 || !radix_fits(&prod, 4, bound)))
        return 0;
    /* the data values span the same values at h and at its twin */
    for (j = 0; j < nd; j++)
    {
        a = data + (long) j * n;
        lo = hi = 0;
        k = 0;
        for (h = 0; h < n; h++)
        {
            if (fp[h] == -1)
                continue;
            if (k == 0 || a[h] < lo)
                lo = a[h];
            if (k == 0 || a[h] > hi)
                hi = a[h];
            k = 1;
        }
        if (!radix_fits(&prod, (long long) hi - lo + 1, bound))
            return 0;
        if (edge && !radix_fits(&prod, (long long) hi - lo + 1, bound))
            return 0;
        dmin[j] = lo;
        drange[j] = (int) ((long long) hi - lo + 1);
    }

    /* every value below is below prod <= bound, which fits in int */
    memset(count, 0, (size_t) prod * sizeof(int));
    for (h = 0; h < n; h++)
    {
        if (fp[h] == -1)
            continue;
        v = (vdeg[h] - vmin) * rf + (fdeg[h] - fmin);
        t = fp[h ^ 1] != -1 ? h ^ 1 : h;
        if (edge)
            v = ((v * rv + (vdeg[t] - vmin)) * rf + (fdeg[t] - fmin)) * 4
                + 2 * (vid[h] == vid[t]) + (fid[h] == fid[t]);
        for (j = 0; j < nd; j++)
            v = v * drange[j]
                + (int) ((long long) data[(long) j * n + h] - dmin[j]);
        if (edge)
            for (j = 0; j < nd; j++)
                v = v * drange[j]
                    + (int) ((long long) data[(long) j * n + t] - dmin[j]);
        pv[h] = v;
        count[v]++;
    }
    c = choose_class(count, (int) prod, rule);
    k = 0;
    for (h = 0; h < n; h++)
        if (fp[h] != -1 && pv[h] == c)
            out[k++] = h;
    return k;
}

/*
 * The class tried among the m classes of root_classes, cnt[c] being the
 * number of half-edges of the class c that may be roots (at least one class
 * has one): under *_MAX the largest class of positive count, under *_FEWEST
 * the class of smallest positive count, ties broken by the largest class.
 * The flag DOMAP_ISO_ROOTS_DATA does not change the choice.
 */
static int
choose_class(const int *cnt, int m, domap_iso_root_rule rule)
{
    int c, best = -1;
    int base_rule = rule & ~DOMAP_ISO_ROOTS_DATA;

    if (base_rule == DOMAP_ISO_ROOTS_HALF_EDGE_MAX
        || base_rule == DOMAP_ISO_ROOTS_EDGE_MAX)
    {
        for (c = m - 1; cnt[c] == 0; c--)
            ;
        return c;
    }
    for (c = 0; c < m; c++)
        if (cnt[c] > 0 && (best == -1 || cnt[c] <= cnt[best]))
            best = c;
    return best;
}

#ifdef OMAP_DEBUG
/*
 * Check the class of packed_class, the half-edges out[0:len], against the
 * class chosen among the classes of root_classes on the same input, in
 * buffers of its own.
 */
static void
debug_packed_class(const int *out, int len, const int *vp, const int *fp,
                   int n, const int *data, int nd, domap_iso_root_rule rule)
{
    int *buf = malloc(((size_t) 8 * n + 1 + (size_t) nd * n) * sizeof(int));
    int *order = buf, *cls = buf + n, *start = buf + 2 * n,
        *tmp = buf + 3 * n + 1;
    int *ranks = buf + 7 * n + 1;
    int m, c, i;

    OMAP_CHECK(buf != NULL);
    m = root_classes(order, cls, start, ranks, vp, fp, n, data, nd, rule,
                     tmp);
    /* the sizes of the classes, in tmp */
    for (c = 0; c < m; c++)
        tmp[c] = start[c + 1] - start[c];
    c = choose_class(tmp, m, rule);
    OMAP_CHECK(len == start[c + 1] - start[c]);
    for (i = 0; i < len; i++)
        OMAP_CHECK(out[i] == order[start[c] + i]);
    free(buf);
}
#endif

/*
 * The walk of omap_relabelling_from from root, on r all -1 on entry, with the
 * face permutation part of its key compared with best as it is built. The
 * entry i of that part is r[fp[rinv[i]]], or -1 at an inactive label. The
 * entries below the label of the next half-edge of the queue are final: an
 * even label is on a half-edge of a face walked, and an odd one is inactive
 * or on a half-edge taken from the queue. The final entries are compared
 * with best after a face walk once their number beyond those compared
 * reaches a step that starts at 2 and doubles after each comparison, and
 * at the end of the walk. A rejected root has thus given about twice
 * the labels needed to make the differing entry final, plus those of the
 * half-edges still waiting in the queue, and a tied root runs O(log n)
 * comparison passes. The doubling step avoids one comparison pass per face
 * on maps where every root ties.
 *
 * At the first entry larger than the one of best, the walk stops and 1 is
 * returned; at the first entry smaller, the comparison stops, the walk
 * goes on to the end and -1 is returned; 0 is returned when the walk ends
 * with every entry equal.
 *
 * The number of labels given is stored into *k; r and rinv hold them, and
 * reset_local sets r back to -1.
 */
static int
walk_compare(int *r, int *rinv, int *k, const int *fp, const int *best,
             int root, int *queue)
{
    int head = 0, tail = 0, used = 0, p = 0, cmp = 0, step = 1, next = 0;
    int h, x, v, lim;

    queue[tail++] = root;
    while (head < tail)
    {
        h = queue[head++];
        x = h;
        do
        {
            if (r[x] == -1)
            {
                r[x] = used;
                r[x ^ 1] = used + 1;
                rinv[used] = x;
                if (fp[x ^ 1] != -1)
                {
                    rinv[used + 1] = x ^ 1;
                    queue[tail++] = x ^ 1;
                }
                else
                {
                    rinv[used + 1] = -1;
                }
                used += 2;
            }
            x = fp[x];
        } while (x != h);
        if (cmp != 0)
            continue;
        lim = head < tail ? r[queue[head]] : used;
        if (lim < next && head < tail)
            continue;
        step *= 2;
        next = lim + step;
        for (; p < lim; p++)
        {
            x = rinv[p];
            v = x == -1 ? -1 : r[fp[x]];
#ifdef OMAP_DEBUG
            OMAP_CHECK(v != -1 || x == -1);
#endif
            if (v != best[p])
            {
                if (v > best[p])
                {
                    *k = used;
                    return 1;
                }
                cmp = -1;
                break;
            }
        }
    }
#ifdef OMAP_DEBUG
    OMAP_CHECK(cmp != 0 || p == used);
#endif
    *k = used;
    return cmp;
}

/*
 * The data and class parts of the key of (r, rinv) using k labels, written
 * into cur and compared with best as they are built, both from the position
 * k of the key. Return 1 at the first entry larger than the one of best,
 * with cur only partly written; otherwise write the parts in full and
 * return -1 or 0 as the comparison. names is all -1 on entry and on return.
 */
static int
tail_compare(int *cur, const domap_t m, const int *rinv, const int *best,
             int k, int *names)
{
    int n = m->map.n, nd = m->nd, nc = m->nc;
    const int *data = m->data, *classes = m->classes;
    int i, j, h, v, pos = k, next, cmp = 0;
    const int *a;

    for (j = 0; j < nd; j++)
    {
        a = data + (long) j * n;
        for (i = 0; i < k; i++)
        {
            h = rinv[i];
            if (h == -1)
                continue;
            v = a[h];
            if (cmp == 0 && v != best[pos])
            {
                if (v > best[pos])
                    return 1;
                cmp = -1;
            }
            cur[pos++] = v;
        }
    }
    for (j = 0; j < nc; j++)
    {
        a = classes + (long) j * n;
        next = 0;
        for (i = 0; i < k; i++)
        {
            h = rinv[i];
            if (h == -1)
                continue;
            if (names[a[h]] == -1)
                names[a[h]] = next++;
            v = names[a[h]];
            if (cmp == 0 && v != best[pos])
            {
                if (v > best[pos])
                {
                    /* the values named so far are on the labels 0, ..., i */
                    for (; i >= 0; i--)
                        if (rinv[i] != -1)
                            names[a[rinv[i]]] = -1;
                    return 1;
                }
                cmp = -1;
            }
            cur[pos++] = v;
        }
        for (i = 0; i < k; i++)
        {
            h = rinv[i];
            if (h != -1)
                names[a[h]] = -1;
        }
    }
    return cmp;
}

static void reset_local(int *lr, const int *lrinv, int k);

#ifdef OMAP_DEBUG
/*
 * Check the outcome cmp of the root against the key computed in full by
 * omap_relabelling_from and _domap_iso_key (len ints): it compares with
 * best as cmp says, and equals cur when cmp < 0. scratch holds 3 n +
 * DOMAP_ISO_KEY_SIZE(n, nd, nc) ints.
 */
static void
debug_root(const domap_t m, const int *best, const int *cur, int len,
           int root, int cmp, int *names, int *scratch)
{
    int n = m->map.n;
    int *r2 = scratch, *rinv2 = scratch + n, *queue2 = scratch + 2 * n,
        *key2 = scratch + 3 * n;
    int k2, len2, i, ref = 0;

    k2 = _omap_relabelling_from(r2, rinv2, m->map.fp, n, root, queue2);
    len2 = _domap_iso_key(key2, m, r2, rinv2, k2, names);
    OMAP_CHECK(len2 == len);
    for (i = 0; i < len; i++)
    {
        if (key2[i] != best[i])
        {
            ref = key2[i] < best[i] ? -1 : 1;
            break;
        }
    }
    OMAP_CHECK(ref == cmp);
    if (cmp < 0)
        for (i = 0; i < len; i++)
            OMAP_CHECK(cur[i] == key2[i]);
}

/*
 * Check that h -> rinv[rbest[h]] on the active half-edges, from the
 * relabellings of two roots of equal keys, is an automorphism of the
 * decorated map: a permutation of the active half-edges that commutes with
 * fp and with the twin, keeps the data arrays and the class arrays up to a
 * bijection of their values. scratch holds 4 n ints.
 */
static void
debug_automorphism(const domap_t m, const int *rbest, const int *rinv,
                   int *scratch)
{
    const int *fp = m->map.fp;
    int n = m->map.n, nd = m->nd, nc = m->nc;
    const int *data = m->data, *classes = m->classes;
    int *aut = scratch, *seen = scratch + n, *fwd = scratch + 2 * n,
        *bwd = scratch + 3 * n;
    int h, j, v, w;
    const int *a;

    for (h = 0; h < n; h++)
        seen[h] = 0;
    for (h = 0; h < n; h++)
    {
        if (fp[h] == -1)
            continue;
        OMAP_CHECK(0 <= rbest[h] && rbest[h] < n);
        aut[h] = rinv[rbest[h]];
        OMAP_CHECK(0 <= aut[h] && aut[h] < n && fp[aut[h]] != -1
                   && !seen[aut[h]]);
        seen[aut[h]] = 1;
    }
    for (h = 0; h < n; h++)
    {
        if (fp[h] == -1)
            continue;
        OMAP_CHECK(aut[fp[h]] == fp[aut[h]]);
        if (fp[h ^ 1] != -1)
            OMAP_CHECK(aut[h ^ 1] == (aut[h] ^ 1));
        else
            OMAP_CHECK(fp[aut[h] ^ 1] == -1);
    }
    for (j = 0; j < nd; j++)
    {
        a = data + (long) j * n;
        for (h = 0; h < n; h++)
            if (fp[h] != -1)
                OMAP_CHECK(a[aut[h]] == a[h]);
    }
    for (j = 0; j < nc; j++)
    {
        a = classes + (long) j * n;
        for (v = 0; v < n; v++)
        {
            fwd[v] = -1;
            bwd[v] = -1;
        }
        for (h = 0; h < n; h++)
        {
            if (fp[h] == -1)
                continue;
            v = a[h];
            w = a[aut[h]];
            OMAP_CHECK((fwd[v] == -1 && bwd[w] == -1)
                       || (fwd[v] == w && bwd[w] == v));
            fwd[v] = w;
            bwd[w] = v;
        }
    }
}
#endif

/*
 * The state of _domap_canonical while it tries the roots after the first one:
 * the buffers, the best key so far and the candidate key (swapped when a
 * root beats the best one), the length of the keys, the best root and its
 * relabelling, and the parents of the union-find of the orbits of the
 * automorphisms found, set up at the first tie only (merged nonzero): until
 * then every orbit is a single half-edge and no root is skipped, so that a
 * map whose roots are all rejected pays nothing for the union-find.
 */
struct canonical_state
{
    domap_struct m;
    int ne, len, best_root, merged;
    int *r, *rinv, *queue, *names, *best, *cur, *rbest, *parent;
#ifdef OMAP_DEBUG
    int *scratch;
#endif
};

/*
 * Try root after the first one: skip it when its orbit has a smaller
 * half-edge, which has the same key and was tried or skipped before;
 * otherwise compare its key with the best one, merge the orbits by the
 * automorphism it gives when it ties, and make it the best root when it
 * beats it. r is -1 on entry and on return.
 *
 * A skipped root has a key larger than the best one, or is in the orbit of
 * the best root: the earlier root of its key was rejected, or tied with or
 * was the best root then, and the later best roots only have smaller keys.
 * Hence a root that ties is outside the orbit of the best root.
 */
static inline void
canonical_root(struct canonical_state *s, int root)
{
    const int *fp = s->m.map.fp;
    int n = s->m.map.n, k, cmp, h;
    int *r = s->r, *rinv = s->rinv, *tmp;

    if (s->merged && uf_find_halving(s->parent, root) < root)
    {
#ifdef OMAP_DEBUG
        cmp = uf_find_halving(s->parent, root)
            == uf_find_halving(s->parent, s->best_root) ? 0 : 1;
        debug_root(&s->m, s->best, s->cur, s->len, root, cmp, s->names,
                   s->scratch);
#endif
        return;
    }

    cmp = walk_compare(r, rinv, &k, fp, s->best, root, s->queue);
    if (cmp < 0)
    {
        _domap_iso_key(s->cur, &s->m, r, rinv, k, s->names);
    }
    else if (cmp == 0)
    {
        cmp = tail_compare(s->cur, &s->m, rinv, s->best, k, s->names);
        if (cmp < 0)
            memcpy(s->cur, s->best, (size_t) k * sizeof(int));
    }
#ifdef OMAP_DEBUG
    OMAP_CHECK(cmp == 1 || k == 2 * s->ne);
    debug_root(&s->m, s->best, s->cur, s->len, root, cmp, s->names,
               s->scratch);
#endif
    if (cmp < 0)
    {
        tmp = s->best;
        s->best = s->cur;
        s->cur = tmp;
        s->best_root = root;
        memcpy(s->rbest, r, (size_t) n * sizeof(int));
    }
    else if (cmp == 0)
    {
        /* the automorphism h -> rinv[rbest[h]] sends the best root to root */
        if (!s->merged)
        {
            for (h = 0; h < n; h++)
                s->parent[h] = h;
            s->merged = 1;
        }
#ifdef OMAP_DEBUG
        OMAP_CHECK(uf_find_halving(s->parent, root)
                   != uf_find_halving(s->parent, s->best_root));
        debug_automorphism(&s->m, s->rbest, rinv, s->scratch);
#endif
        for (h = 0; h < n; h++)
            if (fp[h] != -1)
                uf_link(s->parent, uf_find_halving(s->parent, h),
                        uf_find_halving(s->parent, rinv[s->rbest[h]]));
    }
    /* past n / 2 labels, a sequential clear (all bytes 0xff, the int
     * -1) is cheaper than the scattered writes of reset_local */
    if (2 * k >= n)
        memset(r, -1, (size_t) n * sizeof(int));
    else
        reset_local(r, rinv, k);
#ifdef OMAP_DEBUG
    {
        int i;
        for (i = 0; i < n; i++)
            OMAP_CHECK(r[i] == -1);
        debug_names(s->names, n);
    }
#endif
}

/*
 * Search the minimal key over the roots of the class of the rule. A root is
 * rejected at the first entry of its key larger than the best one, in walks
 * of doubling length; a tie gives an automorphism, whose orbits skip the
 * later roots of a known key, so that there are at most log2 |Aut| ties at
 * O(n (1 + nd + nc)) each.
 * See doc/c/omap_algorithms.rst, section "Canonical key of a connected map".
 */
int
_domap_canonical(int *key, int *relabelling, int *roots, int *nroots,
                 const domap_t m, domap_iso_root_rule rule, int *work)
{
    const int *fp = m->map.fp;
    int n = m->map.n, nd = m->nd, nc = m->nc;
    const int *data = m->data;
    int *r = work;
    int *rinv = work + n;
    int *queue = work + 2 * n;
    int *names = work + 3 * n;
    int *cur = work + 4 * n;
    int *order = cur + DOMAP_ISO_KEY_SIZE(n, nd, nc);
    int *cls = order + n;
    int *start = cls + n;
    int *rbest = start + n + 1;
    int *parent = rbest + n;
    struct canonical_state s;
    int all = (rule & ~DOMAP_ISO_ROOTS_DATA) == DOMAP_ISO_ROOTS_ALL;
    int ne, first, root, k, i, c, ncls, rep, count, lo = 0, hi = 0;

    ne = _omap_num_edges(fp, n);
    for (first = 0; first < n && fp[first] == -1; first++)
        ;
    if (first == n)
    {
        if (relabelling != NULL)
            for (i = 0; i < n; i++)
                relabelling[i] = -1;
        if (nroots != NULL)
            *nroots = 0;
        return 0;
    }

    if (!all)
    {
        /* r, rinv, queue and names (4 n) are free until the first walk, cls
         * and start until the walks (they serve root_classes only), cur (the
         * numbers of the data values of root_classes) until the first key,
         * rbest and parent (2 n, contiguous) until the first key too. The
         * packed values of packed_class go into cls, the ranges of the data
         * into start, its counts into rbest and parent, which bounds the
         * product of the ranges by 2 n; the sizes of the classes of
         * root_classes go into r. */
        hi = packed_class(order, m->map.vp, fp, n, data, nd, rule,
                          2 * (long long) n, work, cls, start, rbest);
        if (hi == 0)
        {
            ncls = root_classes(order, cls, start, cur, m->map.vp, fp, n,
                                data, nd, rule, work);
            for (c = 0; c < ncls; c++)
                r[c] = start[c + 1] - start[c];
            c = choose_class(r, ncls, rule);
            lo = start[c];
            hi = start[c + 1];
        }
#ifdef OMAP_DEBUG
        else
        {
            debug_packed_class(order, hi, m->map.vp, fp, n, data, nd, rule);
        }
#endif
        first = order[lo];
    }

    for (i = 0; i < n; i++)
        names[i] = -1;

    /* the first root, walked in full: connectedness does not depend on the
     * root, so only this walk can fail */
    k = _omap_relabelling_from(r, rinv, fp, n, first, queue);
    if (k != 2 * ne)
        return DOMAP_ISO_DISCONNECTED;
    s.m = *m;
    s.ne = ne;
    s.r = r;
    s.rinv = rinv;
    s.queue = queue;
    s.names = names;
    s.best = key;
    s.cur = cur;
    s.rbest = rbest;
    s.parent = parent;
    s.best_root = first;
    s.merged = 0;
    s.len = _domap_iso_key(s.best, m, r, rinv, k, names);
    memcpy(rbest, r, (size_t) n * sizeof(int));
    reset_local(r, rinv, k);
#ifdef OMAP_DEBUG
    s.scratch = malloc(((size_t) 3 * n
                        + (size_t) DOMAP_ISO_KEY_SIZE(n, nd, nc))
                       * sizeof(int));
    OMAP_CHECK(s.scratch != NULL);
#endif

    /* the other roots: r is -1 on entry of each, and only the labels a root
     * gave are reset after it. The loop is chosen once per call, so that
     * DOMAP_ISO_ROOTS_ALL tests no class per root. */
    if (all)
    {
        for (root = first + 1; root < n; root++)
            if (fp[root] != -1)
                canonical_root(&s, root);
    }
    else
    {
        for (i = lo + 1; i < hi; i++)
            canonical_root(&s, order[i]);
    }
#ifdef OMAP_DEBUG
    free(s.scratch);
#endif
    if (s.best != key)
        memcpy(key, s.best, (size_t) s.len * sizeof(int));
    if (relabelling != NULL)
        memcpy(relabelling, rbest, (size_t) n * sizeof(int));

    /* the roots of minimal key: the orbit of the best root, which is its
     * smallest half-edge (a smaller one would have been tried first) */
    if (!s.merged)
    {
        if (roots != NULL)
            roots[0] = s.best_root;
        count = 1;
    }
    else
    {
        rep = uf_find_halving(parent, s.best_root);
#ifdef OMAP_DEBUG
        OMAP_CHECK(rep == s.best_root);
#endif
        count = 0;
        if (all)
        {
            for (root = rep; root < n; root++)
            {
                if (fp[root] != -1
                    && uf_find_halving(parent, root) == rep)
                {
                    if (roots != NULL)
                        roots[count] = root;
                    count++;
                }
            }
        }
        else
        {
            /* the orbit lies in the class, increasing in order[lo..hi) */
            for (i = lo; i < hi; i++)
            {
                root = order[i];
                if (uf_find_halving(parent, root) == rep)
                {
                    if (roots != NULL)
                        roots[count] = root;
                    count++;
                }
            }
        }
    }
    if (nroots != NULL)
        *nroots = count;
    return s.len;
}

int
domap_canonical(int *key, int *relabelling, int *roots, int *nroots,
                const domap_t m, domap_iso_root_rule rule)
{
    /* allocates the work buffer of the raw form */
    int *work = combisurf_malloc(
        (size_t) DOMAP_CANONICAL_WORK_SIZE(m->map.n, m->nd, m->nc)
        * sizeof(int));
    int ret = _domap_canonical(key, relabelling, roots, nroots, m, rule, work);

    free(work);
    return ret;
}

/*
 * The breadth-first walk of omap_relabelling_from restricted to the
 * component of root: lr must be -1 on the half-edges of that component and
 * is left set there (reset_local undoes it), so that a walk costs the size
 * of the component and not n.
 */
static int
local_walk(int *lr, int *lrinv, const int *fp, int root, int *queue)
{
    int head = 0, tail = 0, k = 0;
    int h, x;

    queue[tail++] = root;
    while (head < tail)
    {
        h = queue[head++];
        x = h;
        do
        {
            if (lr[x] == -1)
            {
                lr[x] = k;
                lr[x ^ 1] = k + 1;
                lrinv[k] = x;
                if (fp[x ^ 1] != -1)
                {
                    lrinv[k + 1] = x ^ 1;
                    queue[tail++] = x ^ 1;
                }
                else
                {
                    lrinv[k + 1] = -1;
                }
                k += 2;
            }
            x = fp[x];
        } while (x != h);
    }
    return k;
}

/* Set lr back to -1 on the component walked by local_walk. */
static void
reset_local(int *lr, const int *lrinv, int k)
{
    int i, h;
    /* an even label is always on an active half-edge, the one met in a face */
    for (i = 0; i < k; i += 2)
    {
        h = lrinv[i];
        lr[h] = -1;
        lr[h ^ 1] = -1;
    }
}

/*
 * Write into block the block of the component walked by local_walk (k
 * labels), the class values numbered from next[j] on in the numbering names
 * (nc arrays of n ints, -1 at the values without number). With commit, the
 * numbering keeps the new values and next is updated; otherwise names is
 * left as on entry. Return the length of the block.
 */
static int
component_block(int *block, const domap_t m, const int *lr, const int *lrinv,
                int k, int commit, int *names, int *next)
{
    const int *fp = m->map.fp;
    int n = m->map.n, nd = m->nd, nc = m->nc;
    const int *data = m->data, *classes = m->classes;
    int i, j, h, pos = 0, count;
    const int *a;
    int *nm;

    block[pos++] = k;
    for (i = 0; i < k; i++)
    {
        h = lrinv[i];
        block[pos++] = h == -1 ? -1 : lr[fp[h]];
    }
    for (j = 0; j < nd; j++)
    {
        a = data + (long) j * n;
        for (i = 0; i < k; i++)
        {
            h = lrinv[i];
            if (h != -1)
                block[pos++] = a[h];
        }
    }
    for (j = 0; j < nc; j++)
    {
        a = classes + (long) j * n;
        nm = names + (long) j * n;
        count = next[j];
        for (i = 0; i < k; i++)
        {
            h = lrinv[i];
            if (h != -1)
            {
#ifdef OMAP_DEBUG
                OMAP_CHECK(0 <= a[h] && a[h] < n);
#endif
                if (nm[a[h]] == -1)
                    nm[a[h]] = count++;
                block[pos++] = nm[a[h]];
            }
        }
        if (commit)
        {
            next[j] = count;
        }
        else
        {
            /* the values numbered by this component are those numbered
             * from next[j] on */
            for (i = 0; i < k; i++)
            {
                h = lrinv[i];
                if (h != -1 && nm[a[h]] >= next[j])
                    nm[a[h]] = -1;
            }
        }
    }
    return pos;
}

/*
 * Compare the flat sequences a (la ints) and b (lb ints)
 * lexicographically: -1, 0 or 1.
 */
static int
compare_sequences(const int *a, int la, const int *b, int lb)
{
    int i, m = la < lb ? la : lb;
    for (i = 0; i < m; i++)
        if (a[i] != b[i])
            return a[i] < b[i] ? -1 : 1;
    return (la > lb) - (la < lb);
}

/*
 * Return whether the sequence a (la ints) exceeds the prefix of b (lb
 * ints) with as many blocks. The blocks of equal prefixes have equal
 * lengths, so the first difference, if any, lies within min(la, lb) ints.
 */
static int
exceeds_prefix(const int *a, int la, const int *b, int lb)
{
    int i, m = la < lb ? la : lb;
    for (i = 0; i < m; i++)
        if (a[i] != b[i])
            return a[i] > b[i];
    return la > lb;
}

/*
 * The block of the component of root, not labelled yet, written into cand
 * (see component_block, without commit); return its length. lr is -1 on
 * entry and on return.
 */
static inline int
root_block(int *cand, const domap_t m, int root, int *lr, int *lrinv,
           int *queue, int *names, int *next)
{
    int k = local_walk(lr, lrinv, m->map.fp, root, queue);
    int len = component_block(cand, m, lr, lrinv, k, 0, names, next);
    reset_local(lr, lrinv, k);
    return len;
}

/*
 * Return whether the block (len ints, written by component_block) numbers a
 * class value from next[j] on in some class array j, that is, whether its
 * component extends the numbering names.
 */
static int
block_extends(const int *block, int len, const int *next, int nd, int nc)
{
    int k = block[0];
    int na = (len - 1 - k) / (nd + nc);
    const int *b = block + 1 + k + (long) nd * na;
    int i, j;

    for (j = 0; j < nc; j++, b += na)
        for (i = 0; i < na; i++)
            if (b[i] >= next[j])
                return 1;
    return 0;
}

/*
 * Number into names, from next[j] on, the class values met first in the
 * component of root (not labelled yet) with set, or set back to -1 those
 * numbered from next[j] on without set; lr is -1 on entry and on return.
 */
static void
number_component(int *names, const domap_t m, const int *next, int root,
                 int set, int *lr, int *lrinv, int *queue)
{
    int n = m->map.n, nc = m->nc;
    const int *classes = m->classes;
    int k = local_walk(lr, lrinv, m->map.fp, root, queue);
    int i, j, h, count;
    const int *a;
    int *nm;

    for (j = 0; j < nc; j++)
    {
        a = classes + (long) j * n;
        nm = names + (long) j * n;
        count = next[j];
        for (i = 0; i < k; i++)
        {
            h = lrinv[i];
            if (h == -1)
                continue;
            if (set && nm[a[h]] == -1)
                nm[a[h]] = count++;
            else if (!set && nm[a[h]] >= next[j])
                nm[a[h]] = -1;
        }
    }
    reset_local(lr, lrinv, k);
}

/*
 * Return whether the roots root1 and root2, not labelled yet, whose blocks
 * are both equal to block (len ints), extend the numbering names (numbers
 * from next[j] on) by the same values in the same order. The numbering of
 * root1 is written into names, and the class values of the component of
 * root2 must then get from it the numbers of block, which are the ones the
 * numbering of root2 gives them: the two numberings agree on every value of
 * the component of root2, hence on the values it numbers, which are all
 * those numbered from next[j] on. names is left as on entry; lr is -1 on
 * entry and on return. Costs three walks of the components.
 */
static int
same_extension(const domap_t m, const int *block, int len, const int *next,
               int root1, int root2, int *lr, int *lrinv, int *queue,
               int *names)
{
    int n = m->map.n, nd = m->nd, nc = m->nc;
    const int *classes = m->classes;
    int k = block[0];
    int na = (len - 1 - k) / (nd + nc);
    const int *b = block + 1 + k + (long) nd * na;
    int i, j, h, idx, same = 1;
    const int *a;
    const int *nm;

    number_component(names, m, next, root1, 1, lr, lrinv, queue);
    i = local_walk(lr, lrinv, m->map.fp, root2, queue);
#ifdef OMAP_DEBUG
    OMAP_CHECK(i == k);
#endif
    for (j = 0; j < nc && same; j++, b += na)
    {
        a = classes + (long) j * n;
        nm = names + (long) j * n;
        for (i = 0, idx = 0; i < k; i++)
        {
            h = lrinv[i];
            if (h == -1)
                continue;
            if (nm[a[h]] != b[idx++])
            {
                same = 0;
                break;
            }
        }
    }
    reset_local(lr, lrinv, k);
    number_component(names, m, next, root1, 0, lr, lrinv, queue);
    return same;
}

/*
 * Search the minimal sequence of blocks depth first, one component per level,
 * following only the roots of minimal block, and cut a state whose sequence
 * exceeds the prefix of the best one. Without all_leaves, a tied root of the
 * same extension of the class values as an earlier one is skipped, and a leaf
 * equal to the best one backjumps to the first level where they differ; the
 * search stays exponential in the worst case.
 * See doc/c/omap_algorithms.rst, section "Search over the components".
 *
 * The levels of the search store the position of their root in the scan of
 * the roots: the half-edge itself under DOMAP_ISO_ROOTS_ALL, its index in
 * order (see root_classes) under the other rules. The scans are written
 * twice, once per kind of position, so that DOMAP_ISO_ROOTS_ALL tests no
 * class per root.
 */
int
_domap_iso_search(int *sequence, int *length, int *relabelling,
                  const domap_t m, domap_iso_root_rule rule, int all_leaves,
                  domap_iso_leaf_callback callback, void *arg, int *work)
{
    const int *fp = m->map.fp;
    int n = m->map.n, nd = m->nd, nc = m->nc;
    const int *data = m->data, *classes = m->classes;
    long levels = n / 2 + 1;
    int *r = work;
    int *rinv = r + n;
    int *lr = rinv + n;
    int *lrinv = lr + n;
    int *queue = lrinv + n;
    int *names = queue + n;
    int *lvl_root = names + (long) nc * n;
    int *lvl_used = lvl_root + levels;
    int *lvl_pos = lvl_used + levels;
    int *lvl_len = lvl_pos + levels;
    int *lvl_best = lvl_len + levels;
    int *lvl_next = lvl_best + levels;
    int *next = lvl_next + (long) nc * levels;
    int *seq = next + nc;
    int *cand = seq + DOMAP_ISO_SEQUENCE_SIZE(n, nd, nc);
    int *order = cand + DOMAP_ISO_KEY_SIZE(n, nd, nc) + 1;
    int *cls = order + n;
    int *start = cls + n;
    int *cnt = start + n + 1;
    int follow_ties = nc > 0 || all_leaves;
    int cut_ties = nc > 0 && !all_leaves;
    int all = (rule & ~DOMAP_ISO_ROOTS_DATA) == DOMAP_ISO_ROOTS_ALL;
    int ne, depth = 0, used = 0, seqlen = 0, bestlen = -1, ncls = 0, target, t;
    int root, best_root, blen, len, k, cmp, ret, i, j, h, c, found;

    ne = _omap_num_edges(fp, n);
    if (!all)
    {
        /* r, rinv, lr and lrinv (4 n) are free until they are cleared below,
         * cand (the numbers of the data values) until the first block;
         * cnt[c] is the number of half-edges of the class c not labelled */
        ncls = root_classes(order, cls, start, cand, m->map.vp, fp, n, data,
                            nd, rule, work);
        for (c = 0; c < ncls; c++)
            cnt[c] = start[c + 1] - start[c];
    }
    for (h = 0; h < n; h++)
    {
        r[h] = -1;
        lr[h] = -1;
    }
    for (i = 0; i < nc * n; i++)
        names[i] = -1;
    for (j = 0; j < nc; j++)
        next[j] = 0;

    for (;;)
    {
        /* enter the state: depth levels labelled, used labels, the blocks
         * in seq[0:seqlen] */
        int cut = bestlen >= 0
            && exceeds_prefix(seq, seqlen, sequence, bestlen);

        /* the deepest level whose further tied roots are tried when the
         * search backtracks from this state */
        target = depth - 1;
        if (!cut && used == 2 * ne)
        {
            cmp = bestlen < 0
                ? -1 : compare_sequences(seq, seqlen, sequence, bestlen);
            if (cmp < 0)
            {
                memcpy(sequence, seq, (size_t) seqlen * sizeof(int));
                bestlen = seqlen;
                memcpy(lvl_best, lvl_root, (size_t) depth * sizeof(int));
                if (relabelling != NULL)
                    memcpy(relabelling, r, (size_t) n * sizeof(int));
                if (callback != NULL && (ret = callback(r, 1, arg)) != 0)
                    return ret;
            }
            else if (cmp == 0 && all_leaves)
            {
                if (callback != NULL && (ret = callback(r, 0, arg)) != 0)
                    return ret;
            }
            else if (cmp == 0)
            {
                /* Backjump. Both leaves have the same sequence, so the
                 * relabelling of the best leaf followed by the inverse of
                 * the one of this leaf is an automorphism of fp and data,
                 * and of classes up to a relabelling of their values.
                 * It fixes the components labelled at the levels before the
                 * first one, t, where the roots of the two leaves differ, and
                 * with them the state at level t. The root of this leaf at t
                 * comes after the one of the best leaf, whose subtree has been
                 * searched; the subtree of the root of this leaf is its image,
                 * with the same sequences, none smaller than the best. */
                for (t = 0; t < depth && lvl_root[t] == lvl_best[t]; t++)
                    ;
#ifdef OMAP_DEBUG
                OMAP_CHECK(t < depth && lvl_root[t] > lvl_best[t]);
#endif
                target = t;
            }
            cut = 1;
        }

        if (!cut)
        {
            /* the minimal block over the roots not labelled yet, written
             * in place in seq, and its first root */
            best_root = -1;
            blen = 0;
            if (all)
            {
                for (root = 0; root < n; root++)
                {
                    if (fp[root] == -1 || r[root] != -1)
                        continue;
                    len = root_block(cand, m, root, lr, lrinv, queue, names,
                                     next);
                    if (best_root == -1
                        || compare_sequences(cand, len, seq + seqlen,
                                             blen) < 0)
                    {
                        memcpy(seq + seqlen, cand, (size_t) len * sizeof(int));
                        blen = len;
                        best_root = root;
                    }
                }
            }
            else
            {
                c = choose_class(cnt, ncls, rule);
                for (i = start[c]; i < start[c + 1]; i++)
                {
                    root = order[i];
                    if (r[root] != -1)
                        continue;
                    len = root_block(cand, m, root, lr, lrinv, queue, names,
                                     next);
                    if (best_root == -1
                        || compare_sequences(cand, len, seq + seqlen,
                                             blen) < 0)
                    {
                        memcpy(seq + seqlen, cand, (size_t) len * sizeof(int));
                        blen = len;
                        best_root = i;
                    }
                }
            }
#ifdef OMAP_DEBUG
            OMAP_CHECK(best_root != -1 && depth < levels);
#endif
            lvl_root[depth] = best_root;
            lvl_used[depth] = used;
            lvl_pos[depth] = seqlen;
            lvl_len[depth] = blen;
            for (j = 0; j < nc; j++)
                lvl_next[(long) depth * nc + j] = next[j];
        }
        else
        {
            /* backtrack to the deepest level with a further tied root */
            found = 0;
            while (!found && depth > 0)
            {
                depth--;
                if (!all)
                    for (i = lvl_used[depth]; i < used; i++)
                        if (rinv[i] != -1)
                            cnt[cls[rinv[i]]]++;
                for (i = lvl_used[depth]; i < used; i++)
                {
                    h = rinv[i];
                    if (h == -1)
                        continue;
                    for (j = 0; j < nc; j++)
                    {
                        int *nm = names + (long) j * n;
                        const int *a = classes + (long) j * n;
                        if (nm[a[h]] >= lvl_next[(long) depth * nc + j])
                            nm[a[h]] = -1;
                    }
                    r[h] = -1;
                    r[h ^ 1] = -1;
                }
                for (j = 0; j < nc; j++)
                    next[j] = lvl_next[(long) depth * nc + j];
                used = lvl_used[depth];
                seqlen = lvl_pos[depth];
                if (!follow_ties || depth > target)
                    continue;
                /* Without all_leaves, a tied root is skipped when an earlier
                 * tied root of the level extends the numbering of the class
                 * values by the same values in the same order. Both
                 * components then carry the same values at the same local
                 * labels, so that exchanging them (or relabelling the
                 * component from one root to the other when it is the same
                 * one), the identity elsewhere, is an automorphism of fp,
                 * data and classes fixing the components labelled before.
                 * It sends the subtree of the later root onto the one of
                 * the earlier root, searched first with the same minimum.
                 * The extensions of tied roots have equal lengths: when
                 * they are empty, every tied root after the first is
                 * skipped. */
                blen = lvl_len[depth];
                if (cut_ties
                    && !block_extends(seq + seqlen, blen, next, nd, nc))
                    continue;
                /* the tied roots after the current one: their block is
                 * the one of the level, already in seq; under a rule, they
                 * are in the class of the root of the level, the one chosen
                 * at the level since the state is back to it */
                if (all)
                {
                    for (root = lvl_root[depth] + 1; root < n && !found;
                         root++)
                    {
                        if (fp[root] == -1 || r[root] != -1)
                            continue;
                        len = root_block(cand, m, root, lr, lrinv, queue,
                                         names, next);
                        if (len != blen
                            || compare_sequences(cand, len, seq + seqlen,
                                                 len) != 0)
                            continue;
                        found = 1;
                        for (t = 0; cut_ties && t < root; t++)
                        {
                            if (fp[t] == -1 || r[t] != -1)
                                continue;
                            len = root_block(cand, m, t, lr, lrinv, queue,
                                             names, next);
                            if (len == blen
                                && compare_sequences(cand, len, seq + seqlen,
                                                     len) == 0
                                && same_extension(m, seq + seqlen, blen, next,
                                                  t, root, lr, lrinv, queue,
                                                  names))
                            {
                                found = 0;
                                break;
                            }
                        }
                        if (found)
                            lvl_root[depth] = root;
                    }
                }
                else
                {
                    c = cls[order[lvl_root[depth]]];
                    for (i = lvl_root[depth] + 1; i < start[c + 1] && !found;
                         i++)
                    {
                        root = order[i];
                        if (r[root] != -1)
                            continue;
                        len = root_block(cand, m, root, lr, lrinv, queue,
                                         names, next);
                        if (len != blen
                            || compare_sequences(cand, len, seq + seqlen,
                                                 len) != 0)
                            continue;
                        found = 1;
                        for (t = start[c]; cut_ties && t < i; t++)
                        {
                            if (r[order[t]] != -1)
                                continue;
                            len = root_block(cand, m, order[t], lr, lrinv,
                                             queue, names, next);
                            if (len == blen
                                && compare_sequences(cand, len, seq + seqlen,
                                                     len) == 0
                                && same_extension(m, seq + seqlen, blen, next,
                                                  order[t], root, lr, lrinv,
                                                  queue, names))
                            {
                                found = 0;
                                break;
                            }
                        }
                        if (found)
                            lvl_root[depth] = i;
                    }
                }
            }
            if (!found)
                break;
        }

        /* label the component of the root of the level depth */
        root = all ? lvl_root[depth] : order[lvl_root[depth]];
        k = local_walk(lr, lrinv, fp, root, queue);
        for (i = 0; i < k; i += 2)
        {
            h = lrinv[i];
            r[h] = used + i;
            r[h ^ 1] = used + i + 1;
            rinv[used + i] = h;
            rinv[used + i + 1] = lrinv[i + 1];
        }
        component_block(cand, m, lr, lrinv, k, 1, names, next);
        if (!all)
            for (i = 0; i < k; i++)
                if (lrinv[i] != -1)
                    cnt[cls[lrinv[i]]]--;
        reset_local(lr, lrinv, k);
        used += k;
        seqlen += lvl_len[depth];
        depth++;
    }

#ifdef OMAP_DEBUG
    /* every level undone */
    OMAP_CHECK(bestlen >= 0 && used == 0 && seqlen == 0);
    for (h = 0; h < n; h++)
        OMAP_CHECK(r[h] == -1 && lr[h] == -1);
    for (i = 0; i < nc * n; i++)
        OMAP_CHECK(names[i] == -1);
    for (c = 0; c < ncls; c++)
        OMAP_CHECK(cnt[c] == start[c + 1] - start[c]);
#endif
    *length = bestlen;
    return 0;
}

int
domap_iso_search(int *sequence, int *length, int *relabelling,
                 const domap_t m, domap_iso_root_rule rule, int all_leaves,
                 domap_iso_leaf_callback callback, void *arg)
{
    /* allocates the work buffer of the raw form */
    int *work = combisurf_malloc(
        (size_t) DOMAP_ISO_SEARCH_WORK_SIZE(m->map.n, m->nd, m->nc)
        * sizeof(int));
    int ret = _domap_iso_search(sequence, length, relabelling, m, rule,
                                all_leaves, callback, arg, work);

    free(work);
    return ret;
}

/*
 * The state of _domap_stabiliser: the map, its layout, the maps of values fw
 * and bw of each class array (-1 where unset, the array q at the offset
 * q * n), the trail of the entries of fw set, in order, and its length
 * ntrail, the copies already images (used), and per level its image target,
 * the index choice of its label among the labels of the roots of its kind,
 * the length mark of the trail before its step and the flag dead, set when
 * a placement of the level adding no entry to the trail had no completion
 * since the level was entered (the later placements adding none are then
 * skipped, see extend). rb, rbinv, queue, perm and ident are buffers of the
 * size of the largest component, ident holding the identity.
 */
struct stabiliser_state
{
    domap_struct m;
    domap_wreath_struct w;
    int ntrail;
    int *fw, *bw, *trail, *used, *target, *choice, *mark, *dead;
    int *rb, *rbinv, *queue, *perm, *ident;
};

/*
 * Write into perm the automorphism of the first copy of the kind t sending
 * the half-edge of label 0 to the one of label l, on the labels: perm[i] is
 * the label of the image of the half-edge of label i.
 */
static void
relabel(struct stabiliser_state *s, int t, int l)
{
    int o = s->w.coff[t], size = s->w.csize[t], h;
    const int *cra = s->w.cra + o;
    int *rb = s->rb, *rbinv = s->rbinv, *perm = s->perm;

    /* _omap_relabelling_from writes -1 into rbinv at the labels of the
     * inactive halves of folded edges, which are rebuilt from rb */
    _omap_relabelling_from(rb, rbinv, s->w.cfp + o, size, s->w.crainv[o + l],
                           s->queue);
    for (h = 0; h < size; h++)
        rbinv[rb[h]] = h;
    for (h = 0; h < size; h++)
        perm[h] = cra[rbinv[h]];
}

/*
 * Extend the maps of values by the isomorphism from the level a to the level
 * b sending the half-edge of label i to the one of label perm[i]; return 0
 * on a conflict (the entries set before it stay on the trail) and 1
 * otherwise.
 */
static int
step(struct stabiliser_state *s, int a, int b, const int *perm)
{
    int n = s->m.map.n, nc = s->m.nc, size = s->w.csize[s->w.kind[a]];
    const int *f = s->m.map.fp, *cls = s->m.classes;
    const int *la = s->w.lab + s->w.off[a], *lb = s->w.lab + s->w.off[b];
    int *fw = s->fw, *bw = s->bw, *trail = s->trail;
    int i, q, h, y, v, w, base;

    for (i = 0; i < size; i++)
    {
        h = la[i];
        if (f[h] == -1)
            continue;
        y = lb[perm[i]];
        for (q = 0; q < nc; q++)
        {
            base = q * n;
            v = cls[base + h];
            w = cls[base + y];
            if (fw[base + v] == -1)
            {
                if (bw[base + w] != -1)
                    return 0;
                fw[base + v] = w;
                bw[base + w] = v;
                trail[s->ntrail++] = base + v;
            }
            else if (fw[base + v] != w)
                return 0;
        }
    }
    return 1;
}

/*
 * Unset the entries of the maps of values set after the first length
 * entries of the trail.
 */
static void
undo(struct stabiliser_state *s, int length)
{
    int n = s->m.map.n, e;
    int *fw = s->fw, *bw = s->bw;

    while (s->ntrail > length)
    {
        e = s->trail[--s->ntrail];
        bw[(e / n) * n + fw[e]] = -1;
        fw[e] = -1;
    }
}

/*
 * Search images for the levels i + 1, ..., ncomp - 1, each a copy of its
 * kind not used yet, through any label of its roots, in order; return 1 at
 * the first leaf, the choices left in target and choice, and 0 if there is
 * none. A placement adding no entry to the trail is skipped when one such
 * placement of its level already failed since the level was entered (dead).
 */
static int
extend(struct stabiliser_state *s, int i)
{
    const int *kind = s->w.kind, *start = s->w.start, *loff = s->w.loff,
        *labels = s->w.labels;
    int *used = s->used, *target = s->target, *choice = s->choice,
        *mark = s->mark, *dead = s->dead;
    int ncomp = s->w.ncomp, j = i + 1, t, b, placed;

    if (j < ncomp)
    {
        target[j] = start[kind[j]];
        choice[j] = 0;
        dead[j] = 0;
    }
    for (;;)
    {
        if (j == ncomp)
            return 1;
        if (j <= i)
            return 0;
        t = kind[j];
        placed = 0;
        while (target[j] < start[t + 1])
        {
            b = target[j];
            if (used[b] || choice[j] >= loff[t + 1] - loff[t])
            {
                target[j]++;
                choice[j] = 0;
                continue;
            }
            mark[j] = s->ntrail;
            relabel(s, t, labels[loff[t] + choice[j]]);
            if (step(s, j, b, s->perm))
            {
                /* Two placements of the level adding no entry to the maps
                 * of values, onto the copies b1 and b2 (or twice onto one
                 * copy), put equal values at the labels they match: the
                 * automorphism of the copy relating their labels is an
                 * isomorphism from b1 to b2 carrying the values, and
                 * exchanging b1 and b2 through it sends a completion of
                 * either placement onto one of the other. */
                if (s->ntrail == mark[j] && dead[j])
                {
                    choice[j]++;
                    continue;
                }
                used[b] = 1;
                placed = 1;
                break;
            }
            undo(s, mark[j]);
            choice[j]++;
        }
        if (placed)
        {
            j++;
            if (j < ncomp)
            {
                target[j] = start[kind[j]];
                choice[j] = 0;
                dead[j] = 0;
            }
        }
        else
        {
            j--;
            if (j > i)
            {
                /* the trail is back to its length after the step of the
                 * level j, whose subtree failed */
                if (s->ntrail == mark[j])
                    dead[j] = 1;
                undo(s, mark[j]);
                used[target[j]] = 0;
                choice[j]++;
            }
        }
    }
}

/*
 * Write into g (n ints) the automorphism of the map given by the choices of
 * the levels i, ..., ncomp - 1 (the identity on the other levels and on the
 * inactive edges).
 */
static void
leaf(int *g, struct stabiliser_state *s, int i)
{
    const int *lab = s->w.lab, *perm = s->perm;
    int n = s->m.map.n, j, t, k, oa, ob, h;

    for (h = 0; h < n; h++)
        g[h] = h;
    for (j = i; j < s->w.ncomp; j++)
    {
        t = s->w.kind[j];
        relabel(s, t, s->w.labels[s->w.loff[t] + s->choice[j]]);
        oa = s->w.off[j];
        ob = s->w.off[s->target[j]];
        for (k = 0; k < s->w.csize[t]; k++)
            g[lab[oa + k]] = lab[ob + perm[k]];
    }
}

#ifdef OMAP_DEBUG
/*
 * Check that the generator g is a permutation of the half-edges, sending the
 * inactive ones to inactive ones, and an automorphism of the decorated map
 * (see debug_automorphism).
 */
static void
debug_generator(const domap_t m, const int *g)
{
    const int *fp = m->map.fp;
    int n = m->map.n, h;
    int *buf = malloc(((size_t) 6 * n + 1) * sizeof(int));
    int *ident = buf, *seen = buf + n, *scratch = buf + 2 * n;

    OMAP_CHECK(buf != NULL);
    for (h = 0; h < n; h++)
    {
        ident[h] = h;
        seen[h] = 0;
    }
    for (h = 0; h < n; h++)
    {
        OMAP_CHECK(0 <= g[h] && g[h] < n && !seen[g[h]]);
        OMAP_CHECK((fp[h] == -1) == (fp[g[h]] == -1));
        seen[g[h]] = 1;
    }
    debug_automorphism(m, g, ident, scratch);
    free(buf);
}
#endif

/*
 * Sims' subgroup search over the base of the canonical labels 0 of the
 * levels, from the last level to the first: each point of the orbit of b_i
 * outside the orbit of the generators found so far is tried by a depth first
 * search for an element sending b_i to it, which skips the placements adding
 * no entry to the bijections of values once one of them failed. Exponential
 * in the worst case.
 * See doc/c/omap_algorithms.rst, section "Stabiliser of the classes".
 */
int
_domap_stabiliser(int *orbits, const domap_t m, const domap_wreath_t w,
                  domap_stabiliser_callback callback, void *arg, int *work)
{
    int n = m->map.n, nc = m->nc, ncomp = w->ncomp;
    int levels = n / 2;
    struct stabiliser_state s;
    int *p, *sz, *g;
    int i, t, b, q, h, x, y, bi, nl, dead, empty, found, ret;
    const int *start = w->start, *loff = w->loff, *labels = w->labels,
        *lab = w->lab, *off = w->off;

    s.m = *m;
    s.w = *w;
    s.ntrail = 0;
    s.fw = work;
    s.bw = s.fw + nc * n;
    s.trail = s.bw + nc * n;
    s.used = s.trail + nc * n;
    s.target = s.used + levels;
    s.choice = s.target + levels;
    s.mark = s.choice + levels;
    s.dead = s.mark + levels;
    s.rb = s.dead + levels;
    s.rbinv = s.rb + n;
    s.queue = s.rbinv + n;
    s.perm = s.queue + n;
    s.ident = s.perm + n;
    p = s.ident + n;
    sz = p + n;
    g = sz + n;

    for (i = 0; i < nc * n; i++)
    {
        s.fw[i] = -1;
        s.bw[i] = -1;
    }
    for (i = 0; i < ncomp; i++)
    {
        s.used[i] = 0;
        s.dead[i] = 0;
    }
    for (h = 0; h < n; h++)
    {
        s.ident[h] = h;
        p[h] = h;
        sz[h] = 1;
    }

    for (i = ncomp - 1; i >= 0; i--)
    {
        t = w->kind[i];
        nl = loff[t + 1] - loff[t];
        bi = lab[off[i]];
        if (start[t + 1] - start[t] == 1 && nl == 1)
        {
            orbits[i] = 1;
            continue;
        }
        /* the levels < i are fixed: their values are mapped to themselves */
        for (b = 0; b < i; b++)
        {
            s.used[b] = 1;
            step(&s, b, b, s.ident);
        }
        dead = 0;
        for (b = start[t] > i ? start[t] : i; b < start[t + 1]; b++)
        {
            for (q = 0; q < nl; q++)
            {
                x = lab[off[b] + labels[loff[t] + q]];
                if (x == bi
                    || uf_find_halving(p, x) == uf_find_halving(p, bi))
                    continue;
                s.target[i] = b;
                s.choice[i] = q;
                s.mark[i] = s.ntrail;
                s.used[b] = 1;
                relabel(&s, t, labels[loff[t] + q]);
                found = 0;
                if (step(&s, i, b, s.perm))
                {
                    /* as in extend: a candidate adding no entry to the maps
                     * of values is the image of an earlier one adding none
                     * by an element of H fixing the levels < i, so it is in
                     * the orbit of b_i exactly when that one is */
                    empty = s.ntrail == s.mark[i];
                    if (!(empty && dead))
                    {
                        found = extend(&s, i);
                        if (empty && !found)
                            dead = 1;
                    }
                }
                if (found)
                {
                    leaf(g, &s, i);
#ifdef OMAP_DEBUG
                    debug_generator(m, g);
#endif
                    if (callback != NULL && (ret = callback(g, arg)) != 0)
                        return ret;
                    for (h = 0; h < n; h++)
                    {
                        x = uf_find_halving(p, h);
                        y = uf_find_halving(p, g[h]);
                        if (x != y)
                        {
                            p[x] = y;
                            sz[y] += sz[x];
                        }
                    }
                }
                undo(&s, s.mark[i]);
                for (h = i; h < ncomp; h++)
                    s.used[h] = 0;
            }
        }
        undo(&s, 0);
        for (b = 0; b < i; b++)
            s.used[b] = 0;
        orbits[i] = sz[uf_find_halving(p, bi)];
    }
    return 0;
}

int
domap_stabiliser(int *orbits, const domap_t m, const domap_wreath_t w,
                 domap_stabiliser_callback callback, void *arg)
{
    /* allocates the work buffer of the raw form */
    int *work = combisurf_malloc(
        (size_t) DOMAP_STABILISER_WORK_SIZE(m->map.n, m->nc) * sizeof(int));
    int ret = _domap_stabiliser(orbits, m, w, callback, arg, work);

    free(work);
    return ret;
}
