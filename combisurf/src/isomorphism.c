/*
 * Canonical keys, canonical labellings and automorphisms of maps. See
 * isomorphism.h for the conventions and the sizes of the buffers.
 *
 * On a connected map (iso_canonical), every active half-edge is tried as a
 * root. The first root is walked and keyed in full; for the others, the
 * walk is compared with the best key so far as it goes and stops at the
 * first entry larger, so that a rejected root costs the labels it gave
 * (they alone are reset afterwards), and a root costs O(n) only when it
 * ties with or beats the best one. A tie gives an automorphism, and the
 * orbits of the automorphisms found are kept in a union-find: a root of an
 * orbit already met costs one find, so that ties are only met outside the
 * orbit of the best root, at most log2 |Aut| of them at O(n) each. A map
 * whose roots all tie (a regular map) costs O(n log2 |Aut|) instead of
 * O(n^2) for n half-edges.
 *
 * On any map (iso_search), the search goes component after component. Its
 * memory does not grow with the depth beyond a few ints per level: the
 * state of a level is its root, its first label, the position and length
 * of its block in the current sequence and the counts of class values
 * before it. A level is undone by walking the labels it gave, and its next
 * tied root is found again by scanning the roots after the current one for
 * a block equal to the one of the level, which is still in the sequence.
 *
 * Under a rule other than ISO_ROOTS_ALL, the classes of the half-edges are
 * computed once per call (root_classes: a radix sort of the active
 * half-edges by invariant, O(n)). Under ISO_ROOTS_DATA, whose invariants
 * end with data values, arbitrary ints, each data array is first replaced
 * by the ranks of its values, O(n log n), so that every entry of the radix
 * sort lies in 0, ..., n - 1. The flag is masked wherever the base rule is
 * tested. iso_canonical then tries the half-edges
 * of the chosen class only. iso_search keeps the number of half-edges not
 * labelled yet in each class, updated as components are labelled and
 * undone, and chooses the class at each step from these numbers (O(number
 * of classes)); its scans run over the half-edges of that class. Under
 * ISO_ROOTS_ALL the loops over the roots are those without classes.
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

/*
 * Write into id the index of the cycle of p (n ints, -1 on the inactive
 * half-edges) of each active half-edge, -1 on the inactive ones, and into
 * size the length of that cycle (not written on the inactive half-edges).
 */
static void cycles(int n, const int *p, int *id, int *size)
{
    int h, x, len, c = 0;

    for (h = 0; h < n; h++)
        id[h] = -1;
    for (h = 0; h < n; h++) {
        if (p[h] == -1 || id[h] != -1)
            continue;
        len = 0;
        x = h;
        do {
            id[x] = c;
            len++;
            x = p[x];
        } while (x != h);
        do {
            size[x] = len;
            x = p[x];
        } while (x != h);
        c++;
    }
}

/*
 * The entry j of the invariant of the active half-edge h (see
 * isomorphism.h), shifted into 0, ..., n - 1: the degrees minus 1, and loop
 * and monogon as 0 or 1. vid and fid are the indices of the vertices and of
 * the faces, vdeg and fdeg their degrees.
 */
static int invariant_entry(int j, int h, const int *fp, const int *vid, const int *vdeg,
                           const int *fid, const int *fdeg)
{
    int t = fp[h ^ 1] != -1 ? h ^ 1 : h;
    switch (j) {
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
 * The entries of the invariants of a rule other than ISO_ROOTS_ALL, each in
 * 0, ..., n - 1: the base entries of invariant_entry, then, under
 * ISO_ROOTS_DATA, the ranks of the data values at h (nd entries), and under
 * the edge rules those at the twin t (nd more). rank holds the nd arrays of
 * ranks, nd being 0 without the flag.
 */
struct invariants {
    const int *fp, *vid, *vdeg, *fid, *fdeg, *rank;
    int n, nd, base, edge;
};

/* The entry j of the invariant of the active half-edge h. */
static inline int entry(const struct invariants *iv, int j, int h)
{
    if (j < iv->base)
        return invariant_entry(j, h, iv->fp, iv->vid, iv->vdeg, iv->fid, iv->fdeg);
    j -= iv->base;
    if (j >= iv->nd) {
        j -= iv->nd;
        if (iv->fp[h ^ 1] != -1)
            h ^= 1;
    }
    return iv->rank[(long) j * iv->n + h];
}

/* Sift down the entry i of the max-heap a[0:len] ordered by the values v[a[.]]. */
static void sift_down(int *a, int i, int len, const int *v)
{
    int c, x = a[i];
    while ((c = 2 * i + 1) < len) {
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
 * Write into rank (n ints) the rank of the value of the array a at each
 * active half-edge of fp among the distinct values of a on the active
 * half-edges, in 0, ..., k - 1 in increasing order of value (not written on
 * the inactive half-edges). tmp holds n ints. The active half-edges are
 * sorted by value by a heap sort, which needs no buffer beyond them and
 * runs in O(n log n) whatever the values.
 */
static void data_ranks(int n, const int *fp, const int *a, int *rank, int *tmp)
{
    int h, i, x, na = 0, k = 0;

    for (h = 0; h < n; h++)
        if (fp[h] != -1)
            tmp[na++] = h;
    for (i = na / 2 - 1; i >= 0; i--)
        sift_down(tmp, i, na, a);
    for (i = na - 1; i > 0; i--) {
        x = tmp[0];
        tmp[0] = tmp[i];
        tmp[i] = x;
        sift_down(tmp, 0, i, a);
    }
    for (i = 0; i < na; i++) {
        if (i > 0 && a[tmp[i]] != a[tmp[i - 1]])
            k++;
        rank[tmp[i]] = k;
    }
}

#ifdef ISO_DEBUG
/*
 * Compare the invariants of the active half-edges h1 and h2 computed
 * directly, from the degrees and the data values themselves rather than
 * their ranks: -1, 0 or 1.
 */
static int debug_compare_invariants(const struct invariants *iv, const int *data, int h1, int h2)
{
    int j, x, y, side, t1, t2;
    for (j = 0; j < iv->base; j++) {
        x = invariant_entry(j, h1, iv->fp, iv->vid, iv->vdeg, iv->fid, iv->fdeg);
        y = invariant_entry(j, h2, iv->fp, iv->vid, iv->vdeg, iv->fid, iv->fdeg);
        if (x != y)
            return x < y ? -1 : 1;
    }
    for (side = 0; side <= iv->edge; side++) {
        t1 = side && iv->fp[h1 ^ 1] != -1 ? h1 ^ 1 : h1;
        t2 = side && iv->fp[h2 ^ 1] != -1 ? h2 ^ 1 : h2;
        for (j = 0; j < iv->nd; j++) {
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
static void debug_classes(const struct invariants *iv, const int *data, int na, const int *order,
                          const int *pos_cls)
{
    int i, cmp;
    ISO_CHECK(pos_cls[0] == 0);
    for (i = 1; i < na; i++) {
        cmp = debug_compare_invariants(iv, data, order[i - 1], order[i]);
        if (pos_cls[i] == pos_cls[i - 1]) {
            ISO_CHECK(cmp == 0 && order[i - 1] < order[i]);
        } else {
            ISO_CHECK(pos_cls[i] == pos_cls[i - 1] + 1 && cmp < 0);
        }
    }
}
#endif

/*
 * The classes of the rule (not ISO_ROOTS_ALL) on the active half-edges of
 * fp, from the invariants of the whole map. Write into order the active
 * half-edges sorted by invariant, the half-edges of equal invariant in
 * increasing order, into cls (n ints) the class of each half-edge, -1 on the
 * inactive ones, the classes being numbered 0, 1, ... in increasing order
 * of invariant, and into start the start of each class in order, followed
 * by the number of active half-edges (the class c is order[start[c]],
 * ..., order[start[c + 1] - 1]). order holds n ints, start n + 1 and tmp
 * 4 n. Under ISO_ROOTS_DATA, the nd data arrays (nd n ints) are read and
 * ranks (nd n ints) receives their ranks. Return the number of classes.
 *
 * The sort is a least significant digit radix sort: one stable counting
 * sort per entry of the invariant (2 or 6 of them, plus nd or 2 nd ranks
 * under ISO_ROOTS_DATA, each in 0, ..., n - 1), from the last to the
 * first, O(n) each.
 */
static int root_classes(int n, const int *fp, iso_root_rule rule, int nd, const int *data,
                        int *order, int *cls, int *start, int *tmp, int *ranks)
{
    int *vdeg = tmp, *fdeg = tmp + n, *aux = tmp + 2 * n, *count = tmp + 3 * n;
    /* the indices of the vertices and of the faces, until the classes are
     * written over them */
    int *vid = cls, *fid = start;
    int base_rule = rule & ~ISO_ROOTS_DATA;
    int edge = base_rule == ISO_ROOTS_EDGE_MAX || base_rule == ISO_ROOTS_EDGE_FEWEST;
    struct invariants iv;
    int entries;
    int h, t, i, j, v, sum, na = 0, m;
    int *src = order, *dst = aux, *swap;

    if (!(rule & ISO_ROOTS_DATA))
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

    /* the ranks first, while aux is free */
    for (j = 0; j < nd; j++)
        data_ranks(n, fp, data + (long) j * n, ranks + (long) j * n, aux);

    /* the vertex permutation, in aux */
    for (h = 0; h < n; h++)
        aux[h] = -1;
    for (h = 0; h < n; h++) {
        if (fp[h] == -1)
            continue;
        t = fp[h ^ 1] != -1 ? h ^ 1 : h;
        aux[fp[t]] = h;
    }
    cycles(n, aux, vid, vdeg);
    cycles(n, fp, fid, fdeg);

    for (h = 0; h < n; h++)
        if (fp[h] != -1)
            order[na++] = h;
    if (na == 0) {
        for (h = 0; h < n; h++)
            cls[h] = -1;
        start[0] = 0;
        return 0;
    }

    for (j = entries - 1; j >= 0; j--) {
        for (v = 0; v < n; v++)
            count[v] = 0;
        for (i = 0; i < na; i++)
            count[entry(&iv, j, src[i])]++;
        sum = 0;
        for (v = 0; v < n; v++) {
            t = count[v];
            count[v] = sum;
            sum += t;
        }
        for (i = 0; i < na; i++) {
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
    for (i = 1; i < na; i++) {
        aux[i] = aux[i - 1];
        for (j = 0; j < entries; j++) {
            if (entry(&iv, j, order[i - 1]) != entry(&iv, j, order[i])) {
                aux[i]++;
                break;
            }
        }
    }
#ifdef ISO_DEBUG
    /* before the classes are written over vid and fid */
    debug_classes(&iv, data, na, order, aux);
#endif
    m = aux[na - 1] + 1;
    for (h = 0; h < n; h++)
        cls[h] = -1;
    for (i = 0; i < na; i++) {
        cls[order[i]] = aux[i];
        if (i == 0 || aux[i] != aux[i - 1])
            start[aux[i]] = i;
    }
    start[m] = na;
    return m;
}

/*
 * The class tried among the m classes of root_classes, cnt[c] being the
 * number of half-edges of the class c that may be roots (at least one class
 * has one): under *_MAX the largest class of positive count, under *_FEWEST
 * the class of smallest positive count, ties broken by the largest class.
 * The flag ISO_ROOTS_DATA does not change the choice.
 */
static int choose_class(iso_root_rule rule, int m, const int *cnt)
{
    int c, best = -1;
    int base_rule = rule & ~ISO_ROOTS_DATA;

    if (base_rule == ISO_ROOTS_HALF_EDGE_MAX || base_rule == ISO_ROOTS_EDGE_MAX) {
        for (c = m - 1; cnt[c] == 0; c--)
            ;
        return c;
    }
    for (c = 0; c < m; c++)
        if (cnt[c] > 0 && (best == -1 || cnt[c] <= cnt[best]))
            best = c;
    return best;
}

/*
 * The walk of iso_relabelling_from from root, on r all -1 on entry, with the
 * face permutation part of its key compared with best as it is built. The
 * entry i of that part is r[fp[rinv[i]]], or -1 at an inactive label. The
 * entries below the label of the next half-edge of the queue are final: an
 * even label is on a half-edge of a face walked, and an odd one is inactive
 * or on a half-edge taken from the queue. The final entries are compared
 * with best after a face walk once their number beyond those compared
 * reaches a step that starts at 2 and doubles after each comparison, and
 * at the end of the walk. A rejected root has thus given about twice
 * the labels needed to make the differing entry final, plus those of the
 * half-edges still waiting in the queue, and a tied root runs O(log n) comparison
 * passes. On regular maps (S_5 to PSL(2, 17), 120 to 2448 half-edges),
 * where every root ties, comparing after every face was 11 to 32 % slower
 * than walking and keying every root in full, and this step 5 to 12 %.
 *
 * At the first entry larger than the one of best, the walk stops and 1 is
 * returned; at the first entry smaller, the comparison stops, the walk
 * goes on to the end and -1 is returned; 0 is returned when the walk ends
 * with every entry equal.
 *
 * The number of labels given is stored into *k; r and rinv hold them, and
 * reset_local sets r back to -1.
 */
static int walk_compare(const int *fp, int root, int *r, int *rinv, int *queue,
                        const int *best, int *k)
{
    int head = 0, tail = 0, used = 0, p = 0, cmp = 0, step = 1, next = 0;
    int h, x, v, lim;

    queue[tail++] = root;
    while (head < tail) {
        h = queue[head++];
        x = h;
        do {
            if (r[x] == -1) {
                r[x] = used;
                r[x ^ 1] = used + 1;
                rinv[used] = x;
                if (fp[x ^ 1] != -1) {
                    rinv[used + 1] = x ^ 1;
                    queue[tail++] = x ^ 1;
                } else {
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
        for (; p < lim; p++) {
            x = rinv[p];
            v = x == -1 ? -1 : r[fp[x]];
#ifdef ISO_DEBUG
            ISO_CHECK(v != -1 || x == -1);
#endif
            if (v != best[p]) {
                if (v > best[p]) {
                    *k = used;
                    return 1;
                }
                cmp = -1;
                break;
            }
        }
    }
#ifdef ISO_DEBUG
    ISO_CHECK(cmp != 0 || p == used);
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
static int tail_compare(int n, int nd, const int *data, int nc, const int *classes,
                        const int *rinv, int k, const int *best, int *cur, int *names)
{
    int i, j, h, v, pos = k, next, cmp = 0;
    const int *a;

    for (j = 0; j < nd; j++) {
        a = data + (long) j * n;
        for (i = 0; i < k; i++) {
            h = rinv[i];
            if (h == -1)
                continue;
            v = a[h];
            if (cmp == 0 && v != best[pos]) {
                if (v > best[pos])
                    return 1;
                cmp = -1;
            }
            cur[pos++] = v;
        }
    }
    for (j = 0; j < nc; j++) {
        a = classes + (long) j * n;
        next = 0;
        for (i = 0; i < k; i++) {
            h = rinv[i];
            if (h == -1)
                continue;
            if (names[a[h]] == -1)
                names[a[h]] = next++;
            v = names[a[h]];
            if (cmp == 0 && v != best[pos]) {
                if (v > best[pos]) {
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
        for (i = 0; i < k; i++) {
            h = rinv[i];
            if (h != -1)
                names[a[h]] = -1;
        }
    }
    return cmp;
}

static void reset_local(int k, const int *lrinv, int *lr);

#ifdef ISO_DEBUG
/*
 * Check the outcome cmp of the root against the key computed in full by
 * iso_relabelling_from and iso_key (len ints): it compares with best as
 * cmp says, and equals cur when cmp < 0. scratch holds 3 n +
 * ISO_KEY_SIZE(n, nd, nc) ints.
 */
static void debug_root(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                       int root, int cmp, const int *best, const int *cur, int len, int *names,
                       int *scratch)
{
    int *r2 = scratch, *rinv2 = scratch + n, *queue2 = scratch + 2 * n, *key2 = scratch + 3 * n;
    int k2, len2, i, ref = 0;

    k2 = iso_relabelling_from(n, fp, root, r2, rinv2, queue2);
    len2 = iso_key(n, fp, nd, data, nc, classes, r2, rinv2, k2, key2, names);
    ISO_CHECK(len2 == len);
    for (i = 0; i < len; i++) {
        if (key2[i] != best[i]) {
            ref = key2[i] < best[i] ? -1 : 1;
            break;
        }
    }
    ISO_CHECK(ref == cmp);
    if (cmp < 0)
        for (i = 0; i < len; i++)
            ISO_CHECK(cur[i] == key2[i]);
}

/*
 * Check that h -> rinv[rbest[h]] on the active half-edges, from the
 * relabellings of two roots of equal keys, is an automorphism of the
 * decorated map: a permutation of the active half-edges that commutes with
 * fp and with the twin, keeps the data arrays and the class arrays up to a
 * bijection of their values. scratch holds 4 n ints.
 */
static void debug_automorphism(int n, const int *fp, int nd, const int *data, int nc,
                               const int *classes, const int *rbest, const int *rinv, int *scratch)
{
    int *aut = scratch, *seen = scratch + n, *fwd = scratch + 2 * n, *bwd = scratch + 3 * n;
    int h, j, v, w;
    const int *a;

    for (h = 0; h < n; h++)
        seen[h] = 0;
    for (h = 0; h < n; h++) {
        if (fp[h] == -1)
            continue;
        ISO_CHECK(0 <= rbest[h] && rbest[h] < n);
        aut[h] = rinv[rbest[h]];
        ISO_CHECK(0 <= aut[h] && aut[h] < n && fp[aut[h]] != -1 && !seen[aut[h]]);
        seen[aut[h]] = 1;
    }
    for (h = 0; h < n; h++) {
        if (fp[h] == -1)
            continue;
        ISO_CHECK(aut[fp[h]] == fp[aut[h]]);
        if (fp[h ^ 1] != -1)
            ISO_CHECK(aut[h ^ 1] == (aut[h] ^ 1));
        else
            ISO_CHECK(fp[aut[h] ^ 1] == -1);
    }
    for (j = 0; j < nd; j++) {
        a = data + (long) j * n;
        for (h = 0; h < n; h++)
            if (fp[h] != -1)
                ISO_CHECK(a[aut[h]] == a[h]);
    }
    for (j = 0; j < nc; j++) {
        a = classes + (long) j * n;
        for (v = 0; v < n; v++) {
            fwd[v] = -1;
            bwd[v] = -1;
        }
        for (h = 0; h < n; h++) {
            if (fp[h] == -1)
                continue;
            v = a[h];
            w = a[aut[h]];
            ISO_CHECK((fwd[v] == -1 && bwd[w] == -1) || (fwd[v] == w && bwd[w] == v));
            fwd[v] = w;
            bwd[w] = v;
        }
    }
}
#endif

/*
 * The representative of the orbit of h in the union-find parent, its
 * smallest half-edge, with path halving.
 */
static inline int orbit_find(int *parent, int h)
{
    while (parent[h] != h) {
        parent[h] = parent[parent[h]];
        h = parent[h];
    }
    return h;
}

/* Merge the orbits of a and b, the larger representative under the smaller. */
static inline void orbit_union(int *parent, int a, int b)
{
    a = orbit_find(parent, a);
    b = orbit_find(parent, b);
    if (a < b)
        parent[b] = a;
    else if (b < a)
        parent[a] = b;
}

/*
 * The state of iso_canonical while it tries the roots after the first one:
 * the buffers, the best key so far and the candidate key (swapped when a
 * root beats the best one), the length of the keys, the best root and its
 * relabelling, and the parents of the union-find of the orbits of the
 * automorphisms found, set up at the first tie only (merged nonzero): until
 * then every orbit is a single half-edge and no root is skipped, so that a
 * map whose roots are all rejected pays nothing for the union-find.
 */
struct canonical_state {
    int n, nd, nc, ne, len, best_root, merged;
    const int *fp, *data, *classes;
    int *r, *rinv, *queue, *names, *best, *cur, *rbest, *parent;
#ifdef ISO_DEBUG
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
static inline void canonical_root(struct canonical_state *s, int root)
{
    int n = s->n, k, cmp, h;
    int *r = s->r, *rinv = s->rinv, *tmp;

    if (s->merged && orbit_find(s->parent, root) < root) {
#ifdef ISO_DEBUG
        cmp = orbit_find(s->parent, root) == orbit_find(s->parent, s->best_root) ? 0 : 1;
        debug_root(n, s->fp, s->nd, s->data, s->nc, s->classes, root, cmp, s->best, s->cur,
                   s->len, s->names, s->scratch);
#endif
        return;
    }

    cmp = walk_compare(s->fp, root, r, rinv, s->queue, s->best, &k);
    if (cmp < 0) {
        iso_key(n, s->fp, s->nd, s->data, s->nc, s->classes, r, rinv, k, s->cur, s->names);
    } else if (cmp == 0) {
        cmp = tail_compare(n, s->nd, s->data, s->nc, s->classes, rinv, k, s->best, s->cur, s->names);
        if (cmp < 0)
            memcpy(s->cur, s->best, (size_t) k * sizeof(int));
    }
#ifdef ISO_DEBUG
    ISO_CHECK(cmp == 1 || k == 2 * s->ne);
    debug_root(n, s->fp, s->nd, s->data, s->nc, s->classes, root, cmp, s->best, s->cur, s->len,
               s->names, s->scratch);
#endif
    if (cmp < 0) {
        tmp = s->best;
        s->best = s->cur;
        s->cur = tmp;
        s->best_root = root;
        memcpy(s->rbest, r, (size_t) n * sizeof(int));
    } else if (cmp == 0) {
        /* the automorphism h -> rinv[rbest[h]] sends the best root to root */
        if (!s->merged) {
            for (h = 0; h < n; h++)
                s->parent[h] = h;
            s->merged = 1;
        }
#ifdef ISO_DEBUG
        ISO_CHECK(orbit_find(s->parent, root) != orbit_find(s->parent, s->best_root));
        debug_automorphism(n, s->fp, s->nd, s->data, s->nc, s->classes, s->rbest, rinv,
                           s->scratch);
#endif
        for (h = 0; h < n; h++)
            if (s->fp[h] != -1)
                orbit_union(s->parent, h, rinv[s->rbest[h]]);
    }
    /* past n / 2 labels, a sequential clear (all bytes 0xff, the int
     * -1) is cheaper than the scattered writes of reset_local: 6 to
     * 13 % of the time of a regular map, where every root walks all */
    if (2 * k >= n)
        memset(r, -1, (size_t) n * sizeof(int));
    else
        reset_local(k, rinv, r);
#ifdef ISO_DEBUG
    {
        int i;
        for (i = 0; i < n; i++)
            ISO_CHECK(r[i] == -1);
        debug_names(n, s->names);
    }
#endif
}

int iso_canonical(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                  iso_root_rule rule, int *key, int *relabelling, int *roots, int *nroots,
                  int *work)
{
    int *r = work;
    int *rinv = work + n;
    int *queue = work + 2 * n;
    int *names = work + 3 * n;
    int *cur = work + 4 * n;
    int *order = cur + ISO_KEY_SIZE(n, nd, nc);
    int *cls = order + n;
    int *start = cls + n;
    int *rbest = start + n + 1;
    int *parent = rbest + n;
    struct canonical_state s;
    int all = (rule & ~ISO_ROOTS_DATA) == ISO_ROOTS_ALL;
    int ne, first, root, k, i, c, m, rep, count, lo = 0, hi = 0;

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

    if (!all) {
        /* r, rinv, queue and names (4 n) are free until the first walk, cur
         * (the ranks of the data values) until the first key; the sizes of
         * the classes go into r */
        m = root_classes(n, fp, rule, nd, data, order, cls, start, work, cur);
        for (c = 0; c < m; c++)
            r[c] = start[c + 1] - start[c];
        c = choose_class(rule, m, r);
        lo = start[c];
        hi = start[c + 1];
        first = order[lo];
    }

    for (i = 0; i < n; i++)
        names[i] = -1;

    /* the first root, walked in full: connectedness does not depend on the
     * root, so only this walk can fail */
    k = iso_relabelling_from(n, fp, first, r, rinv, queue);
    if (k != 2 * ne)
        return ISO_DISCONNECTED;
    s.n = n;
    s.nd = nd;
    s.nc = nc;
    s.ne = ne;
    s.fp = fp;
    s.data = data;
    s.classes = classes;
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
    s.len = iso_key(n, fp, nd, data, nc, classes, r, rinv, k, s.best, names);
    memcpy(rbest, r, (size_t) n * sizeof(int));
    reset_local(k, rinv, r);
#ifdef ISO_DEBUG
    s.scratch = malloc(((size_t) 3 * n + (size_t) ISO_KEY_SIZE(n, nd, nc)) * sizeof(int));
    ISO_CHECK(s.scratch != NULL);
#endif

    /* the other roots: r is -1 on entry of each, and only the labels a root
     * gave are reset after it. The loop is chosen once per call, so that
     * ISO_ROOTS_ALL tests no class per root. */
    if (all) {
        for (root = first + 1; root < n; root++)
            if (fp[root] != -1)
                canonical_root(&s, root);
    } else {
        for (i = lo + 1; i < hi; i++)
            canonical_root(&s, order[i]);
    }
#ifdef ISO_DEBUG
    free(s.scratch);
#endif
    if (s.best != key)
        memcpy(key, s.best, (size_t) s.len * sizeof(int));
    if (relabelling != NULL)
        memcpy(relabelling, rbest, (size_t) n * sizeof(int));

    /* the roots of minimal key: the orbit of the best root, which is its
     * smallest half-edge (a smaller one would have been tried first) */
    if (!s.merged) {
        if (roots != NULL)
            roots[0] = s.best_root;
        count = 1;
    } else {
        rep = orbit_find(parent, s.best_root);
#ifdef ISO_DEBUG
        ISO_CHECK(rep == s.best_root);
#endif
        count = 0;
        if (all) {
            for (root = rep; root < n; root++) {
                if (fp[root] != -1 && orbit_find(parent, root) == rep) {
                    if (roots != NULL)
                        roots[count] = root;
                    count++;
                }
            }
        } else {
            /* the orbit lies in the class, increasing in order[lo..hi) */
            for (i = lo; i < hi; i++) {
                root = order[i];
                if (orbit_find(parent, root) == rep) {
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

/*
 * The breadth-first walk of iso_relabelling_from restricted to the
 * component of root: lr must be -1 on the half-edges of that component and
 * is left set there (reset_local undoes it), so that a walk costs the size
 * of the component and not n.
 */
static int local_walk(const int *fp, int root, int *lr, int *lrinv, int *queue)
{
    int head = 0, tail = 0, k = 0;
    int h, x;

    queue[tail++] = root;
    while (head < tail) {
        h = queue[head++];
        x = h;
        do {
            if (lr[x] == -1) {
                lr[x] = k;
                lr[x ^ 1] = k + 1;
                lrinv[k] = x;
                if (fp[x ^ 1] != -1) {
                    lrinv[k + 1] = x ^ 1;
                    queue[tail++] = x ^ 1;
                } else {
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
static void reset_local(int k, const int *lrinv, int *lr)
{
    int i, h;
    /* an even label is always on an active half-edge, the one met in a face */
    for (i = 0; i < k; i += 2) {
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
static int component_block(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                           const int *lr, const int *lrinv, int k, int *names, int *next,
                           int commit, int *block)
{
    int i, j, h, pos = 0, count;
    const int *a;
    int *nm;

    block[pos++] = k;
    for (i = 0; i < k; i++) {
        h = lrinv[i];
        block[pos++] = h == -1 ? -1 : lr[fp[h]];
    }
    for (j = 0; j < nd; j++) {
        a = data + (long) j * n;
        for (i = 0; i < k; i++) {
            h = lrinv[i];
            if (h != -1)
                block[pos++] = a[h];
        }
    }
    for (j = 0; j < nc; j++) {
        a = classes + (long) j * n;
        nm = names + (long) j * n;
        count = next[j];
        for (i = 0; i < k; i++) {
            h = lrinv[i];
            if (h != -1) {
#ifdef ISO_DEBUG
                ISO_CHECK(0 <= a[h] && a[h] < n);
#endif
                if (nm[a[h]] == -1)
                    nm[a[h]] = count++;
                block[pos++] = nm[a[h]];
            }
        }
        if (commit) {
            next[j] = count;
        } else {
            /* the values numbered by this component are those numbered
             * from next[j] on */
            for (i = 0; i < k; i++) {
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
static int compare_sequences(const int *a, int la, const int *b, int lb)
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
static int exceeds_prefix(const int *a, int la, const int *b, int lb)
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
static inline int root_block(int n, const int *fp, int nd, const int *data, int nc,
                             const int *classes, int root, int *lr, int *lrinv, int *queue,
                             int *names, int *next, int *cand)
{
    int k = local_walk(fp, root, lr, lrinv, queue);
    int len = component_block(n, fp, nd, data, nc, classes, lr, lrinv, k, names, next, 0, cand);
    reset_local(k, lrinv, lr);
    return len;
}

/*
 * The levels of the search store the position of their root in the scan of
 * the roots: the half-edge itself under ISO_ROOTS_ALL, its index in order
 * (see root_classes) under the other rules. The scans are written twice, once
 * per kind of position, so that ISO_ROOTS_ALL tests no class per root.
 */
int iso_search(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
               iso_root_rule rule, int all_leaves, int *sequence, int *length, int *relabelling,
               iso_leaf_callback callback, void *arg, int *work)
{
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
    int *lvl_next = lvl_len + levels;
    int *next = lvl_next + (long) nc * levels;
    int *seq = next + nc;
    int *cand = seq + ISO_SEQUENCE_SIZE(n, nd, nc);
    int *order = cand + ISO_KEY_SIZE(n, nd, nc) + 1;
    int *cls = order + n;
    int *start = cls + n;
    int *cnt = start + n + 1;
    int follow_ties = nc > 0 || all_leaves;
    int all = (rule & ~ISO_ROOTS_DATA) == ISO_ROOTS_ALL;
    int ne, depth = 0, used = 0, seqlen = 0, bestlen = -1, m = 0;
    int root, best_root, blen, len, k, cmp, ret, i, j, h, c, found;

    ne = iso_num_edges(n, fp);
    if (!all) {
        /* r, rinv, lr and lrinv (4 n) are free until they are cleared below,
         * cand (the ranks of the data values) until the first block;
         * cnt[c] is the number of half-edges of the class c not labelled */
        m = root_classes(n, fp, rule, nd, data, order, cls, start, work, cand);
        for (c = 0; c < m; c++)
            cnt[c] = start[c + 1] - start[c];
    }
    for (h = 0; h < n; h++) {
        r[h] = -1;
        lr[h] = -1;
    }
    for (i = 0; i < nc * n; i++)
        names[i] = -1;
    for (j = 0; j < nc; j++)
        next[j] = 0;

    for (;;) {
        /* enter the state: depth levels labelled, used labels, the blocks
         * in seq[0:seqlen] */
        int cut = bestlen >= 0 && exceeds_prefix(seq, seqlen, sequence, bestlen);

        if (!cut && used == 2 * ne) {
            cmp = bestlen < 0 ? -1 : compare_sequences(seq, seqlen, sequence, bestlen);
            if (cmp < 0) {
                memcpy(sequence, seq, (size_t) seqlen * sizeof(int));
                bestlen = seqlen;
                if (relabelling != NULL)
                    memcpy(relabelling, r, (size_t) n * sizeof(int));
                if (callback != NULL && (ret = callback(r, 1, arg)) != 0)
                    return ret;
            } else if (cmp == 0 && all_leaves && callback != NULL) {
                if ((ret = callback(r, 0, arg)) != 0)
                    return ret;
            }
            cut = 1;
        }

        if (!cut) {
            /* the minimal block over the roots not labelled yet, written
             * in place in seq, and its first root */
            best_root = -1;
            blen = 0;
            if (all) {
                for (root = 0; root < n; root++) {
                    if (fp[root] == -1 || r[root] != -1)
                        continue;
                    len = root_block(n, fp, nd, data, nc, classes, root, lr, lrinv, queue,
                                     names, next, cand);
                    if (best_root == -1 || compare_sequences(cand, len, seq + seqlen, blen) < 0) {
                        memcpy(seq + seqlen, cand, (size_t) len * sizeof(int));
                        blen = len;
                        best_root = root;
                    }
                }
            } else {
                c = choose_class(rule, m, cnt);
                for (i = start[c]; i < start[c + 1]; i++) {
                    root = order[i];
                    if (r[root] != -1)
                        continue;
                    len = root_block(n, fp, nd, data, nc, classes, root, lr, lrinv, queue,
                                     names, next, cand);
                    if (best_root == -1 || compare_sequences(cand, len, seq + seqlen, blen) < 0) {
                        memcpy(seq + seqlen, cand, (size_t) len * sizeof(int));
                        blen = len;
                        best_root = i;
                    }
                }
            }
#ifdef ISO_DEBUG
            ISO_CHECK(best_root != -1 && depth < levels);
#endif
            lvl_root[depth] = best_root;
            lvl_used[depth] = used;
            lvl_pos[depth] = seqlen;
            lvl_len[depth] = blen;
            for (j = 0; j < nc; j++)
                lvl_next[(long) depth * nc + j] = next[j];
        } else {
            /* backtrack to the deepest level with a further tied root */
            found = 0;
            while (!found && depth > 0) {
                depth--;
                if (!all)
                    for (i = lvl_used[depth]; i < used; i++)
                        if (rinv[i] != -1)
                            cnt[cls[rinv[i]]]++;
                for (i = lvl_used[depth]; i < used; i++) {
                    h = rinv[i];
                    if (h == -1)
                        continue;
                    for (j = 0; j < nc; j++) {
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
                if (!follow_ties)
                    continue;
                /* the tied roots after the current one: their block is
                 * the one of the level, already in seq; under a rule, they
                 * are in the class of the root of the level, the one chosen
                 * at the level since the state is back to it */
                if (all) {
                    for (root = lvl_root[depth] + 1; root < n; root++) {
                        if (fp[root] == -1 || r[root] != -1)
                            continue;
                        len = root_block(n, fp, nd, data, nc, classes, root, lr, lrinv, queue,
                                         names, next, cand);
                        if (len == lvl_len[depth] && compare_sequences(cand, len, seq + seqlen, len) == 0) {
                            lvl_root[depth] = root;
                            found = 1;
                            break;
                        }
                    }
                } else {
                    c = cls[order[lvl_root[depth]]];
                    for (i = lvl_root[depth] + 1; i < start[c + 1]; i++) {
                        root = order[i];
                        if (r[root] != -1)
                            continue;
                        len = root_block(n, fp, nd, data, nc, classes, root, lr, lrinv, queue,
                                         names, next, cand);
                        if (len == lvl_len[depth] && compare_sequences(cand, len, seq + seqlen, len) == 0) {
                            lvl_root[depth] = i;
                            found = 1;
                            break;
                        }
                    }
                }
            }
            if (!found)
                break;
        }

        /* label the component of the root of the level depth */
        root = all ? lvl_root[depth] : order[lvl_root[depth]];
        k = local_walk(fp, root, lr, lrinv, queue);
        for (i = 0; i < k; i += 2) {
            h = lrinv[i];
            r[h] = used + i;
            r[h ^ 1] = used + i + 1;
            rinv[used + i] = h;
            rinv[used + i + 1] = lrinv[i + 1];
        }
        component_block(n, fp, nd, data, nc, classes, lr, lrinv, k, names, next, 1, cand);
        if (!all)
            for (i = 0; i < k; i++)
                if (lrinv[i] != -1)
                    cnt[cls[lrinv[i]]]--;
        reset_local(k, lrinv, lr);
        used += k;
        seqlen += lvl_len[depth];
        depth++;
    }

#ifdef ISO_DEBUG
    /* every level undone */
    ISO_CHECK(bestlen >= 0 && used == 0 && seqlen == 0);
    for (h = 0; h < n; h++)
        ISO_CHECK(r[h] == -1 && lr[h] == -1);
    for (i = 0; i < nc * n; i++)
        ISO_CHECK(names[i] == -1);
    for (c = 0; c < m; c++)
        ISO_CHECK(cnt[c] == start[c + 1] - start[c]);
#endif
    *length = bestlen;
    return 0;
}
