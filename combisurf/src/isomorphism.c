/*
 * Canonical keys, canonical labellings and automorphisms of maps. See
 * isomorphism.h for the conventions and the sizes of the buffers.
 *
 * On a connected map (iso_canonical), every active half-edge is tried as a
 * root: the relabelling from it and its key are computed in full and
 * compared with the best key so far, which costs O(n) per root and O(n^2)
 * in all for n half-edges.
 *
 * On any map (iso_search), the search goes component after component. Its
 * memory does not grow with the depth beyond a few ints per level: the
 * state of a level is its root, its first label, the position and length
 * of its block in the current sequence and the counts of class values
 * before it. A level is undone by walking the labels it gave, and its next
 * tied root is found again by scanning the roots after the current one for
 * a block equal to the one of the level, which is still in the sequence.
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

int iso_search(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
               int all_leaves, int *sequence, int *length, int *relabelling,
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
    int follow_ties = nc > 0 || all_leaves;
    int ne, depth = 0, used = 0, seqlen = 0, bestlen = -1;
    int root, best_root, blen, len, k, cmp, ret, i, j, h, found;

    ne = iso_num_edges(n, fp);
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
            for (root = 0; root < n; root++) {
                if (fp[root] == -1 || r[root] != -1)
                    continue;
                k = local_walk(fp, root, lr, lrinv, queue);
                len = component_block(n, fp, nd, data, nc, classes, lr, lrinv, k, names, next, 0, cand);
                reset_local(k, lrinv, lr);
                if (best_root == -1 || compare_sequences(cand, len, seq + seqlen, blen) < 0) {
                    memcpy(seq + seqlen, cand, (size_t) len * sizeof(int));
                    blen = len;
                    best_root = root;
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
                 * the one of the level, already in seq */
                for (root = lvl_root[depth] + 1; root < n; root++) {
                    if (fp[root] == -1 || r[root] != -1)
                        continue;
                    k = local_walk(fp, root, lr, lrinv, queue);
                    len = component_block(n, fp, nd, data, nc, classes, lr, lrinv, k, names, next, 0, cand);
                    reset_local(k, lrinv, lr);
                    if (len == lvl_len[depth] && compare_sequences(cand, len, seq + seqlen, len) == 0) {
                        lvl_root[depth] = root;
                        found = 1;
                        break;
                    }
                }
            }
            if (!found)
                break;
        }

        /* label the component of the root of the level depth */
        root = lvl_root[depth];
        k = local_walk(fp, root, lr, lrinv, queue);
        for (i = 0; i < k; i += 2) {
            h = lrinv[i];
            r[h] = used + i;
            r[h ^ 1] = used + i + 1;
            rinv[used + i] = h;
            rinv[used + i + 1] = lrinv[i + 1];
        }
        component_block(n, fp, nd, data, nc, classes, lr, lrinv, k, names, next, 1, cand);
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
#endif
    *length = bestlen;
    return 0;
}
