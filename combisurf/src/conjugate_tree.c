/*
 * Conjugate trees: suffix trees of the conjugates of a finite set of words.
 * See conjugate_tree.h for the interface.
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
#include "conjugate_tree.h"

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*
 * The layouts of the children (see conjugate_tree.h).
 *
 * A dense table costs 4 * alphabet bytes per node, where the sibling lists
 * cost 8 and the rest of a node 24, and every new node clears its row and
 * every listing of the leaves scans it. A single word of length 100000 over
 * 128 letters has 225000 nodes, that is 110 MB of dense table against 1.7 MB
 * of sibling lists: past a small alphabet, the table is only affordable for
 * the few nodes that have many children.
 *
 * Those are the nodes near the root. On two random cyclically reduced words
 * over n = 1024 letters and their inverses, the root has 1024 children, the
 * nodes of depth 1 have 226 on average at 2^16 letters of each word (646 at
 * 2^18), and the deeper ones about 2. Walking the siblings then takes 510
 * steps per lookup at the root and 101 at depth 1, all but 0.01 % of the steps,
 * each a cache miss into the word buffer for the first letter. The rows
 * layout reads that letter from flet instead, and gives the nodes with at
 * least CT_PROMOTE children a dense row. Times to build the tree of the
 * words and their inverses, rows against sparse:
 *
 *     n        letters per word      8      1000     2^16     2^18
 *     1024     sparse           1.61 us   8.06 ms   1.19 s   12.2 s
 *              rows             1.25 us    686 us   31.5 ms   708 ms
 *     128      sparse           1.61 us   1.29 ms   225 ms   1.64 s
 *              rows             1.16 us    217 us   35.4 ms   186 ms
 *
 * With words drawn with Zipf's law of exponent 1, the nodes of depth 2 have
 * up to 958 children (mean 7.5) at 2^18 letters, and rows are 25 to 27 times
 * faster than the sparse layout at n = 1024 and 2^16 or 2^18 letters, 13
 * times at 1000 letters. On words with few factors (Fibonacci, the ruler
 * sequence, random Fibonacci), whose nodes have 2 to 16 children, rows take
 * 0.8 to 1.1 times the sparse time.
 *
 * Giving a row at 8, 16 or 32 children makes no more than 7 % of difference
 * on random words at n = 1024; 16 is the fastest on the Zipf words at 2^18
 * letters, where 8 takes 222 MB against 155 MB for 16 (71 MB for sparse
 * lists); at n = 128 and 2^18 letters, 32 and 64 are 1.6 and 3 times slower
 * than 16. Hence CT_PROMOTE.
 *
 * A hash table (node, letter) -> child with linear probing was at most 27 %
 * slower than the rows on random and Zipf words, but 1.3 to 2.8 times slower
 * than the sparse lists on the words with few factors, and took more memory
 * (130 MB against 96 MB at n = 1024 and 2^18 letters).
 *
 * Up to 32 letters the default is the dense table: it costs at most 128
 * bytes per node, and at 32 letters the rows are 3 to 50 % slower than it
 * (1000 letters: 147 us dense, 221 us rows).
 *
 * Large alphabets. A row costs O(alphabet) time and memory, and with rows
 * from CT_PROMOTE children on there may be one per CT_PROMOTE leaves; with
 * n = 16384 letters, two random words of 2^16 letters and their inverses
 * took 1.1 GB in 8739 rows. A node therefore gets its row from
 * max(CT_PROMOTE, n / CT_ROW_RATIO) children on, so that the rows cost at
 * most CT_ROW_RATIO ints per child; below that, from crit_min children on,
 * the children of a node are also held in a crit-bit tree, a binary trie on
 * the bits of their first letters with one branching crit node per child
 * but one: a lookup reads at most one crit node per bit of the letters. Up
 * to CT_PROMOTE * CT_ROW_RATIO = 1024 letters the rows come at CT_PROMOTE
 * children as before, and crit_min is CT_PROMOTE, so that no crit-bit tree
 * is built: one built at 8 children, then dropped at 16 for the row, made
 * the build 7 to 11 % slower at n = 128 and 1024 (Markov and Zipf words of
 * 2^16 letters: 36.3 against 33.1 ms, 43.7 against 39.3 ms). Past 1024
 * letters, a crit-bit tree from 8 children on beats 16 and 32 (n = 16384,
 * 2^16 letters per word, nodes of depth 1 with 16 children on average: 30.5,
 * 36.4 and 37.7 ms). CT_ROW_RATIO = 32, 64 or 128 made no difference
 * there (29.6, 30.3, 30.7 ms), nor rows at 16 or 64 children at n = 4096
 * and 3000 letters per word (709 and 711 us).
 *
 * The worst case, with b the number of bits of the letters (log2 of the
 * alphabet, rounded up, or 31 without alphabet) and d the number of children
 * of the node at hand, in the rows layout:
 * - a lookup costs O(1) with a row, O(b) with a crit-bit tree, and fewer
 *   than crit_min steps along the siblings otherwise;
 * - adding a child costs the same, plus: counting the siblings, fewer than
 *   crit_min; building the crit-bit tree, O(d b) once per node at d =
 *   crit_min children; giving the row, O(n + d) once per node at d >=
 *   n / CT_ROW_RATIO children, that is O(CT_ROW_RATIO) per child;
 * - replacing a child by the node that splits its edge costs the same as a
 *   lookup: the row and the crit-bit tree are walked by letter, the siblings
 *   of a node with a row or a crit-bit tree are linked both ways (psib), and
 *   the other nodes have fewer than crit_min children to walk.
 * Ukkonen's construction makes O(l) of these operations on words of l
 * letters in all, hence O(l b) for the build. These bounds assume that the
 * allocations succeed: a node whose row or crit-bit tree cannot be allocated
 * walks its siblings. The leaves are listed in O(n + l log n): the children
 * of a node are sorted by an insertion sort up to 16 of them, by placing
 * them in an array of n slots from n / 8 of them on (O(n), at most 8 per
 * child), and by a merge sort in between. What remains in O(n) per call is
 * the reading and checking of order and pivot in ct_sorted_leaves, with its
 * arrays of n entries; the dense layout costs O(n) per node by design, and
 * the sparse layout O(d) per lookup.
 *
 * Measured, before against after (best of 5, two random words and their
 * inverses, so 4 l letters in the tree; build and listing of the leaves):
 *
 *     n        l          before               after
 *     1024     3000     1.58 ms  335 us     838 us  193 us
 *     4096     3000     8.17 ms  2.14 ms    714 us  155 us
 *     16384    3000     16.2 ms  10.4 ms    600 us  164 us
 *     16384    2^16      320 ms  44.3 ms   30.3 ms  7.83 ms   (1.1 GB -> 31 MB)
 *     1024     2^16     30.7 ms  13.0 ms   19.1 ms  7.69 ms
 *
 * At 1024 letters, where no crit-bit tree is built, the gain of the build is
 * the replacement of a child, which walked the siblings before, up to 1024
 * steps at the root, and that of the listing the sort of the children of the
 * root. On the words whose nodes have 2 to 8 children (Fibonacci over 1024
 * letters, 2^16 letters per word) the build takes 17.4 ms against 15.6 ms:
 * psib costs 4 bytes more per node. The first version, with psib kept for
 * every node and separate arrays for the number of children and the root of
 * the crit-bit tree, took 20.7 ms.
 */

/* ------------------------------------------------------------------ */
/* allocation                                                          */
/* ------------------------------------------------------------------ */

/* Resize *p to n ints; *p is left alone on failure. */
static int resize_ints(int **p, int n)
{
    int *q;
    if ((size_t) n > SIZE_MAX / sizeof(int))
        return CT_ETOOLARGE;
    q = (int *) realloc(*p, (size_t) n * sizeof(int));
    if (q == NULL && n != 0)
        return CT_ENOMEM;
    *p = q;
    return CT_OK;
}

/* The capacity to grow from cap to hold at least need, doubling, at most max. */
static int grown(int cap, int need, int max, int start)
{
    int c = cap ? cap : start;
    if (c > max)
        c = max;
    while (c < need)
        c = c > max / 2 ? max : 2 * c;
    return c;
}

/* The largest number of nodes the node indices and the dense table allow. */
static int max_nodes(const ct_tree *T)
{
    return T->layout == CT_LAYOUT_DENSE ? INT_MAX / T->alphabet_size : INT_MAX;
}

/* ------------------------------------------------------------------ */
/* the rows of CT_LAYOUT_ROWS                                          */
/* ------------------------------------------------------------------ */

/*
 * A node gets its row when it gets its promote-th child. The number of such
 * nodes is only bounded by nstates / promote, which is far more than the
 * trees met in practice have (at n = 1024 letters and 2^19 letters of words,
 * 1025 nodes out of 1.1 million), so the rows are not reserved ahead:
 * promote_node allocates, and when it cannot the node keeps its sibling
 * list alone, which gives the same answers, and tries again at its next
 * child. An insertion therefore never fails because of the rows.
 */
static void crit_release(ct_tree *T, int p);

/* Give s a row; return 1 when done, 0 when it could not be allocated. */
static int promote_node(ct_tree *T, int s)
{
    size_t a = (size_t) T->alphabet_size;
    int r, t, c, *row;
    if (T->nrows == T->rows_capacity) {
        int cap = T->rows_capacity;
        int *p;
        if (cap > INT_MAX / 2)
            return 0;
        cap = cap ? 2 * cap : 4;
        if ((size_t) cap > SIZE_MAX / sizeof(int) / a)
            return 0;
        p = (int *) realloc(T->rows, (size_t) cap * a * sizeof(int));
        if (p == NULL)
            return 0;
        T->rows = p;
        T->rows_capacity = cap;
    }
    r = T->nrows++;
    row = T->rows + (size_t) r * a;
    for (c = 0; c < T->alphabet_size; c++)
        row[c] = -1;
    /* the children of a node with a row or a crit-bit tree are linked both
     * ways (see replace_child) */
    for (c = -1, t = T->fchild[s]; t != -1; c = t, t = T->nsib[t]) {
        row[T->flet[t]] = t;
        T->psib[t] = c;
    }
    /* the row answers every lookup from now on */
    crit_release(T, T->look[s]);
    T->look[s] = r;
    return 1;
}

/* ------------------------------------------------------------------ */
/* the crit-bit trees of CT_LAYOUT_ROWS                                */
/* ------------------------------------------------------------------ */

/*
 * The crit nodes are allocated from one pool, crit nodes released by a node
 * that gets its row going to a free list. The pool keeps the invariant
 *
 *     crit_capacity - ncrit >= crit_promised - nstates,
 *
 * where crit_promised is the number of nodes that the last reservations
 * (ct_process, ct_reserve) promised to make room for. Inserting a child into
 * an existing crit-bit tree takes one crit node, and comes with a new node of
 * the tree (each call to add_child pairs with a call to add_node), so these
 * insertions never run out of crit nodes. Building the crit-bit tree of a node
 * with d children takes d crit nodes at once (its header and d - 1 branching
 * ones) that no reservation accounted for: crit_build grows the pool itself,
 * and when it cannot the node keeps its sibling list alone and tries again at
 * its next child.
 */

/* Grow the pool so that it holds extra crit nodes on top of the room that
 * the invariant above asks for with crit_promised = promised. */
static int crit_grow(ct_tree *T, int extra, int promised)
{
    int64_t need = (int64_t) T->ncrit + extra + promised - T->nstates;
    int cap, err;
    if (need <= T->crit_capacity)
        return CT_OK;
    if (need > INT_MAX / 3)
        return CT_ETOOLARGE;
    cap = grown(T->crit_capacity, (int) need, INT_MAX / 3, 8);
    if ((err = resize_ints(&T->crit, 3 * cap)))
        return err;
    T->crit_capacity = cap;
    return CT_OK;
}

static inline int high_bit(unsigned x)
{
#if defined(__GNUC__)
    return 31 - __builtin_clz(x);
#else
    int b = 0;
    while (x >>= 1)
        b++;
    return b;
#endif
}

/* The header of the crit-bit tree whose reference is look. */
static inline int *crit_header(const ct_tree *T, int look)
{
    return T->crit + 3 * (size_t) (-2 - look);
}

/* The child reached from the reference p by the bits of letter: the only one
 * that may start with letter. */
static inline int crit_find(const ct_tree *T, int p, int letter)
{
    while (p < -1) {
        const int *c = T->crit + 3 * (size_t) (-2 - p);
        p = c[1 + ((letter >> c[0]) & 1)];
    }
    return p;
}

/* A crit node out of the free list or the pool, whose room the caller
 * checked; -1 (and the tree marked broken) when there is none. */
static inline int crit_new(ct_tree *T)
{
    int j = T->crit_free;
    if (j != -1) {
        T->crit_free = T->crit[3 * (size_t) j + 1];
        return j;
    }
    if (T->ncrit == T->crit_capacity) {
        T->broken = CT_EINTERNAL;
        return -1;
    }
    return T->ncrit++;
}

/*
 * Insert the child t into the crit-bit tree of header h, which holds no child
 * with the letter of t yet. O(log alphabet): two descents of at most one crit
 * node per bit.
 */
static int crit_insert(ct_tree *T, int *h, int t)
{
    int key = T->flet[t];
    int p, b, j, side, *slot, *c;
    p = crit_find(T, h[1], key);
    b = high_bit((unsigned) (T->flet[p] ^ key));
    if ((j = crit_new(T)) < 0)
        return -1;
    /* the crit nodes above b agree with key on their bit: key goes where
     * the first crit node with a lower bit, or the child, was */
    slot = h + 1;
    while (*slot < -1) {
        c = T->crit + 3 * (size_t) (-2 - *slot);
        if (c[0] < b)
            break;
        slot = c + 1 + ((key >> c[0]) & 1);
    }
    c = T->crit + 3 * (size_t) j;
    side = (key >> b) & 1;
    c[0] = b;
    c[1 + side] = t;
    c[2 - side] = *slot;
    *slot = -2 - j;
    return 0;
}

/*
 * Give s, which has neither row nor crit-bit tree, a crit-bit tree over its d
 * children, if the pool can make room for it. O(d log alphabet). Return -1
 * when the tree is broken, else 0.
 */
static int crit_build(ct_tree *T, int s, int d)
{
    int t = T->fchild[s];
    int j, prev, *h;
    if (crit_grow(T, d, T->crit_promised))
        return 0;
    if ((j = crit_new(T)) < 0)
        return -1;
    h = T->crit + 3 * (size_t) j;
    h[0] = -2;
    h[1] = t;
    h[2] = d;
    /* the children of a node with a row or a crit-bit tree are linked both
     * ways (see replace_child) */
    T->psib[t] = -1;
    for (prev = t, t = T->nsib[t]; t != -1; prev = t, t = T->nsib[t]) {
        T->psib[t] = prev;
        if (crit_insert(T, h, t))
            return -1;
    }
    T->look[s] = -2 - j;
    return 0;
}

/* Put the crit nodes under the reference p (a header or a branching crit
 * node) on the free list. */
static void crit_release(ct_tree *T, int p)
{
    int *c;
    int j;
    if (p >= -1)
        return;
    j = -2 - p;
    c = T->crit + 3 * (size_t) j;
    crit_release(T, c[1]);
    if (c[0] != -2)
        crit_release(T, c[2]);
    c[0] = -1;
    c[1] = T->crit_free;
    c[2] = -1;
    T->crit_free = j;
}

/* ------------------------------------------------------------------ */
/* nodes                                                               */
/* ------------------------------------------------------------------ */

/*
 * Make room for need nodes. A failure leaves the tree as it was: arrays that
 * were already resized are merely larger than the capacity says.
 */
static int reserve_nodes(ct_tree *T, int need)
{
    int cap, err;
    if (need <= T->capacity)
        return CT_OK;
    if (need > max_nodes(T))
        return CT_ETOOLARGE;
    cap = grown(T->capacity, need, max_nodes(T), 8);

    if ((err = resize_ints(&T->dep, cap)) ||
        (err = resize_ints(&T->sl, cap)) ||
        (err = resize_ints(&T->parent, cap)) ||
        (err = resize_ints(&T->tword, cap)) ||
        (err = resize_ints(&T->tstart, cap)) ||
        (err = resize_ints(&T->tend, cap)))
        return err;
    if (T->layout == CT_LAYOUT_DENSE) {
        /* cap * alphabet_size <= INT_MAX by max_nodes */
        if ((err = resize_ints(&T->trans, cap * T->alphabet_size)))
            return err;
    } else {
        if ((err = resize_ints(&T->fchild, cap)) ||
            (err = resize_ints(&T->nsib, cap)))
            return err;
    }
    if (T->layout == CT_LAYOUT_ROWS &&
        ((err = resize_ints(&T->flet, cap)) || (err = resize_ints(&T->look, cap)) ||
         (err = resize_ints(&T->psib, cap))))
        return err;
    T->capacity = cap;
    return CT_OK;
}

/*
 * Make room for num_words more words of num_letters letters in total, and
 * for num_nodes more nodes. Called before any change to the tree, so that an
 * insertion never fails half way: a word of length l creates at most l
 * leaves and at most one internal node per leaf.
 */
static int reserve(ct_tree *T, int num_words, int num_letters, int num_nodes)
{
    int err, cap;
    if (num_words > INT_MAX - T->nwords ||
        num_letters > INT_MAX - T->wbuf_size ||
        num_nodes > INT_MAX - T->nstates)
        return CT_ETOOLARGE;
    if (T->nwords + num_words > T->words_capacity) {
        cap = grown(T->words_capacity, T->nwords + num_words, INT_MAX, 4);
        if ((err = resize_ints(&T->wstart, cap)) ||
            (err = resize_ints(&T->wlen, cap)))
            return err;
        T->words_capacity = cap;
    }
    if (T->wbuf_size + num_letters > T->wbuf_capacity) {
        cap = grown(T->wbuf_capacity, T->wbuf_size + num_letters, INT_MAX, 4);
        if ((err = resize_ints(&T->wbuf, cap)))
            return err;
        T->wbuf_capacity = cap;
    }
    if ((err = reserve_nodes(T, T->nstates + num_nodes)))
        return err;
    if (T->layout == CT_LAYOUT_ROWS) {
        int promised = T->nstates + num_nodes;
        if (promised < T->crit_promised)
            promised = T->crit_promised;
        if ((err = crit_grow(T, 0, promised)))
            return err;
        T->crit_promised = promised;
    }
    return CT_OK;
}

/*
 * Append a node with nothing but its children initialized and return its
 * index, or -2 when there is no room left (the reservation made by the caller
 * should prevent it; the tree is then marked broken).
 */
static inline int add_node(ct_tree *T)
{
    int n = T->nstates;
    int c;
    if (n == T->capacity) {
        T->broken = CT_EINTERNAL;
        return -2;
    }
    T->nstates = n + 1;
    /* NOTE: -2 is the code for uninitialized, as in ConjugateTreeNaive */
    T->dep[n] = -2;
    T->sl[n] = -2;
    T->parent[n] = -2;
    T->tword[n] = -2;
    T->tstart[n] = -2;
    T->tend[n] = -2;
    if (T->layout == CT_LAYOUT_DENSE) {
        int *row = T->trans + (size_t) n * (size_t) T->alphabet_size;
        for (c = 0; c < T->alphabet_size; c++)
            row[c] = -1;
    } else {
        T->fchild[n] = -1;
        T->nsib[n] = -1;
    }
    if (T->layout == CT_LAYOUT_ROWS) {
        T->flet[n] = -1;
        T->look[n] = -1;
        T->psib[n] = -1;
    }
    return n;
}

int ct_init(ct_tree *T, int alphabet, int reserve, int layout)
{
    int err;
    memset(T, 0, sizeof(ct_tree));
    T->max_letter = -1;
    if (alphabet < 0 || reserve < 0 || layout > CT_LAYOUT_ROWS ||
        (layout == CT_LAYOUT_DENSE && alphabet == 0))
        return CT_EINVALID;
    T->alphabet_size = alphabet;
    if (layout < 0)
        layout = 0 < alphabet && alphabet <= CT_DENSE_MAX_ALPHABET ? CT_LAYOUT_DENSE : CT_LAYOUT_ROWS;
    T->layout = layout;
    T->promote = CT_PROMOTE;
    if (alphabet / CT_ROW_RATIO >= CT_PROMOTE)
        T->promote = alphabet / CT_ROW_RATIO + (alphabet % CT_ROW_RATIO != 0);
    /* a node that gets its row at CT_PROMOTE children does without a
     * crit-bit tree before (see the comment at the top) */
    T->crit_min = alphabet && T->promote == CT_PROMOTE ? CT_PROMOTE : CT_CRIT_MIN;
    T->crit_free = -1;

    if ((err = reserve_nodes(T, reserve > 1 ? reserve : 1)))
        return err;
    add_node(T);
    /* NOTE: the root doubles as the transition of length one out of the
     * imaginary node -1 that canonize starts from, hence its (-4, -3) */
    T->dep[0] = 0;
    T->sl[0] = -1;
    T->parent[0] = -1;
    T->tword[0] = 0;
    T->tstart[0] = -4;
    T->tend[0] = -3;
    return CT_OK;
}

void ct_free(ct_tree *T)
{
    free(T->wbuf);
    free(T->wstart);
    free(T->wlen);
    free(T->dep);
    free(T->sl);
    free(T->parent);
    free(T->tword);
    free(T->tstart);
    free(T->tend);
    free(T->trans);
    free(T->fchild);
    free(T->nsib);
    free(T->flet);
    free(T->psib);
    free(T->look);
    free(T->rows);
    free(T->crit);
    memset(T, 0, sizeof(ct_tree));
}

/* ------------------------------------------------------------------ */
/* words                                                               */
/* ------------------------------------------------------------------ */

static inline int letter_of(const ct_tree *T, int i, int k)
{
    int l = T->wlen[i];
    k %= l;
    if (k < 0)
        k += l;
    return T->wbuf[T->wstart[i] + k];
}

int ct_letter(const ct_tree *T, int i, int k)
{
    if (i < 0 || i >= T->nwords)
        return CT_EINVALID;
    return letter_of(T, i, k);
}

/* Check the letters of w and write the largest one to *max. */
static int check_word(const ct_tree *T, const int *w, int len, int *max)
{
    int j, m = -1;
    if (len < 0)
        return CT_EINVALID;
    if (len == 0)
        return CT_EEMPTY;
    for (j = 0; j < len; j++) {
        if (w[j] < 0)
            return CT_ENEGATIVE;
        if (T->alphabet_size && w[j] >= T->alphabet_size)
            return CT_EALPHABET;
        if (w[j] > m)
            m = w[j];
    }
    *max = m;
    return CT_OK;
}

/* Append w, whose letters were checked and whose room was reserved. */
static void append_word(ct_tree *T, const int *w, int len, int max)
{
    memcpy(T->wbuf + T->wbuf_size, w, (size_t) len * sizeof(int));
    if (max > T->max_letter)
        T->max_letter = max;
    T->wstart[T->nwords] = T->wbuf_size;
    T->wlen[T->nwords] = len;
    T->wbuf_size += len;
    T->nwords++;
}

/* Undo the last append_word. */
static void pop_word(ct_tree *T)
{
    T->nwords--;
    T->wbuf_size = T->wstart[T->nwords];
}

/* Keep only the first size letters of the last word. */
static void truncate_last_word(ct_tree *T, int size)
{
    int i = T->nwords - 1;
    T->wlen[i] = size;
    T->wbuf_size = T->wstart[i] + size;
}

/* ------------------------------------------------------------------ */
/* children                                                            */
/* ------------------------------------------------------------------ */

/* The first letter of the label into the node t. */
static inline int first_letter(const ct_tree *T, int t)
{
    if (T->layout == CT_LAYOUT_ROWS)
        return T->flet[t];
    return letter_of(T, T->tword[t], T->tstart[t]);
}

/*
 * The child of s by letter in the layouts with sibling lists. Kept out of
 * child(), which is inlined everywhere: inlining this too made the dense
 * layout 5 to 10 % slower at 32 letters, for the same instructions.
 */
static int child_in_lists(const ct_tree *T, int s, int letter)
{
    int t;
    if (T->layout == CT_LAYOUT_SPARSE) {
        for (t = T->fchild[s]; t != -1; t = T->nsib[t])
            if (letter_of(T, T->tword[t], T->tstart[t]) == letter)
                return t;
        return -1;
    }
    t = T->look[s];
    if (t >= 0)
        return T->rows[(size_t) t * (size_t) T->alphabet_size + (size_t) letter];
    if (t < -1) {
        t = crit_find(T, crit_header(T, t)[1], letter);
        return T->flet[t] == letter ? t : -1;
    }
    for (t = T->fchild[s]; t != -1; t = T->nsib[t])
        if (T->flet[t] == letter)
            return t;
    return -1;
}

static inline int child(const ct_tree *T, int s, int letter)
{
    if (T->layout == CT_LAYOUT_DENSE)
        return T->trans[s * T->alphabet_size + letter];
    return child_in_lists(T, s, letter);
}

int ct_child(const ct_tree *T, int s, int letter)
{
    if (s < 0 || s >= T->nstates || letter < 0 ||
        (T->alphabet_size && letter >= T->alphabet_size))
        return CT_EINVALID;
    return child(T, s, letter);
}

/*
 * Make t, whose label starts with letter, a child of s. Return -1 when the
 * tree is broken (a crit node missing in spite of the reservation), else 0.
 */
static inline int add_child(ct_tree *T, int s, int letter, int t)
{
    int d, l, u, lim;
    if (T->layout == CT_LAYOUT_DENSE) {
        T->trans[s * T->alphabet_size + letter] = t;
        return 0;
    }
    u = T->fchild[s];
    T->nsib[t] = u;
    T->fchild[s] = t;
    if (T->layout != CT_LAYOUT_ROWS)
        return 0;
    T->flet[t] = letter;
    l = T->look[s];
    if (l != -1) {
        T->psib[t] = -1;
        if (u != -1)
            T->psib[u] = t;
    }
    if (l >= 0) {
        T->rows[(size_t) l * (size_t) T->alphabet_size + (size_t) letter] = t;
        return 0;
    }
    if (l < -1) {
        int *h = crit_header(T, l);
        d = ++h[2];
        if (crit_insert(T, h, t))
            return -1;
        if (d >= T->promote && T->alphabet_size)
            promote_node(T, s);
        return 0;
    }
    /* a node with neither row nor crit-bit tree has fewer than lim children
     * but in the rare case where an allocation failed, so that counting them
     * (t included) costs as much as the lookup that preceded the insertion */
    lim = T->alphabet_size && T->promote < T->crit_min ? T->promote : T->crit_min;
    for (d = 1; u != -1 && d < lim; u = T->nsib[u])
        d++;
    if (d < lim)
        return 0;
    if (T->alphabet_size && d >= T->promote && promote_node(T, s))
        return 0;
    /* the exact count, past lim only after a failed allocation */
    for (; u != -1; u = T->nsib[u])
        d++;
    if (d >= T->crit_min)
        return crit_build(T, s, d);
    return 0;
}

/*
 * Put new where old, whose label starts with letter, sits among the children
 * of s.
 *
 * In the sibling lists old is found by its index and not by its letter, so
 * that the caller may change the label of old before or after the call. When
 * s has neither row nor crit-bit tree (and always in the sparse layout), the
 * siblings before old are walked, which costs no more than the lookup that
 * found old. The children of a node with a row or a crit-bit tree are linked
 * both ways instead, through psib, which is left stale elsewhere: keeping it
 * up to date in every node made the build of the tree of the Fibonacci words
 * over 1024 letters (2^16 letters per word, nodes with 2 children) 5 %
 * slower, 17.8 against 17.0 ms.
 */
static inline void replace_child(ct_tree *T, int s, int letter, int old, int new_)
{
    int prev, next, c, l = -1;
    if (T->layout == CT_LAYOUT_DENSE) {
        T->trans[s * T->alphabet_size + letter] = new_;
        return;
    }
    if (T->layout == CT_LAYOUT_ROWS)
        l = T->look[s];
    if (l == -1) {
        prev = -1;
        c = T->fchild[s];
        while (c != old) {
            prev = c;
            c = T->nsib[c];
        }
    } else {
        prev = T->psib[old];
        T->psib[new_] = prev;
    }
    next = T->nsib[old];
    T->nsib[new_] = next;
    if (prev == -1)
        T->fchild[s] = new_;
    else
        T->nsib[prev] = new_;
    T->nsib[old] = -1;
    if (T->layout == CT_LAYOUT_ROWS) {
        if (l != -1 && next != -1)
            T->psib[next] = new_;
        T->flet[new_] = letter;
        c = l;
        if (c >= 0) {
            T->rows[(size_t) c * (size_t) T->alphabet_size + (size_t) letter] = new_;
        } else if (c < -1) {
            int *slot = crit_header(T, c) + 1;
            while (*slot < -1) {
                int *cn = T->crit + 3 * (size_t) (-2 - *slot);
                slot = cn + 1 + ((letter >> cn[0]) & 1);
            }
            *slot = new_;
        }
    }
}

/* ------------------------------------------------------------------ */
/* construction                                                        */
/* ------------------------------------------------------------------ */

#ifdef CT_DEBUG
/*
 * Local checks, run after each step of an insertion when CT_DEBUG is
 * defined. ct_check does not hold in the middle of an insertion (the suffix
 * link of the last internal node is set one step late, and the leaves of the
 * word being inserted are not all there yet), but these do, and each costs
 * O(1) except debug_active, which costs the length of the active edge. On a
 * failure they mark the tree broken.
 */

static int debug_fail(ct_tree *T)
{
    T->broken = CT_EINTERNAL;
    return 0;
}

/* ss was just made explicit between s and t: the three nodes are linked both
 * ways and the depths agree with the labels. */
static int debug_split(ct_tree *T, int s, int ss, int t)
{
    if (T->parent[ss] != s || T->parent[t] != ss ||
        T->tend[ss] <= T->tstart[ss] || T->tstart[t] != T->tend[ss] ||
        T->tword[t] != T->tword[ss] ||
        child(T, s, letter_of(T, T->tword[ss], T->tstart[ss])) != ss ||
        child(T, ss, letter_of(T, T->tword[t], T->tstart[t])) != t ||
        T->dep[ss] != T->dep[s] + T->tend[ss] - T->tstart[ss] ||
        (T->tend[t] != -1 && T->dep[t] != T->dep[ss] + T->tend[t] - T->tstart[t]))
        return debug_fail(T);
    return 1;
}

/* rr was just added as a leaf under r, reading letter first. */
static int debug_leaf(ct_tree *T, int r, int rr, int letter)
{
    if (T->parent[rr] != r || T->tend[rr] != -1 ||
        letter_of(T, T->tword[rr], T->tstart[rr]) != letter ||
        child(T, r, letter) != rr)
        return debug_fail(T);
    return 1;
}

/* The suffix link of the internal node r was just set. */
static int debug_link(ct_tree *T, int r)
{
    int t = T->sl[r];
    if (t < 0 || t >= T->nstates || (t != 0 && T->tend[t] == -1) ||
        T->dep[t] != T->dep[r] - 1)
        return debug_fail(T);
    return 1;
}

/* (s, i, k, p) is canonical and reads the letters k, ..., p - 1 of the i-th
 * word along the tree. */
static int debug_active(ct_tree *T, int s, int i, int k, int p)
{
    int t, j;
    if (k > p || s < -1 || s >= T->nstates || (s > 0 && T->tend[s] == -1))
        return debug_fail(T);
    if (k == p)
        return 1;
    if (s == -1)
        return debug_fail(T);
    t = child(T, s, letter_of(T, i, k));
    if (t == -1 || (T->tend[t] != -1 && T->tend[t] - T->tstart[t] <= p - k))
        return debug_fail(T);
    for (j = 0; j < p - k; j++)
        if (letter_of(T, i, k + j) != letter_of(T, T->tword[t], T->tstart[t] + j))
            return debug_fail(T);
    return 1;
}

/* The whole tree, after an insertion. A CT_ENOMEM of ct_check says nothing
 * about the tree, whose insertion went through: only an inconsistency is an
 * error. */
static int debug_tree(ct_tree *T)
{
    if (ct_check(T) == CT_EINTERNAL) {
        T->broken = CT_EINTERNAL;
        return CT_EINTERNAL;
    }
    return CT_OK;
}
#endif

/*
 * Given the canonical reference (s, i, k, p), test whether reading letter
 * from it creates a branching. Return -1 if the transition exists, otherwise
 * the node from which the new transition starts, which is made explicit if
 * it was not. Return -2 on failure.
 */
static inline int test_and_split(ct_tree *T, int s, int i, int k, int p, int letter)
{
    int t, ii, kk, index, lletter, ss, first;
    if (k < p) {
        /* implicit state: get the transition from s starting with
         * word[i][k] and test whether its (p - k)-th letter is letter */
        t = child(T, s, letter_of(T, i, k));
        ii = T->tword[t];
        kk = T->tstart[t];
        index = kk + p - k;
        lletter = letter_of(T, ii, index);
        if (letter == lletter)
            /* the node already exists */
            return -1;
        /* make the node explicit: the new node ss splits s ---> t into
         * s --> ss --> t */
        first = letter_of(T, ii, kk);
        ss = add_node(T);
        if (ss < 0)
            return -2;

        T->tword[ss] = ii;
        T->tstart[ss] = kk;
        T->tend[ss] = index;
        T->parent[ss] = s;
        T->dep[ss] = T->dep[s] + index - kk;

        replace_child(T, s, first, t, ss);
        T->tstart[t] = index;
        T->parent[t] = ss;
        if (add_child(T, ss, lletter, t))
            return -2;
#ifdef CT_DEBUG
        if (!debug_split(T, s, ss, t))
            return -2;
#endif
        return ss;
    }
    /* explicit state */
    if (s == -1 || child(T, s, letter) != -1)
        /* the node already exists */
        return -1;
    return s;
}

/*
 * Canonize the reference (s, i, k, p) of the (explicit or implicit) state
 * obtained after reading word[i][k:p] from s, writing the answer back into
 * sp and kp; i and p do not change.
 */
static inline void canonize(const ct_tree *T, int *sp, int i, int *kp, int p)
{
    int s = *sp;
    int k = *kp;
    int ss, kk, pp;
    if (k >= p) {
        /* already explicit */
        *kp = p;
        return;
    }
    ss = s == -1 ? 0 : child(T, s, letter_of(T, i, k));
    kk = T->tstart[ss];
    pp = T->tend[ss];
    while (pp != -1 && pp - kk < p - k) {
        k += pp - kk;
        s = ss;
        ss = child(T, s, letter_of(T, i, k));
        kk = T->tstart[ss];
        pp = T->tend[ss];
    }
    if (pp != -1 && pp - kk == p - k) {
        /* explicit */
        *sp = ss;
        *kp = p;
    } else {
        /* implicit */
        *sp = s;
        *kp = k;
    }
}

/*
 * The same as canonize, for a reference that may not describe a path of the
 * tree: canonize reads tstart[-1] when a child is missing.
 */
int ct_canonize(const ct_tree *T, int *sp, int i, int *kp, int p)
{
    int s = *sp;
    int k = *kp;
    int ss;
    if (s < -1 || s >= T->nstates || i < 0 || i >= T->nwords || k < 0 || k > p)
        return CT_EINVALID;
    if (k == p)
        return CT_OK;
    ss = s == -1 ? 0 : child(T, s, letter_of(T, i, k));
    while (ss != -1 && T->tend[ss] != -1 && T->tend[ss] - T->tstart[ss] < p - k) {
        k += T->tend[ss] - T->tstart[ss];
        s = ss;
        ss = child(T, s, letter_of(T, i, k));
    }
    if (ss == -1)
        return CT_EINVALID;
    if (T->tend[ss] != -1 && T->tend[ss] - T->tstart[ss] == p - k) {
        *sp = ss;
        *kp = p;
    } else {
        *sp = s;
        *kp = k;
    }
    return CT_OK;
}

/*
 * Read the letter p of the i-th word, where (*sp, i, *kp, p) is the canonical
 * reference of the active state. Return the number of leaves created, or -1
 * on failure.
 */
static int update(ct_tree *T, int *sp, int i, int *kp, int p)
{
    /* (s, k, p): active state, the first state along the boundary path
     * which is not an active leaf; r: closest branching from s (either s or
     * one of its ancestors) */
    int s = *sp;
    int k = *kp;
    int letter = letter_of(T, i, p);
    int old_r = 0;
    int created = 0;
    int r, rr;
    r = test_and_split(T, s, i, k, p, letter);
    while (r != -1) {
        if (r < 0)
            return -1;
        /* create a leaf */
        rr = add_node(T);
        if (rr < 0)
            return -1;
        created++;
        if (add_child(T, r, letter, rr))
            return -1;
        T->parent[rr] = r;
        T->tword[rr] = i;
        T->tstart[rr] = p;
        T->tend[rr] = -1;
#ifdef CT_DEBUG
        if (!debug_leaf(T, r, rr, letter))
            return -1;
#endif
        if (old_r != 0) {
            T->sl[old_r] = r;
#ifdef CT_DEBUG
            if (!debug_link(T, old_r))
                return -1;
#endif
        }
        old_r = r;
        /* NOTE: s is never -1 here; test_and_split returns -1 on the
         * imaginary node, which ends the loop */
        s = T->sl[s];
        canonize(T, &s, i, &k, p);
        r = test_and_split(T, s, i, k, p, letter);
    }

    if (old_r != 0) {
        T->sl[old_r] = s;
#ifdef CT_DEBUG
        if (!debug_link(T, old_r))
            return -1;
#endif
    }

    *sp = s;
    *kp = k;
    return created;
}

/*
 * Insert the last word of the word buffer, for which room was reserved, and
 * write what ct_process returns to *result.
 */
static int insert_last(ct_tree *T, int *result)
{
    int i = T->nwords - 1;
    int l = T->wlen[i];
    int s = 0;
    int k = 0;
    int p = 0;
    int num_leaves = 0;
    int ss, ii, pp, created;

    /* To ensure that we find all conjugates we must create as many leaves as
     * the length of w (assuming it is primitive) */
    for (;;) {
        if (p != k) {
            ss = child(T, s, letter_of(T, i, k));
            ii = T->tword[ss];
            pp = T->tend[ss];
        } else {
            ii = -1;
            pp = -2;
        }
        created = update(T, &s, i, &k, p);
        if (created < 0)
            return T->broken;
        num_leaves += created;
        canonize(T, &s, i, &k, p + 1);
#ifdef CT_DEBUG
        if (!debug_active(T, s, i, k, p + 1))
            return T->broken;
#endif

        /* halt condition */
        if (num_leaves == l)
            /* w is primitive */
            break;
        else if (ii == i && p >= 2 * l)
            /* w is non primitive */
            break;
        else if (p >= l && num_leaves == 0 && ii != -1 && pp == -1 && l % T->wlen[ii] == 0)
            /* w is conjugate to a power of the ii-th word */
            break;

        p++;
    }

    if (num_leaves == 0) {
        pop_word(T);
        *result = -ii;
#ifdef CT_DEBUG
        return debug_tree(T);
#else
        return CT_OK;
#endif
    }
    if (l % num_leaves) {
        /* the length is not a multiple of the number of new leaves */
        T->broken = CT_EINTERNAL;
        return CT_EINTERNAL;
    }
    if (num_leaves != l)
        /* NOTE: only store primitive words */
        truncate_last_word(T, num_leaves);
    *result = l / num_leaves;
#ifdef CT_DEBUG
    return debug_tree(T);
#else
    return CT_OK;
#endif
}

int ct_process(ct_tree *T, const int *w, int len, int *result)
{
    int err, max;
    if (T->broken)
        return CT_EBROKEN;
    if (T->nstates == 0)
        return CT_EINVALID;
    if ((err = check_word(T, w, len, &max)))
        return err;
    if (len > (INT_MAX - T->nstates) / 2)
        return CT_ETOOLARGE;
    if ((err = reserve(T, 1, len, 2 * len)))
        return err;
    append_word(T, w, len, max);
    return insert_last(T, result);
}

int ct_reserve(ct_tree *T, int words, int letters)
{
    if (T->broken)
        return CT_EBROKEN;
    if (T->nstates == 0 || words < 0 || letters < 0)
        return CT_EINVALID;
    /* the same bound as in ct_process, for the letters of all the words */
    if (letters > (INT_MAX - T->nstates) / 2)
        return CT_ETOOLARGE;
    return reserve(T, words, letters, 2 * letters);
}

/* ------------------------------------------------------------------ */
/* inspection                                                          */
/* ------------------------------------------------------------------ */

int ct_leaf_as_conjugate(const ct_tree *T, int s, int *i, int *k)
{
    int l, ans;
    if (s <= 0 || s >= T->nstates || T->tend[s] != -1)
        return CT_EINVALID;
    *i = T->tword[s];
    l = T->wlen[*i];
    ans = (T->tstart[s] - T->dep[T->parent[s]]) % l;
    if (ans < 0)
        ans += l;
    *k = ans;
    return CT_OK;
}

int64_t ct_size(const ct_tree *T)
{
    int64_t ans = 0;
    int s;
    if (T->nstates == 0)
        return CT_EINVALID;
    for (s = 0; s < T->nstates; s++) {
        if (T->tend[s] == -1)
            ans += 1;
        else
            ans += T->tend[s] - T->tstart[s];
    }
    return ans;
}

/*
 * The children of a node are sorted in three ways: by an insertion sort up
 * to CT_SORT_RUN of them, the usual case (most internal nodes have 2 to 10
 * children whatever the alphabet); by placing them at their key in an array
 * of n slots and reading it back when there are at least n / CT_SORT_BUCKET
 * of them, which costs O(n), at most CT_SORT_BUCKET per child; and by a merge
 * sort, O(d log d), in between.
 */
#define CT_SORT_RUN 16
#define CT_SORT_BUCKET 8

/*
 * Sort a[:d] by decreasing value with tmp as scratch: insertion sorts of runs
 * of CT_SORT_RUN, then bottom-up merges.
 */
static void sort_decreasing(uint64_t *a, uint64_t *tmp, int d)
{
    uint64_t *src = a, *dst = tmp, *sw, x;
    size_t n = (size_t) d, w, lo, mid, hi, i, j, k, r;
    for (lo = 0; lo < n; lo += CT_SORT_RUN) {
        hi = lo + CT_SORT_RUN < n ? lo + CT_SORT_RUN : n;
        for (i = lo + 1; i < hi; i++) {
            x = a[i];
            for (j = i; j > lo && a[j - 1] < x; j--)
                a[j] = a[j - 1];
            a[j] = x;
        }
    }
    for (w = CT_SORT_RUN; w < n; w *= 2) {
        for (lo = 0; lo < n; lo += 2 * w) {
            mid = lo + w < n ? lo + w : n;
            hi = lo + 2 * w < n ? lo + 2 * w : n;
            i = lo;
            j = mid;
            for (k = lo; k < hi; k++)
                dst[k] = j == hi || (i < mid && src[i] > src[j]) ? src[i++] : src[j++];
        }
        sw = src;
        src = dst;
        dst = sw;
    }
    if (src != a)
        for (r = 0; r < n; r++)
            a[r] = src[r];
}

/*
 * Push kids[:d] onto stack by decreasing keys (distinct, smaller than n),
 * so that popping them takes them by increasing key. slot is an array of n
 * entries -1, which it is again on return; a and tmp have room for n values.
 *
 * The insertion sort works on the two int arrays: packing each pair into one
 * uint64_t made the listing of the leaves of the Fibonacci words over 1024
 * letters (2^16 letters per word, almost every internal node with 2 children)
 * 50 % slower, 11.8 ms against 7.8 ms.
 */
static int push_sorted(int *stack, int top, int *kids, int *keys, uint64_t *a,
                       uint64_t *tmp, int *slot, int d, int n)
{
    int j, jj, kid, key, c;
    if (d <= CT_SORT_RUN) {
        for (j = 1; j < d; j++) {
            kid = kids[j];
            key = keys[j];
            jj = j - 1;
            while (jj >= 0 && keys[jj] < key) {
                kids[jj + 1] = kids[jj];
                keys[jj + 1] = keys[jj];
                jj--;
            }
            kids[jj + 1] = kid;
            keys[jj + 1] = key;
        }
        for (j = 0; j < d; j++)
            stack[top + j] = kids[j];
        return top + d;
    }
    if ((int64_t) d * CT_SORT_BUCKET >= n) {
        for (j = 0; j < d; j++)
            slot[keys[j]] = kids[j];
        for (c = n - 1; c >= 0; c--) {
            if (slot[c] != -1) {
                stack[top++] = slot[c];
                slot[c] = -1;
            }
        }
        return top;
    }
    for (j = 0; j < d; j++)
        a[j] = (uint64_t) keys[j] << 32 | (uint64_t) kids[j];
    sort_decreasing(a, tmp, d);
    for (j = 0; j < d; j++)
        stack[top + j] = (int) (a[j] & 0xffffffffu);
    return top + d;
}

int ct_sorted_leaves(const ct_tree *T, const int *order, const int *pivot, int n,
                     int *out, int *num)
{
    int *stack, *kids, *keys, *slot;
    uint64_t *a, *tmp;
    int top = 0;
    int count = 0;
    int c, s, t, d, key, base;

    if (T->nstates == 0 || n <= 0 || T->max_letter >= n)
        return CT_EINVALID;
    for (c = 0; c < n; c++)
        if (order[c] < 0 || order[c] >= n || pivot[c] < 0 || pivot[c] >= n)
            return CT_EINVALID;

    /* each node is pushed at most once, and a node has at most n children
     * since its children start with distinct letters, all smaller than n */
    stack = (int *) malloc((size_t) T->nstates * sizeof(int));
    kids = (int *) malloc((size_t) n * sizeof(int));
    keys = (int *) malloc((size_t) n * sizeof(int));
    slot = (int *) malloc((size_t) n * sizeof(int));
    a = (uint64_t *) malloc((size_t) n * sizeof(uint64_t));
    tmp = (uint64_t *) malloc((size_t) n * sizeof(uint64_t));
    if (stack == NULL || kids == NULL || keys == NULL || slot == NULL || a == NULL ||
        tmp == NULL) {
        free(stack);
        free(kids);
        free(keys);
        free(slot);
        free(a);
        free(tmp);
        return CT_ENOMEM;
    }

    /* order is a permutation: slot[] marks the values seen, which also
     * makes the keys of the children of a node distinct */
    for (c = 0; c < n; c++)
        slot[c] = -1;
    for (c = 0; c < n; c++) {
        if (slot[order[c]] != -1) {
            free(stack);
            free(kids);
            free(keys);
            free(slot);
            free(a);
            free(tmp);
            return CT_EINVALID;
        }
        slot[order[c]] = 0;
    }
    for (c = 0; c < n; c++)
        slot[c] = -1;

    /* the children of the root, ordered by order[] of their first letter */
    d = 0;
    if (T->layout == CT_LAYOUT_DENSE) {
        for (c = 0; c < T->alphabet_size; c++) {
            t = T->trans[c];
            if (t != -1) {
                kids[d] = t;
                keys[d] = order[c];
                d++;
            }
        }
    } else {
        for (t = T->fchild[0]; t != -1; t = T->nsib[t]) {
            kids[d] = t;
            keys[d] = order[first_letter(T, t)];
            d++;
        }
    }
    top = push_sorted(stack, top, kids, keys, a, tmp, slot, d, n);

    while (top) {
        s = stack[--top];
        if (T->tend[s] == -1) {
            out[count++] = s;
            continue;
        }
        /* further down, order[] is read from the pivot of the last letter
         * of the label */
        base = pivot[letter_of(T, T->tword[s], T->tend[s] - 1)];
        d = 0;
        if (T->layout == CT_LAYOUT_DENSE) {
            const int *row = T->trans + (size_t) s * (size_t) T->alphabet_size;
            for (c = 0; c < T->alphabet_size; c++) {
                t = row[c];
                if (t != -1) {
                    key = order[c] - base;
                    if (key < 0)
                        key += n;
                    kids[d] = t;
                    keys[d] = key;
                    d++;
                }
            }
        } else {
            for (t = T->fchild[s]; t != -1; t = T->nsib[t]) {
                key = order[first_letter(T, t)] - base;
                if (key < 0)
                    key += n;
                kids[d] = t;
                keys[d] = key;
                d++;
            }
        }
        top = push_sorted(stack, top, kids, keys, a, tmp, slot, d, n);
    }

    free(stack);
    free(kids);
    free(keys);
    free(slot);
    free(a);
    free(tmp);
    *num = count;
    return CT_OK;
}

const char *ct_strerror(int code)
{
    switch (code) {
    case CT_OK: return "no error";
    case CT_ENOMEM: return "out of memory";
    case CT_EEMPTY: return "empty word";
    case CT_ENEGATIVE: return "negative letter";
    case CT_EALPHABET: return "letter outside the alphabet";
    case CT_ETOOLARGE: return "tree too large for int indices or for a dense table";
    case CT_EINVALID: return "invalid argument";
    case CT_EINTERNAL: return "internal inconsistency";
    case CT_EBROKEN: return "tree broken by an earlier internal inconsistency";
    default: return "unknown error";
    }
}

/* ------------------------------------------------------------------ */
/* self-checks                                                         */
/* ------------------------------------------------------------------ */

/* Write the word of the internal node (or root) s to buf[0..dep[s]-1]. */
static void node_word(const ct_tree *T, int s, int *buf)
{
    int pos = T->dep[s];
    int j, len;
    while (s != 0) {
        len = T->tend[s] - T->tstart[s];
        pos -= len;
        for (j = 0; j < len; j++)
            buf[pos + j] = letter_of(T, T->tword[s], T->tstart[s] + j);
        s = T->parent[s];
    }
}

/*
 * Check the part under the reference p of the crit-bit tree of s, whose crit
 * nodes have bits below above: every crit node is in no other tree and has
 * both sides, its children are children of s seen nowhere else, and on each
 * side of a crit node of bit b the first letters agree above b and have the
 * bit b of the side. Accumulate the number of children in *count and the AND
 * and OR of their first letters in *all and *any. The bits decrease strictly
 * downwards, so that the recursion is at most 32 deep. Return 1 when it holds.
 */
static int check_crit(const ct_tree *T, int s, int p, int above, char *owned, char *seen,
                      int *count, unsigned *all, unsigned *any)
{
    int j, b, c0 = 0, c1 = 0;
    unsigned all0, any0, all1, any1;
    const int *c;
    if (p >= 0) {
        if (p >= T->nstates || T->parent[p] != s || seen[p])
            return 0;
        seen[p] = 1;
        *count = 1;
        *all = *any = (unsigned) T->flet[p];
        return 1;
    }
    if (p == -1 || -(p + 2) >= T->ncrit)
        return 0;
    j = -(p + 2);
    if (owned[j])
        return 0;
    owned[j] = 1;
    c = T->crit + 3 * (size_t) j;
    b = c[0];
    if (b < 0 || b >= above)
        return 0;
    if (!check_crit(T, s, c[1], b, owned, seen, &c0, &all0, &any0) ||
        !check_crit(T, s, c[2], b, owned, seen, &c1, &all1, &any1))
        return 0;
    if ((any0 >> b) & 1 || !((all1 >> b) & 1))
        return 0;
    *all = all0 & all1;
    *any = any0 | any1;
    if ((*all ^ *any) >> b >> 1)
        return 0;
    *count = c0 + c1;
    return 1;
}

#define CHECK(cond) do { if (!(cond)) { err = CT_EINTERNAL; goto done; } } while (0)

int ct_check(const ct_tree *T)
{
    int n = T->nstates;
    int err = CT_OK;
    int *buf = NULL, *buf2 = NULL, *owned = NULL;
    char *crit_owned = NULL, *seen = NULL;
    int i, j, s, t, r, c, k, count, steps, max_dep, ei, ek;

    CHECK(n >= 1 && n <= T->capacity);
    CHECK(T->broken == 0);
    CHECK(T->alphabet_size >= 0);
    CHECK(T->layout >= CT_LAYOUT_SPARSE && T->layout <= CT_LAYOUT_ROWS);
    CHECK(T->layout != CT_LAYOUT_DENSE || T->alphabet_size > 0);
    if (T->layout == CT_LAYOUT_ROWS) {
        CHECK(T->promote >= 1);
        CHECK(T->nrows >= 0 && T->nrows <= T->rows_capacity);
        /* without an alphabet, rows cannot be allocated */
        CHECK(T->alphabet_size > 0 || T->nrows == 0);
        CHECK(T->crit_min >= 1);
        CHECK(T->ncrit >= 0 && T->ncrit <= T->crit_capacity);
        CHECK(T->crit_free >= -1 && T->crit_free < T->ncrit);
        /* the room promised to the insertions into the crit-bit trees */
        CHECK((int64_t) T->crit_capacity - T->ncrit >= (int64_t) T->crit_promised - n);
    }

    /* the words: laid end to end, non-empty, letters in range */
    CHECK(T->nwords >= 0 && T->nwords <= T->words_capacity);
    CHECK(T->wbuf_size >= 0 && T->wbuf_size <= T->wbuf_capacity);
    for (i = 0; i < T->nwords; i++) {
        CHECK(T->wlen[i] > 0);
        CHECK(T->wstart[i] == (i ? T->wstart[i - 1] + T->wlen[i - 1] : 0));
    }
    CHECK(T->wbuf_size == (T->nwords ? T->wstart[T->nwords - 1] + T->wlen[T->nwords - 1] : 0));
    for (j = 0; j < T->wbuf_size; j++) {
        CHECK(T->wbuf[j] >= 0 && T->wbuf[j] <= T->max_letter);
        CHECK(!T->alphabet_size || T->wbuf[j] < T->alphabet_size);
    }

    /* the root */
    CHECK(T->parent[0] == -1 && T->sl[0] == -1 && T->dep[0] == 0);
    CHECK(T->tword[0] == 0 && T->tstart[0] == -4 && T->tend[0] == -3);

    /* the labels, and the depths of the internal nodes; since a label is
     * non-empty, dep increases strictly from a node to its children, which
     * rules out a cycle of parents */
    max_dep = 0;
    for (s = 1; s < n; s++) {
        CHECK(T->parent[s] >= 0 && T->parent[s] < n);
        CHECK(T->parent[s] == 0 || T->tend[T->parent[s]] != -1);
        CHECK(T->tword[s] >= 0 && T->tword[s] < T->nwords);
        CHECK(T->tstart[s] >= 0);
        CHECK(T->tend[s] == -1 || T->tend[s] > T->tstart[s]);
        if (T->tend[s] != -1) {
            /* in 64 bits, since a corrupted depth may be close to INT_MAX */
            CHECK(T->dep[s] == (int64_t) T->dep[T->parent[s]] + T->tend[s] - T->tstart[s]);
            if (T->dep[s] > max_dep)
                max_dep = T->dep[s];
        }
    }

    /* the first letters */
    if (T->layout == CT_LAYOUT_ROWS) {
        CHECK(T->flet[0] == -1);
        for (s = 1; s < n; s++)
            CHECK(T->flet[s] == letter_of(T, T->tword[s], T->tstart[s]));
    }

    /* the children: a child of s has parent s and its slot is the first
     * letter of its label; a node has children if and only if it is not a
     * leaf; every node but the root is found under its parent by its first
     * letter, and there are n - 1 children in all, so each one exactly once */
    if (T->layout == CT_LAYOUT_ROWS) {
        owned = (int *) calloc((size_t) T->nrows + 1, sizeof(int));
        crit_owned = (char *) calloc((size_t) T->ncrit + 1, 1);
        seen = (char *) calloc((size_t) n, 1);
        if (owned == NULL || crit_owned == NULL || seen == NULL) {
            err = CT_ENOMEM;
            goto done;
        }
    }
    count = 0;
    for (s = 0; s < n; s++) {
        int d = 0;
        if (T->layout == CT_LAYOUT_DENSE) {
            for (c = 0; c < T->alphabet_size; c++) {
                t = T->trans[s * T->alphabet_size + c];
                if (t == -1)
                    continue;
                CHECK(t >= 1 && t < n && T->parent[t] == s);
                CHECK(letter_of(T, T->tword[t], T->tstart[t]) == c);
                d++;
            }
        } else {
            int prev = -1;
            steps = 0;
            for (t = T->fchild[s]; t != -1; t = T->nsib[t]) {
                CHECK(++steps < n);
                CHECK(t >= 1 && t < n && T->parent[t] == s);
                /* the rows layout links the children of a node with a row
                 * or a crit-bit tree both ways */
                CHECK(T->layout != CT_LAYOUT_ROWS || T->look[s] == -1 || T->psib[t] == prev);
                prev = t;
                d++;
            }
        }
        CHECK(s == 0 || (T->tend[s] == -1) == (d == 0));
        count += d;
        /* a crit-bit tree belongs to one node with at least crit_min
         * children, counts them in its header, and holds exactly them; a node
         * with that many children may lack it (see crit_build) */
        if (T->layout == CT_LAYOUT_ROWS && T->look[s] < -1) {
            int e = 0;
            unsigned all, any;
            const int *h;
            j = -(T->look[s] + 2);
            CHECK(j < T->ncrit && !crit_owned[j]);
            crit_owned[j] = 1;
            h = T->crit + 3 * (size_t) j;
            CHECK(h[0] == -2 && h[2] == d && d >= T->crit_min);
            CHECK(check_crit(T, s, h[1], 32, crit_owned, seen, &e, &all, &any));
            CHECK(e == d);
        }
        /* a row belongs to one node with at least promote children, and
         * holds exactly its children, each at its first letter; a node with
         * that many children may lack a row (see promote_node) */
        if (T->layout == CT_LAYOUT_ROWS && T->look[s] >= 0) {
            const int *row;
            int e = 0;
            r = T->look[s];
            CHECK(r < T->nrows && !owned[r]);
            owned[r] = 1;
            CHECK(d >= T->promote);
            row = T->rows + (size_t) r * (size_t) T->alphabet_size;
            for (c = 0; c < T->alphabet_size; c++) {
                t = row[c];
                if (t == -1)
                    continue;
                CHECK(t >= 1 && t < n && T->parent[t] == s && T->flet[t] == c);
                e++;
            }
            CHECK(e == d);
        }
    }
    CHECK(count == n - 1);
    if (T->layout == CT_LAYOUT_ROWS) {
        CHECK(T->psib[0] == -1);
        for (r = 0; r < T->nrows; r++)
            CHECK(owned[r]);
        /* every crit node is either in a crit-bit tree or free */
        for (j = T->crit_free; j != -1; j = T->crit[3 * (size_t) j + 1]) {
            CHECK(j >= 0 && j < T->ncrit && !crit_owned[j]);
            CHECK(T->crit[3 * (size_t) j] == -1);
            crit_owned[j] = 1;
        }
        for (j = 0; j < T->ncrit; j++)
            CHECK(crit_owned[j]);
    }
    for (t = 1; t < n; t++)
        CHECK(child(T, T->parent[t], letter_of(T, T->tword[t], T->tstart[t])) == t);

    buf = (int *) malloc(((size_t) max_dep + 1) * sizeof(int));
    buf2 = (int *) malloc(((size_t) max_dep + 1) * sizeof(int));
    if (buf == NULL || buf2 == NULL) {
        err = CT_ENOMEM;
        goto done;
    }

    /* the suffix link of an internal node reads its word minus the first
     * letter */
    for (s = 1; s < n; s++) {
        if (T->tend[s] == -1)
            continue;
        t = T->sl[s];
        CHECK(t >= 0 && t < n && (t == 0 || T->tend[t] != -1));
        CHECK(T->dep[t] == T->dep[s] - 1);
        node_word(T, s, buf);
        node_word(T, t, buf2);
        for (j = 0; j + 1 < T->dep[s]; j++)
            CHECK(buf[j + 1] == buf2[j]);
    }

    /* the word of the parent of a node is read in the word of its label,
     * right before the label */
    for (s = 1; s < n; s++) {
        r = T->parent[s];
        i = T->tword[s];
        node_word(T, r, buf);
        for (j = 0; j < T->dep[r]; j++)
            CHECK(buf[j] == letter_of(T, i, T->tstart[s] - T->dep[r] + j));
    }

    /* the leaves, by increasing index, are the conjugates (i, k) in
     * lexicographic order */
    ei = 0;
    ek = 0;
    for (s = 1; s < n; s++) {
        if (T->tend[s] != -1)
            continue;
        CHECK(ei < T->nwords);
        CHECK(ct_leaf_as_conjugate(T, s, &i, &k) == CT_OK);
        CHECK(i == ei && k == ek);
        if (++ek == T->wlen[ei]) {
            ei++;
            ek = 0;
        }
    }
    CHECK(ei == T->nwords && ek == 0);

done:
    free(buf);
    free(buf2);
    free(owned);
    free(crit_owned);
    free(seen);
    return err;
}
