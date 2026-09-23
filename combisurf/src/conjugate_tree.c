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
 * A dense table costs 4 * alphabet bytes per node, where the sibling lists
 * cost 8 and the rest of a node 24, and every new node clears its row and
 * every listing of the leaves scans it. Measured on the pair of curves u, v
 * of the one-vertex map with n half-edges (so an alphabet of n letters), the
 * dense table makes the pairing slower on short curves and faster on long
 * ones, by these ratios of dense over sparse time:
 *
 *     n = alphabet     32     64    128    256
 *     length 8       1.03   1.06   1.11   1.20
 *     length 100     0.91   0.87   0.86   0.88
 *     length 1000    0.88   0.85   0.77   0.69
 *
 * A single word of length 100000 over 128 letters has 225000 nodes, that is
 * 110 MB of dense table against 1.7 MB of sibling lists. Up to 32 letters the
 * dense table costs at most 3 % on short curves and 128 bytes per node; past
 * it the loss on short curves, which are the common case, grows with the
 * alphabet, and so does the memory. So past it we walk the children of a node
 * instead of indexing them (CT_DENSE_MAX_ALPHABET).
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
    return T->dense ? INT_MAX / T->alphabet_size : INT_MAX;
}

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
    if (T->dense) {
        /* cap * alphabet_size <= INT_MAX by max_nodes */
        if ((err = resize_ints(&T->trans, cap * T->alphabet_size)))
            return err;
    } else {
        if ((err = resize_ints(&T->fchild, cap)) ||
            (err = resize_ints(&T->nsib, cap)))
            return err;
    }
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
    return reserve_nodes(T, T->nstates + num_nodes);
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
    if (T->dense) {
        int *row = T->trans + (size_t) n * (size_t) T->alphabet_size;
        for (c = 0; c < T->alphabet_size; c++)
            row[c] = -1;
    } else {
        T->fchild[n] = -1;
        T->nsib[n] = -1;
    }
    return n;
}

int ct_init(ct_tree *T, int alphabet, int reserve, int dense)
{
    int err;
    memset(T, 0, sizeof(ct_tree));
    T->max_letter = -1;
    if (alphabet < 0 || reserve < 0 || (dense > 0 && alphabet == 0))
        return CT_EINVALID;
    T->alphabet_size = alphabet;
    if (dense < 0)
        T->dense = 0 < alphabet && alphabet <= CT_DENSE_MAX_ALPHABET;
    else
        T->dense = dense != 0;

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

static inline int child(const ct_tree *T, int s, int letter)
{
    int t;
    if (T->dense)
        return T->trans[s * T->alphabet_size + letter];
    t = T->fchild[s];
    while (t != -1) {
        if (letter_of(T, T->tword[t], T->tstart[t]) == letter)
            return t;
        t = T->nsib[t];
    }
    return -1;
}

int ct_child(const ct_tree *T, int s, int letter)
{
    return child(T, s, letter);
}

static inline void add_child(ct_tree *T, int s, int letter, int t)
{
    if (T->dense) {
        T->trans[s * T->alphabet_size + letter] = t;
    } else {
        T->nsib[t] = T->fchild[s];
        T->fchild[s] = t;
    }
}

/*
 * Put new where old sits among the children of s.
 *
 * In the sparse representation old is found by its index and not by its
 * letter, since the letter of a node is read off the label of the edge into
 * it and a caller splitting that edge is about to move it.
 */
static void replace_child(ct_tree *T, int s, int letter, int old, int new_)
{
    int prev, c;
    if (T->dense) {
        T->trans[s * T->alphabet_size + letter] = new_;
        return;
    }
    prev = -1;
    c = T->fchild[s];
    while (c != old) {
        prev = c;
        c = T->nsib[c];
    }
    T->nsib[new_] = T->nsib[old];
    if (prev == -1)
        T->fchild[s] = new_;
    else
        T->nsib[prev] = new_;
    T->nsib[old] = -1;
}

/* ------------------------------------------------------------------ */
/* construction                                                        */
/* ------------------------------------------------------------------ */

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

        /* NOTE: ss takes the place of t under s before the label of t is
         * shortened, since that label is where its first letter is read */
        replace_child(T, s, first, t, ss);
        T->tstart[t] = index;
        T->parent[t] = ss;
        add_child(T, ss, lletter, t);

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

void ct_canonize(const ct_tree *T, int *s, int i, int *k, int p)
{
    canonize(T, s, i, k, p);
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
        add_child(T, r, letter, rr);
        T->parent[rr] = r;
        T->tword[rr] = i;
        T->tstart[rr] = p;
        T->tend[rr] = -1;
        if (old_r != 0)
            T->sl[old_r] = r;
        old_r = r;
        /* NOTE: s is never -1 here; test_and_split returns -1 on the
         * imaginary node, which ends the loop */
        s = T->sl[s];
        canonize(T, &s, i, &k, p);
        r = test_and_split(T, s, i, k, p, letter);
    }

    if (old_r != 0)
        T->sl[old_r] = s;

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
        return CT_OK;
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
    return CT_OK;
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
    for (s = 0; s < T->nstates; s++) {
        if (T->tend[s] == -1)
            ans += 1;
        else
            ans += T->tend[s] - T->tstart[s];
    }
    return ans;
}

/*
 * Push kids[:d] onto stack by decreasing keys, so that popping them takes
 * them by increasing key.
 *
 * An insertion sort is the right one here: the number of children of an
 * internal node is 2 to 10 on average whatever the alphabet and whatever the
 * length of the words.
 */
static int push_sorted(int *stack, int top, int *kids, int *keys, int d)
{
    int j, jj, kid, key;
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

int ct_sorted_leaves(const ct_tree *T, const int *order, const int *pivot, int n,
                     int *out, int *num)
{
    int *stack, *kids, *keys;
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
    if (stack == NULL || kids == NULL || keys == NULL) {
        free(stack);
        free(kids);
        free(keys);
        return CT_ENOMEM;
    }

    /* order is a permutation: keys[] marks the values seen, which also
     * makes the keys of the children of a node distinct */
    for (c = 0; c < n; c++)
        keys[c] = 0;
    for (c = 0; c < n; c++) {
        if (keys[order[c]]) {
            free(stack);
            free(kids);
            free(keys);
            return CT_EINVALID;
        }
        keys[order[c]] = 1;
    }

    /* the children of the root, ordered by order[] of their first letter */
    d = 0;
    if (T->dense) {
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
            keys[d] = order[letter_of(T, T->tword[t], T->tstart[t])];
            d++;
        }
    }
    top = push_sorted(stack, top, kids, keys, d);

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
        if (T->dense) {
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
                key = order[letter_of(T, T->tword[t], T->tstart[t])] - base;
                if (key < 0)
                    key += n;
                kids[d] = t;
                keys[d] = key;
                d++;
            }
        }
        top = push_sorted(stack, top, kids, keys, d);
    }

    free(stack);
    free(kids);
    free(keys);
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

#define CHECK(cond) do { if (!(cond)) { err = CT_EINTERNAL; goto done; } } while (0)

int ct_check(const ct_tree *T)
{
    int n = T->nstates;
    int err = CT_OK;
    int *buf = NULL, *buf2 = NULL;
    int i, j, s, t, r, c, k, count, steps, max_dep, ei, ek;

    CHECK(n >= 1 && n <= T->capacity);
    CHECK(T->broken == 0);
    CHECK(T->alphabet_size >= 0);
    CHECK(!T->dense || T->alphabet_size > 0);

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
            CHECK(T->dep[s] == T->dep[T->parent[s]] + T->tend[s] - T->tstart[s]);
            if (T->dep[s] > max_dep)
                max_dep = T->dep[s];
        }
    }

    /* the children: a child of s has parent s and its slot is the first
     * letter of its label; a node has children if and only if it is not a
     * leaf; every node but the root is found under its parent by its first
     * letter, and there are n - 1 children in all, so each one exactly once */
    count = 0;
    for (s = 0; s < n; s++) {
        int d = 0;
        if (T->dense) {
            for (c = 0; c < T->alphabet_size; c++) {
                t = T->trans[s * T->alphabet_size + c];
                if (t == -1)
                    continue;
                CHECK(t >= 1 && t < n && T->parent[t] == s);
                CHECK(letter_of(T, T->tword[t], T->tstart[t]) == c);
                d++;
            }
        } else {
            steps = 0;
            for (t = T->fchild[s]; t != -1; t = T->nsib[t]) {
                CHECK(++steps < n);
                CHECK(t >= 1 && t < n && T->parent[t] == s);
                d++;
            }
        }
        CHECK(s == 0 || (T->tend[s] == -1) == (d == 0));
        count += d;
    }
    CHECK(count == n - 1);
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
    return err;
}
