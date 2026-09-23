/*
 * Tests of the conjugate trees of combisurf/src/conjugate_tree.c, in C.
 *
 * Usage: test_conjugate_tree [seed]
 *
 * Exits 0 and prints one summary line when every case passes; otherwise
 * prints the failing case and exits 1. The random cases are reproducible from
 * the seed (1 by default).
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
#include "../../combisurf/src/conjugate_tree.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* harness                                                             */
/* ------------------------------------------------------------------ */

static const char *current_case = "";
static long num_checks = 0;

#define CHECK(cond) do { \
        num_checks++; \
        if (!(cond)) { \
            fprintf(stderr, "FAIL [%s] %s:%d: %s\n", current_case, __FILE__, __LINE__, #cond); \
            exit(1); \
        } \
    } while (0)

/* CHECK that f returns the code expected, printing both when not */
#define CHECK_CODE(f, expected) do { \
        int code_ = (f); \
        num_checks++; \
        if (code_ != (expected)) { \
            fprintf(stderr, "FAIL [%s] %s:%d: %s returned %d (%s), expected %d\n", \
                    current_case, __FILE__, __LINE__, #f, code_, ct_strerror(code_), (int) (expected)); \
            exit(1); \
        } \
    } while (0)

static const char *layout_names[] = {"sparse", "dense", "rows"};

/* The promotion threshold of the rows layout in the tests: small, so that
 * small trees have rows, and changed by the cases that vary it. */
static int test_promote = 2;

/* ct_init, with the promotion threshold of the tests for the rows layout. */
static int init_tree(ct_tree *T, int alphabet, int reserve, int layout)
{
    int err = ct_init(T, alphabet, reserve, layout);
    if (err == CT_OK && T->layout == CT_LAYOUT_ROWS)
        T->promote = test_promote;
    return err;
}

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (p == NULL) {
        fprintf(stderr, "out of memory\n");
        exit(2);
    }
    return p;
}

/* splitmix64 */
static uint64_t rng_state;

static uint64_t rng_next(void)
{
    uint64_t z = (rng_state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* uniform in {0, ..., n - 1}, n > 0 */
static int rnd(int n)
{
    return (int) (rng_next() % (uint64_t) n);
}

/* ------------------------------------------------------------------ */
/* brute force model                                                   */
/* ------------------------------------------------------------------ */

/* The primitive roots stored in the tree, in the order of insertion. */
typedef struct {
    int n;
    int cap;
    int **w;
    int *len;
} model;

static void model_init(model *M)
{
    M->n = 0;
    M->cap = 0;
    M->w = NULL;
    M->len = NULL;
}

static void model_free(model *M)
{
    int j;
    for (j = 0; j < M->n; j++)
        free(M->w[j]);
    free(M->w);
    free(M->len);
    model_init(M);
}

static void model_push(model *M, const int *w, int len)
{
    if (M->n == M->cap) {
        M->cap = M->cap ? 2 * M->cap : 8;
        M->w = (int **) realloc(M->w, (size_t) M->cap * sizeof(int *));
        M->len = (int *) realloc(M->len, (size_t) M->cap * sizeof(int));
        if (M->w == NULL || M->len == NULL) {
            fprintf(stderr, "out of memory\n");
            exit(2);
        }
    }
    M->w[M->n] = (int *) xmalloc((size_t) len * sizeof(int));
    memcpy(M->w[M->n], w, (size_t) len * sizeof(int));
    M->len[M->n] = len;
    M->n++;
}

/* The smallest period of w that divides its length. */
static int primitive_period(const int *w, int len)
{
    int p, j;
    for (p = 1; p < len; p++) {
        if (len % p)
            continue;
        for (j = p; j < len && w[j] == w[j - p]; j++)
            ;
        if (j == len)
            return p;
    }
    return len;
}

/* Whether u and v, of the same length, are conjugate. */
static int conjugate(const int *u, const int *v, int len)
{
    int k, j;
    for (k = 0; k < len; k++) {
        for (j = 0; j < len && u[(k + j) % len] == v[j]; j++)
            ;
        if (j == len)
            return 1;
    }
    return 0;
}

/*
 * What ct_process must return on w: the exponent of w if its primitive root
 * is not conjugate to a stored word (and then the root is stored), minus the
 * index of that stored word otherwise.
 */
static int model_process(model *M, const int *w, int len)
{
    int p = primitive_period(w, len);
    int j;
    for (j = 0; j < M->n; j++)
        if (M->len[j] == p && conjugate(w, M->w[j], p))
            return -j;
    model_push(M, w, p);
    return len / p;
}

static int model_letters(const model *M)
{
    int j, total = 0;
    for (j = 0; j < M->n; j++)
        total += M->len[j];
    return total;
}

/* ------------------------------------------------------------------ */
/* helpers on trees                                                    */
/* ------------------------------------------------------------------ */

static int num_leaves(const ct_tree *T)
{
    int s, count = 0;
    for (s = 1; s < T->nstates; s++)
        count += T->tend[s] == -1;
    return count;
}

/* The words of T are those of M. */
static void check_words(const ct_tree *T, const model *M)
{
    int i;
    CHECK(T->nwords == M->n);
    for (i = 0; i < M->n; i++) {
        CHECK(T->wlen[i] == M->len[i]);
        CHECK(memcmp(T->wbuf + T->wstart[i], M->w[i], (size_t) M->len[i] * sizeof(int)) == 0);
    }
    CHECK(num_leaves(T) == model_letters(M));
}

/* The fields that a failed call must leave alone. */
typedef struct {
    int nstates, nwords, wbuf_size, max_letter, broken;
    int capacity, wbuf_capacity, words_capacity;
} snapshot;

static snapshot snap(const ct_tree *T)
{
    snapshot S;
    S.nstates = T->nstates;
    S.nwords = T->nwords;
    S.wbuf_size = T->wbuf_size;
    S.max_letter = T->max_letter;
    S.broken = T->broken;
    S.capacity = T->capacity;
    S.wbuf_capacity = T->wbuf_capacity;
    S.words_capacity = T->words_capacity;
    return S;
}

static int same_content(snapshot a, snapshot b)
{
    return a.nstates == b.nstates && a.nwords == b.nwords &&
           a.wbuf_size == b.wbuf_size && a.max_letter == b.max_letter &&
           a.broken == b.broken;
}

static int same_snapshot(snapshot a, snapshot b)
{
    return same_content(a, b) && a.capacity == b.capacity &&
           a.wbuf_capacity == b.wbuf_capacity && a.words_capacity == b.words_capacity;
}

/* Whether A and B hold the same tree, with the same node indices. */
static int trees_equal(const ct_tree *A, const ct_tree *B)
{
    size_t n = (size_t) A->nstates, sz = sizeof(int);
    if (A->alphabet_size != B->alphabet_size || A->layout != B->layout ||
        A->max_letter != B->max_letter || A->broken != B->broken ||
        A->nwords != B->nwords || A->wbuf_size != B->wbuf_size ||
        A->nstates != B->nstates)
        return 0;
    /* no word, no word buffer (memcmp must not see NULL) */
    if (A->nwords && (memcmp(A->wbuf, B->wbuf, (size_t) A->wbuf_size * sz) ||
                      memcmp(A->wstart, B->wstart, (size_t) A->nwords * sz) ||
                      memcmp(A->wlen, B->wlen, (size_t) A->nwords * sz)))
        return 0;
    if (memcmp(A->dep, B->dep, n * sz) || memcmp(A->sl, B->sl, n * sz) ||
        memcmp(A->parent, B->parent, n * sz) || memcmp(A->tword, B->tword, n * sz) ||
        memcmp(A->tstart, B->tstart, n * sz) || memcmp(A->tend, B->tend, n * sz))
        return 0;
    if (A->layout == CT_LAYOUT_DENSE)
        return memcmp(A->trans, B->trans, n * (size_t) A->alphabet_size * sz) == 0;
    if (memcmp(A->fchild, B->fchild, n * sz) || memcmp(A->nsib, B->nsib, n * sz))
        return 0;
    if (A->layout == CT_LAYOUT_ROWS)
        return A->promote == B->promote && A->nrows == B->nrows &&
               memcmp(A->flet, B->flet, n * sz) == 0 && memcmp(A->row, B->row, n * sz) == 0 &&
               (A->nrows == 0 || memcmp(A->rows, B->rows, (size_t) A->nrows * (size_t) A->alphabet_size * sz) == 0);
    return 1;
}

/* Insert w, compare with the model and check the tree. */
static void process_checked(ct_tree *T, model *M, const int *w, int len)
{
    int result = 0x7fff, expected;
    expected = model_process(M, w, len);
    CHECK_CODE(ct_process(T, w, len, &result), CT_OK);
    CHECK(result == expected);
    CHECK_CODE(ct_check(T), CT_OK);
    check_words(T, M);
}

/* ------------------------------------------------------------------ */
/* ct_sorted_leaves against a brute force                              */
/* ------------------------------------------------------------------ */

static const ct_tree *cmp_tree;
static const int *cmp_order, *cmp_pivot;
static int cmp_n;

/*
 * Compare two leaves by their infinite periodic words u and v: at the first
 * position m where they differ, the branching node reads u[0..m-1], and the
 * children are ordered by order[] of u[m], shifted by pivot[u[m - 1]] below
 * the root. Two leaves differ within the sum of the two periods (Fine and
 * Wilf), since the stored words are primitive and pairwise not conjugate.
 */
static int cmp_leaves(const void *a, const void *b)
{
    const ct_tree *T = cmp_tree;
    int s = *(const int *) a, t = *(const int *) b;
    int i, k, j, l, m, bound, x = 0, y = 0, kx, ky;
    if (s == t)
        return 0;
    if (ct_leaf_as_conjugate(T, s, &i, &k) || ct_leaf_as_conjugate(T, t, &j, &l)) {
        fprintf(stderr, "FAIL [%s]: not a leaf\n", current_case);
        exit(1);
    }
    bound = T->wlen[i] + T->wlen[j];
    for (m = 0; m < bound; m++) {
        x = ct_letter(T, i, k + m);
        y = ct_letter(T, j, l + m);
        if (x != y)
            break;
    }
    if (m == bound) {
        fprintf(stderr, "FAIL [%s]: leaves %d and %d read the same word\n", current_case, s, t);
        exit(1);
    }
    kx = cmp_order[x];
    ky = cmp_order[y];
    if (m > 0) {
        int base = cmp_pivot[ct_letter(T, i, k + m - 1)];
        kx = (kx - base + cmp_n) % cmp_n;
        ky = (ky - base + cmp_n) % cmp_n;
    }
    return kx < ky ? -1 : 1;
}

static void random_permutation(int *p, int n)
{
    int j, r, tmp;
    for (j = 0; j < n; j++)
        p[j] = j;
    for (j = n - 1; j > 0; j--) {
        r = rnd(j + 1);
        tmp = p[j];
        p[j] = p[r];
        p[r] = tmp;
    }
}

/* ct_sorted_leaves with a random order and pivot, against the brute force,
 * and its refusal of an order that is not a permutation. */
static void check_sorted_leaves(const ct_tree *T)
{
    int n = (T->max_letter + 1 > T->alphabet_size ? T->max_letter + 1 : T->alphabet_size) + rnd(3);
    int *order, *pivot, *out, *expected, *seen;
    int j, num, leaves, s;
    if (n == 0)
        n = 1;
    order = (int *) xmalloc((size_t) n * sizeof(int));
    pivot = (int *) xmalloc((size_t) n * sizeof(int));
    out = (int *) xmalloc((size_t) T->nstates * sizeof(int));
    expected = (int *) xmalloc((size_t) T->nstates * sizeof(int));
    seen = (int *) xmalloc((size_t) T->nstates * sizeof(int));
    random_permutation(order, n);
    for (j = 0; j < n; j++)
        pivot[j] = rnd(n);

    num = -12345;
    CHECK_CODE(ct_sorted_leaves(T, order, pivot, n, out, &num), CT_OK);
    leaves = 0;
    for (s = 1; s < T->nstates; s++)
        if (T->tend[s] == -1)
            expected[leaves++] = s;
    CHECK(num == leaves);
    for (s = 0; s < T->nstates; s++)
        seen[s] = 0;
    for (j = 0; j < num; j++) {
        CHECK(out[j] > 0 && out[j] < T->nstates && T->tend[out[j]] == -1);
        CHECK(!seen[out[j]]);
        seen[out[j]] = 1;
    }
    cmp_tree = T;
    cmp_order = order;
    cmp_pivot = pivot;
    cmp_n = n;
    qsort(expected, (size_t) leaves, sizeof(int), cmp_leaves);
    CHECK(memcmp(out, expected, (size_t) leaves * sizeof(int)) == 0);

    /* not a permutation: a repeated value, then a value out of range */
    if (n >= 2) {
        int a = rnd(n), b = (a + 1 + rnd(n - 1)) % n, saved = order[b];
        order[b] = order[a];
        num = -12345;
        CHECK_CODE(ct_sorted_leaves(T, order, pivot, n, out, &num), CT_EINVALID);
        CHECK(num == -12345);
        order[b] = saved;
    }
    j = rnd(n);
    order[j] = rnd(2) ? n : -1;
    num = -12345;
    CHECK_CODE(ct_sorted_leaves(T, order, pivot, n, out, &num), CT_EINVALID);
    CHECK(num == -12345);
    /* a pivot out of range, and n not above every letter */
    random_permutation(order, n);
    pivot[rnd(n)] = n;
    CHECK_CODE(ct_sorted_leaves(T, order, pivot, n, out, &num), CT_EINVALID);
    CHECK(num == -12345);
    if (T->max_letter >= 0) {
        for (j = 0; j < n; j++)
            pivot[j] = 0;
        random_permutation(order, T->max_letter);
        CHECK_CODE(ct_sorted_leaves(T, order, pivot, T->max_letter, out, &num), CT_EINVALID);
        CHECK(num == -12345);
    }

    free(order);
    free(pivot);
    free(out);
    free(expected);
    free(seen);
}

/* ------------------------------------------------------------------ */
/* 1 and 3: random trees against the brute force                       */
/* ------------------------------------------------------------------ */

/*
 * Write to w a random word of length at most maxlen: a fresh word, a power of
 * one, or a conjugate of a power of a stored word (so that ct_process finds
 * it). The letters are drawn from a random sub-alphabet of {0, ..., A - 1} of
 * size 1 to A, small ones being more likely, so that some trees are deep.
 */
static int random_word(int *w, int maxlen, int A, const model *M)
{
    int len, base, j, rep, i, k, sub;
    int kind = rnd(10);
    if (kind < 3 && M->n > 0) {
        /* a conjugate of a power of a stored word */
        i = rnd(M->n);
        rep = 1 + rnd(3);
        if (M->len[i] * rep > maxlen)
            rep = 1;
        if (M->len[i] <= maxlen) {
            len = M->len[i] * rep;
            k = rnd(M->len[i]);
            for (j = 0; j < len; j++)
                w[j] = M->w[i][(k + j) % M->len[i]];
            return len;
        }
    }
    sub = rnd(3) ? 1 + rnd(A) : 1 + rnd(A < 3 ? A : 3);
    base = 1 + rnd(rnd(4) ? (maxlen < 12 ? maxlen : 12) : maxlen);
    for (j = 0; j < base; j++)
        w[j] = rnd(sub);
    if (sub < A) {
        /* move the sub-alphabet to a random place of the alphabet */
        int shift = rnd(A - sub + 1);
        for (j = 0; j < base; j++)
            w[j] += shift;
    }
    rep = kind < 6 ? 1 : 1 + rnd(4);
    len = base * rep;
    if (len > maxlen)
        len = base;
    for (j = base; j < len; j++)
        w[j] = w[j - base];
    return len;
}

/*
 * One random tree: alphabet (0 for unknown), letter bound A (A = alphabet
 * when it is known), layout (negative for the default), reserve, and up to
 * nw words of length at most maxlen.
 */
static void random_tree(int alphabet, int A, int layout, int reserve, int nw, int maxlen)
{
    ct_tree T;
    model M;
    int *w = (int *) xmalloc((size_t) maxlen * sizeof(int));
    int j, len;
    CHECK_CODE(init_tree(&T, alphabet, reserve, layout), CT_OK);
    model_init(&M);
    for (j = 0; j < nw; j++) {
        len = random_word(w, maxlen, A, &M);
        process_checked(&T, &M, w, len);
    }
    check_sorted_leaves(&T);
    ct_free(&T);
    model_free(&M);
    free(w);
}

static void test_random(void)
{
    int trial, alphabet, A, layout;
    current_case = "random small trees";
    for (trial = 0; trial < 5000; trial++) {
        /* the alphabet: unknown, small, and past the dense threshold */
        alphabet = trial % 3 == 0 ? 0 : 1 + rnd(trial % 3 == 1 ? 32 : 200);
        A = alphabet ? alphabet : 1 + rnd(10);
        /* the layout: any, dense only when the alphabet is known */
        layout = rnd(4) - 1;
        if (layout == CT_LAYOUT_DENSE && !alphabet)
            layout = CT_LAYOUT_DEFAULT;
        test_promote = 1 + rnd(4);
        random_tree(alphabet, A, layout, rnd(2) ? 0 : rnd(64), 1 + rnd(15), 1 + rnd(12));
    }

    /* long words past the dense threshold, in every layout but dense */
    for (layout = CT_LAYOUT_SPARSE; layout <= CT_LAYOUT_ROWS; layout++) {
        if (layout == CT_LAYOUT_DENSE)
            continue;
        current_case = layout_names[layout];
        for (trial = 0; trial < 60; trial++) {
            alphabet = trial % 2 ? 0 : 33 + rnd(168);
            A = alphabet ? alphabet : 2 + rnd(199);
            test_promote = 1 + rnd(16);
            random_tree(alphabet, A, layout, 0, 1 + rnd(8), 300);
        }
    }
    test_promote = 2;

    current_case = "random long words, dense past 32 letters";
    for (trial = 0; trial < 40; trial++) {
        alphabet = 33 + rnd(80);
        random_tree(alphabet, alphabet, CT_LAYOUT_DENSE, 0, 1 + rnd(8), 300);
    }

    current_case = "random long words, small alphabets";
    for (trial = 0; trial < 60; trial++) {
        alphabet = trial % 3 ? 1 + rnd(4) : 0;
        A = alphabet ? alphabet : 1 + rnd(4);
        layout = trial % 4 - 1;
        if (layout == CT_LAYOUT_DENSE && !alphabet)
            layout = CT_LAYOUT_DEFAULT;
        random_tree(alphabet, A, layout, 0, 1 + rnd(8), 300);
    }
}

/* ------------------------------------------------------------------ */
/* 2: hand-picked words                                                */
/* ------------------------------------------------------------------ */

/* Insert w into T expecting result, keeping the model in step. */
static void expect(ct_tree *T, model *M, const int *w, int len, int result)
{
    int r = 0x7fff;
    CHECK(model_process(M, w, len) == result);
    CHECK_CODE(ct_process(T, w, len, &r), CT_OK);
    CHECK(r == result);
    CHECK_CODE(ct_check(T), CT_OK);
    check_words(T, M);
}

/* ct_process on w fails with code and leaves T alone. */
static void expect_error(ct_tree *T, const int *w, int len, int code)
{
    snapshot before = snap(T);
    int r = 0x7fff;
    CHECK_CODE(ct_process(T, w, len, &r), code);
    CHECK(r == 0x7fff);
    CHECK(same_snapshot(before, snap(T)));
    CHECK_CODE(ct_check(T), CT_OK);
}

/* Write the Lyndon words of length at most maxlen over k letters to out, one
 * after the other, their lengths to lens; return their number (Duval). */
static int lyndon_words(int k, int maxlen, int *out, int *lens)
{
    int w[64];
    int len = 1, num = 0, pos = 0, j;
    w[0] = 0;
    for (;;) {
        memcpy(out + pos, w, (size_t) len * sizeof(int));
        lens[num++] = len;
        pos += len;
        /* extend to maxlen by repetition, then strip trailing maximal
         * letters and increment */
        for (j = len; j < maxlen; j++)
            w[j] = w[j - len];
        len = maxlen;
        while (len > 0 && w[len - 1] == k - 1)
            len--;
        if (len == 0)
            return num;
        w[len - 1]++;
    }
}

/* The hand-picked cases, in the given layout and over the given alphabet
 * (0, or at least 16). */
static void test_hand_picked(int layout, int alphabet)
{
    ct_tree T;
    model M;
    int w[64], lyn[1024], lens[64];
    int k, j, i, num, pos;

    current_case = layout_names[layout];

    /* a single letter */
    CHECK_CODE(init_tree(&T, alphabet, 0, layout), CT_OK);
    model_init(&M);
    w[0] = 3;
    expect(&T, &M, w, 1, 1);
    CHECK(num_leaves(&T) == 1 && T.max_letter == 3);
    expect(&T, &M, w, 1, 0);
    ct_free(&T);
    model_free(&M);

    /* a^k, each in a fresh tree and then all in one */
    for (k = 1; k <= 6; k++) {
        CHECK_CODE(init_tree(&T, alphabet, 0, layout), CT_OK);
        model_init(&M);
        for (j = 0; j < k; j++)
            w[j] = 5;
        expect(&T, &M, w, k, k);
        CHECK(T.nwords == 1 && T.wlen[0] == 1 && num_leaves(&T) == 1);
        ct_free(&T);
        model_free(&M);
    }
    CHECK_CODE(init_tree(&T, alphabet, 0, layout), CT_OK);
    model_init(&M);
    for (k = 1; k <= 6; k++)
        expect(&T, &M, w, k, k == 1 ? 1 : 0);
    ct_free(&T);
    model_free(&M);

    /* (ab)^3 */
    CHECK_CODE(init_tree(&T, alphabet, 0, layout), CT_OK);
    model_init(&M);
    for (j = 0; j < 6; j++)
        w[j] = j % 2 ? 7 : 2;
    expect(&T, &M, w, 6, 3);
    CHECK(T.nwords == 1 && T.wlen[0] == 2 && num_leaves(&T) == 2);
    ct_free(&T);
    model_free(&M);

    /* a word and all its conjugates; a word then one of its conjugates,
     * after another word so that the index is not 0 */
    {
        int u[] = {0, 1, 1, 0, 2, 1, 0};
        int l = 7;
        CHECK_CODE(init_tree(&T, alphabet, 0, layout), CT_OK);
        model_init(&M);
        w[0] = 4;
        expect(&T, &M, w, 1, 1);
        expect(&T, &M, u, l, 1);
        for (k = 0; k < l; k++) {
            for (j = 0; j < l; j++)
                w[j] = u[(k + j) % l];
            expect(&T, &M, w, l, -1);
            /* and a power of it */
            for (j = l; j < 2 * l; j++)
                w[j] = w[j - l];
            expect(&T, &M, w, 2 * l, -1);
        }
        CHECK(num_leaves(&T) == 1 + l);
        ct_free(&T);
        model_free(&M);
    }

    /* the Lyndon words of length at most 6 over 2 letters: each is new and
     * primitive, and each conjugate of each is found */
    num = lyndon_words(2, 6, lyn, lens);
    CHECK(num == 23);
    CHECK_CODE(init_tree(&T, alphabet, 0, layout), CT_OK);
    model_init(&M);
    for (i = 0, pos = 0; i < num; pos += lens[i], i++)
        expect(&T, &M, lyn + pos, lens[i], 1);
    for (i = 0, pos = 0; i < num; pos += lens[i], i++) {
        for (k = 0; k < lens[i]; k++) {
            for (j = 0; j < lens[i]; j++)
                w[j] = lyn[pos + (k + j) % lens[i]];
            expect(&T, &M, w, lens[i], -i);
        }
    }
    /* 2 + 1 + 2 + 3 + 6 + 9 Lyndon words of lengths 1 to 6 */
    CHECK(num_leaves(&T) == 2 + 2 + 6 + 12 + 30 + 54);
    ct_free(&T);
    model_free(&M);

    /* every letter distinct */
    CHECK_CODE(init_tree(&T, alphabet, 0, layout), CT_OK);
    model_init(&M);
    for (j = 0; j < 16; j++)
        w[j] = (7 * j) % 16;
    expect(&T, &M, w, 16, 1);
    CHECK(T.nstates == 17);
    ct_free(&T);
    model_free(&M);

    /* errors, on an empty tree and on a tree with words */
    CHECK_CODE(init_tree(&T, alphabet, 0, layout), CT_OK);
    model_init(&M);
    for (k = 0; k < 2; k++) {
        w[0] = 1;
        w[1] = -1;
        expect_error(&T, w, 0, CT_EEMPTY);
        expect_error(&T, w, -1, CT_EINVALID);
        expect_error(&T, w, 2, CT_ENEGATIVE);
        w[0] = -3;
        expect_error(&T, w, 1, CT_ENEGATIVE);
        if (alphabet) {
            w[0] = 1;
            w[1] = alphabet;
            expect_error(&T, w, 2, CT_EALPHABET);
            w[1] = INT_MAX;
            expect_error(&T, w, 2, CT_EALPHABET);
        }
        w[0] = 0;
        w[1] = 1;
        w[2] = 0;
        expect(&T, &M, w, 3, k ? 0 : 1);
    }

    /* words that are not reduced in a free group are ordinary words */
    w[0] = 0;
    w[1] = 1;
    expect(&T, &M, w, 2, 1);
    w[0] = 0;
    w[1] = 0;
    expect(&T, &M, w, 2, 2);
    ct_free(&T);
    model_free(&M);
}

/* ------------------------------------------------------------------ */
/* 4: corruption                                                       */
/* ------------------------------------------------------------------ */

/* What a node reads: its label for an internal node, its conjugate for a
 * leaf. Two trees that agree on this for every node, with the same links,
 * are the same tree written differently, and ct_check cannot tell them apart. */
static int *label_letters(const ct_tree *T, int *len)
{
    int s, j, total = 0, pos = 0, i, k;
    int *out;
    for (s = 1; s < T->nstates; s++)
        total += T->tend[s] == -1 ? 2 : 1 + T->tend[s] - T->tstart[s];
    out = (int *) xmalloc((size_t) total * sizeof(int));
    for (s = 1; s < T->nstates; s++) {
        if (T->tend[s] == -1) {
            if (ct_leaf_as_conjugate(T, s, &i, &k))
                i = k = -1;
            out[pos++] = i;
            out[pos++] = k;
        } else {
            out[pos++] = T->tend[s] - T->tstart[s];
            for (j = T->tstart[s]; j < T->tend[s]; j++)
                out[pos++] = ct_letter(T, T->tword[s], j);
        }
    }
    *len = total;
    return out;
}

static int corrupt_tries, corrupt_equivalent;

/*
 * Set *field to value, expect ct_check to fail, and restore it. For the
 * fields of a label, a value that describes the same labels (another
 * occurrence of the same letters) is not a corruption and is skipped.
 */
static void corrupt(ct_tree *T, int *field, int value, int label, const char *name, int s)
{
    int saved = *field;
    int code, len0 = 0, len1 = 0, *before = NULL, *after = NULL, same = 0;
    if (value == saved)
        return;
    if (label)
        before = label_letters(T, &len0);
    *field = value;
    code = ct_check(T);
    if (code != CT_EINTERNAL && label) {
        /* reading the corrupted labels is only safe once they are within
         * range, which a CT_OK of ct_check guarantees */
        after = label_letters(T, &len1);
        same = len0 == len1 && memcmp(before, after, (size_t) len0 * sizeof(int)) == 0;
    }
    *field = saved;
    free(before);
    free(after);
    if (same) {
        corrupt_equivalent++;
        return;
    }
    corrupt_tries++;
    if (code != CT_EINTERNAL) {
        fprintf(stderr, "FAIL [%s]: %s[%d] = %d (was %d) not detected: %d\n",
                current_case, name, s, value, saved, code);
        exit(1);
    }
}

/* Values to try for a field that holds x, in a tree of n nodes. */
static void values_for(int x, int n, int wl, int *v)
{
    v[0] = x + 1;
    v[1] = x - 1;
    v[2] = -1;
    v[3] = -2;
    v[4] = n;
    v[5] = x + wl;
    v[6] = 0;
    v[7] = INT_MAX;
}

#define NVALUES 8

static void corrupt_node(ct_tree *T, int s)
{
    int v[NVALUES];
    int j, c, n = T->nstates, wl = T->wlen[T->tword[s]];
    int leaf = s > 0 && T->tend[s] == -1;
    int *fields[6];
    const char *names[6] = {"dep", "sl", "parent", "tword", "tstart", "tend"};
    fields[0] = T->dep;
    fields[1] = T->sl;
    fields[2] = T->parent;
    fields[3] = T->tword;
    fields[4] = T->tstart;
    fields[5] = T->tend;
    for (j = 0; j < 6; j++) {
        /* the depth and the suffix link of a leaf are not used */
        if (leaf && j < 2)
            continue;
        values_for(fields[j][s], n, wl, v);
        for (c = 0; c < NVALUES; c++)
            corrupt(T, fields[j] + s, v[c], j >= 3, names[j], s);
        if (j == 3)
            corrupt(T, fields[j] + s, (fields[j][s] + 1) % T->nwords, 1, names[j], s);
    }
    if (T->layout == CT_LAYOUT_DENSE) {
        int *row = T->trans + (size_t) s * (size_t) T->alphabet_size;
        for (c = 0; c < T->alphabet_size; c++) {
            values_for(row[c], n, 0, v);
            for (j = 0; j < 5; j++)
                corrupt(T, row + c, v[j], 0, "trans", s);
            corrupt(T, row + c, s, 0, "trans", s);
            corrupt(T, row + c, 1 + rnd(n - 1), 0, "trans", s);
        }
    } else {
        values_for(T->fchild[s], n, 0, v);
        for (j = 0; j < 5; j++)
            corrupt(T, T->fchild + s, v[j], 0, "fchild", s);
        corrupt(T, T->fchild + s, s, 0, "fchild", s);
        corrupt(T, T->fchild + s, 1 + rnd(n - 1), 0, "fchild", s);
        if (s > 0) {
            values_for(T->nsib[s], n, 0, v);
            for (j = 0; j < 5; j++)
                corrupt(T, T->nsib + s, v[j], 0, "nsib", s);
            corrupt(T, T->nsib + s, s, 0, "nsib", s);
            corrupt(T, T->nsib + s, 1 + rnd(n - 1), 0, "nsib", s);
        }
    }
    if (T->flet) {
        values_for(T->flet[s], n, 0, v);
        for (j = 0; j < 5; j++)
            corrupt(T, T->flet + s, v[j], 0, "flet", s);
        corrupt(T, T->flet + s, rnd(T->alphabet_size), 0, "flet", s);
    }
    if (T->layout == CT_LAYOUT_ROWS) {
        values_for(T->row[s], T->nrows, 0, v);
        for (j = 0; j < 7; j++)
            corrupt(T, T->row + s, v[j], 0, "row", s);
        if (T->row[s] != -1) {
            int *row = T->rows + (size_t) T->row[s] * (size_t) T->alphabet_size;
            for (c = 0; c < T->alphabet_size; c++) {
                values_for(row[c], n, 0, v);
                for (j = 0; j < 5; j++)
                    corrupt(T, row + c, v[j], 0, "rows", s);
                corrupt(T, row + c, s, 0, "rows", s);
                corrupt(T, row + c, 1 + rnd(n - 1), 0, "rows", s);
            }
        }
    }
}

static void test_corruption(int layout)
{
    ct_tree T;
    model M;
    int w1[] = {0, 2, 0, 3, 1, 0, 2}, w2[] = {2, 0, 2, 0, 1}, w3[] = {3, 3, 1};
    int s, j, x;
    current_case = layout_names[layout];
    CHECK_CODE(init_tree(&T, 4, 0, layout), CT_OK);
    model_init(&M);
    process_checked(&T, &M, w1, 7);
    process_checked(&T, &M, w2, 5);
    process_checked(&T, &M, w3, 3);

    /* every node, the root included */
    for (s = 0; s < T.nstates; s++)
        corrupt_node(&T, s);
    if (T.layout == CT_LAYOUT_ROWS) {
        /* the root has 4 children, so it has a row */
        CHECK(T.nrows > 0 && T.row[0] != -1);
        corrupt(&T, &T.nrows, T.nrows - 1, 0, "nrows", 0);
        if (T.nrows < T.rows_capacity)
            corrupt(&T, &T.nrows, T.nrows + 1, 0, "nrows", 0);
        corrupt(&T, &T.promote, 0, 0, "promote", 0);
        corrupt(&T, &T.promote, 5, 0, "promote", 0);
    }

    /* the words */
    for (j = 0; j < T.wbuf_size; j++) {
        x = T.wbuf[j];
        corrupt(&T, T.wbuf + j, (x + 1) % 4, 0, "wbuf", j);
        corrupt(&T, T.wbuf + j, -1, 0, "wbuf", j);
        corrupt(&T, T.wbuf + j, 4, 0, "wbuf", j);
    }
    for (j = 0; j < T.nwords; j++) {
        corrupt(&T, T.wlen + j, T.wlen[j] + 1, 0, "wlen", j);
        corrupt(&T, T.wlen + j, 0, 0, "wlen", j);
        corrupt(&T, T.wstart + j, T.wstart[j] + 1, 0, "wstart", j);
    }
    corrupt(&T, &T.nwords, T.nwords - 1, 0, "nwords", 0);
    corrupt(&T, &T.wbuf_size, T.wbuf_size - 1, 0, "wbuf_size", 0);
    corrupt(&T, &T.max_letter, 2, 0, "max_letter", 0);
    corrupt(&T, &T.nstates, T.nstates - 1, 0, "nstates", 0);
    corrupt(&T, &T.broken, CT_EINTERNAL, 0, "broken", 0);
    corrupt(&T, &T.layout, CT_LAYOUT_ROWS + 1, 0, "layout", 0);
    corrupt(&T, &T.layout, -1, 0, "layout", 0);

    /* nothing was left corrupted */
    CHECK_CODE(ct_check(&T), CT_OK);
    ct_free(&T);
    model_free(&M);
}

/* ------------------------------------------------------------------ */
/* 5: odd states                                                       */
/* ------------------------------------------------------------------ */

/* Every call on a tree that holds no node fails cleanly. */
static void check_unusable(ct_tree *T)
{
    int w[] = {0, 1}, order[] = {0, 1}, pivot[] = {1, 0}, out[4];
    int r = 0x7fff, num = -12345, i = -1, k = -1, s = 0;
    CHECK_CODE(ct_process(T, w, 2, &r), CT_EINVALID);
    CHECK(r == 0x7fff);
    CHECK_CODE(ct_reserve(T, 1, 2), CT_EINVALID);
    CHECK_CODE(ct_sorted_leaves(T, order, pivot, 2, out, &num), CT_EINVALID);
    CHECK(num == -12345);
    CHECK_CODE(ct_check(T), CT_EINTERNAL);
    CHECK_CODE(ct_leaf_as_conjugate(T, 0, &i, &k), CT_EINVALID);
    CHECK_CODE(ct_leaf_as_conjugate(T, 1, &i, &k), CT_EINVALID);
    CHECK(i == -1 && k == -1);
    CHECK(ct_size(T) == CT_EINVALID);
    CHECK_CODE(ct_child(T, 0, 0), CT_EINVALID);
    CHECK_CODE(ct_letter(T, 0, 0), CT_EINVALID);
    k = 0;
    CHECK_CODE(ct_canonize(T, &s, 0, &k, 1), CT_EINVALID);
    s = -1;
    CHECK_CODE(ct_canonize(T, &s, 0, &k, 1), CT_EINVALID);
    CHECK(s == -1 && k == 0);
    ct_free(T);
    ct_free(T);
}

static void test_odd_states(void)
{
    ct_tree T;
    current_case = "odd states";

    memset(&T, 0, sizeof(ct_tree));
    ct_free(&T);
    check_unusable(&T);

    /* failed ct_init: invalid arguments, and a dense table too large */
    CHECK_CODE(ct_init(&T, -1, 0, -1), CT_EINVALID);
    check_unusable(&T);
    CHECK_CODE(ct_init(&T, 4, -1, -1), CT_EINVALID);
    check_unusable(&T);
    CHECK_CODE(ct_init(&T, 0, 0, 1), CT_EINVALID);
    check_unusable(&T);
    CHECK_CODE(ct_init(&T, 4, 0, CT_LAYOUT_ROWS + 1), CT_EINVALID);
    check_unusable(&T);
    CHECK_CODE(ct_init(&T, INT_MAX / 2, 10, 1), CT_ETOOLARGE);
    check_unusable(&T);

    /* ct_free twice on a used tree */
    CHECK_CODE(ct_init(&T, 3, 0, -1), CT_OK);
    ct_free(&T);
    ct_free(&T);

    /* the layouts that ct_init chooses */
    CHECK_CODE(ct_init(&T, 0, 0, CT_LAYOUT_DEFAULT), CT_OK);
    CHECK(T.layout == CT_LAYOUT_ROWS);
    ct_free(&T);
    CHECK_CODE(ct_init(&T, CT_DENSE_MAX_ALPHABET, 0, CT_LAYOUT_DEFAULT), CT_OK);
    CHECK(T.layout == CT_LAYOUT_DENSE);
    ct_free(&T);
    CHECK_CODE(ct_init(&T, CT_DENSE_MAX_ALPHABET + 1, 0, CT_LAYOUT_DEFAULT), CT_OK);
    CHECK(T.layout == CT_LAYOUT_ROWS && T.promote == CT_PROMOTE);
    ct_free(&T);

    /* without an alphabet, the rows layout never promotes a node */
    {
        int w[] = {0, 1, 2, 3, 4, 5, 6, 7}, r;
        CHECK_CODE(init_tree(&T, 0, 0, CT_LAYOUT_ROWS), CT_OK);
        CHECK(T.layout == CT_LAYOUT_ROWS && T.promote == test_promote);
        CHECK_CODE(ct_process(&T, w, 8, &r), CT_OK);
        CHECK(T.nrows == 0 && T.row[0] == -1);
        CHECK_CODE(ct_check(&T), CT_OK);
        ct_free(&T);
    }

    /* the checked queries on a tree that holds words */
    {
        int w[] = {0, 1, 2, 0, 1, 1};
        int r, s, k, i;
        CHECK_CODE(ct_init(&T, 3, 0, -1), CT_OK);
        CHECK_CODE(ct_process(&T, w, 6, &r), CT_OK);
        CHECK_CODE(ct_letter(&T, 1, 0), CT_EINVALID);
        CHECK_CODE(ct_letter(&T, -1, 0), CT_EINVALID);
        CHECK(ct_letter(&T, 0, -1) == 1 && ct_letter(&T, 0, 8) == 2);
        CHECK_CODE(ct_child(&T, T.nstates, 0), CT_EINVALID);
        CHECK_CODE(ct_child(&T, -1, 0), CT_EINVALID);
        CHECK_CODE(ct_child(&T, 0, 3), CT_EINVALID);
        CHECK_CODE(ct_child(&T, 0, -1), CT_EINVALID);
        CHECK(ct_child(&T, 0, 2) > 0);
        CHECK_CODE(ct_leaf_as_conjugate(&T, T.nstates, &i, &k), CT_EINVALID);
        CHECK(ct_size(&T) > 0);
        /* ct_canonize: out of range, and a path that leaves the tree (a
         * leaf has no child, and the letters 1 2 then 2 are not read) */
        s = 0; k = 0;
        CHECK_CODE(ct_canonize(&T, &s, 0, &k, -1), CT_EINVALID);
        s = 0; k = -1;
        CHECK_CODE(ct_canonize(&T, &s, 0, &k, 2), CT_EINVALID);
        s = 0; k = 3;
        CHECK_CODE(ct_canonize(&T, &s, 0, &k, 2), CT_EINVALID);
        s = T.nstates; k = 0;
        CHECK_CODE(ct_canonize(&T, &s, 0, &k, 2), CT_EINVALID);
        s = -2; k = 0;
        CHECK_CODE(ct_canonize(&T, &s, 0, &k, 2), CT_EINVALID);
        s = 0; k = 0;
        CHECK_CODE(ct_canonize(&T, &s, 1, &k, 2), CT_EINVALID);
        CHECK(s == 0 && k == 0);
        for (s = 1; s < T.nstates; s++) {
            if (T.tend[s] == -1) {
                int ss = s, kk = 0;
                CHECK_CODE(ct_canonize(&T, &ss, 0, &kk, 1), CT_EINVALID);
                CHECK(ss == s && kk == 0);
            }
        }
        /* reading the whole word from the root reaches the leaf of the
         * conjugate 0, implicitly */
        s = 0; k = 0;
        CHECK_CODE(ct_canonize(&T, &s, 0, &k, 6), CT_OK);
        CHECK(s == 0 || T.tend[s] != -1);
        s = -1; k = 0;
        CHECK_CODE(ct_canonize(&T, &s, 0, &k, 1), CT_OK);
        CHECK(s == 0 && k == 1);
        ct_free(&T);
    }
}

/* ------------------------------------------------------------------ */
/* 6: ct_reserve                                                       */
/* ------------------------------------------------------------------ */

/* Build over (alphabet, layout) the words of pre, reserve for the words of
 * post, insert them, and compare with the tree built without the
 * reservation. */
static void check_reserve(int alphabet, int layout, int npre, int npost, int maxlen)
{
    ct_tree A, B;
    int *words = (int *) xmalloc((size_t) (npre + npost) * (size_t) maxlen * sizeof(int));
    int *lens = (int *) xmalloc((size_t) (npre + npost) * sizeof(int));
    int j, total = 0, r1, r2, A_ = alphabet ? alphabet : 6;
    snapshot S;
    model M;
    model_init(&M);
    for (j = 0; j < npre + npost; j++) {
        lens[j] = random_word(words + (size_t) j * (size_t) maxlen, maxlen, A_, &M);
        model_process(&M, words + (size_t) j * (size_t) maxlen, lens[j]);
        if (j >= npre)
            total += lens[j];
    }
    model_free(&M);

    CHECK_CODE(init_tree(&A, alphabet, 0, layout), CT_OK);
    CHECK_CODE(init_tree(&B, alphabet, 0, layout), CT_OK);
    for (j = 0; j < npre; j++) {
        CHECK_CODE(ct_process(&A, words + (size_t) j * (size_t) maxlen, lens[j], &r1), CT_OK);
        CHECK_CODE(ct_process(&B, words + (size_t) j * (size_t) maxlen, lens[j], &r2), CT_OK);
    }
    CHECK_CODE(ct_reserve(&A, npost, total), CT_OK);
    CHECK_CODE(ct_check(&A), CT_OK);
    CHECK(trees_equal(&A, &B));
    S = snap(&A);
    for (j = npre; j < npre + npost; j++) {
        CHECK_CODE(ct_process(&A, words + (size_t) j * (size_t) maxlen, lens[j], &r1), CT_OK);
        CHECK_CODE(ct_process(&B, words + (size_t) j * (size_t) maxlen, lens[j], &r2), CT_OK);
        CHECK(r1 == r2);
        CHECK(A.capacity == S.capacity && A.wbuf_capacity == S.wbuf_capacity &&
              A.words_capacity == S.words_capacity);
    }
    CHECK_CODE(ct_check(&A), CT_OK);
    CHECK(trees_equal(&A, &B));
    ct_free(&A);
    ct_free(&B);
    free(words);
    free(lens);
}

/* ct_reserve fails with code and leaves the tree alone. */
static void expect_reserve_error(ct_tree *T, int words, int letters, int code)
{
    snapshot before = snap(T);
    CHECK_CODE(ct_reserve(T, words, letters), code);
    CHECK(same_snapshot(before, snap(T)));
    CHECK_CODE(ct_check(T), CT_OK);
}

static void test_reserve(void)
{
    ct_tree T, Z;
    int w[] = {0, 1, 1}, r, trial;
    current_case = "ct_reserve";
    for (trial = 0; trial < 300; trial++) {
        int alphabet = trial % 3 == 0 ? 0 : 1 + rnd(trial % 3 == 1 ? 32 : 100);
        int layout = rnd(4) - 1;
        if (layout == CT_LAYOUT_DENSE && !alphabet)
            layout = CT_LAYOUT_DEFAULT;
        test_promote = 1 + rnd(4);
        check_reserve(alphabet, layout, trial % 2 ? 0 : 1 + rnd(5), 1 + rnd(10), 1 + rnd(40));
    }
    test_promote = 2;

    for (trial = 0; trial < 6; trial++) {
        CHECK_CODE(init_tree(&T, 4, 0, trial % 3), CT_OK);
        if (trial >= 3)
            CHECK_CODE(ct_process(&T, w, 3, &r), CT_OK);
        expect_reserve_error(&T, 1, INT_MAX, CT_ETOOLARGE);
        expect_reserve_error(&T, 1, (INT_MAX - T.nstates) / 2 + 1, CT_ETOOLARGE);
        expect_reserve_error(&T, INT_MAX, INT_MAX, CT_ETOOLARGE);
        if (T.nwords)
            expect_reserve_error(&T, INT_MAX, 1, CT_ETOOLARGE);
        expect_reserve_error(&T, -1, 1, CT_EINVALID);
        expect_reserve_error(&T, 1, -1, CT_EINVALID);
        expect_reserve_error(&T, INT_MIN, INT_MIN, CT_EINVALID);
        CHECK_CODE(ct_reserve(&T, 0, 0), CT_OK);
        ct_free(&T);
    }

    /* past the dense table: nothing is inserted, whatever grew */
    CHECK_CODE(ct_init(&T, 1000, 0, 1), CT_OK);
    CHECK_CODE(ct_process(&T, w, 3, &r), CT_OK);
    {
        snapshot before = snap(&T);
        CHECK_CODE(ct_reserve(&T, 1, INT_MAX / 1000), CT_ETOOLARGE);
        CHECK(same_content(before, snap(&T)));
        CHECK_CODE(ct_check(&T), CT_OK);
        CHECK_CODE(ct_process(&T, w + 1, 2, &r), CT_OK);
        CHECK_CODE(ct_check(&T), CT_OK);
    }
    ct_free(&T);

    memset(&Z, 0, sizeof(ct_tree));
    CHECK_CODE(ct_reserve(&Z, 1, 1), CT_EINVALID);
    CHECK_CODE(ct_reserve(&Z, 0, 0), CT_EINVALID);
    ct_free(&Z);
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    uint64_t seed = 1;
    if (argc > 2) {
        fprintf(stderr, "usage: %s [seed]\n", argv[0]);
        return 2;
    }
    if (argc == 2)
        seed = strtoull(argv[1], NULL, 10);
    rng_state = seed;

    {
        int layout;
        for (layout = CT_LAYOUT_SPARSE; layout <= CT_LAYOUT_ROWS; layout++) {
            test_hand_picked(layout, 16);
            if (layout != CT_LAYOUT_DENSE)
                test_hand_picked(layout, 0);
            test_corruption(layout);
        }
    }
    test_odd_states();
    test_reserve();
    test_random();

    printf("test_conjugate_tree: seed %llu, %ld checks, %d corruptions caught "
           "(%d equivalent rewritings skipped): all passed\n",
           (unsigned long long) seed, num_checks, corrupt_tries, corrupt_equivalent);
    return 0;
}
