/*
 * Benchmark of the conjugate trees of combisurf/src/conjugate_tree.c.
 *
 * Usage: bench <word file> [rounds] [--with-inverse] [--layout <layout>]
 *              [--promote <P>]
 *
 * Reads a word file (see extra/word_files.py), builds a fresh tree over its
 * alphabet with the given layout (default, sparse, dense or rows; with
 * --promote, the rows layout gives a row to a node at P children)
 * and inserts every word with ct_process;
 * with --with-inverse, it inserts each new word followed by the inverse of
 * its primitive root (reversed, each letter h replaced by h ^ 1), as
 * geometric_intersection does. It prints the best time over the rounds of
 * one build, which a round repeats until it has lasted 10 ms, the shape of
 * the tree, and the memory it holds (from the capacities of its arrays).
 *
 * It also times, on the same sequence of sizes:
 * - the copy of the words into the word buffer of the tree, and the growth of
 *   that buffer, alone: an upper bound on what the tree would save by keeping
 *   pointers to the words of the caller instead of copying them;
 * - with --with-inverse, ct_sorted_leaves with order the identity and
 *   pivot[b] = b ^ 1, and the check that order is a permutation alone.
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
#define _POSIX_C_SOURCE 199309L

#include "../../combisurf/src/conjugate_tree.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void die(const char *msg)
{
    fprintf(stderr, "bench: %s\n", msg);
    exit(1);
}

static void *xrealloc(void *p, size_t n)
{
    p = realloc(p, n ? n : 1);
    if (p == NULL)
        die("out of memory");
    return p;
}

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) ts.tv_sec + 1e-9 * (double) ts.tv_nsec;
}

/* ------------------------------------------------------------------ */
/* word files                                                          */
/* ------------------------------------------------------------------ */

static int alphabet;
static int nwords;
static int **words;
static int *lens;
static long total_letters;

static void read_file(const char *path)
{
    FILE *f = fopen(path, "r");
    size_t cap = 0, lcap = 0;
    int c, len, *w;
    if (f == NULL)
        die("cannot open the word file");
    if (fscanf(f, "# alphabet %d", &alphabet) != 1 || alphabet <= 0)
        die("the first line must be '# alphabet <n>'");
    while ((c = fgetc(f)) != '\n' && c != EOF)
        ;
    for (;;) {
        c = fgetc(f);
        if (c == EOF)
            break;
        if (c == '#') {
            while ((c = fgetc(f)) != '\n' && c != EOF)
                ;
            continue;
        }
        if (c == '\n')
            continue;
        ungetc(c, f);
        len = 0;
        w = NULL;
        lcap = 0;
        for (;;) {
            int x;
            if (fscanf(f, "%d", &x) != 1)
                die("malformed word");
            if (x < 0 || x >= alphabet)
                die("letter outside the alphabet");
            if ((size_t) len == lcap) {
                lcap = lcap ? 2 * lcap : 64;
                w = (int *) xrealloc(w, lcap * sizeof(int));
            }
            w[len++] = x;
            while ((c = fgetc(f)) == ' ')
                ;
            if (c == '\n' || c == EOF)
                break;
            ungetc(c, f);
        }
        if ((size_t) nwords == cap) {
            cap = cap ? 2 * cap : 16;
            words = (int **) xrealloc(words, cap * sizeof(int *));
            lens = (int *) xrealloc(lens, cap * sizeof(int));
        }
        words[nwords] = w;
        lens[nwords] = len;
        nwords++;
        total_letters += len;
    }
    fclose(f);
    if (nwords == 0)
        die("no word in the file");
}

/* ------------------------------------------------------------------ */
/* the build, recording what the word buffer sees                      */
/* ------------------------------------------------------------------ */

/*
 * The calls that touch the word buffer, in order: a reservation of need
 * letters in all, then (len > 0) the copy of len letters of the word src at
 * the end of the buffer, after which the buffer holds size letters.
 */
typedef struct {
    int need;
    int len;
    int size;
    const int *src;
} wevent;

static wevent *events;
static int nevents, events_cap;

static void record(int need, int len, int size, const int *src)
{
    if (nevents == events_cap) {
        events_cap = events_cap ? 2 * events_cap : 64;
        events = (wevent *) xrealloc(events, (size_t) events_cap * sizeof(wevent));
    }
    events[nevents].need = need;
    events[nevents].len = len;
    events[nevents].size = size;
    events[nevents].src = src;
    nevents++;
}

static int *invbuf;
static int opt_layout = CT_LAYOUT_DEFAULT, opt_promote = 0;

static void process(ct_tree *T, const int *w, int len, int *result, int rec)
{
    int before = T->wbuf_size;
    if (ct_process(T, w, len, result))
        die("ct_process failed");
    if (rec)
        record(before + len, len, T->wbuf_size, w);
}

/* Build the tree of the words of the file in T. */
static void build(ct_tree *T, int with_inverse, int rec)
{
    int j, k, r, status, check, i;
    if (ct_init(T, alphabet, 0, opt_layout))
        die("ct_init failed");
    if (opt_promote)
        T->promote = opt_promote;
    for (j = 0; j < nwords; j++) {
        if (!with_inverse) {
            process(T, words[j], lens[j], &status, rec);
            continue;
        }
        /* as _tree_add_with_inverse of crossing_arcs.pyx */
        if (ct_reserve(T, 2, 2 * lens[j]))
            die("ct_reserve failed");
        if (rec)
            record(T->wbuf_size + 2 * lens[j], 0, T->wbuf_size, NULL);
        process(T, words[j], lens[j], &status, rec);
        if (status <= 0)
            continue;
        i = T->nwords - 1;
        r = T->wlen[i];
        for (k = 0; k < r; k++)
            invbuf[k] = T->wbuf[T->wstart[i] + r - 1 - k] ^ 1;
        process(T, invbuf, r, &check, rec);
        if (check != 1)
            die("the inverse of a new word is not new and primitive");
    }
}

/* ------------------------------------------------------------------ */
/* the copy alone                                                      */
/* ------------------------------------------------------------------ */

/* grown() of conjugate_tree.c */
static int grown(int cap, int need, int max, int start)
{
    int c = cap ? cap : start;
    if (c > max)
        c = max;
    while (c < need)
        c = c > max / 2 ? max : 2 * c;
    return c;
}

static volatile int sink;

/* Replay the recorded reservations and copies on a buffer of its own. */
static void copy_alone(void)
{
    int *buf = NULL;
    int cap = 0, size = 0, e;
    for (e = 0; e < nevents; e++) {
        if (events[e].need > cap) {
            cap = grown(cap, events[e].need, INT_MAX, 4);
            buf = (int *) xrealloc(buf, (size_t) cap * sizeof(int));
        }
        if (events[e].len) {
            memcpy(buf + size, events[e].src, (size_t) events[e].len * sizeof(int));
            size = events[e].size;
        }
    }
    sink = buf[0];
    free(buf);
}

/* ------------------------------------------------------------------ */
/* timing                                                              */
/* ------------------------------------------------------------------ */

static int opt_inverse;
static ct_tree sort_tree;
static int *order, *pivot, *out;

static void run_build(void)
{
    ct_tree T;
    build(&T, opt_inverse, 0);
    ct_free(&T);
}

static void run_sort(void)
{
    int num;
    if (ct_sorted_leaves(&sort_tree, order, pivot, alphabet, out, &num))
        die("ct_sorted_leaves failed");
    sink = num;
}

/* The check of ct_sorted_leaves that order is a permutation, copied, with
 * the allocation of its array (which ct_sorted_leaves shares with the sort). */
static void run_check(void)
{
    int c, n = alphabet;
    int *keys = (int *) malloc((size_t) n * sizeof(int));
    if (keys == NULL)
        die("out of memory");
    for (c = 0; c < n; c++)
        keys[c] = 0;
    for (c = 0; c < n; c++) {
        if (keys[order[c]])
            die("not a permutation");
        keys[order[c]] = 1;
    }
    sink = keys[n - 1];
    free(keys);
}

/* The best over rounds of the mean time of f, a round lasting 10 ms. */
static double best_time(void (*f)(void), int rounds)
{
    double best = 1e300, t0, t;
    long reps;
    int r;
    for (r = 0; r < rounds; r++) {
        reps = 0;
        t0 = now();
        do {
            f();
            reps++;
            t = now() - t0;
        } while (t < 0.01);
        if (t / (double) reps < best)
            best = t / (double) reps;
    }
    return best;
}

static void print_time(const char *what, double t)
{
    if (t < 1e-3)
        printf("%s %.3f us", what, t * 1e6);
    else if (t < 1)
        printf("%s %.3f ms", what, t * 1e3);
    else
        printf("%s %.3f s", what, t);
}

/* ------------------------------------------------------------------ */
/* the shape of a tree                                                 */
/* ------------------------------------------------------------------ */

static void print_shape(const ct_tree *T)
{
    int n = T->nstates, s, t, d, leaves = 0;
    int *depth = (int *) xrealloc(NULL, (size_t) n * sizeof(int));
    int *kids = (int *) xrealloc(NULL, (size_t) n * sizeof(int));
    long num[4] = {0, 0, 0, 0}, sum[4] = {0, 0, 0, 0}, max[4] = {0, 0, 0, 0};
    int *stack = (int *) xrealloc(NULL, (size_t) n * sizeof(int));
    int top;
    /* a node may be created after its children (a split), so the depths
     * (in nodes) are computed from the root down */
    for (s = 0; s < n; s++)
        depth[s] = -1;
    depth[0] = 0;
    for (s = 1; s < n; s++) {
        top = 0;
        for (t = s; depth[t] == -1; t = T->parent[t])
            stack[top++] = t;
        while (top) {
            int u = stack[--top];
            depth[u] = depth[T->parent[u]] + 1;
        }
    }
    for (s = 0; s < n; s++)
        kids[s] = 0;
    for (s = 1; s < n; s++) {
        kids[T->parent[s]]++;
        leaves += T->tend[s] == -1;
    }
    for (s = 0; s < n; s++) {
        if (s > 0 && T->tend[s] == -1)
            continue;
        d = depth[s] < 3 ? depth[s] : 3;
        num[d]++;
        sum[d] += kids[s];
        if (kids[s] > max[d])
            max[d] = kids[s];
    }
    printf("shape: %d nodes, %d leaves; children of internal nodes (max / mean) by depth:", n, leaves);
    for (d = 0; d < 4; d++) {
        if (num[d])
            printf(" %s%d: %ld / %.2f (%ld nodes)", d == 3 ? ">=" : "", d, max[d],
                   (double) sum[d] / (double) num[d], num[d]);
        else
            printf(" %s%d: -", d == 3 ? ">=" : "", d);
    }
    printf("\n");
    free(depth);
    free(kids);
    free(stack);
}

/* The bytes allocated for T, from the capacities. */
static double tree_bytes(const ct_tree *T)
{
    double per_node = 6, b;
    if (T->layout == CT_LAYOUT_DENSE)
        per_node += T->alphabet_size;
    else
        per_node += 2;
    if (T->layout == CT_LAYOUT_ROWS)
        per_node += 2;
    b = 4 * (per_node * T->capacity + T->wbuf_capacity + 2.0 * T->words_capacity);
    b += 4.0 * T->rows_capacity * T->alphabet_size;
    return b;
}

static const char *layout_name(int layout)
{
    switch (layout) {
    case CT_LAYOUT_SPARSE: return "sparse";
    case CT_LAYOUT_DENSE: return "dense";
    case CT_LAYOUT_ROWS: return "rows";
    default: return "?";
    }
}

/* ------------------------------------------------------------------ */

int main(int argc, char **argv)
{
    const char *path = NULL;
    int rounds = 5, j, maxlen = 0, positional = 0;
    double t_build, t_copy;
    ct_tree T;

    for (j = 1; j < argc; j++) {
        if (strcmp(argv[j], "--with-inverse") == 0)
            opt_inverse = 1;
        else if (strcmp(argv[j], "--layout") == 0 && j + 1 < argc) {
            const char *l = argv[++j];
            opt_layout = !strcmp(l, "default") ? CT_LAYOUT_DEFAULT :
                         !strcmp(l, "sparse") ? CT_LAYOUT_SPARSE :
                         !strcmp(l, "dense") ? CT_LAYOUT_DENSE :
                         !strcmp(l, "rows") ? CT_LAYOUT_ROWS : -2;
            if (opt_layout == -2)
                die("unknown layout");
        } else if (strcmp(argv[j], "--promote") == 0 && j + 1 < argc) {
            opt_promote = atoi(argv[++j]);
            if (opt_promote <= 0)
                die("--promote needs a positive integer");
        } else if (argv[j][0] == '-')
            die("usage: bench <word file> [rounds] [--with-inverse] [--layout <layout>] [--promote <P>]");
        else if (positional++ == 0)
            path = argv[j];
        else
            rounds = atoi(argv[j]);
    }
    if (path == NULL || rounds <= 0)
        die("usage: bench <word file> [rounds] [--with-inverse] [--layout <layout>] [--promote <P>]");
    read_file(path);
    if (opt_inverse && alphabet % 2)
        die("an alphabet of odd size has no inverses");
    for (j = 0; j < nwords; j++)
        if (lens[j] > maxlen)
            maxlen = lens[j];
    invbuf = (int *) xrealloc(NULL, (size_t) maxlen * sizeof(int));

    t_build = best_time(run_build, rounds);

    build(&T, opt_inverse, 1);
    t_copy = best_time(copy_alone, rounds);

    printf("%s: n = %d, %d words, %ld letters, %s, layout %s", path, alphabet, nwords,
           total_letters, opt_inverse ? "with inverse" : "plain", layout_name(T.layout));
    if (T.layout == CT_LAYOUT_ROWS)
        printf(" (promote %d, %d rows)", T.promote, T.nrows);
    printf("\n");
    print_time("build:", t_build);
    printf(" per tree, %.2f ns per letter\n", t_build / (double) total_letters * 1e9);
    print_time("copy alone:", t_copy);
    printf(" per tree, %.1f %% of the build\n", 100 * t_copy / t_build);
    print_shape(&T);
    printf("memory: %.0f bytes\n", tree_bytes(&T));

    if (opt_inverse) {
        double t_sort, t_check;
        order = (int *) xrealloc(NULL, (size_t) alphabet * sizeof(int));
        pivot = (int *) xrealloc(NULL, (size_t) alphabet * sizeof(int));
        out = (int *) xrealloc(NULL, (size_t) T.nstates * sizeof(int));
        for (j = 0; j < alphabet; j++) {
            order[j] = j;
            pivot[j] = j ^ 1;
        }
        sort_tree = T;
        t_sort = best_time(run_sort, rounds);
        t_check = best_time(run_check, rounds);
        print_time("sorted leaves:", t_sort);
        printf(", ");
        print_time("permutation check:", t_check);
        printf(", %.2f %% of the sort, %.2f %% of build + sort\n",
               100 * t_check / t_sort, 100 * t_check / (t_build + t_sort));
        free(order);
        free(pivot);
        free(out);
    }

    ct_free(&T);
    for (j = 0; j < nwords; j++)
        free(words[j]);
    free(words);
    free(lens);
    free(events);
    free(invbuf);
    return 0;
}
