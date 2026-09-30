/*
 * Tests of the canonical keys of combisurf/src/isomorphism.c, in C.
 *
 * Usage: test_isomorphism [seed]
 *
 * Exits 0 and prints one summary line when every case passes; otherwise
 * prints the failing case and exits 1. The fixed cases are the examples of
 * combisurf/isomorphism_naive.py; the random cases, reproducible from the
 * seed (1 by default), check that the key and the number of automorphisms
 * do not change under a random relabelling of the map.
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
#include "../../combisurf/src/isomorphism.h"

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

/* the largest map of the random cases, in half-edges */
#define MAX_N 256
#define MAX_ARRAYS 2
#define KEY_SIZE ISO_KEY_SIZE(MAX_N, MAX_ARRAYS, MAX_ARRAYS)
#define WORK_SIZE ISO_WORK_SIZE(MAX_N, MAX_ARRAYS, MAX_ARRAYS)

static int key[KEY_SIZE];
static int key2[KEY_SIZE];
static int work[WORK_SIZE];
static int relabelling[MAX_N];
static int roots[MAX_N];

/* splitmix64 */
static uint64_t rng_state;

static uint64_t rng_next(void)
{
    uint64_t z = (rng_state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* a random integer in 0, ..., m - 1 */
static int rng_below(int m)
{
    return (int) (rng_next() % (uint64_t) m);
}

/* the relabelling r applied to fp, into out (-1 on the inactive images) */
static void relabel_fp(int n, const int *fp, const int *r, int *out)
{
    int h;
    for (h = 0; h < n; h++)
        out[h] = -1;
    for (h = 0; h < n; h++)
        if (fp[h] != -1)
            out[r[h]] = r[fp[h]];
}

/* the relabelling r applied to the m arrays of a (each of n ints) */
static void relabel_arrays(int n, int m, const int *a, const int *r, int *out)
{
    int j, h;
    for (j = 0; j < m; j++)
        for (h = 0; h < n; h++)
            out[j * n + r[h]] = a[j * n + h];
}

/* ------------------------------------------------------------------ */
/* fixed cases                                                         */
/* ------------------------------------------------------------------ */

/*
 * Check iso_canonical on fp against the expected key, relabelling (-1 on
 * the inactive edges) and number of automorphisms.
 */
static void check_case(const char *name, int n, const int *fp, int nd, const int *data,
                       int nc, const int *classes, int len, const int *expected_key,
                       const int *expected_relabelling, int expected_count)
{
    int count = -1, i;
    current_case = name;
    CHECK(iso_canonical(n, fp, nd, data, nc, classes, key, relabelling, roots, &count, work) == len);
    for (i = 0; i < len; i++)
        CHECK(key[i] == expected_key[i]);
    for (i = 0; i < n; i++)
        CHECK(relabelling[i] == expected_relabelling[i]);
    CHECK(count == expected_count);
    for (i = 1; i < count; i++)
        CHECK(roots[i - 1] < roots[i]);
    /* the key and the count do not depend on the optional outputs */
    CHECK(iso_canonical(n, fp, nd, data, nc, classes, key2, NULL, NULL, NULL, work) == len);
    CHECK(memcmp(key, key2, (size_t) len * sizeof(int)) == 0);
}

static void test_fixed(void)
{
    /* the faces (0,1,2)(~0,~2,~1) */
    static const int triangle[] = {2, 5, 4, 1, 0, 3};
    static const int identity6[] = {0, 1, 2, 3, 4, 5};
    /* one marked half-edge */
    static const int marked[] = {0, 0, 1, 0, 0, 0};
    static const int marked_key[] = {2, 5, 4, 1, 0, 3, 0, 0, 0, 0, 0, 1};
    static const int marked_relabelling[] = {1, 0, 5, 4, 3, 2};
    /* the faces (0,1,2)(~0,~1,~2), with the two faces marked */
    static const int two_faces[] = {2, 3, 4, 5, 0, 1};
    static const int faces[] = {0, 1, 0, 1, 0, 1};
    static const int two_faces_key[] = {2, 3, 4, 5, 0, 1, 0, 1, 0, 1, 0, 1};
    /* a folded edge: the faces (0,1)(~1) */
    static const int folded[] = {2, -1, 0, 3};
    static const int folded_key[] = {0, 2, 1, -1};
    static const int folded_relabelling[] = {2, 3, 1, 0};
    /* a label gap: the face (0,2,~2), the edge 1 inactive */
    static const int gap[] = {4, -1, -1, -1, 5, 0};
    static const int gap_key[] = {1, 2, 0, -1};
    static const int gap_relabelling[] = {2, 3, -1, -1, 0, 1};
    /* a single half-edge */
    static const int single[] = {0, -1};
    static const int single_key[] = {0, -1};
    static const int single_relabelling[] = {0, 1};
    /* no active half-edge */
    static const int empty[] = {-1, -1};
    static const int empty_relabelling[] = {-1, -1};
    /* two components (0)(~0) and (1)(~1) */
    static const int disconnected[] = {0, 1, 2, 3};
    int count = -1;

    check_case("triangle", 6, triangle, 0, NULL, 0, NULL, 6, triangle, identity6, 6);
    check_case("triangle, one marked half-edge", 6, triangle, 1, marked, 0, NULL,
               12, marked_key, marked_relabelling, 1);
    check_case("triangle, one marked half-edge as a class", 6, triangle, 0, NULL, 1, marked,
               12, marked_key, marked_relabelling, 1);
    check_case("two faces, data", 6, two_faces, 1, faces, 0, NULL, 12, two_faces_key, identity6, 3);
    check_case("two faces, classes", 6, two_faces, 0, NULL, 1, faces, 12, two_faces_key, identity6, 6);
    check_case("folded edge", 4, folded, 0, NULL, 0, NULL, 4, folded_key, folded_relabelling, 1);
    check_case("label gap", 6, gap, 0, NULL, 0, NULL, 4, gap_key, gap_relabelling, 1);
    check_case("single half-edge", 2, single, 0, NULL, 0, NULL, 2, single_key, single_relabelling, 1);
    check_case("no active half-edge", 2, empty, 0, NULL, 0, NULL, 0, NULL, empty_relabelling, 0);

    current_case = "disconnected";
    CHECK(iso_canonical(4, disconnected, 0, NULL, 0, NULL, key, relabelling, roots, &count, work)
          == ISO_DISCONNECTED);
    CHECK(count == -1);
}

/* ------------------------------------------------------------------ */
/* random cases                                                        */
/* ------------------------------------------------------------------ */

/*
 * A random connected face permutation on n = 2 * ne half-edges, into fp:
 * the active half-edges (the odd one of an edge missing with probability
 * folded / 4) in a random cyclic order of cycles, retried until connected.
 * Returns 0 if no connected one was found in 1000 attempts.
 */
static int random_connected(int ne, int folded, int *fp)
{
    int n = 2 * ne, active[MAX_N], images[MAX_N];
    int attempt, na, i, j, t;
    for (attempt = 0; attempt < 1000; attempt++) {
        na = 0;
        for (i = 0; i < ne; i++) {
            active[na++] = 2 * i;
            if (rng_below(4) >= folded)
                active[na++] = 2 * i + 1;
        }
        memcpy(images, active, (size_t) na * sizeof(int));
        for (i = na - 1; i > 0; i--) {
            j = rng_below(i + 1);
            t = images[i];
            images[i] = images[j];
            images[j] = t;
        }
        for (i = 0; i < n; i++)
            fp[i] = -1;
        for (i = 0; i < na; i++)
            fp[active[i]] = images[i];
        if (iso_canonical(n, fp, 0, NULL, 0, NULL, key, NULL, NULL, NULL, work) != ISO_DISCONNECTED)
            return 1;
    }
    return 0;
}

/*
 * A random relabelling of fp that preserves the edges: a random permutation
 * of the edges, each non-folded edge flipped at random. The active
 * half-edge of a folded edge may land on either half of its new edge.
 */
static void random_relabelling(int n, int *r)
{
    int ne = n / 2, sigma[MAX_N / 2], i, j, t, swap;
    for (i = 0; i < ne; i++)
        sigma[i] = i;
    for (i = ne - 1; i > 0; i--) {
        j = rng_below(i + 1);
        t = sigma[i];
        sigma[i] = sigma[j];
        sigma[j] = t;
    }
    for (i = 0; i < ne; i++) {
        swap = rng_below(2);
        r[2 * i] = 2 * sigma[i] + swap;
        r[2 * i + 1] = 2 * sigma[i] + 1 - swap;
    }
}

static void test_random(void)
{
    int fp[MAX_N], fp2[MAX_N], fp3[MAX_N], r[MAX_N];
    int dec[2 * MAX_ARRAYS * MAX_N], dec2[2 * MAX_ARRAYS * MAX_N];
    int trial, ne, n, nd, nc, len, len2, count, count2, i, h;

    current_case = "random";
    for (trial = 0; trial < 300; trial++) {
        ne = 1 + rng_below(trial < 250 ? 20 : MAX_N / 2);
        n = 2 * ne;
        CHECK(random_connected(ne, rng_below(3), fp));
        nd = rng_below(MAX_ARRAYS + 1);
        nc = rng_below(MAX_ARRAYS + 1);
        /* few values, so that automorphisms may survive */
        for (i = 0; i < (nd + nc) * n; i++)
            dec[i] = rng_below(2);
        len = iso_canonical(n, fp, nd, dec, nc, dec + nd * n, key, relabelling, roots, &count, work);
        CHECK(len >= n);
        CHECK(count >= 1);
        for (i = 0; i < count; i++)
            CHECK(fp[roots[i]] != -1 && (i == 0 || roots[i - 1] < roots[i]));

        /* the canonical relabelling sends fp to the face permutation of the key */
        relabel_fp(n, fp, relabelling, fp3);
        for (h = 0; h < n; h++)
            CHECK(fp3[h] == key[h]);

        /* a relabelled copy, with the class values renamed, has the same key */
        random_relabelling(n, r);
        relabel_fp(n, fp, r, fp2);
        relabel_arrays(n, nd + nc, dec, r, dec2);
        for (i = nd * n; i < (nd + nc) * n; i++)
            dec2[i] = 1 - dec2[i];
        len2 = iso_canonical(n, fp2, nd, dec2, nc, dec2 + nd * n, key2, NULL, NULL, &count2, work);
        CHECK(len2 == len);
        CHECK(count2 == count);
        CHECK(memcmp(key, key2, (size_t) len * sizeof(int)) == 0);

        /* a new value on an active half-edge changes the multiset of the
         * values of the first data array, hence the key */
        if (nd > 0) {
            for (h = 0; fp2[h] == -1; h++)
                ;
            dec2[h] = 2;
            len2 = iso_canonical(n, fp2, nd, dec2, nc, dec2 + nd * n, key2, NULL, NULL, NULL, work);
            CHECK(len2 == len);
            CHECK(memcmp(key, key2, (size_t) len * sizeof(int)) != 0);
        }
    }
}

int main(int argc, char **argv)
{
    rng_state = argc > 1 ? strtoull(argv[1], NULL, 10) : 1;
    test_fixed();
    test_random();
    printf("test_isomorphism: %ld checks passed\n", num_checks);
    return 0;
}
