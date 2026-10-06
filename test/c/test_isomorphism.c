/*
 * Tests of the canonical keys of combisurf/src/isomorphism.c, in C.
 *
 * Usage: test_isomorphism [seed]
 *
 * Exits 0 and prints one summary line when every case passes; otherwise
 * prints the failing case and exits 1. The fixed cases are the examples of
 * combisurf/isomorphism_naive.py; the random cases, reproducible from the
 * seed (1 by default), check that the key and the number of automorphisms
 * do not change under a random relabelling of the map, for connected maps
 * (iso_canonical) and for disjoint unions (iso_search, where the minimal
 * sequence of blocks and the number of leaves reaching it are checked).
 * iso_canonical, which rejects a root as soon as its key exceeds the best
 * one, is also compared with a reference that walks and keys every root in
 * full, on random maps, on square-tiled surfaces decorated by h % 4 and on
 * relabelled copies of the tetrahedron and of the cube, whose roots all tie.
 * Under the rules restricting the roots, both are compared with a
 * reference of the classes that computes the invariants by walking the
 * cycles and compares every pair of half-edges: iso_canonical tries the
 * roots of the class of the reference, and each block of iso_search starts
 * at a half-edge of the class of the reference among the half-edges not
 * labelled before; the number of automorphisms does not depend on the
 * rule. The rules run with and without the flag ISO_ROOTS_DATA, whose
 * reference invariants read the data values themselves; half of the random
 * cases take their data values among -7, n + 11, INT_MIN and INT_MAX, far
 * from the ranks 0, ..., n - 1 the C code sorts by.
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

/* the rules, ISO_ROOTS_ALL first, then the base rules with the flag */
static const iso_root_rule rules[] = {
    ISO_ROOTS_ALL, ISO_ROOTS_HALF_EDGE_MAX, ISO_ROOTS_HALF_EDGE_FEWEST,
    ISO_ROOTS_EDGE_MAX, ISO_ROOTS_EDGE_FEWEST,
    ISO_ROOTS_HALF_EDGE_MAX | ISO_ROOTS_DATA, ISO_ROOTS_HALF_EDGE_FEWEST | ISO_ROOTS_DATA,
    ISO_ROOTS_EDGE_MAX | ISO_ROOTS_DATA, ISO_ROOTS_EDGE_FEWEST | ISO_ROOTS_DATA
};
#define NUM_RULES 9

/* the base rule of a rule, without the flag */
static int base_rule(iso_root_rule rule)
{
    return rule & ~ISO_ROOTS_DATA;
}

/*
 * Data values far from 0, ..., n - 1, negative and larger than n, among
 * which the random cases draw: index v of the table replaces the value v.
 */
static void spread_values(int n, int count, int *a)
{
    const int wide[] = {-7, n + 11, INT_MIN, INT_MAX};
    int i;
    for (i = 0; i < count; i++)
        a[i] = wide[a[i]];
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
    CHECK(iso_canonical(n, fp, nd, data, nc, classes, ISO_ROOTS_ALL, key, relabelling, roots, &count,
                        work) == len);
    for (i = 0; i < len; i++)
        CHECK(key[i] == expected_key[i]);
    for (i = 0; i < n; i++)
        CHECK(relabelling[i] == expected_relabelling[i]);
    CHECK(count == expected_count);
    for (i = 1; i < count; i++)
        CHECK(roots[i - 1] < roots[i]);
    /* the key and the count do not depend on the optional outputs */
    CHECK(iso_canonical(n, fp, nd, data, nc, classes, ISO_ROOTS_ALL, key2, NULL, NULL, NULL, work)
          == len);
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
    CHECK(iso_canonical(4, disconnected, 0, NULL, 0, NULL, ISO_ROOTS_ALL, key, relabelling, roots,
                        &count, work) == ISO_DISCONNECTED);
    CHECK(count == -1);
}

/* ------------------------------------------------------------------ */
/* search over the components                                          */
/* ------------------------------------------------------------------ */

#define SEQUENCE_SIZE ISO_SEQUENCE_SIZE(MAX_N, MAX_ARRAYS, MAX_ARRAYS)
#define SEARCH_WORK_SIZE ISO_SEARCH_WORK_SIZE(MAX_N, MAX_ARRAYS, MAX_ARRAYS)

static int sequence[SEQUENCE_SIZE];
static int sequence2[SEQUENCE_SIZE];
static int search_work[SEARCH_WORK_SIZE];

/*
 * The state of the leaf callback: the leaves since the last one reported
 * with first = 1 (count) and in all (calls). Every leaf reported after a
 * first one must send the map and its data arrays to the same arrays as
 * the first one (image, of (1 + nd) * n ints).
 */
struct leaves {
    int n;
    const int *fp;
    int nd;
    const int *data;
    long count;
    long calls;
    long stop_after;
    int image[(1 + MAX_ARRAYS) * MAX_N];
};

static void leaf_image(const struct leaves *s, const int *leaf, int *image)
{
    int n = s->n, h, j;
    for (h = 0; h < (1 + s->nd) * n; h++)
        image[h] = -2;
    for (h = 0; h < n; h++) {
        if (s->fp[h] == -1)
            continue;
        CHECK(0 <= leaf[h] && leaf[h] < n && leaf[h ^ 1] == (leaf[h] ^ 1));
        image[leaf[h]] = leaf[s->fp[h]];
        for (j = 0; j < s->nd; j++)
            image[(1 + j) * n + leaf[h]] = s->data[j * n + h];
    }
}

static int count_leaf(const int *leaf, int first, void *arg)
{
    struct leaves *s = arg;
    int image[(1 + MAX_ARRAYS) * MAX_N];
    s->calls++;
    if (first) {
        s->count = 1;
        leaf_image(s, leaf, s->image);
    } else {
        CHECK(s->count >= 1);
        s->count++;
        leaf_image(s, leaf, image);
        CHECK(memcmp(image, s->image, (size_t) (1 + s->nd) * s->n * sizeof(int)) == 0);
    }
    return s->stop_after > 0 && s->calls >= s->stop_after ? 7 : 0;
}

/*
 * Run iso_search on fp with all leaves, into sequence and relabelling;
 * store its length into *len and return the number of leaves reaching the
 * minimal sequence. The search without all leaves, which backjumps from
 * the leaves equal to the best one, must give the same sequence and
 * relabelling.
 */
static long search(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                   iso_root_rule rule, int *len, long *calls)
{
    static int first_sequence[SEQUENCE_SIZE];
    static int first_relabelling[MAX_N];
    struct leaves s;
    int len2;
    s.n = n;
    s.fp = fp;
    s.nd = nd;
    s.data = data;
    s.count = 0;
    s.calls = 0;
    s.stop_after = 0;
    CHECK(iso_search(n, fp, nd, data, nc, classes, rule, 1, sequence, len, relabelling,
                     count_leaf, &s, search_work) == 0);
    CHECK(s.count >= 1);
    if (calls != NULL)
        *calls = s.calls;
    CHECK(iso_search(n, fp, nd, data, nc, classes, rule, 0, first_sequence, &len2, first_relabelling,
                     NULL, NULL, search_work) == 0);
    CHECK(len2 == *len);
    CHECK(memcmp(first_sequence, sequence, (size_t) len2 * sizeof(int)) == 0);
    CHECK(memcmp(first_relabelling, relabelling, (size_t) n * sizeof(int)) == 0);
    return s.count;
}

static void test_search_fixed(void)
{
    static const int two_edges[] = {0, 1, 2, 3};
    static const int two_edges_sequence[] = {2, 0, 1, 2, 0, 1};
    static const int shared[] = {0, 1, 1, 2};
    static const int shared_sequence[] = {2, 0, 1, 0, 1, 2, 0, 1, 0, 2};
    static const int shared_relabelling[] = {1, 0, 2, 3};
    static const int four_edges[] = {0, 1, 2, 3, 4, 5, 6, 7};
    static const int pairs[] = {0, 1, 0, 1, 2, 3, 2, 3};
    static const int triangle_and_folded[] = {2, 5, 4, 1, 0, 3, 6, -1};
    static const int triangle_and_folded_relabelling[] = {2, 3, 4, 5, 6, 7, 0, 1};
    static const int folded_first_sequence[] = {2, 0, -1, 6, 2, 5, 4, 1, 0, 3};
    static const int triangle_first_sequence[] = {6, 2, 5, 4, 1, 0, 3, 2, 0, -1};
    static const int identity8[] = {0, 1, 2, 3, 4, 5, 6, 7};
    static const int empty[] = {-1, -1};
    struct leaves s;
    int len, i, count;
    long calls;

    current_case = "search, two edges";
    CHECK(search(4, two_edges, 0, NULL, 0, NULL, ISO_ROOTS_ALL, &len, NULL) == 8);
    CHECK(len == 6);
    for (i = 0; i < len; i++)
        CHECK(sequence[i] == two_edges_sequence[i]);
    for (i = 0; i < 4; i++)
        CHECK(relabelling[i] == i);

    current_case = "search, two edges, shared classes";
    /* the first leaf, from the roots 0 then 2, is replaced by the one from
     * the roots 1 then 2, followed by the one from 2 then 1 */
    CHECK(search(4, two_edges, 0, NULL, 1, shared, ISO_ROOTS_ALL, &len, &calls) == 2);
    CHECK(calls == 3);
    CHECK(len == 10);
    for (i = 0; i < len; i++)
        CHECK(sequence[i] == shared_sequence[i]);
    for (i = 0; i < 4; i++)
        CHECK(relabelling[i] == shared_relabelling[i]);

    current_case = "search, four edges";
    CHECK(search(8, four_edges, 0, NULL, 0, NULL, ISO_ROOTS_ALL, &len, NULL) == 384);
    CHECK(len == 12);
    CHECK(search(8, four_edges, 0, NULL, 1, pairs, ISO_ROOTS_ALL, &len, NULL) == 32);
    /* the same classes kept exactly */
    CHECK(search(8, four_edges, 1, pairs, 0, NULL, ISO_ROOTS_ALL, &len, NULL) == 4);

    current_case = "search, triangle and folded edge";
    CHECK(search(8, triangle_and_folded, 0, NULL, 0, NULL, ISO_ROOTS_ALL, &len, NULL) == 6);
    CHECK(len == 3 + 7);
    for (i = 0; i < 8; i++)
        CHECK(relabelling[i] == triangle_and_folded_relabelling[i]);

    /* the outputs of isomorphism_naive: the class changes from the first
     * step to the second, the folded edge (its own class) then the
     * triangle under *_FEWEST, the triangle (vertices of degree 2, faces of
     * degree 3) then the folded edge under *_MAX */
    for (i = 1; i < NUM_RULES; i++) {
        const int *expected_sequence, *expected_relabelling;
        int j;
        current_case = "search, triangle and folded edge, under a rule";
        if (base_rule(rules[i]) == ISO_ROOTS_HALF_EDGE_FEWEST || base_rule(rules[i]) == ISO_ROOTS_EDGE_FEWEST) {
            expected_sequence = folded_first_sequence;
            expected_relabelling = triangle_and_folded_relabelling;
        } else {
            expected_sequence = triangle_first_sequence;
            expected_relabelling = identity8;
        }
        CHECK(search(8, triangle_and_folded, 0, NULL, 0, NULL, rules[i], &len, NULL) == 6);
        CHECK(len == 10);
        for (j = 0; j < len; j++)
            CHECK(sequence[j] == expected_sequence[j]);
        for (j = 0; j < 8; j++)
            CHECK(relabelling[j] == expected_relabelling[j]);
        current_case = "connected triangle and folded edge, under a rule";
        CHECK(iso_canonical(8, triangle_and_folded, 0, NULL, 0, NULL, rules[i], key, relabelling,
                            roots, &count, work) == ISO_DISCONNECTED);
    }

    current_case = "search, no active half-edge";
    CHECK(search(2, empty, 0, NULL, 0, NULL, ISO_ROOTS_ALL, &len, NULL) == 1);
    CHECK(len == 0);
    CHECK(relabelling[0] == -1 && relabelling[1] == -1);

    current_case = "search, stopped by the callback";
    s.n = 8;
    s.fp = four_edges;
    s.nd = 0;
    s.data = NULL;
    s.count = 0;
    s.calls = 0;
    s.stop_after = 5;
    CHECK(iso_search(8, four_edges, 0, NULL, 0, NULL, ISO_ROOTS_ALL, 1, sequence, &len, NULL,
                     count_leaf, &s, search_work) == 7);
    CHECK(s.calls == 5);
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
        if (iso_canonical(n, fp, 0, NULL, 0, NULL, ISO_ROOTS_ALL, key, NULL, NULL, NULL, work)
                != ISO_DISCONNECTED)
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
        len = iso_canonical(n, fp, nd, dec, nc, dec + nd * n, ISO_ROOTS_ALL, key, relabelling, roots,
                            &count, work);
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
        len2 = iso_canonical(n, fp2, nd, dec2, nc, dec2 + nd * n, ISO_ROOTS_ALL, key2, NULL, NULL,
                             &count2, work);
        CHECK(len2 == len);
        CHECK(count2 == count);
        CHECK(memcmp(key, key2, (size_t) len * sizeof(int)) == 0);

        /* a new value on an active half-edge changes the multiset of the
         * values of the first data array, hence the key */
        if (nd > 0) {
            for (h = 0; fp2[h] == -1; h++)
                ;
            dec2[h] = 2;
            len2 = iso_canonical(n, fp2, nd, dec2, nc, dec2 + nd * n, ISO_ROOTS_ALL, key2, NULL, NULL,
                                 NULL, work);
            CHECK(len2 == len);
            CHECK(memcmp(key, key2, (size_t) len * sizeof(int)) != 0);
        }
    }
}

/* ------------------------------------------------------------------ */
/* iso_canonical against every root walked and keyed in full           */
/* ------------------------------------------------------------------ */

/*
 * The invariant of the active half-edge h under rule, into inv (INV_SIZE
 * ints, the unused entries 0), computed by walking the cycles from h and
 * reading the data values themselves, not ranks of them: the vertex of h
 * is h, vp[h], ... with vp[fp[t]] = h for t the twin of h (h itself on a
 * folded edge), that is vp[x] = twin(fp^-1(x)).
 */
static int ref_twin(const int *fp, int h)
{
    return fp[h ^ 1] != -1 ? h ^ 1 : h;
}

/*
 * The length of the cycle of h in the vertex permutation (vertex nonzero)
 * or in fp, and whether g is on it; finv is the inverse of fp.
 */
static int ref_cycle(const int *fp, const int *finv, int vertex, int h, int g, int *meets)
{
    int x = h, len = 0;
    *meets = 0;
    do {
        if (x == g)
            *meets = 1;
        len++;
        x = vertex ? ref_twin(fp, finv[x]) : fp[x];
    } while (x != h);
    return len;
}

/* the length of a reference invariant: 6 base entries, 2 per data array */
#define INV_SIZE (6 + 2 * MAX_ARRAYS)

/*
 * Under ISO_ROOTS_DATA, the base entries are followed by the nd data values
 * at h, then under the edge rules by those at t, in the entries 6 on.
 */
static void ref_invariant(const int *fp, const int *finv, iso_root_rule rule, int n, int nd,
                          const int *data, int h, int *inv)
{
    int t = ref_twin(fp, h), loop, monogon, i, base = base_rule(rule);
    int edge = base == ISO_ROOTS_EDGE_MAX || base == ISO_ROOTS_EDGE_FEWEST;
    for (i = 0; i < INV_SIZE; i++)
        inv[i] = 0;
    if (base == ISO_ROOTS_ALL)
        return;
    inv[0] = ref_cycle(fp, finv, 1, h, t, &loop);
    inv[1] = ref_cycle(fp, finv, 0, h, t, &monogon);
    if (edge) {
        inv[2] = ref_cycle(fp, finv, 1, t, t, &i);
        inv[3] = ref_cycle(fp, finv, 0, t, t, &i);
        inv[4] = loop;
        inv[5] = monogon;
    }
    if (rule & ISO_ROOTS_DATA) {
        for (i = 0; i < nd; i++) {
            inv[6 + i] = data[i * n + h];
            if (edge)
                inv[6 + nd + i] = data[i * n + t];
        }
    }
}

/*
 * The reference of the class of a rule: into mask (n ints), 1 on the
 * half-edges among those with avail[h] nonzero (active ones) whose
 * invariant is the one chosen by the rule among them, 0 elsewhere. Every
 * pair of half-edges is compared, O(n^2) invariants.
 */
static void ref_class(int n, const int *fp, iso_root_rule rule, int nd, const int *data,
                      const int *avail, int *mask)
{
    static int inv[MAX_N][INV_SIZE], finv[MAX_N];
    int h, g, cnt, best = -1, best_cnt = 0, cmp;
    int fewest = base_rule(rule) == ISO_ROOTS_HALF_EDGE_FEWEST || base_rule(rule) == ISO_ROOTS_EDGE_FEWEST;

    for (h = 0; h < n; h++)
        if (fp[h] != -1)
            finv[fp[h]] = h;
    for (h = 0; h < n; h++)
        if (avail[h])
            ref_invariant(fp, finv, rule, n, nd, data, h, inv[h]);
    for (h = 0; h < n; h++) {
        if (!avail[h])
            continue;
        cnt = 0;
        for (g = 0; g < n; g++)
            if (avail[g] && memcmp(inv[g], inv[h], sizeof(inv[h])) == 0)
                cnt++;
        cmp = best == -1 ? 1 : 0;
        if (cmp == 0 && fewest && cnt != best_cnt)
            cmp = cnt < best_cnt ? 1 : -1;
        if (cmp == 0) {
            for (g = 0; g < INV_SIZE && inv[h][g] == inv[best][g]; g++)
                ;
            cmp = g == INV_SIZE ? 0 : (inv[h][g] > inv[best][g] ? 1 : -1);
        }
        if (cmp > 0) {
            best = h;
            best_cnt = cnt;
        }
    }
    for (h = 0; h < n; h++)
        mask[h] = best != -1 && avail[h] && memcmp(inv[h], inv[best], sizeof(inv[h])) == 0;
}

/*
 * The reference of iso_canonical: every active half-edge of the class of the
 * rule (ref_class) as a root, its relabelling by iso_relabelling_from and
 * its key by iso_key in full, then compared with the best key so far.
 */
static int ref_canonical(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                         iso_root_rule rule, int *best, int *best_relabelling, int *best_roots,
                         int *nroots)
{
    static int r[MAX_N], rinv[MAX_N], queue[MAX_N], names[MAX_N], cur[KEY_SIZE];
    static int avail[MAX_N], mask[MAX_N];
    int ne = iso_num_edges(n, fp), root, k, len = 0, l, count = 0, cmp, i;

    for (i = 0; i < n; i++) {
        names[i] = -1;
        best_relabelling[i] = -1;
        avail[i] = fp[i] != -1;
    }
    ref_class(n, fp, rule, nd, data, avail, mask);
    for (root = 0; root < n; root++) {
        if (!mask[root])
            continue;
        k = iso_relabelling_from(n, fp, root, r, rinv, queue);
        if (k != 2 * ne)
            return ISO_DISCONNECTED;
        l = iso_key(n, fp, nd, data, nc, classes, r, rinv, k, cur, names);
        if (count == 0) {
            cmp = -1;
        } else {
            CHECK(l == len);
            cmp = 0;
            for (i = 0; i < len && cmp == 0; i++)
                if (cur[i] != best[i])
                    cmp = cur[i] < best[i] ? -1 : 1;
        }
        if (cmp < 0) {
            len = l;
            memcpy(best, cur, (size_t) l * sizeof(int));
            memcpy(best_relabelling, r, (size_t) n * sizeof(int));
            count = 0;
        }
        if (cmp <= 0)
            best_roots[count++] = root;
    }
    *nroots = count;
    return len;
}

/*
 * Check iso_canonical against ref_canonical on fp under every rule. The
 * number of roots of minimal key, the order of the group, does not depend
 * on the rule.
 */
static void check_reference(int n, const int *fp, int nd, const int *data, int nc, const int *classes)
{
    static int ref_key[KEY_SIZE], ref_relabelling[MAX_N], ref_roots[MAX_N];
    int len, ref_len, count = -1, ref_count = -1, all_count = -1, i;

    for (i = 0; i < NUM_RULES; i++) {
        ref_len = ref_canonical(n, fp, nd, data, nc, classes, rules[i], ref_key, ref_relabelling,
                                ref_roots, &ref_count);
        len = iso_canonical(n, fp, nd, data, nc, classes, rules[i], key, relabelling, roots, &count,
                            work);
        CHECK(len == ref_len);
        if (len == ISO_DISCONNECTED)
            return;
        CHECK(count == ref_count);
        CHECK(memcmp(key, ref_key, (size_t) len * sizeof(int)) == 0);
        CHECK(memcmp(relabelling, ref_relabelling, (size_t) n * sizeof(int)) == 0);
        CHECK(memcmp(roots, ref_roots, (size_t) count * sizeof(int)) == 0);
        if (i == 0)
            all_count = count;
        CHECK(count == all_count);
    }
}

/*
 * The square-tiled surface of the permutations right and up of the squares
 * 0, ..., m - 1 as a quadrangulation on 4 m half-edges, into fp: the edge
 * 2 i glues the right side 4 i of the square i to the left side 4 i + 1 of
 * right[i], the edge 2 i + 1 its top side 4 i + 2 to the bottom side
 * 4 i + 3 of up[i], so that h % 4 tells the side of the half-edge h.
 */
static void origami(int m, const int *right, const int *up, int *fp)
{
    int left_of[MAX_N / 4], below[MAX_N / 4], i, bottom, left;
    for (i = 0; i < m; i++) {
        left_of[right[i]] = i;
        below[up[i]] = i;
    }
    for (i = 0; i < m; i++) {
        bottom = 4 * below[i] + 3;
        left = 4 * left_of[i] + 1;
        fp[bottom] = 4 * i;
        fp[4 * i] = 4 * i + 2;
        fp[4 * i + 2] = left;
        fp[left] = bottom;
    }
}

/* A random permutation of 0, ..., m - 1, into p. */
static void random_permutation(int m, int *p)
{
    int i, j, t;
    for (i = 0; i < m; i++)
        p[i] = i;
    for (i = m - 1; i > 0; i--) {
        j = rng_below(i + 1);
        t = p[i];
        p[i] = p[j];
        p[j] = t;
    }
}

/*
 * The map of the nf faces of degree deg of a polyhedron on the vertices 0,
 * ..., nv - 1, each face given by its vertices in counterclockwise order
 * (faces holds nf * deg ints), into fp; return its number of half-edges.
 * The edges are numbered by first appearance, the half-edge 2 e going in
 * the direction of that appearance.
 */
static int polyhedron(int nv, int nf, int deg, const int *faces, int *fp)
{
    int he[16][16], f, i, a, b, c, n = 0;
    CHECK(nv <= 16);
    for (a = 0; a < nv; a++)
        for (b = 0; b < nv; b++)
            he[a][b] = -1;
    for (f = 0; f < nf; f++) {
        for (i = 0; i < deg; i++) {
            a = faces[f * deg + i];
            b = faces[f * deg + (i + 1) % deg];
            if (he[a][b] == -1) {
                CHECK(he[b][a] == -1 && n + 2 <= MAX_N);
                he[a][b] = n;
                he[b][a] = n + 1;
                n += 2;
            }
        }
    }
    for (f = 0; f < nf; f++) {
        for (i = 0; i < deg; i++) {
            a = faces[f * deg + i];
            b = faces[f * deg + (i + 1) % deg];
            c = faces[f * deg + (i + 2) % deg];
            fp[he[a][b]] = he[b][c];
        }
    }
    return n;
}

/*
 * Regular maps, where every root ties in the face permutation: the
 * tetrahedron (12 rotations) and the cube (24 rotations), relabelled at
 * random so that the roots of minimal key come in various orders, without
 * decoration and with a constant data array, under every rule.
 */
static void test_reference_regular(void)
{
    static const int tetrahedron[] = {0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2};
    static const int cube[] = {0, 3, 2, 1, 4, 5, 6, 7, 0, 1, 5, 4, 1, 2, 6, 5, 2, 3, 7, 6, 3, 0, 4, 7};
    int fp[MAX_N], fp2[MAX_N], r[MAX_N], zero[MAX_N];
    int trial, n, i, len, count, solid;

    for (i = 0; i < MAX_N; i++)
        zero[i] = 0;
    for (solid = 0; solid < 2; solid++) {
        current_case = solid == 0 ? "reference, tetrahedron" : "reference, cube";
        n = solid == 0 ? polyhedron(4, 4, 3, tetrahedron, fp) : polyhedron(8, 6, 4, cube, fp);
        for (trial = 0; trial < 20; trial++) {
            if (trial == 0) {
                memcpy(fp2, fp, (size_t) n * sizeof(int));
            } else {
                random_relabelling(n, r);
                relabel_fp(n, fp, r, fp2);
            }
            check_reference(n, fp2, 0, NULL, 0, NULL);
            check_reference(n, fp2, 1, zero, 0, NULL);
            check_reference(n, fp2, 0, NULL, 1, zero);
            for (i = 0; i < NUM_RULES; i++) {
                len = iso_canonical(n, fp2, 0, NULL, 0, NULL, rules[i], key, relabelling, roots,
                                    &count, work);
                CHECK(len == n && count == (solid == 0 ? 12 : 24));
            }
        }
    }
}

static void test_reference(void)
{
    int fp[MAX_N], dec[2 * MAX_ARRAYS * MAX_N], right[MAX_N / 4], up[MAX_N / 4];
    int trial, ne, n, m, nd, nc, i, j, h, s;

    current_case = "reference, random maps";
    for (trial = 0; trial < 300; trial++) {
        ne = 1 + rng_below(trial < 250 ? 20 : MAX_N / 2);
        n = 2 * ne;
        CHECK(random_connected(ne, rng_below(4), fp));
        nd = rng_below(MAX_ARRAYS + 1);
        nc = rng_below(MAX_ARRAYS + 1);
        /* two or three values: ties in the face permutation part are
         * broken late, or not at all (class values stay below n) */
        s = n > 2 ? 2 + rng_below(2) : 2;
        for (i = 0; i < (nd + nc) * n; i++)
            dec[i] = rng_below(s);
        if (trial % 2)
            spread_values(n, nd * n, dec);
        check_reference(n, fp, nd, dec, nc, dec + nd * n);
        check_reference(n, fp, 0, NULL, 0, NULL);
    }

    /* data arrays whose values span exactly n - 1 (the offsets of the
     * values by their minimum) or exactly n (the ranks), from a minimum
     * that is negative, 0, or near INT_MAX */
    current_case = "reference, data range at the threshold";
    for (trial = 0; trial < 300; trial++) {
        const int low[] = {-5, 0, 3, INT_MIN, INT_MAX};
        int lo, span, first, last;
        ne = 2 + rng_below(trial < 250 ? 20 : MAX_N / 2 - 2);
        n = 2 * ne;
        CHECK(random_connected(ne, rng_below(4), fp));
        nd = 1 + rng_below(MAX_ARRAYS);
        nc = rng_below(MAX_ARRAYS + 1);
        for (i = 0; i < nc * n; i++)
            dec[nd * n + i] = rng_below(2);
        for (j = 0; j < nd; j++) {
            span = n - 1 + rng_below(2);
            lo = low[rng_below(5)];
            /* the values lo ... lo + span, without overflow */
            if (lo > INT_MAX - span)
                lo = INT_MAX - span;
            first = last = -1;
            for (h = 0; h < n; h++) {
                if (fp[h] == -1)
                    continue;
                if (first == -1)
                    first = h;
                last = h;
                dec[j * n + h] = lo + (rng_below(3) == 0 ? rng_below(span + 1) : rng_below(2) * span);
            }
            dec[j * n + first] = lo;
            dec[j * n + last] = lo + span;
            for (h = 0; h < n; h++)
                if (fp[h] == -1)
                    dec[j * n + h] = 0;
        }
        check_reference(n, fp, nd, dec, nc, dec + nd * n);
    }

    current_case = "reference, origamis";
    for (trial = 0; trial < 300; trial++) {
        m = 1 + rng_below(trial < 250 ? 16 : MAX_N / 4);
        n = 4 * m;
        if (trial % 3 == 0) {
            /* a cyclic origami: every root ties in the face permutation */
            s = rng_below(m);
            for (i = 0; i < m; i++) {
                right[i] = (i + 1) % m;
                up[i] = (i + s) % m;
            }
        } else {
            random_permutation(m, right);
            random_permutation(m, up);
        }
        origami(m, right, up, fp);
        for (h = 0; h < n; h++) {
            dec[h] = h % 4;
            dec[n + h] = h % 4;
        }
        check_reference(n, fp, 1, dec, 0, NULL);
        check_reference(n, fp, 0, NULL, 1, dec);
        check_reference(n, fp, 1, dec, 1, dec + n);
        check_reference(n, fp, 0, NULL, 0, NULL);
    }
}

/*
 * Whether iso_canonical takes its counting path under rule (not
 * ISO_ROOTS_ALL) rather than the sort of all the classes: the product of
 * the ranges (max - min + 1 over the active half-edges) of the entries of
 * the reference invariant, loop and monogon counting 2 each, is at most
 * 2 n, and 2 nd <= n + 1. This mirrors the bound of isomorphism.c, so that
 * the cases below are known to cover both paths.
 */
static int counting_path(int n, const int *fp, iso_root_rule rule, int nd, const int *data)
{
    static int inv[MAX_N][INV_SIZE], finv[MAX_N];
    int h, i, lo, hi, first = 1, base = base_rule(rule);
    int edge = base == ISO_ROOTS_EDGE_MAX || base == ISO_ROOTS_EDGE_FEWEST;
    int entries = (edge ? 6 : 2) + ((rule & ISO_ROOTS_DATA) ? (edge ? 2 : 1) * nd : 0);
    long long prod = 1;

    if ((rule & ISO_ROOTS_DATA) && 2 * nd > n + 1)
        return 0;
    for (h = 0; h < n; h++)
        if (fp[h] != -1)
            finv[fp[h]] = h;
    for (h = 0; h < n; h++)
        if (fp[h] != -1)
            ref_invariant(fp, finv, rule, n, nd, data, h, inv[h]);
    for (i = 0; i < entries; i++) {
        /* the data entries start at 6 */
        int e = i < (edge ? 6 : 2) ? i : 6 + i - (edge ? 6 : 2);
        if (e == 4 || e == 5) {
            prod *= 2;
        } else {
            first = 1;
            lo = hi = 0;
            for (h = 0; h < n; h++) {
                if (fp[h] == -1)
                    continue;
                if (first || inv[h][e] < lo)
                    lo = inv[h][e];
                if (first || inv[h][e] > hi)
                    hi = inv[h][e];
                first = 0;
            }
            prod *= (long long) hi - lo + 1;
        }
        if (prod > 2 * (long long) n)
            return 0;
    }
    return 1;
}

/*
 * The two faces (0,1,...,k-1)(~(k-1),...,~1,~0) of a polygon on k edges,
 * into fp: every vertex of degree 2 and both faces of degree k, so that the
 * degrees span one value each.
 */
static void polygon(int k, int *fp)
{
    int i;
    for (i = 0; i < k; i++) {
        fp[2 * i] = 2 * ((i + 1) % k);
        fp[2 * i + 1] = 2 * ((i + k - 1) % k) + 1;
    }
}

/*
 * iso_canonical against the reference under every rule, on inputs that take
 * its counting path and inputs that take the sort of all the classes:
 * relabelled polygons (one degree each for the vertices and the faces)
 * with data of two values, random maps with data of two values, of values
 * spread over n + 1 (wide ranges) or of the values far apart of
 * spread_values, and random maps without data, whose degrees spread over
 * more values as the maps grow. The cases run until each rule other than
 * ISO_ROOTS_ALL has taken each path at least 20 times.
 */
static void test_reference_paths(void)
{
    int fp[MAX_N], fp2[MAX_N], r[MAX_N], dec[2 * MAX_ARRAYS * MAX_N];
    int fast[NUM_RULES], slow[NUM_RULES];
    int trial, done = 0, kind, ne, n, nd, nc, i, s;

    current_case = "reference, both paths of the classes";
    for (i = 0; i < NUM_RULES; i++)
        fast[i] = slow[i] = 0;
    for (trial = 0; !done; trial++) {
        CHECK(trial < 20000);
        kind = trial % 4;
        if (kind == 0) {
            ne = 1 + rng_below(20);
            n = 2 * ne;
            polygon(ne, fp);
            random_relabelling(n, r);
            relabel_fp(n, fp, r, fp2);
            memcpy(fp, fp2, (size_t) n * sizeof(int));
        } else {
            ne = 1 + rng_below(kind == 3 ? 40 : 20);
            n = 2 * ne;
            CHECK(random_connected(ne, rng_below(4), fp));
        }
        nd = kind == 3 ? 0 : rng_below(MAX_ARRAYS + 1);
        nc = rng_below(MAX_ARRAYS + 1);
        s = kind == 2 ? n + 1 : 2;
        for (i = 0; i < nd * n; i++)
            dec[i] = rng_below(s);
        if (kind == 1 && trial % 8 == 5)
            spread_values(n, nd * n, dec);
        for (i = 0; i < nc * n; i++)
            dec[nd * n + i] = rng_below(2);
        check_reference(n, fp, nd, dec, nc, dec + nd * n);
        done = 1;
        for (i = 1; i < NUM_RULES; i++) {
            if (counting_path(n, fp, rules[i], nd, dec))
                fast[i]++;
            else
                slow[i]++;
            if (fast[i] < 20 || slow[i] < 20)
                done = 0;
        }
    }
}

/*
 * A random disjoint union of one to three connected maps of at most 12
 * half-edges, each a relabelled copy of an earlier one with probability
 * 1/2, into fp; return its number of half-edges.
 */
static int random_union(int *fp)
{
    int comp[3][MAX_N], size[3], r[MAX_N], ncomp, c, n = 0, h;
    ncomp = 1 + rng_below(3);
    for (c = 0; c < ncomp; c++) {
        if (c > 0 && rng_below(2)) {
            int src = rng_below(c);
            size[c] = size[src];
            random_relabelling(size[c], r);
            relabel_fp(size[c], comp[src], r, comp[c]);
        } else {
            size[c] = 2 * (1 + rng_below(6));
            CHECK(random_connected(size[c] / 2, rng_below(3), comp[c]));
        }
        for (h = 0; h < size[c]; h++)
            fp[n + h] = comp[c][h] == -1 ? -1 : comp[c][h] + n;
        n += size[c];
    }
    return n;
}

/*
 * Check that the root of each block of the leaf rel of fp (the half-edge of
 * the first label of its component), reaching the sequence seq of len ints
 * with nd + nc decoration arrays, lies in the class of rule (ref_class)
 * among the active half-edges not labelled before it.
 */
static void check_block_roots(int n, const int *fp, int nd, const int *data, int narrays,
                              iso_root_rule rule, const int *rel, const int *seq, int len)
{
    static int avail[MAX_N], mask[MAX_N];
    int pos = 0, first = 0, k, na, root, h;

    while (pos < len) {
        k = seq[pos];
        na = 0;
        root = -1;
        for (h = 0; h < n; h++) {
            avail[h] = fp[h] != -1 && rel[h] >= first;
            if (fp[h] != -1 && rel[h] >= first && rel[h] < first + k)
                na++;
            if (rel[h] == first)
                root = h;
        }
        CHECK(root != -1 && fp[root] != -1);
        ref_class(n, fp, rule, nd, data, avail, mask);
        CHECK(mask[root]);
        pos += 1 + k + narrays * na;
        first += k;
    }
    CHECK(pos == len);
}

static void test_search_random(void)
{
    int fp[MAX_N], fp2[MAX_N], r[MAX_N], relabelling1[MAX_N];
    int dec[2 * MAX_ARRAYS * MAX_N], dec2[2 * MAX_ARRAYS * MAX_N];
    int trial, n, nd, nc, len, len2, count, i, h, j;
    long leaves, leaves2, all_leaves = 0;
    iso_root_rule rule;

    current_case = "search, random unions";
    for (trial = 0; trial < 300; trial++) {
        n = random_union(fp);
        nd = rng_below(MAX_ARRAYS + 1);
        nc = rng_below(MAX_ARRAYS + 1);
        for (i = 0; i < (nd + nc) * n; i++)
            dec[i] = rng_below(2);
        if (trial % 2)
            spread_values(n, nd * n, dec);
        for (j = 0; j < NUM_RULES; j++) {
            rule = rules[j];
            leaves = search(n, fp, nd, dec, nc, dec + nd * n, rule, &len, NULL);
            memcpy(relabelling1, relabelling, (size_t) n * sizeof(int));
            /* the leaves reaching the minimum are the automorphisms, whose
             * number does not depend on the rule */
            if (j == 0)
                all_leaves = leaves;
            CHECK(leaves == all_leaves);
            check_block_roots(n, fp, nd, dec, nd + nc, rule, relabelling1, sequence, len);

            /* the first leaf sends fp to the face permutation of the
             * canonical form, which starts with the first block */
            relabel_fp(n, fp, relabelling1, fp2);
            CHECK(sequence[0] <= n);
            for (h = 0; h < sequence[0]; h++)
                CHECK(fp2[h] == sequence[1 + h]);

            /* a connected map gives one block: k, then its key */
            len2 = iso_canonical(n, fp, nd, dec, nc, dec + nd * n, rule, key, relabelling, NULL,
                                 &count, work);
            if (len2 != ISO_DISCONNECTED) {
                CHECK(len == 1 + len2 && sequence[0] == n);
                CHECK(memcmp(sequence + 1, key, (size_t) len2 * sizeof(int)) == 0);
                CHECK(memcmp(relabelling, relabelling1, (size_t) n * sizeof(int)) == 0);
                CHECK(leaves == count);
            }

            /* a relabelled copy, with the class values renamed, has the
             * same sequence and as many leaves */
            random_relabelling(n, r);
            relabel_fp(n, fp, r, fp2);
            relabel_arrays(n, nd + nc, dec, r, dec2);
            for (i = nd * n; i < (nd + nc) * n; i++)
                dec2[i] = 1 - dec2[i];
            memcpy(sequence2, sequence, (size_t) len * sizeof(int));
            leaves2 = search(n, fp2, nd, dec2, nc, dec2 + nd * n, rule, &len2, NULL);
            CHECK(len2 == len);
            CHECK(leaves2 == leaves);
            CHECK(memcmp(sequence, sequence2, (size_t) len * sizeof(int)) == 0);
        }
    }
}

/*
 * Isolated edges whose halves carry the class values of pairs (2 k ints),
 * into fp and dec; return the number of half-edges.
 */
static int isolated_edges(int k, const int *pairs, int *fp, int *dec)
{
    int h;
    for (h = 0; h < 2 * k; h++) {
        fp[h] = h;
        dec[h] = pairs[h];
    }
    return 2 * k;
}

/*
 * Disjoint unions whose components tie, with classes. With all leaves every
 * root tied in a component is followed; without, only the first one of each
 * extension of the numbering of the class values, and the search backjumps
 * from each leaf equal to the best one: search() checks that the sequence
 * and the relabelling are those of the first leaf of the search with all
 * leaves. The searches with every leaf, whose number of leaves is the
 * number of automorphisms, are kept small.
 */
static void test_search_tied(void)
{
    static const int edge_block[] = {2, 0, 1, 0, 0};
    static const int edge_block_second[] = {2, 0, 1, 1, 1};
    int fp[MAX_N], fp2[MAX_N], r[MAX_N], comp[MAX_N], rc[MAX_N], pairs[MAX_N];
    int dec[2 * MAX_N], dec2[2 * MAX_N], cdec[2 * MAX_N], cdec2[2 * MAX_N];
    int trial, k, n, size, nd, nvals, len, len2, i, j, h, c, b, p, q;
    long leaves, leaves2, expected, all_leaves = 0;

    /* b + 1 isolated edges with the class values 0, 0 and b with 1, 1:
     * (b + 1)! 2^(b + 1) b! 2^b automorphisms; the edges 0, 0 come first */
    current_case = "search, isolated edges 0, 0 beside edges 1, 1";
    for (b = 1; b <= 3; b++) {
        for (h = 0; h < 2 * (2 * b + 1); h++)
            pairs[h] = h < 2 * (b + 1) ? 0 : 1;
        n = isolated_edges(2 * b + 1, pairs, fp, dec);
        expected = 1;
        for (i = 1; i <= b + 1; i++)
            expected *= 2 * i;
        for (i = 1; i <= b; i++)
            expected *= 2 * i;
        for (j = 0; j < NUM_RULES; j++) {
            CHECK(search(n, fp, 0, NULL, 1, dec, rules[j], &len, NULL) == expected);
            CHECK(len == 5 * (2 * b + 1));
            for (i = 0; i < len; i++)
                CHECK(sequence[i] == (i < 5 * (b + 1) ? edge_block[i % 5] : edge_block_second[i % 5]));
        }
    }

    /* the same with 41 edges (b = 20), without all leaves only */
    current_case = "search, 21 isolated edges 0, 0 beside 20 edges 1, 1";
    for (h = 0; h < 82; h++)
        pairs[h] = h < 42 ? 0 : 1;
    n = isolated_edges(41, pairs, fp, dec);
    for (j = 0; j < NUM_RULES; j++) {
        CHECK(iso_search(n, fp, 0, NULL, 1, dec, rules[j], 0, sequence, &len, relabelling, NULL, NULL,
                         search_work) == 0);
        CHECK(len == 5 * 41);
        for (i = 0; i < len; i++)
            CHECK(sequence[i] == (i < 5 * 21 ? edge_block[i % 5] : edge_block_second[i % 5]));
        if (j == 0)
            for (h = 0; h < n; h++)
                CHECK(relabelling[h] == h);
    }

    /* p isolated edges 1, 0, q edges 0, 1 and k edges 1, 1, in an order
     * that puts the edges 1, 1 first: (p + q)! k! 2^k automorphisms, since
     * the edges 1, 1 forbid the exchange of the two values */
    current_case = "search, mixed isolated edges";
    for (p = 1; p <= 3; p++) {
        for (q = 0; q <= 2; q++) {
            for (k = 1; k <= 2; k++) {
                for (h = 0; h < 2 * k; h++)
                    pairs[h] = 1;
                for (i = 0; i < p; i++) {
                    pairs[2 * (k + i)] = 1;
                    pairs[2 * (k + i) + 1] = 0;
                }
                for (i = 0; i < q; i++) {
                    pairs[2 * (k + p + i)] = 0;
                    pairs[2 * (k + p + i) + 1] = 1;
                }
                n = isolated_edges(k + p + q, pairs, fp, dec);
                expected = 1;
                for (i = 1; i <= p + q; i++)
                    expected *= i;
                for (i = 1; i <= k; i++)
                    expected *= 2 * i;
                for (j = 0; j < NUM_RULES; j++)
                    CHECK(search(n, fp, 0, NULL, 1, dec, rules[j], &len, NULL) == expected);
            }
        }
    }

    /* isolated edges with random pairs of values among two to four: the
     * graphs on the values given by the edges, whose tied roots extend the
     * numbering in several ways */
    current_case = "search, isolated edges with random pairs of values";
    for (trial = 0; trial < 300; trial++) {
        k = 1 + rng_below(7);
        /* class values lie below the number of half-edges */
        nvals = 2 + rng_below(3);
        if (nvals > 2 * k)
            nvals = 2 * k;
        for (h = 0; h < 2 * k; h++)
            pairs[h] = rng_below(nvals);
        n = isolated_edges(k, pairs, fp, dec);
        for (j = 0; j < NUM_RULES; j++) {
            leaves = search(n, fp, 0, NULL, 1, dec, rules[j], &len, NULL);
            if (j == 0)
                all_leaves = leaves;
            CHECK(leaves == all_leaves);
        }
    }

    /* k isolated edges with one class value: 2^k k! automorphisms, and k
     * blocks (2, 0, 1, 0, 0) */
    current_case = "search, isolated edges with one class value";
    for (k = 1; k <= 7; k++) {
        n = 2 * k;
        for (h = 0; h < n; h++) {
            fp[h] = h;
            dec[h] = 0;
        }
        expected = 1;
        for (i = 1; i <= k; i++)
            expected *= 2 * i;
        for (j = 0; j < NUM_RULES; j++) {
            CHECK(search(n, fp, 0, NULL, 1, dec, rules[j], &len, NULL) == expected);
            CHECK(len == 5 * k);
            for (i = 0; i < len; i++)
                CHECK(sequence[i] == edge_block[i % 5]);
        }
    }

    /* the same with 20 edges, without all leaves only: 2^20 20! leaves */
    current_case = "search, 20 isolated edges with one class value";
    n = 40;
    for (h = 0; h < n; h++) {
        fp[h] = h;
        dec[h] = 0;
    }
    for (j = 0; j < NUM_RULES; j++) {
        CHECK(iso_search(n, fp, 0, NULL, 1, dec, rules[j], 0, sequence, &len, relabelling, NULL, NULL,
                         search_work) == 0);
        CHECK(len == 100);
        for (i = 0; i < len; i++)
            CHECK(sequence[i] == edge_block[i % 5]);
        if (j == 0)
            for (h = 0; h < n; h++)
                CHECK(relabelling[h] == h);
    }

    /* two to four relabelled copies of a random component, its data and
     * its classes copied along, with one or two class values */
    current_case = "search, copies of a component with classes";
    for (trial = 0; trial < 200; trial++) {
        size = 2 * (1 + rng_below(3));
        k = 2 + rng_below(3);
        nd = rng_below(2);
        nvals = 1 + rng_below(2);
        CHECK(random_connected(size / 2, rng_below(3), comp));
        for (h = 0; h < nd * size; h++)
            cdec[h] = rng_below(2);
        for (h = 0; h < size; h++)
            cdec[nd * size + h] = rng_below(nvals);
        n = k * size;
        for (c = 0; c < k; c++) {
            random_relabelling(size, rc);
            relabel_fp(size, comp, rc, fp2);
            relabel_arrays(size, nd + 1, cdec, rc, cdec2);
            for (h = 0; h < size; h++) {
                fp[c * size + h] = fp2[h] == -1 ? -1 : fp2[h] + c * size;
                for (i = 0; i < nd + 1; i++)
                    dec[i * n + c * size + h] = cdec2[i * size + h];
            }
        }
        for (j = 0; j < NUM_RULES; j++) {
            leaves = search(n, fp, nd, dec, 1, dec + nd * n, rules[j], &len, NULL);
            if (j == 0)
                all_leaves = leaves;
            CHECK(leaves == all_leaves);
            CHECK(sequence[0] == size);

            random_relabelling(n, r);
            relabel_fp(n, fp, r, fp2);
            relabel_arrays(n, nd + 1, dec, r, dec2);
            for (i = nd * n; i < (nd + 1) * n; i++)
                dec2[i] = 1 - dec2[i];
            memcpy(sequence2, sequence, (size_t) len * sizeof(int));
            leaves2 = search(n, fp2, nd, dec2, 1, dec2 + nd * n, rules[j], &len2, NULL);
            CHECK(len2 == len);
            CHECK(leaves2 == leaves);
            CHECK(memcmp(sequence, sequence2, (size_t) len * sizeof(int)) == 0);
        }
    }
}

/*
 * iso_invariants_differ: the multisets of the data values and the multisets
 * of the sizes of the classes, on the active half-edges only.
 */
static void test_invariants(void)
{
    /* three isolated edges, then two edges in labels 0, 1 and 4, 5 */
    static const int fp6[] = {0, 1, 2, 3, 4, 5};
    static const int fp4[] = {1, 0, 3, 2};
    static const int fp8[] = {1, 0, -1, -1, 3, 2, -1, -1};
    static const int fp_junk[] = {1, 0, -1, -1};
    int inv_work[ISO_INVARIANTS_WORK_SIZE(8)];
    int d1[12], d2[12], c1[12], c2[12];
    int i;

    current_case = "invariants, equal multisets";
    {
        /* values renamed and moved: the same multisets */
        static const int data1[] = {0, 1, 0, 1, 0, 1}, data2[] = {1, 0, 1, 1, 0, 0};
        static const int classes1[] = {0, 0, 1, 1, 2, 2}, classes2[] = {5, 3, 3, 0, 0, 5};
        CHECK(iso_invariants_differ(6, fp6, data1, classes1, 6, fp6, data2, classes2, 1, 1, inv_work) == 0);
        CHECK(iso_invariants_differ(6, fp6, data2, classes2, 6, fp6, data1, classes1, 1, 1, inv_work) == 0);
        CHECK(iso_invariants_differ(6, fp6, data1, classes1, 6, fp6, data1, classes1, 1, 1, inv_work) == 0);
    }

    current_case = "invariants, data multisets differ";
    {
        static const int data1[] = {0, 1, 0, 1, 0, 1}, data2[] = {0, 0, 0, 0, 1, 1};
        static const int classes[] = {0, 0, 1, 1, 2, 2};
        CHECK(iso_invariants_differ(6, fp6, data1, classes, 6, fp6, data2, classes, 1, 1, inv_work) == 1);
        CHECK(iso_invariants_differ(6, fp6, data2, classes, 6, fp6, data1, classes, 1, 1, inv_work) == 1);
        /* one value of the extremes of an int */
        for (i = 0; i < 6; i++) {
            d1[i] = data1[i];
            d2[i] = data1[i];
        }
        d1[0] = INT_MIN;
        d2[0] = INT_MAX;
        CHECK(iso_invariants_differ(6, fp6, d1, classes, 6, fp6, d2, classes, 1, 1, inv_work) == 1);
        CHECK(iso_invariants_differ(6, fp6, d1, classes, 6, fp6, d1, classes, 1, 1, inv_work) == 0);
    }

    current_case = "invariants, class-size multisets differ";
    {
        /* the sizes 2, 2, 2 against 3, 2, 1 and against 6 */
        static const int data[] = {0, 1, 0, 1, 0, 1};
        static const int classes1[] = {0, 0, 1, 1, 2, 2}, classes2[] = {0, 0, 0, 1, 1, 2};
        static const int classes3[] = {4, 4, 4, 4, 4, 4};
        CHECK(iso_invariants_differ(6, fp6, data, classes1, 6, fp6, data, classes2, 1, 1, inv_work) == 1);
        CHECK(iso_invariants_differ(6, fp6, data, classes2, 6, fp6, data, classes1, 1, 1, inv_work) == 1);
        CHECK(iso_invariants_differ(6, fp6, data, classes1, 6, fp6, data, classes3, 1, 1, inv_work) == 1);
        CHECK(iso_invariants_differ(6, fp6, data, classes3, 6, fp6, data, classes1, 1, 1, inv_work) == 1);
        /* the buffer is left zero after a failure on a previous call: the
         * same inputs give the same verdict again */
        CHECK(iso_invariants_differ(6, fp6, data, classes1, 6, fp6, data, classes2, 1, 1, inv_work) == 1);
        CHECK(iso_invariants_differ(6, fp6, data, classes1, 6, fp6, data, classes1, 1, 1, inv_work) == 0);
    }

    current_case = "invariants, class values renamed";
    {
        static const int classes1[] = {0, 0, 1, 1, 2, 3}, classes2[] = {5, 4, 2, 2, 4, 0};
        static const int classes3[] = {0, 0, 1, 2, 2, 3};
        CHECK(iso_invariants_differ(6, fp6, NULL, classes1, 6, fp6, NULL, classes2, 0, 1, inv_work) == 0);
        /* the sizes 2, 2, 1, 1 in both */
        CHECK(iso_invariants_differ(6, fp6, NULL, classes1, 6, fp6, NULL, classes3, 0, 1, inv_work) == 0);
    }

    current_case = "invariants, several arrays";
    {
        /* the second data array and the second class array decide */
        static const int data1[] = {0, 1, 0, 1, 0, 1, 7, 7, 7, 8, 8, 8};
        static const int data2[] = {1, 1, 0, 0, 0, 1, 8, 8, 7, 7, 7, 8};
        static const int data3[] = {1, 1, 0, 0, 0, 1, 8, 8, 8, 7, 7, 8};
        static const int classes1[] = {0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 2, 2};
        static const int classes2[] = {2, 2, 2, 0, 0, 0, 4, 4, 3, 3, 3, 3};
        CHECK(iso_invariants_differ(6, fp6, data1, classes1, 6, fp6, data2, classes1, 2, 2, inv_work) == 0);
        CHECK(iso_invariants_differ(6, fp6, data1, classes1, 6, fp6, data3, classes1, 2, 2, inv_work) == 1);
        CHECK(iso_invariants_differ(6, fp6, data1, classes1, 6, fp6, data2, classes2, 2, 2, inv_work) == 1);
        CHECK(iso_invariants_differ(6, fp6, data1, classes1, 6, fp6, data2, classes2, 2, 1, inv_work) == 0);
    }

    current_case = "invariants, different label lengths";
    {
        /* the four active half-edges 0, 1, 4, 5 of fp8 against the four of
         * fp4; the values on the half-edges 2, 3, 6, 7 are junk */
        static const int data4[] = {0, 0, 1, 1};
        static const int data8[] = {0, 0, 9, 9, 1, 1, 9, 9};
        static const int data8_bad[] = {0, 0, 9, 9, 1, 0, 9, 9};
        static const int classes4[] = {0, 0, 1, 1};
        static const int classes8[] = {3, 3, 0, 0, 6, 6, 0, 0};
        static const int classes8_bad[] = {1, 1, 0, 0, 1, 2, 0, 0};
        CHECK(iso_invariants_differ(4, fp4, data4, classes4, 8, fp8, data8, classes8, 1, 1, inv_work) == 0);
        CHECK(iso_invariants_differ(8, fp8, data8, classes8, 4, fp4, data4, classes4, 1, 1, inv_work) == 0);
        CHECK(iso_invariants_differ(4, fp4, data4, classes4, 8, fp8, data8_bad, classes8, 1, 1, inv_work) == 1);
        CHECK(iso_invariants_differ(8, fp8, data8_bad, classes8, 4, fp4, data4, classes4, 1, 1, inv_work) == 1);
        CHECK(iso_invariants_differ(4, fp4, data4, classes4, 8, fp8, data8, classes8_bad, 1, 1, inv_work) == 1);
        CHECK(iso_invariants_differ(8, fp8, data8, classes8_bad, 4, fp4, data4, classes4, 1, 1, inv_work) == 1);
    }

    current_case = "invariants, inactive half-edges";
    {
        /* the inactive half-edges 2 and 3 carry values that differ between
         * the maps, and class values out of range */
        static const int data1[] = {5, 5, 100, 200}, data2[] = {5, 5, -3, INT_MAX};
        static const int classes1[] = {0, 0, 3, 2}, classes2[] = {1, 1, 100000, -4};
        static const int data3[] = {5, 6, 100, 200};
        CHECK(iso_invariants_differ(4, fp_junk, data1, classes1, 4, fp_junk, data2, classes2, 1, 1, inv_work)
              == 0);
        CHECK(iso_invariants_differ(4, fp_junk, data1, classes1, 4, fp_junk, data3, classes2, 1, 1, inv_work)
              == 1);
    }

    current_case = "invariants, no class array";
    {
        static const int data1[] = {0, 1, 0, 1, 0, 1}, data2[] = {1, 0, 1, 1, 0, 0}, data3[] = {1, 1, 1, 1, 0, 0};
        /* nc == 0: work is not read */
        CHECK(iso_invariants_differ(6, fp6, data1, NULL, 6, fp6, data2, NULL, 1, 0, NULL) == 0);
        CHECK(iso_invariants_differ(6, fp6, data1, NULL, 6, fp6, data3, NULL, 1, 0, NULL) == 1);
        CHECK(iso_invariants_differ(6, fp6, NULL, NULL, 6, fp6, NULL, NULL, 0, 0, NULL) == 0);
    }

    /* random: renamed and moved values of a relabelled copy never differ,
     * and a value moved to another class does */
    current_case = "invariants, random";
    for (i = 0; i < 200; i++) {
        int n = 2 * (1 + rng_below(4)), h, g, nvals = 1 + rng_below(3);
        static const int fp_id[] = {0, 1, 2, 3, 4, 5, 6, 7};
        int r[8], fp_copy[8];
        for (h = 0; h < n; h++) {
            d1[h] = rng_below(nvals);
            c1[h] = rng_below(nvals);
        }
        random_relabelling(n, r);
        relabel_fp(n, fp_id, r, fp_copy);
        relabel_arrays(n, 1, d1, r, d2);
        relabel_arrays(n, 1, c1, r, c2);
        for (h = 0; h < n; h++)
            c2[h] = nvals - 1 - c2[h];
        CHECK(iso_invariants_differ(n, fp_id, d1, c1, n, fp_copy, d2, c2, 1, 1, inv_work) == 0);
        /* one data value changed: the multiset changes */
        g = rng_below(n);
        d2[g] = nvals + 5;
        CHECK(iso_invariants_differ(n, fp_id, d1, c1, n, fp_copy, d2, c2, 1, 1, inv_work) == 1);
    }
}

int main(int argc, char **argv)
{
    rng_state = argc > 1 ? strtoull(argv[1], NULL, 10) : 1;
    test_fixed();
    test_random();
    test_reference();
    test_reference_regular();
    test_reference_paths();
    test_search_fixed();
    test_search_random();
    test_search_tied();
    test_invariants();
    printf("test_isomorphism: %ld checks passed\n", num_checks);
    return 0;
}
