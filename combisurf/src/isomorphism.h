/*
 * Canonical keys, canonical labellings and automorphisms of connected maps
 * given by their face permutation, possibly decorated by integer arrays on
 * the half-edges.
 *
 * This file is part of combisurf. It is plain C99 and does not depend on
 * Python, so that it can be tested from standalone C programs. The Cython
 * module combisurf/isomorphism.pyx wraps it; the pure Python module
 * combisurf/isomorphism_naive.py is its reference.
 *
 * Conventions. A map on n half-edges (n even) is given by its face
 * permutation fp, an array of n ints with -1 on the inactive half-edges.
 * The twin of the half-edge h is h ^ 1. An edge is active when at least one
 * of its half-edges is; a folded edge has one active half-edge. fp must be
 * a permutation of the active half-edges; the functions below do not check
 * it (the Cython layer does).
 *
 * Decorations are given as one array of nd * n ints for the arrays kept
 * exactly ("data", the array i at data + i * n) and one array of nc * n
 * ints for the arrays kept up to a relabelling of their values ("classes",
 * the array i at classes + i * n). The values on the inactive half-edges
 * are never read. The values of a class array must lie in 0, ..., n - 1
 * (renumber them before the call, for instance by first appearance).
 *
 * The relabelling from a root is the breadth-first walk of
 * isomorphism_naive.relabelling_from: the half-edges are taken from a
 * queue that initially holds the root; the face of the half-edge taken is
 * walked, each half-edge met on an edge without label gets the next even
 * label and its twin the next odd one, and the twin is pushed when active.
 *
 * The key of a relabelling using k labels, na of which are on active
 * half-edges, is made of k + (nd + nc) * na ints: the face permutation in
 * the new labels (-1 at the labels of inactive half-edges), then each data
 * array read in the order of the new labels of the active half-edges, then
 * each class array read in the same order and renumbered 0, 1, 2, ... by
 * first appearance. Keys are compared as sequences of signed ints.
 *
 * No function allocates memory: every buffer is handed in by the caller,
 * with the sizes documented below.
 *
 *       Copyright (C) 2026 Vincent Delecroix
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 */
#ifndef COMBISURF_ISOMORPHISM_H
#define COMBISURF_ISOMORPHISM_H

/* Return value of iso_canonical on a map that is not connected. */
enum {
    ISO_DISCONNECTED = -1
};

/*
 * Upper bound on the length of a key, in ints: n for the face permutation
 * and n for each decoration array.
 */
#define ISO_KEY_SIZE(n, nd, nc) ((n) * (1 + (nd) + (nc)))

/* Size of the work buffer of iso_canonical, in ints. */
#define ISO_WORK_SIZE(n, nd, nc) (4 * (n) + ISO_KEY_SIZE(n, nd, nc))

/* Return the number of active edges of fp. */
int iso_num_edges(int n, const int *fp);

/*
 * Write into r (n ints) the relabelling of fp by the breadth-first walk
 * from the active half-edge root, -1 on the half-edges of the edges not
 * reached, and into rinv (n ints) its inverse on the labels used:
 * rinv[i] is the half-edge of label i, or -1 when that half-edge is
 * inactive (the twin of a folded edge); the entries of rinv past the
 * labels used are not written. queue is a work buffer of n ints.
 *
 * Return the number k of labels used; the map is connected if and only if
 * k is twice its number of active edges.
 */
int iso_relabelling_from(int n, const int *fp, int root, int *r, int *rinv, int *queue);

/*
 * Write into key the key of the relabelling (r, rinv) using k labels, as
 * computed by iso_relabelling_from. key has room for ISO_KEY_SIZE(n, nd,
 * nc) ints. names is a work buffer of n ints, all -1 on entry, and all -1
 * again on return.
 *
 * Return the length of the key.
 */
int iso_key(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
            const int *r, const int *rinv, int k, int *key, int *names);

/*
 * Try every active half-edge of fp as a root, in increasing order, and
 * write into key (ISO_KEY_SIZE(n, nd, nc) ints) the minimal key.
 *
 * When relabelling is not NULL (n ints), write there the relabelling from
 * the first root of minimal key (-1 on the half-edges of the inactive
 * edges). When roots is not NULL (n ints), write there the roots of
 * minimal key in increasing order. When nroots is not NULL, store there
 * their number, which is the number of automorphisms of the decorated map
 * (0 when the map has no active half-edge). work is a buffer of
 * ISO_WORK_SIZE(n, nd, nc) ints.
 *
 * Return the length of the key (0 when the map has no active half-edge),
 * or ISO_DISCONNECTED when the walk from the first active half-edge does
 * not reach every active edge; nothing is written into the output buffers
 * in that case.
 */
int iso_canonical(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                  int *key, int *relabelling, int *roots, int *nroots, int *work);

#endif
