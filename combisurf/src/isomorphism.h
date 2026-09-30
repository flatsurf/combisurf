/*
 * Canonical keys, canonical labellings and automorphisms of maps given by
 * their face permutation, possibly decorated by integer arrays on the
 * half-edges.
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

/*
 * Maps that may be disconnected: the search over the labellings component
 * after component of isomorphism_naive._search.
 *
 * The block of a component labelled from a root is made of k (its number of
 * labels), its face permutation in its local labels 0, ..., k - 1 (-1 at
 * the labels of inactive half-edges), each data array read in the order of
 * the local labels of its active half-edges, then each class array read in
 * the same order and renumbered by first appearance, continuing the
 * numbering of the components labelled before. A sequence of blocks is
 * stored flat, the blocks one after the other; since equal k and equal face
 * permutations give blocks of equal lengths, the lexicographic order of the
 * flat sequences of ints is the order of the sequences of blocks.
 *
 * A sequence of blocks has at most ne + 2 * ne + (nd + nc) * na ints (one k
 * per component, at most ne components) for ne active edges and na active
 * half-edges, hence at most ISO_SEQUENCE_SIZE(n, nd, nc) ints.
 */
#define ISO_SEQUENCE_SIZE(n, nd, nc) (ISO_KEY_SIZE(n, nd, nc) + (n) / 2)

/*
 * Size of the work buffer of iso_search, in ints: the global and local
 * relabellings with their inverses and the queue (5 n), the numberings of
 * the class values (nc n), the state of each of the at most n / 2 + 1
 * levels of the search (its root, its first label, the position and the
 * length of its block, and the counts of class values before it: (4 + nc)
 * (n / 2 + 1)), the counts of class values (nc), the current sequence
 * (ISO_SEQUENCE_SIZE) and one candidate block (ISO_KEY_SIZE + 1).
 *
 * These sizes are computed in int and overflow for about
 * 2^31 / (9.5 + 2 nd + 3.5 nc) half-edges; the caller checks them first.
 */
#define ISO_SEARCH_WORK_SIZE(n, nd, nc) \
    ((5 + (nc)) * (n) + (4 + (nc)) * ((n) / 2 + 1) + (nc) \
     + ISO_SEQUENCE_SIZE(n, nd, nc) + ISO_KEY_SIZE(n, nd, nc) + 1)

/*
 * Callback of iso_search, called on a leaf: leaf is the relabelling of the
 * leaf (n ints, -1 on the half-edges of the inactive edges), valid during
 * the call only. first is 1 when the leaf is smaller than every leaf
 * reported before, which the caller discards, and 0 when it is equal to
 * the last leaf reported with first = 1. A nonzero return value stops the
 * search, which returns it; the Cython layer returns -1 with a Python
 * exception set.
 */
typedef int (*iso_leaf_callback)(const int *leaf, int first, void *arg);

/*
 * Search the lexicographically minimal sequence of blocks over the
 * labellings of fp component after component, exactly as
 * isomorphism_naive._search: from a state (the components labelled so
 * far), every active half-edge not labelled yet is tried as a root, in
 * increasing order, and only the roots of minimal block are followed,
 * depth first in increasing root order; the component gets the labels
 * following those used. When several roots tie, every one of them is
 * followed if nc > 0 or all_leaves is nonzero, only the first one
 * otherwise. A state whose sequence exceeds the prefix of the best
 * sequence found is cut.
 *
 * Any map is accepted, connected or not; a connected map gives a single
 * block, which is its key preceded by k. The empty map gives the empty
 * sequence and one leaf.
 *
 * The minimal sequence is written into sequence (ISO_SEQUENCE_SIZE(n, nd,
 * nc) ints) and its length into *length. When relabelling is not NULL (n
 * ints), the relabelling of the first leaf reaching the minimal sequence
 * is written there (-1 on the half-edges of the inactive edges). When
 * callback is not NULL, it is called with arg on every leaf that is
 * smaller than the best one so far (first = 1) and, when all_leaves is
 * nonzero, on every leaf equal to the best one so far (first = 0). The
 * leaves reported since the last call with first = 1 are then the leaves
 * reaching the minimal sequence, in the order of the search. work is a
 * buffer of ISO_SEARCH_WORK_SIZE(n, nd, nc) ints.
 *
 * Return 0, or the nonzero value returned by the callback, which stops the
 * search and leaves the outputs undefined.
 */
int iso_search(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
               int all_leaves, int *sequence, int *length, int *relabelling,
               iso_leaf_callback callback, void *arg, int *work);

#endif
