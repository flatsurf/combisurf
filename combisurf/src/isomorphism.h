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
 * The roots tried by iso_canonical and iso_search are restricted by a rule
 * (iso_root_rule), as the keyword root_rule of isomorphism_naive. For an
 * active half-edge h, vdeg(h) and fdeg(h) are the degrees of its vertex and
 * of its face, and its twin t is h ^ 1, or h itself on a folded edge; the
 * vertex permutation is vp[fp[t]] = h. A rule gives each active half-edge
 * an invariant: the pair (vdeg(h), fdeg(h)) under the half-edge rules, the
 * tuple (vdeg(h), fdeg(h), vdeg(t), fdeg(t), loop, monogon) under the edge
 * rules, loop (monogon) being 1 when h and t are at the same vertex (in the
 * same face) and 0 otherwise. The half-edges of equal invariant form a
 * class, and only the half-edges of one class are tried: under *_MAX the
 * class of the largest invariant (tuples compared lexicographically), under
 * *_FEWEST the class with the fewest half-edges, ties broken by the largest
 * invariant. The invariants are computed once per call on the whole map.
 * The base rules read fp alone. The flag ISO_ROOTS_DATA, combined with a
 * rule other than ISO_ROOTS_ALL by |, appends to the invariant the values
 * of the data arrays (signed ints, never the class arrays, whose values are
 * renumbered and so not invariant): data_0[h], ..., data_{nd-1}[h] under the
 * half-edge rules, followed by data_0[t], ..., data_{nd-1}[t] under the edge
 * rules; with nd = 0 it changes nothing. Under ISO_ROOTS_ALL every active
 * half-edge is a root.
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
 * The rules restricting the roots, with the values of the constants
 * ROOTS_ALL, ..., ROOTS_EDGE_FEWEST and ROOTS_DATA of isomorphism_naive:
 * a base rule (0 to 4), possibly combined with the flag ISO_ROOTS_DATA by |.
 * The functions below read the base rule as rule & ~ISO_ROOTS_DATA.
 */
typedef enum {
    ISO_ROOTS_ALL = 0,
    ISO_ROOTS_HALF_EDGE_MAX = 1,
    ISO_ROOTS_HALF_EDGE_FEWEST = 2,
    ISO_ROOTS_EDGE_MAX = 3,
    ISO_ROOTS_EDGE_FEWEST = 4,
    ISO_ROOTS_DATA = 8
} iso_root_rule;

/*
 * Upper bound on the length of a key, in ints: n for the face permutation
 * and n for each decoration array.
 */
#define ISO_KEY_SIZE(n, nd, nc) ((n) * (1 + (nd) + (nc)))

/*
 * Size of the work buffer of iso_canonical, in ints: the relabelling, its
 * inverse, the queue and the numbering of the class values (4 n), a
 * candidate key (ISO_KEY_SIZE), under a rule other than ISO_ROOTS_ALL the
 * active half-edges sorted by invariant, the class of each half-edge and
 * the start of each class in that order (3 n + 1), then the relabelling of
 * the best root and the parents of the union-find of the orbits of the
 * automorphisms found (2 n). Under ISO_ROOTS_DATA the order-preserving
 * numbers of the data values (nd n) live in the candidate key, which is not used before the
 * classes are computed, so that the flag needs no more room.
 */
#define ISO_WORK_SIZE(n, nd, nc) (9 * (n) + 1 + ISO_KEY_SIZE(n, nd, nc))

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
 * Try every active half-edge of fp of the class of the rule as a root (every
 * active half-edge under ISO_ROOTS_ALL), in increasing order, and write
 * into key (ISO_KEY_SIZE(n, nd, nc) ints) the minimal key.
 *
 * When relabelling is not NULL (n ints), write there the relabelling from
 * the first root of minimal key (-1 on the half-edges of the inactive
 * edges). When roots is not NULL (n ints), write there the roots of
 * minimal key in increasing order. When nroots is not NULL, store there
 * their number, which is the number of automorphisms of the decorated map
 * under every rule (0 when the map has no active half-edge): an
 * automorphism of a connected map is determined by the image of one
 * half-edge and preserves the invariants, so the roots of minimal key form
 * one orbit of the group. work is a buffer of ISO_WORK_SIZE(n, nd, nc)
 * ints.
 *
 * Return the length of the key (0 when the map has no active half-edge),
 * or ISO_DISCONNECTED when the walk from the first active half-edge does
 * not reach every active edge; nothing is written into the output buffers
 * in that case.
 *
 * The first root is walked and keyed in full. Every other root is walked
 * with its key compared with the best one as the entries of its face
 * permutation part become final (those below the label of the next
 * half-edge to walk), in passes of doubling length: the root is rejected
 * at the first entry larger than the one of the best key, which costs
 * about twice the labels it needed to differ, not n. At the first entry
 * smaller, it becomes the best root and is walked and keyed in full. When
 * the face permutation parts are equal, the data then the class parts are
 * compared as they are built, with the same rejection. A root that ties
 * with or beats the best one costs O(n (1 + nd + nc)).
 *
 * A root that ties with the best one gives an automorphism of the decorated
 * map, which sends the best root to it; the orbits of the group generated
 * by the automorphisms found are kept in a union-find over the half-edges
 * (O(n) to initialise), represented by their smallest half-edge. Every
 * half-edge of an orbit has the same key, so a root whose representative is
 * smaller than itself has the key of an earlier root and is skipped at the
 * cost of one find: either that key was rejected, or the root is in the
 * orbit of the best root. A tie is thus only met outside the orbit of the
 * best root and strictly enlarges the group found: there are at most
 * log2 |Aut| of them, each costing O(n (1 + nd + nc)) for the key and the
 * automorphism. A map whose roots all tie (a regular map) costs
 * O(n log2 |Aut| (1 + nd + nc)) for them instead of O(n^2 (1 + nd + nc)).
 * The roots of minimal key are the orbit of the best root, written by a
 * final scan. Under a rule other than ISO_ROOTS_ALL, the classes cost O(n)
 * more: the half-edges are sorted by invariant by counting sorts on its 2
 * or 6 entries. Under ISO_ROOTS_DATA each data array is first replaced by
 * order-preserving numbers of 0, ..., n - 1: the values minus their minimum
 * when they span fewer than n values (O(n)), else their ranks among the
 * active half-edges (a heap sort, O(n log n) per array). They add nd or
 * 2 nd entries.
 */
int iso_canonical(int n, const int *fp, int nd, const int *data, int nc, const int *classes,
                  iso_root_rule rule, int *key, int *relabelling, int *roots, int *nroots,
                  int *work);

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
 * (ISO_SEQUENCE_SIZE), one candidate block (ISO_KEY_SIZE + 1), and for the
 * rules other than ISO_ROOTS_ALL the active half-edges sorted by invariant,
 * the class of each half-edge, the start of each class in that order and
 * the number of half-edges not labelled yet in each class (4 n + 1). Under
 * ISO_ROOTS_DATA the order-preserving numbers of the data values (nd n) live
 * in the candidate block, which is not used before the classes are computed, so that the
 * flag needs no more room.
 *
 * These sizes are computed in int and overflow for about
 * 2^31 / (13.5 + 2 nd + 3.5 nc) half-edges; the caller checks them first.
 */
#define ISO_SEARCH_WORK_SIZE(n, nd, nc) \
    ((5 + (nc)) * (n) + (4 + (nc)) * ((n) / 2 + 1) + (nc) \
     + ISO_SEQUENCE_SIZE(n, nd, nc) + ISO_KEY_SIZE(n, nd, nc) + 1 + 4 * (n) + 1)

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
 * far), the class of the rule is chosen among the active half-edges not
 * labelled yet (by the invariants of the whole map, and by the numbers of
 * half-edges not labelled yet in each class under *_FEWEST), the
 * half-edges of that class not labelled yet are tried as roots (every
 * active half-edge not labelled yet under ISO_ROOTS_ALL), in increasing
 * order, and only the roots of minimal block are followed,
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
               iso_root_rule rule, int all_leaves, int *sequence, int *length, int *relabelling,
               iso_leaf_callback callback, void *arg, int *work);

#endif
