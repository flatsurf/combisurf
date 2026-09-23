/*
 * Conjugate trees: suffix trees of the conjugates of a finite set of words.
 *
 * This file is part of combisurf. It is plain C99 and does not depend on
 * Python, so that the tree can be tested and benchmarked from standalone C
 * programs. The Cython module combisurf/conjugate_tree.pyx wraps it.
 *
 *       Copyright (C) 2026 Vincent Delecroix
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 */
#ifndef COMBISURF_CONJUGATE_TREE_H
#define COMBISURF_CONJUGATE_TREE_H

#include <stdint.h>

/*
 * Error codes. Every function that can fail returns CT_OK or one of the
 * negative codes below; results go through out-parameters. On a failure the
 * tree is left exactly as it was before the call, except for CT_EINTERNAL,
 * after which the tree records that it is broken and every later call that
 * changes it returns CT_EBROKEN.
 */
enum {
    CT_OK = 0,
    CT_ENOMEM = -1,        /* out of memory */
    CT_EEMPTY = -2,        /* empty word */
    CT_ENEGATIVE = -3,     /* negative letter */
    CT_EALPHABET = -4,     /* letter outside the alphabet */
    CT_ETOOLARGE = -5,     /* tree too large for int indices or for a dense table */
    CT_EINVALID = -6,      /* invalid argument */
    CT_EINTERNAL = -7,     /* internal inconsistency */
    CT_EBROKEN = -8        /* the tree was broken by an earlier CT_EINTERNAL */
};

/*
 * The tree. The nodes are the integers 0, 1, ..., nstates - 1, the root is
 * 0, and a new node takes the first free index (nodes are never removed).
 * All per-node arrays have room for capacity nodes.
 *
 * The edge into a node s reads the letters tstart[s], ..., tend[s] - 1 of
 * the word tword[s], read cyclically, and tend[s] == -1 when s is a leaf.
 * The word of a node is the concatenation of the labels from the root, and
 * dep[s] is its length (unset on leaves).
 */
typedef struct ct_tree {
    int alphabet_size;   /* size of the alphabet, 0 when it is not known */
    int dense;           /* whether children are a flat alphabet-indexed table */
    int max_letter;      /* largest letter seen so far, -1 when none */
    int broken;          /* 0, or the error code that left the tree inconsistent */

    /* the words, laid end to end in one buffer addressed by offsets */
    int *wbuf;
    int wbuf_size;       /* letters in use */
    int wbuf_capacity;   /* letters allocated */
    int *wstart;         /* word -> its offset in wbuf */
    int *wlen;           /* word -> its length */
    int nwords;
    int words_capacity;

    /* one entry per node in each of these */
    int nstates;
    int capacity;
    int *dep;            /* node -> length of the word it reads */
    int *sl;             /* node -> suffix link */
    int *parent;         /* node -> parent */
    int *tword;          /* node -> word of the label of the edge into it */
    int *tstart;         /* node -> start of that label */
    int *tend;           /* node -> end of that label, -1 for a leaf */

    /* children, dense or sibling lists; only one of the two is allocated */
    int *trans;          /* node * alphabet_size + letter -> child, -1 when absent */
    int *fchild;         /* node -> first child, -1 when none */
    int *nsib;           /* node -> next sibling, -1 when last */
} ct_tree;

/* The largest alphabet for which ct_init picks the dense layout by default. */
#define CT_DENSE_MAX_ALPHABET 32

/*
 * Set up an empty tree over an alphabet of the given size (0 when unknown)
 * with room for reserve nodes. dense is 1 for the dense layout, 0 for the
 * sibling lists, and negative for the default choice. On failure *T is left
 * in a state that ct_free accepts.
 */
int ct_init(ct_tree *T, int alphabet, int reserve, int dense);

/* Release the memory of T. Safe on a zero-filled struct; T is zero-filled after. */
void ct_free(ct_tree *T);

/*
 * Add the word w[0..len-1]. *result is the exponent of w when it was not
 * present (1 if and only if w is primitive; only its primitive root is
 * stored), and otherwise minus the index of the word that w is conjugate to a
 * power of.
 */
int ct_process(ct_tree *T, const int *w, int len, int *result);

/*
 * Make room for the given number of words more, of the given number of
 * letters in all, so that adding them needs no allocation: after it,
 * ct_process on valid words of that total length cannot fail. A caller that
 * must add several words or none calls it first.
 */
int ct_reserve(ct_tree *T, int words, int letters);

/*
 * Write to out, which has room for T->nstates entries, the leaves of T in
 * the order below, and their number to *num.
 *
 * order is a permutation of {0, ..., n - 1} indexed by the letters, pivot
 * any map from the letters to {0, ..., n - 1}, and n must exceed every letter
 * of T. The leaves are listed depth first. The children of the root are
 * visited by increasing order[c], where c is the first letter of their label;
 * the children of an internal node whose label ends with the letter b are
 * visited by increasing (order[c] - pivot[b]) mod n.
 */
int ct_sorted_leaves(const ct_tree *T, const int *order, const int *pivot, int n,
                     int *out, int *num);

/* Check the invariants of T; CT_OK or CT_EINTERNAL (or CT_ENOMEM). */
int ct_check(const ct_tree *T);

/* The k-th letter of the i-th word, read cyclically; CT_EINVALID when there
 * is no i-th word. */
int ct_letter(const ct_tree *T, int i, int k);

/* The child of the node s whose label starts with letter, -1 when none;
 * CT_EINVALID when s is not a node or letter not a letter. */
int ct_child(const ct_tree *T, int s, int letter);

/*
 * Canonize the reference (*s, i, *k, p) of the state reached by reading the
 * letters k, ..., p - 1 of the i-th word from the node s (-1 for the node
 * below the root). CT_EINVALID, with *s and *k unchanged, when s is not a
 * node (or -1), i not a word, when not 0 <= k <= p, or when the letters do not
 * lead along the tree.
 */
int ct_canonize(const ct_tree *T, int *s, int i, int *k, int p);

/* The conjugate (*i, *k) of the leaf s: the i-th word shifted by k. */
int ct_leaf_as_conjugate(const ct_tree *T, int s, int *i, int *k);

/* The number of implicit nodes, where each leaf counts for 1; CT_EINVALID on
 * a tree that was not set up. */
int64_t ct_size(const ct_tree *T);

/* A short description of an error code. */
const char *ct_strerror(int code);

#endif
