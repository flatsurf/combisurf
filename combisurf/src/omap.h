/*
 * Canonical keys, labellings and automorphisms of decorated maps.
 * Documented in doc/c/omap.rst.
 *
 *       Copyright (C) 2026 Vincent Delecroix
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 */
#ifndef COMBISURF_OMAP_H
#define COMBISURF_OMAP_H

#define OMAP_INLINE static inline

/* Types *********************************************************************/

/* A map: its vertex permutation vp and face permutation fp on n half-edges. */
typedef struct
{
    int *vp;
    int *fp;
    int n;
}
omap_struct;

typedef omap_struct omap_t[1];

/* A map with nd data and nc class arrays; casts to an omap_struct *. */
typedef struct
{
    omap_struct map;
    const int *data;
    const int *classes;
    int nd;
    int nc;
}
domap_struct;

typedef domap_struct domap_t[1];

/* Return value of _domap_canonical on a map that is not connected. */
enum
{
    DOMAP_ISO_DISCONNECTED = -1
};

/* The rules restricting the roots: a base rule, possibly with the flag. */
typedef enum
{
    DOMAP_ISO_ROOTS_ALL = 0,
    DOMAP_ISO_ROOTS_HALF_EDGE_MAX = 1,
    DOMAP_ISO_ROOTS_HALF_EDGE_FEWEST = 2,
    DOMAP_ISO_ROOTS_EDGE_MAX = 3,
    DOMAP_ISO_ROOTS_EDGE_FEWEST = 4,
    DOMAP_ISO_ROOTS_DATA = 8
} domap_iso_root_rule;

/* Connected maps ************************************************************/

/* Upper bound on the length of a key, in ints. */
#define DOMAP_ISO_KEY_SIZE(n, nd, nc) ((n) * (1 + (nd) + (nc)))

/* Size of the work buffer of _domap_canonical, in ints. */
#define DOMAP_CANONICAL_WORK_SIZE(n, nd, nc) \
    (9 * (n) + 1 + DOMAP_ISO_KEY_SIZE(n, nd, nc))

/* Write into vp the vertex permutation of the map of face permutation fp. */
void _omap_vp_from_fp(int *vp, const int *fp, int n);

/* Return the number of active edges of p, the vertex or face permutation. */
int _omap_num_edges(const int *p, int n);

OMAP_INLINE int
omap_num_edges(const omap_t m)
{
    return _omap_num_edges(m->fp, m->n);
}

/* Return the half-edge after h around its vertex, or -1 on an invalid h. */
OMAP_INLINE int
omap_next_at_vertex(const omap_t m, int h)
{
    if (h < 0 || h >= m->n || m->fp[h] == -1)
        return -1;
    return m->vp[h];
}

/* Return the half-edge before h around its vertex, or -1 on an invalid h. */
OMAP_INLINE int
omap_previous_at_vertex(const omap_t m, int h)
{
    if (h < 0 || h >= m->n || m->fp[h] == -1)
        return -1;
    return m->fp[m->fp[h ^ 1] != -1 ? h ^ 1 : h];
}

/* Return the half-edge after h in its face, or -1 on an invalid h. */
OMAP_INLINE int
omap_next_in_face(const omap_t m, int h)
{
    if (h < 0 || h >= m->n || m->fp[h] == -1)
        return -1;
    return m->fp[h];
}

/* Return the half-edge before h in its face, or -1 on an invalid h. */
OMAP_INLINE int
omap_previous_in_face(const omap_t m, int h)
{
    int x;

    if (h < 0 || h >= m->n || m->fp[h] == -1)
        return -1;
    x = m->vp[h];
    return m->fp[x ^ 1] != -1 ? x ^ 1 : x;
}

/* Size of the work buffer of _omap_relabelling_from, in ints. */
#define OMAP_RELABELLING_WORK_SIZE(n) (n)

/* Relabel m from root into r and rinv, walking the cycles of p, the vertex or
   face permutation; return the number of labels used. */
int _omap_relabelling_from(int *r, int *rinv, const int *p, int n, int root,
                           int *queue);
int omap_relabelling_from(int *r, int *rinv, const omap_t m, int root);

/* Size of the work buffer of _domap_iso_key, in ints. */
#define DOMAP_ISO_KEY_WORK_SIZE(n) (n)

/* Write the key of m under the relabelling (r, rinv); return its length. */
int _domap_iso_key(int *key, const domap_t m, const int *r, const int *rinv,
                   int k, int *names);
int domap_iso_key(int *key, const domap_t m, const int *r, const int *rinv,
                  int k);

/* Size of the work buffer of _domap_invariants_differ, in ints. */
#define DOMAP_INVARIANTS_WORK_SIZE(n) (2 * (n) + 1)

/* Return 1 when an invariant of m1 and m2 differs, 0 otherwise. */
int _domap_invariants_differ(const domap_t m1, const domap_t m2, int *work);
int domap_invariants_differ(const domap_t m1, const domap_t m2);

/* Write the minimal key of m over its roots; return its length, or
   DOMAP_ISO_DISCONNECTED. */
int _domap_canonical(int *key, int *relabelling, int *roots, int *nroots,
                     const domap_t m, domap_iso_root_rule rule, int *work);
int domap_canonical(int *key, int *relabelling, int *roots, int *nroots,
                    const domap_t m, domap_iso_root_rule rule);

/* Maps that may be disconnected *********************************************/

/* Upper bound on the length of a sequence of blocks, in ints. */
#define DOMAP_ISO_SEQUENCE_SIZE(n, nd, nc) \
    (DOMAP_ISO_KEY_SIZE(n, nd, nc) + (n) / 2)

/* Size of the work buffer of _domap_iso_search, in ints. */
#define DOMAP_ISO_SEARCH_WORK_SIZE(n, nd, nc)                               \
    ((5 + (nc)) * (n) + (5 + (nc)) * ((n) / 2 + 1) + (nc)                   \
     + DOMAP_ISO_SEQUENCE_SIZE(n, nd, nc) + DOMAP_ISO_KEY_SIZE(n, nd, nc)   \
     + 1 + 4 * (n) + 1)

/* Callback of _domap_iso_search, called on a leaf. */
typedef int (*domap_iso_leaf_callback)(const int *leaf, int first, void *arg);

/* Write the minimal sequence of blocks of m; return 0, or the nonzero value
   of the callback. */
int _domap_iso_search(int *sequence, int *length, int *relabelling,
                      const domap_t m, domap_iso_root_rule rule,
                      int all_leaves, domap_iso_leaf_callback callback,
                      void *arg, int *work);
int domap_iso_search(int *sequence, int *length, int *relabelling,
                     const domap_t m, domap_iso_root_rule rule,
                     int all_leaves, domap_iso_leaf_callback callback,
                     void *arg);

/* Maps with classes *********************************************************/

/* The components of a map, grouped into kinds of isomorphic components. */
typedef struct
{
    const int *kind;
    const int *off;
    const int *lab;
    const int *start;
    const int *coff;
    const int *csize;
    const int *cfp;
    const int *cra;
    const int *crainv;
    const int *loff;
    const int *labels;
    int ncomp;
}
domap_wreath_struct;

typedef domap_wreath_struct domap_wreath_t[1];

/* Callback of _domap_stabiliser, called on each generator. */
typedef int (*domap_stabiliser_callback)(const int *g, void *arg);

/* Size of the work buffer of _domap_stabiliser, in ints. */
#define DOMAP_STABILISER_WORK_SIZE(n, nc) \
    (3 * (nc) * (n) + 5 * ((n) / 2) + 8 * (n))

/* Write the orbit sizes of the stabiliser of the classes; return 0, or the
   nonzero value of the callback. */
int _domap_stabiliser(int *orbits, const domap_t m, const domap_wreath_t w,
                      domap_stabiliser_callback callback, void *arg,
                      int *work);
int domap_stabiliser(int *orbits, const domap_t m, const domap_wreath_t w,
                     domap_stabiliser_callback callback, void *arg);

#endif
