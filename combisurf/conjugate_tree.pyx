# distutils: sources = combisurf/src/conjugate_tree.c
# distutils: include_dirs = combisurf/src
r"""
Conjugate trees

A conjugate tree is a generalization of suffix trees from word combinatorics.
It is a compact data-structure that contains all factors of all conjugates
of a finite list of words.

The data structures and algorithms in this module are crucial to compute the
geometric intersection numbers of curves and multicurves on surfaces.

More precisely, the leaves of a conjugate tree are in bijection with the
conjugates. And each edge is labeled with a finite word that makes it a
deterministic automaton (where all vertices of degree 2 have been
removed).

EXAMPLES:

The main class from this module is :class:`ConjugateTree` which can be
initialized with no argument::

    sage: from combisurf.conjugate_tree import ConjugateTree
    sage: T = ConjugateTree()

To populate a conjugate tree one uses the function :meth:`~ConjugateTree.process` that
takes as argument a word on non-negative integers (given as a list)::

    sage: T.process([0])
    1
    sage: T.process([0, 1, 0, 0, 1])
    1
    sage: T.process([1, 0, 1, 0])
    2
    sage: T.process([0, 1])
    -2

The output value of :meth:`~ConjugateTree.process` is an integer: the exponent
of the word when it is new, and otherwise minus the index of the word of the
tree it is conjugate to a power of.

To get a hand on the structure of the tree, one can use the following functions
(the root is always index 0 and is omitted in the output)::

    sage: T.leaves()
    [1, 3, 4, 6, 8, 10, 12, 14]
    sage: T.internal_states()
    [2, 5, 7, 9, 11, 13]

To convert the index of a leaf to a conjugate of the words used to populate
one uses::

    sage: T.leaf_as_conjugate(6)
    (1, 2)

Which means that the leaf index ``6`` coressponds to the word number ``1``
(ie ``[0, 1, 0, 0, 1]``) shifted twice.

Telling the tree how large the alphabet is lets it index the children of a
node by letter instead of walking through them, and telling it how many nodes
to expect lets it allocate them all at once::

    sage: T = ConjugateTree(4, reserve=2 * 9 + 1)
    sage: T.process([0, 2, 0, 0, 3])
    1
    sage: T.process([1, 3, 1, 2])
    1
    sage: T.algorithm()
    'dense'

Neither is required and neither changes an answer; see :class:`ConjugateTree`.

The tree itself is implemented in plain C, in ``combisurf/src/conjugate_tree.c``,
so that it can be tested and benchmarked from C programs without Python; this
module wraps it.

.. SEEALSO::

    :mod:`combisurf.conjugate_tree_naive` holds
    :class:`~combisurf.conjugate_tree_naive.ConjugateTreeNaive`, a pure
    Python conjugate tree answering exactly the same questions with the same
    node numbering. It is the reference this one is tested against.
"""
# ****************************************************************************
#  This file is part of combisurf
#
#       Copyright (C) 2026 Vincent Delecroix
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# ****************************************************************************

from cpython cimport array
from libc.string cimport memcpy, memset

from combisurf.word import word_init


cdef array.array _int_array = array.array('i', [])


cdef array.array _as_int_array(w):
    r"""
    Return ``w`` if it is an array of typecode ``'i'`` and a copy of it as
    one otherwise.
    """
    if type(w) is array.array and (<array.array> w).ob_descr.typecode == b'i':
        return <array.array> w
    return array.array('i', w)


cdef class ConjugateTree:
    r"""
    Tree structure to store all conjugates of a finite set of primitive words.

    The data structure works with words over non-negative integers.  The nodes
    are encoded with integers from 0 to the number of nodes minus one. The root
    always get the index ``0`` and created nodes gets the first available index
    (nodes are never deleted). In all algorithms, a node index is often denoted
    by a variable ``s``.

    INPUT:

    - ``alphabet`` -- (default: ``0``) the number of letters, when it is
      known. The words handed to :meth:`process` are then required to be
      words on ``{0, 1, ..., alphabet - 1}``. The default ``0`` means that
      nothing is known about the alphabet, and any non-negative letter is
      accepted.

    - ``reserve`` -- (default: ``0``) the number of nodes to allocate up
      front. This is a hint and nothing more: the tree grows on demand
      whatever it is given. A tree over words of total length ``T`` has at
      most ``2 T + 1`` nodes, so that value makes a reallocation of the nodes
      impossible (with ``algorithm='rows'`` a node that reaches 16 children
      may still allocate its row).

    - ``algorithm`` -- (default: ``None``) how to hold the children of a node,
      either ``'dense'`` (one slot per letter, which needs ``alphabet``),
      ``'sparse'`` (a linked list of siblings) or ``'rows'`` (a linked list
      of siblings, and one slot per letter for the nodes with at least 16
      children when ``alphabet`` is given). The default takes ``'dense'``
      for an alphabet of at most 32 letters and ``'rows'`` otherwise. It is
      a time against memory trade-off and changes no answer.

    EXAMPLES::

        sage: from combisurf.conjugate_tree import ConjugateTree

        sage: T = ConjugateTree()
        sage: T.process([0, 1, 0, 0, 1])
        1
        sage: T
        ConjugateTree with 9 states, 5 leaves and 11 implicit nodes

    The alphabet is checked when it is declared::

        sage: T = ConjugateTree(2)
        sage: T.process([0, 1, 2])
        Traceback (most recent call last):
        ...
        ValueError: invalid word: letter 2 not in the alphabet {0, 1, ..., 1}
    """
    def __cinit__(self, *args, **kwds):
        r"""
        Set up an empty tree; ``__init__`` reads the arguments.

        TESTS::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: ConjugateTree(2).num_states()
            1
        """
        memset(&self.T, 0, sizeof(ct_tree))
        # no letter yet, as after ct_init
        self.T.max_letter = -1

    def __dealloc__(self):
        ct_free(&self.T)

    def __init__(self, alphabet=0, reserve=0, algorithm=None):
        r"""
        TESTS::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: ConjugateTree(-1)
            Traceback (most recent call last):
            ...
            ValueError: alphabet (=-1) must be a non-negative integer
            sage: ConjugateTree(reserve=-1)
            Traceback (most recent call last):
            ...
            ValueError: reserve (=-1) must be a non-negative integer
            sage: ConjugateTree(algorithm='dense')
            Traceback (most recent call last):
            ...
            ValueError: the 'dense' algorithm needs the alphabet
            sage: ConjugateTree(2, algorithm='bogus')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be None, 'dense', 'sparse' or 'rows'
        """
        cdef int n = alphabet
        cdef int r = reserve
        cdef int d
        if n < 0:
            raise ValueError(f"alphabet (={alphabet}) must be a non-negative integer")
        if r < 0:
            raise ValueError(f"reserve (={reserve}) must be a non-negative integer")

        if algorithm is None:
            d = CT_LAYOUT_DEFAULT
        elif algorithm == 'dense':
            if n == 0:
                raise ValueError("the 'dense' algorithm needs the alphabet")
            d = CT_LAYOUT_DENSE
        elif algorithm == 'sparse':
            d = CT_LAYOUT_SPARSE
        elif algorithm == 'rows':
            d = CT_LAYOUT_ROWS
        else:
            raise ValueError("algorithm must be None, 'dense', 'sparse' or 'rows'")

        ct_free(&self.T)
        cdef int err = ct_init(&self.T, n, r, d)
        if err:
            self._raise(err, None)

    def __repr__(self):
        return "ConjugateTree with {} states, {} leaves and {} implicit nodes".format(self.num_states(), len(self.leaves()), self.size())

    cdef int _raise(self, int err, array.array w) except -1:
        r"""
        Raise the exception for the error code ``err`` of the C library,
        returned on the word ``w`` (or ``None``).
        """
        cdef int j, n = self.T.alphabet_size
        if err == CT_ENOMEM:
            raise MemoryError
        if err == CT_EEMPTY:
            raise ValueError("empty word in input")
        if err == CT_ENEGATIVE:
            raise ValueError("invalid word: must be made of non-negative integers")
        if err == CT_EALPHABET:
            letter = next(a for a in w if a >= n)
            raise ValueError(f"invalid word: letter {letter} not in the alphabet "
                             f"{{0, 1, ..., {n - 1}}}")
        if err == CT_ETOOLARGE:
            raise OverflowError("conjugate tree too large for int indices or for a dense transition table")
        if err == CT_EINVALID:
            raise ValueError("invalid argument")
        raise RuntimeError(f"conjugate tree: {ct_strerror(err).decode()} (error code {err})")

    cdef int _reserve(self, int words, int letters) except -1:
        r"""
        Make room so that adding ``words`` more words of ``letters`` letters
        in all cannot fail, as ``ct_reserve``.

        The C functions are only linked into this module, so another
        extension reaches them through this method.
        """
        cdef int err = ct_reserve(&self.T, words, letters)
        if err:
            self._raise(err, None)
        return 0

    cdef int _process(self, const int *w, int length, int *result) except -1:
        r"""
        Add the word ``w[:length]`` and write to ``result`` what
        :meth:`process` returns, as ``ct_process``.

        The C functions are only linked into this module, so another
        extension reaches them through this method.
        """
        cdef int err = ct_process(&self.T, w, length, result)
        cdef int j
        if err:
            self._raise(err, array.array('i', [w[j] for j in range(max(length, 0))]))
        return 0

    # ------------------------------------------------------------------
    # words
    # ------------------------------------------------------------------

    def alphabet(self):
        r"""
        Return the size of the alphabet this tree was given, or ``0`` if it
        was not given one.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: ConjugateTree().alphabet()
            0
            sage: ConjugateTree(8).alphabet()
            8
        """
        return self.T.alphabet_size

    def algorithm(self):
        r"""
        Return how the children of a node are held, ``'dense'``,
        ``'sparse'`` or ``'rows'``.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: ConjugateTree().algorithm()
            'rows'
            sage: ConjugateTree(8).algorithm()
            'dense'
            sage: ConjugateTree(256).algorithm()
            'rows'
            sage: ConjugateTree(256, algorithm='dense').algorithm()
            'dense'
            sage: ConjugateTree(256, algorithm='sparse').algorithm()
            'sparse'
        """
        if self.T.layout == CT_LAYOUT_DENSE:
            return 'dense'
        if self.T.layout == CT_LAYOUT_SPARSE:
            return 'sparse'
        return 'rows'

    def num_words(self):
        r"""
        Return the number of primitive words that define this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,0,1])
            2
            sage: T.num_words()
            1
        """
        return self.T.nwords

    def word(self, i):
        r"""
        Return the ``i``-th primitive word of this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,0,0,1])
            1
            sage: T.word(0)
            array('i', [0, 1, 0, 0, 1])
        """
        cdef int j = i
        if j < 0 or j >= self.T.nwords:
            raise ValueError(f"i (={i}) must be the index of a word")
        cdef int l = self.T.wlen[j]
        cdef array.array ans = array.clone(_int_array, l, False)
        memcpy(ans.data.as_ints, self.T.wbuf + self.T.wstart[j], l * sizeof(int))
        return ans

    def word_length(self, i):
        r"""
        Return the length of the ``i``-th primitive word of this conjugate
        tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,0,0,1])
            1
            sage: T.word_length(0)
            5
        """
        cdef int j = i
        if j < 0 or j >= self.T.nwords:
            raise ValueError(f"i (={i}) must be the index of a word")
        return self.T.wlen[j]

    def words(self):
        r"""
        Return the list of primitive words that define this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0,1,0,1])
            2
            sage: T.process([0,1,0,0,1])
            1
            sage: T.process([0,1])
            0
            sage: T.process([0,1,0,1,1])
            1
            sage: T.words()
            [array('i', [0, 1]), array('i', [0, 1, 0, 0, 1]), array('i', [0, 1, 0, 1, 1])]
        """
        return [self.word(i) for i in range(self.T.nwords)]

    def _letter_at(self, i, k):
        r"""
        Return the ``k``-th letter of the ``i``-th word, read cyclically.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0, 4, 2, 3])
            1
            sage: T._letter_at(0, 1)
            4
            sage: T._letter_at(0, 19)
            3
            sage: T._letter_at(0, -1)
            3
        """
        cdef int j = i
        if j < 0 or j >= self.T.nwords:
            raise ValueError(f"i (={i}) must be the index of a word")
        return ct_letter(&self.T, j, k)

    def _slice(self, i, k, p):
        r"""
        Return the slice from ``k`` to ``p`` of the ``i``-th word.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0,1,2])
            1
            sage: T._slice(0, 15, 19)
            [0, 1, 2, 0]
            sage: T._slice(0, 16, 23)
            [1, 2, 0, 1, 2, 0, 1]
            sage: T._slice(0, 17, 21)
            [2, 0, 1, 2]
        """
        if k < 0 or p < 0:
            raise ValueError(f"k(={k}) and p=({p}) must be non-negative integers")
        cdef int j = i
        if j < 0 or j >= self.T.nwords:
            raise ValueError(f"i (={i}) must be the index of a word")
        cdef int a = k, b = p, t
        return [ct_letter(&self.T, j, t) for t in range(a, b)]

    # ------------------------------------------------------------------
    # transitions
    # ------------------------------------------------------------------

    def transitions(self, s):
        r"""
        Return the children of the node ``s`` as a dictionary mapping the
        first letter of a transition to its target, ordered by letter.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0,1,1])
            1
            sage: T.transitions(0)
            {0: 1, 1: 3}
            sage: T.transitions(1)
            {}
        """
        cdef int node = s
        if node < 0 or node >= self.T.nstates:
            raise ValueError(f"s (={s}) must be a node")
        cdef int c, t
        ans = {}
        if self.T.layout == CT_LAYOUT_DENSE:
            for c in range(self.T.alphabet_size):
                t = self.T.trans[node * self.T.alphabet_size + c]
                if t != -1:
                    ans[c] = t
        else:
            # NOTE: sorted, so that the answer does not depend on the order
            # the siblings happen to sit in
            pairs = []
            t = self.T.fchild[node]
            while t != -1:
                pairs.append((ct_letter(&self.T, self.T.tword[t], self.T.tstart[t]), t))
                t = self.T.nsib[t]
            pairs.sort()
            ans = dict(pairs)
        return ans

    # ------------------------------------------------------------------
    # the tree
    # ------------------------------------------------------------------

    def num_states(self):
        r"""
        Return the number of states.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,0,0,1])
            1
            sage: T.num_states()
            7
        """
        return self.T.nstates

    def size(self):
        r"""
        Return the size of this conjugate tree.

        The size is the number of implicit states where each leaf accounts for
        1.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0,0,0,0,1])
            1
            sage: T.size()
            9
            sage: T.process([0,0,0,1])
            1
            sage: T.size()
            29

        The total size is the same if we input the two words in the opposite order::

            sage: T = ConjugateTree()
            sage: T.process([0,0,0,1])
            1
            sage: T.size()
            7
            sage: T.process([0,0,0,0,1])
            1
            sage: T.size()
            29
        """
        return ct_size(&self.T)

    def internal_states(self):
        r"""
        Return the internal states in this conjugate tree.

        Note that any further call to :meth:`process` might change the
        structure of the tree but not the word encoded by a given state. In
        particular, internal states remain internal states.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0,0,0,1])
            1
            sage: T.internal_states()
            [2, 4]
        """
        cdef int s
        return [s for s in range(1, self.T.nstates) if self.T.tend[s] != -1]

    def leaves(self):
        r"""
        Return the leaves in this conjugate tree.

        Note that any further call to :meth:`process` might change the
        structure of the tree but not the word encoded by a given state. In
        particular, leaves remain leaves after an update.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0,1,1])
            1
            sage: T.process([0,1])
            1
            sage: len(T.leaves())
            5
        """
        cdef int s
        return [s for s in range(1, self.T.nstates) if self.T.tend[s] == -1]

    def leaf_as_conjugate(self, s):
        r"""
        Return the pair ``(i, k)`` such that the leaf ``s`` corresponds to the
        ``i``-th word shifted by ``k``.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0,0,0,1,0,2])
            1
            sage: T.process([1,2,1,2,1,1,2])
            1
            sage: T.process([0])
            1
            sage: [T.leaf_as_conjugate(s) for s in T.leaves()] == [(i, k) for i, w in enumerate(T.words()) for k in range(len(w))]
            True

        TESTS::

            sage: T.internal_states()
            [2, 4, 9, 11, 13, 15, 18, 20, 22]
            sage: T.leaf_as_conjugate(2)
            Traceback (most recent call last):
            ...
            ValueError: s (=2) must be a leaf
            sage: T.leaf_as_conjugate(0)
            Traceback (most recent call last):
            ...
            ValueError: s (=0) must be a leaf
            sage: T.leaf_as_conjugate(T.num_states())
            Traceback (most recent call last):
            ...
            ValueError: s (=24) must be a leaf
        """
        cdef int node = s
        cdef int i, k
        if ct_leaf_as_conjugate(&self.T, node, &i, &k):
            raise ValueError(f"s (={s}) must be a leaf")
        return (i, k)

    def _leaf_shift(self, s):
        r"""
        Given a leaf with index ``s`` return the leaf corresponding to its shifted word.

        The function ``leaf_shift`` is a permutation of the leaves of this conjugate
        tree whose orbits represent conjugate words. There is no need for this function
        as each processing of a word provides a cycle of the created leaves (by
        increasing order).

        This function is kept only for testing purposes.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()

        Get the orbit `(1,2,4,5,6,8,10)` from the first word::

            sage: T.process([0, 1, 0, 2, 0, 0, 1])
            1
            sage: T.leaves()
            [1, 2, 4, 5, 6, 8, 10]
            sage: for s in T.leaves():
            ....:     print(f"{s:2} -> {T._leaf_shift(s)}")
             1 -> 2
             2 -> 4
             4 -> 5
             5 -> 6
             6 -> 8
             8 -> 10
            10 -> 1

        Adding a word of length 2 creates a new orbit of length two `(12,14)`::

            sage: T.process([1, 2])
            1
            sage: T.leaves()
            [1, 2, 4, 5, 6, 8, 10, 12, 14]
            sage: for s in T.leaves():
            ....:     print(f"{s:2} -> {T._leaf_shift(s)}")
             1 -> 2
             2 -> 4
             4 -> 5
             5 -> 6
             6 -> 8
             8 -> 10
            10 -> 1
            12 -> 14
            14 -> 12
        """
        # NOTE: in the case the transition to s is made of a single letter
        # we have to go through the tree
        cdef int node = s
        if node < 0 or node >= self.T.nstates:
            raise ValueError(f"s (={s}) must be a node")
        if self.T.tend[node] != -1:
            raise ValueError(f"s(={s}) not a leaf")
        cdef int i = self.T.tword[node]
        cdef int k = self.T.tstart[node]
        cdef int letter = ct_letter(&self.T, i, k)
        cdef int ss = self.T.sl[self.T.parent[node]]

        if ss == -1:
            ss = 0
        else:
            ss = ct_child(&self.T, ss, letter)
        while self.T.tend[ss] != -1:
            k += self.T.tend[ss] - self.T.tstart[ss]
            letter = ct_letter(&self.T, i, k)
            ss = ct_child(&self.T, ss, letter)
        return ss

    def internal_state_word(self, s):
        r"""
        Return the word corresponding to the explicit state ``s``

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0,1,0,0,1])
            1
            sage: T.process([0])
            1
            sage: for s in T.internal_states():
            ....:     print(s, T.internal_state_word(s))
            3 [0]
            5 [0, 1, 0]
            7 [1, 0]
            9 [0, 0]
        """
        cdef int node = s
        if node < 0 or node >= self.T.nstates:
            raise ValueError("s must be a node")
        # NOTE: no path[-1] here; this module is compiled with
        # wraparound=False, under which a negative index on a list is not
        # caught but read out of bounds
        path = []
        while node != 0:
            path.append(node)
            node = self.T.parent[node]
        ans = []
        cdef int i, k, p
        for node in reversed(path):
            i = self.T.tword[node]
            k = self.T.tstart[node]
            p = self.T.tend[node]
            if p == -1:
                p = self.T.wlen[i]
            ans.extend(self._slice(i, k, p))
        return ans

    def cyclically_sorted_leaves(self, order, pivot):
        r"""
        Return the leaves sorted by the order of the letters ``order``, the
        order below a node being turned by ``pivot``.

        The leaves are listed depth first. The children of the root are
        visited by increasing ``order[c]``, where ``c`` is the first letter of
        their label, and the children of an internal node whose label ends
        with the letter ``b`` by increasing ``(order[c] - pivot[b]) % n``.
        So the leaf of a conjugate comes before the leaf of another when
        ``order`` of its first letter is the smaller, ties being broken, at the
        first letter ``c`` where they differ, by ``(order[c] - pivot[b]) % n``
        where ``b`` is the letter before ``c``.

        INPUT:

        - ``order`` -- a permutation of ``{0, 1, ..., n - 1}``, indexed by the
          letters, where ``n`` is larger than every letter of this tree

        - ``pivot`` -- a sequence of ``n`` values in ``{0, 1, ..., n - 1}``,
          indexed by the letters

        OUTPUT: a list of leaves

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0,1,1])
            1
            sage: T.process([0,1])
            1
            sage: T.cyclically_sorted_leaves([0, 1], [1, 0])
            [6, 1, 8, 4, 2]

        With a constant ``pivot``, the children of every node are ordered by
        ``order`` as the ones of the root, which is the lexicographic order of
        the conjugates when ``order`` is the identity; other pivots turn the
        order below each node::

            sage: T = ConjugateTree()
            sage: T.process([0, 1, 2])
            1
            sage: T.process([0, 2, 1, 1])
            1
            sage: leaves = T.cyclically_sorted_leaves([0, 1, 2], [0, 0, 0])
            sage: [T.leaf_as_conjugate(s) for s in leaves]
            [(0, 0), (1, 0), (1, 3), (1, 2), (0, 1), (0, 2), (1, 1)]
            sage: leaves = T.cyclically_sorted_leaves([0, 1, 2], [2, 0, 1])
            sage: [T.leaf_as_conjugate(s) for s in leaves]
            [(1, 0), (0, 0), (1, 3), (1, 2), (0, 1), (1, 1), (0, 2)]

        TESTS::

            sage: T.cyclically_sorted_leaves([], [])
            Traceback (most recent call last):
            ...
            ValueError: order must be non-empty
            sage: T.cyclically_sorted_leaves([0, 2], [0, 1])
            Traceback (most recent call last):
            ...
            ValueError: order must be a permutation of {0, 1}
            sage: T.cyclically_sorted_leaves([1, 1], [0, 1])
            Traceback (most recent call last):
            ...
            ValueError: order must be a permutation of {0, 1}
            sage: T.cyclically_sorted_leaves([0], [0])
            Traceback (most recent call last):
            ...
            ValueError: the letters of this tree do not fit in an alphabet of size 1
            sage: T.cyclically_sorted_leaves([0, 1], [0, 1])
            Traceback (most recent call last):
            ...
            ValueError: the letters of this tree do not fit in an alphabet of size 2
            sage: T.cyclically_sorted_leaves([0, 1, 2], [0, 3, 0])
            Traceback (most recent call last):
            ...
            ValueError: pivot must be a sequence of 3 values in {0, 1, 2}
            sage: T.cyclically_sorted_leaves([0, 1, 2], [0, -1, 0])
            Traceback (most recent call last):
            ...
            ValueError: pivot must be a sequence of 3 values in {0, 1, 2}
            sage: T.cyclically_sorted_leaves([0, 1, 2], [0])
            Traceback (most recent call last):
            ...
            ValueError: pivot must be a sequence of 3 values in {0, 1, 2}
            sage: ConjugateTree.__new__(ConjugateTree).cyclically_sorted_leaves([0, 1], [0, 1])
            Traceback (most recent call last):
            ...
            ValueError: invalid argument
        """
        cdef array.array a_order = self._order_array(order)
        cdef array.array a_pivot = self._pivot_array(pivot, len(a_order))
        cdef array.array a_leaves = array.clone(_int_array, self.T.nstates, False)
        cdef int *out = a_leaves.data.as_ints
        cdef int num = self._sorted_leaves(a_order, a_pivot, out)
        cdef int j
        return [out[j] for j in range(num)]

    def sorted_leaves_as_conjugates(self, order, pivot):
        r"""
        Return the conjugates of the leaves in the order of
        :meth:`cyclically_sorted_leaves`.

        INPUT:

        - ``order``, ``pivot`` -- as in :meth:`cyclically_sorted_leaves`

        OUTPUT: two arrays of typecode ``'i'``, with one entry per leaf: the
        index ``i`` of its word and its shift ``k``, as in
        :meth:`leaf_as_conjugate`

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0, 1, 1])
            1
            sage: T.process([0, 1])
            1
            sage: T.sorted_leaves_as_conjugates([0, 1], [1, 0])
            (array('i', [1, 0, 1, 0, 0]), array('i', [0, 0, 1, 2, 1]))
            sage: [T.leaf_as_conjugate(s) for s in T.cyclically_sorted_leaves([0, 1], [1, 0])]
            [(1, 0), (0, 0), (1, 1), (0, 2), (0, 1)]

        TESTS::

            sage: T.sorted_leaves_as_conjugates([0, 1], [0, 2])
            Traceback (most recent call last):
            ...
            ValueError: pivot must be a sequence of 2 values in {0, 1}
        """
        cdef array.array a_order = self._order_array(order)
        cdef array.array a_pivot = self._pivot_array(pivot, len(a_order))
        cdef array.array a_leaves = array.clone(_int_array, self.T.nstates, False)
        cdef int *out = a_leaves.data.as_ints
        cdef int num = self._sorted_leaves(a_order, a_pivot, out)
        cdef array.array a_word = array.clone(_int_array, num, False)
        cdef array.array a_shift = array.clone(_int_array, num, False)
        cdef int j
        for j in range(num):
            ct_leaf_as_conjugate(&self.T, out[j], a_word.data.as_ints + j, a_shift.data.as_ints + j)
        return (a_word, a_shift)

    cdef array.array _order_array(self, order):
        r"""
        Return ``order`` as an array of typecode ``'i'``.

        Its values are checked by ``ct_sorted_leaves``, and on a failure
        :meth:`_sorted_leaves` says which condition does not hold.
        """
        if not len(order):
            raise ValueError("order must be non-empty")
        return _as_int_array(order)

    cdef array.array _pivot_array(self, pivot, int n):
        r"""
        Return ``pivot`` as an array of typecode ``'i'`` after checking that
        it has ``n`` entries.

        Its values are checked by ``ct_sorted_leaves``, and on a failure
        :meth:`_sorted_leaves` says which condition does not hold.
        """
        cdef array.array a_pivot = _as_int_array(pivot)
        if len(a_pivot) != n:
            raise ValueError(f"pivot must be a sequence of {n} values in "
                             "{%s}" % ", ".join(str(j) for j in range(n)))
        return a_pivot

    cdef int _sorted_leaves(self, array.array order, array.array pivot, int *out) except -1:
        r"""
        Write the leaves in the order of :meth:`cyclically_sorted_leaves` to
        ``out``, which has room for ``self.num_states()`` entries, and return
        their number.

        ``order`` and ``pivot`` must be arrays of typecode ``'i'`` of the
        same length, as returned by :meth:`_order_array` and
        :meth:`_pivot_array`.
        """
        cdef int n = len(order)
        cdef int num
        cdef int err = ct_sorted_leaves(&self.T, order.data.as_ints, pivot.data.as_ints,
                                        n, out, &num)
        if err == CT_EINVALID:
            # NOTE: the arguments are only looked at again on a failure, to
            # say which condition does not hold
            if sorted(order) != list(range(n)):
                raise ValueError("order must be a permutation of {%s}" % ", ".join(str(j) for j in range(n)))
            if self.T.max_letter >= n:
                raise ValueError(f"the letters of this tree do not fit in an alphabet of size {n}")
            if any(p < 0 or p >= n for p in pivot):
                raise ValueError(f"pivot must be a sequence of {n} values in "
                                 "{%s}" % ", ".join(str(j) for j in range(n)))
        if err:
            self._raise(err, None)
        return num

    def graph(self):
        r"""
        Return this conjugate tree as a directed graph, each edge labelled by
        the letters it reads, or by its first letter only if it ends at a
        leaf (such an edge reads an infinite periodic word).

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,1])
            1
            sage: T.graph()
            Digraph on 5 vertices
            sage: sorted(T.graph().edges())
            [(0, 1, '0'), (0, 3, '1'), (3, 2, '1'), (3, 4, '0')]
        """
        from sage.graphs.digraph import DiGraph
        G = DiGraph(self.num_states(), loops=False, multiedges=False)
        cdef int s
        for s in range(self.T.nstates):
            for t in self.transitions(s).values():
                G.add_edge(s, t, self._edge_label(t))
        return G

    def _edge_label(self, t):
        r"""
        Return the label of the edge ending at the node ``t``: the letters it
        reads, or only its first letter if ``t`` is a leaf (such an edge
        reads an infinite periodic word).

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,1])
            1
            sage: T._edge_label(1)
            '0'
            sage: T._edge_label(3)
            '1'

        TESTS::

            sage: T._edge_label(0)
            Traceback (most recent call last):
            ...
            ValueError: t (=0) must be a node
            sage: T._edge_label(-1)
            Traceback (most recent call last):
            ...
            ValueError: t (=-1) must be a node
            sage: T._edge_label(T.num_states())
            Traceback (most recent call last):
            ...
            ValueError: t (=5) must be a node
        """
        cdef int node = t
        if node <= 0 or node >= self.T.nstates:
            raise ValueError(f"t (={t}) must be a node")
        if self.T.tend[node] == -1:
            return str(ct_letter(&self.T, self.T.tword[node], self.T.tstart[node]))
        return ''.join(map(str, self._slice(self.T.tword[node], self.T.tstart[node], self.T.tend[node])))

    def _pprint(self):
        r"""
        Print the transitions of this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,1])
            1
            sage: T._pprint()
             0 --0(i=0, k=0)-->  1
             0 --1(array('i', [1]))-->  3
             3 --0(i=0, k=3)-->  4
             3 --1(i=0, k=2)-->  2
        """
        ans = []
        cdef int s, i, k, p
        for s in range(self.T.nstates):
            transitions = self.transitions(s)
            for letter in sorted(transitions):
                ss = transitions[letter]
                i = self.T.tword[ss]
                k = self.T.tstart[ss]
                p = self.T.tend[ss]
                if p != -1:
                    ans.append(f"{s:2} --{letter}({self.word(i)[k:p]})--> {ss:2}")
                else:
                    ans.append(f"{s:2} --{letter}(i={i}, k={k})--> {ss:2}")
        print("\n".join(ans))

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------

    def _canonize_state(self, s, i, k, p):
        r"""
        Canonize the quadruple ``(s, i, k, p)`` representing
        the (explicit or implicit) state obtained after reading word[i][k:p]
        from s.

        Return a pair ``(s, k)`` (as ``i`` and ``p`` do not change).

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,0,0,1])
            1
            sage: T._canonize_state(0, 0, 0, 0)
            (0, 0)
            sage: T._canonize_state(0, 0, 0, 2)
            (3, 1)
            sage: T._canonize_state(-1, 0, 0, 3)
            (7, 3)

        TESTS::

            sage: T._canonize_state(0, 0, 3, 2)
            Traceback (most recent call last):
            ...
            ValueError: (s, i, k, p) = (0, 0, 3, 2) is not a reference of a state of the tree
            sage: T._canonize_state(0, 0, -1, 2)
            Traceback (most recent call last):
            ...
            ValueError: (s, i, k, p) = (0, 0, -1, 2) is not a reference of a state of the tree
            sage: T = ConjugateTree(3)
            sage: T.process([0, 1])
            1
            sage: T._canonize_state(0, 0, 0, 1)
            (0, 0)
            sage: T.process([2])
            1
            sage: T._canonize_state(1, 1, 0, 1)
            Traceback (most recent call last):
            ...
            ValueError: (s, i, k, p) = (1, 1, 0, 1) is not a reference of a state of the tree
        """
        cdef int ss = s
        cdef int kk = k
        cdef int ii = i
        cdef int pp = p
        if ss < -1 or ss >= self.T.nstates:
            raise ValueError(f"s (={s}) must be a node")
        if ii < 0 or ii >= self.T.nwords:
            raise ValueError(f"i (={i}) must be the index of a word")
        if ct_canonize(&self.T, &ss, ii, &kk, pp):
            raise ValueError(f"(s, i, k, p) = {(s, i, k, p)} is not a reference of a state of the tree")
        return (ss, kk)

    def process(self, w, check=True):
        r"""
        Add the word ``w`` in this conjugate tree.

        The output value is an integer. Depending on its sign it encodes
        different information.

        - a positive ``exponent`` if the word ``w`` is not already present
          and ``exponent`` is the exponent of ``w`` (which is ``1`` if and only
          if ``w`` is primitive). In that case, the number of leaves in the tree
          increases by the period of ``w`` (which is its length divided by the
          exponent).

        - a non-negative ``-index`` if the word ``w`` is already present, that
          is, if it is conjugate to a power of the word of index ``index`` of
          this conjugate tree

        INPUT:

        - ``w`` -- a non-empty word

        - ``check`` -- boolean (default: ``True``); whether to convert ``w``
          with :func:`~combisurf.word.word_init`, so that ``w`` can be any
          input that :func:`~combisurf.word.word_init` accepts, including a
          string. With ``check=False``, ``w`` is used as it is if it is an
          ``array('i')`` and copied with ``array('i', w)`` otherwise, so that
          it can be any iterable of integers. In both cases the letters are
          checked by the tree: a negative letter, or a letter outside of the
          alphabet when it is declared, raises a ``ValueError``

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0, 1, 0, 0, 1])
            1
            sage: T.process([0, 1, 0, 1])
            2
            sage: T.process([1, 0, 0, 1, 0])
            0
            sage: T.process([0, 1, 0, 1, 0, 1])
            -1

        TESTS::

            sage: T.process([])
            Traceback (most recent call last):
            ...
            ValueError: empty word in input
            sage: T.process([0, -1])
            Traceback (most recent call last):
            ...
            ValueError: invalid word: must be made of non-negative integers
        """
        if not w:
            raise ValueError("empty word in input")
        if check:
            w = word_init(w)
        cdef array.array a = _as_int_array(w)
        cdef int result
        cdef int err = ct_process(&self.T, a.data.as_ints, len(a), &result)
        if err:
            self._raise(err, a)
        return result

    # ------------------------------------------------------------------
    # self-checks
    # ------------------------------------------------------------------

    def _check_structural(self):
        r"""
        Check the per-node invariants of this conjugate tree.

        Unlike :meth:`_check_bijection`, these invariants hold at every
        intermediate step of an insertion, not only once it returns, since
        they say nothing about the leaves of the word currently being added.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,0,0,1])
            1
            sage: T._check_structural()
        """
        cdef int s, ss, i, k, p, n
        n = self.T.nstates

        for i in range(self.T.nwords):
            assert self.T.wlen[i] > 0, i

        assert self.T.parent[0] == -1
        assert self.T.tstart[0] == -4, self.T.tstart[0]
        assert self.T.tend[0] == -3, self.T.tend[0]

        for s in range(1, n):
            assert self.T.parent[s] >= 0, s
            assert self.T.tstart[s] >= 0, s
            assert self.T.tend[s] > -2, s
            assert self.T.tword[s] != -2, s

            if self.T.tend[s] != -1:
                # suffix link are only for internal nodes different from the root
                assert self.T.sl[s] != -2, s

        for s in range(n):
            transitions = self.transitions(s)
            for letter, ss in transitions.items():
                assert letter == ct_letter(&self.T, self.T.tword[ss], self.T.tstart[ss])
                assert self.T.parent[ss] == s

            # the leaves should correspond to the -1 states
            if s != 0:
                k = self.T.tstart[s]
                p = self.T.tend[s]
                assert p == -1 or p - k > 0, (s, k, p)

                assert (p == -1) == (not transitions)

                if p != -1:
                    # branching
                    ss = self.T.sl[s]
                    assert s != ss
                    w0 = self.internal_state_word(ss)
                    w1 = self.internal_state_word(s)
                    assert w0 == w1[1:], (ss, w0, s, w1)

    def _check_bijection(self):
        r"""
        Check that the leaves of this conjugate tree are in bijection with
        the conjugates of the words in :meth:`words`.

        This only holds once :meth:`process` has returned; while it is
        running, the word being added does not yet have all of its leaves.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,0,0,1])
            1
            sage: T._check_bijection()
        """
        assert [self.leaf_as_conjugate(s) for s in self.leaves()] == [(i, k) for i, w in enumerate(self.words()) for k in range(len(w))]

    def _check(self):
        r"""
        Check all invariants of this conjugate tree.

        Runs the checks of the C library and, as an independent second
        implementation of the same invariants, :meth:`_check_structural` and
        :meth:`_check_bijection`.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,0,0,1])
            1
            sage: T._check()
        """
        cdef int err = ct_check(&self.T)
        if err:
            raise AssertionError(f"ct_check failed with error code {err}")
        self._check_structural()
        self._check_bijection()

    # ------------------------------------------------------------------
    # plotting
    # ------------------------------------------------------------------

    def plot(self, state_size=.25, xscale=1, yscale=1, reverse=False):
        r"""
        Return a plot of this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([1,0,0,1,0,1,1,0])
            1
            sage: T.plot()
            Graphics object consisting of ... graphics primitives
            sage: T.plot(xscale=1.2, yscale=0.6)
            Graphics object consisting of ... graphics primitives
            sage: T.process([0,0,1,1])
            1
            sage: T.plot()
            Graphics object consisting of ... graphics primitives
        """
        children = [self.transitions(s) for s in range(self.T.nstates)]

        # compute lexicographically sorted leaves
        leaves = []
        queue = [children[0][letter] for letter in sorted(children[0], reverse=reverse)]
        while queue:
            s = queue.pop()
            if self.T.tend[s] == -1:
                leaves.append(s)
            else:
                queue.extend(children[s][letter] for letter in sorted(children[s], reverse=reverse))

        pos = {}
        for i, s in enumerate(leaves):
            pos[s] = (xscale * self.T.dep[self.T.parent[s]] + 1, yscale * i)

        queue = set(range(self.T.nstates))
        queue.difference_update(leaves)
        while queue:
            treated = []
            for s in queue:
                assert children[s], "got a leaf!"
                if any(ss not in pos for ss in children[s].values()):
                    continue
                x = xscale * self.T.dep[s]
                y = sum(pos[ss][1] for ss in children[s].values()) / len(children[s])
                pos[s] = (x, y)
                treated.append(s)
            assert treated
            queue.difference_update(treated)

        from sage.plot.graphics import Graphics
        from sage.plot.circle import circle
        from sage.plot.text import text
        from sage.plot.line import line2d
        import matplotlib as mpl

        cmap = None
        colors = None
        if self.T.nwords == 1:
            colors = ["gainsboro"]
        if self.T.nwords <= 10:
            cmap = mpl.cm.tab10
        elif self.T.nwords <= 20:
            cmap = mpl.cm.tab20
        else:
            raise NotImplementedError
        if colors is None and cmap is not None:
            colors = [tuple(row[:3]) for row in cmap(range(self.T.nwords))]
        G = Graphics()
        for s in range(self.T.nstates):
            if self.T.tend[s] == -1:
                # leaf
                G += circle(pos[s], state_size, color=colors[self.T.tword[s]], fill=True, zorder=1)
            else:
                G += circle(pos[s], state_size, color="silver", fill=True, zorder=1)
                G += circle(pos[s], state_size, color="black", fill=False, zorder=2)
            G += text(str(s), pos[s], color="black", zorder=3)
            for ss in children[s].values():
                G += line2d([pos[s], pos[ss]], color="grey", zorder=0)
                mid = ((pos[s][0]+pos[ss][0])/2, (pos[s][1]+pos[ss][1])/2)
                G += text(self._edge_label(ss), mid, color="blue")
        G.axes(False)
        return G
