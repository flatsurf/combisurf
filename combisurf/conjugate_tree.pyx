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

The output value of :meth:`~ConjugateTree.process` is either a pair ``(False,
exponent)`` if the word is not present or ``(True, position)``.

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
from libc.limits cimport INT_MAX

from combisurf.word import word_init


cdef enum:
    # Above this many letters a dense transition table costs more than all the
    # rest of a node put together, while the nodes stay just as sparse: the
    # mean degree of an internal node is between 2 and 10 for every alphabet
    # from 4 to 256 letters and every word length. So past it we walk the
    # children of a node instead of indexing them.
    DENSE_MAX_ALPHABET = 32


cdef int _push_sorted(int *stack, int top, int *kids, int *keys, int d) noexcept nogil:
    r"""
    Push ``kids[:d]`` onto ``stack`` by decreasing ``keys``, so that popping
    them takes them by increasing key.

    An insertion sort is the right one here: the number of children of an
    internal node is 2 to 10 on average whatever the alphabet and whatever the
    length of the words.
    """
    cdef int j, jj, kid, key
    for j in range(1, d):
        kid = kids[j]
        key = keys[j]
        jj = j - 1
        while jj >= 0 and keys[jj] < key:
            kids[jj + 1] = kids[jj]
            keys[jj + 1] = keys[jj]
            jj -= 1
        kids[jj + 1] = kid
        keys[jj + 1] = key
    for j in range(d):
        stack[top + j] = kids[j]
    return top + d


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
      most ``2 T + 1`` nodes, so that value makes a reallocation impossible.

    - ``algorithm`` -- (default: ``None``) how to hold the children of a node,
      either ``'dense'`` (one slot per letter, which needs ``alphabet``) or
      ``'sparse'`` (a linked list of siblings). The default takes ``'dense'``
      for an alphabet of at most 32 letters and ``'sparse'`` otherwise. It is
      a time against memory trade-off and changes no answer.

    EXAMPLES::

        sage: from combisurf.conjugate_tree import ConjugateTree

        sage: T = ConjugateTree()
        sage: T.process([0, 1, 0, 0, 1])
        1
        sage: T
        SuffixTree with 9 states, 5 leaves and 11 implicit nodes

    The alphabet is checked when it is declared::

        sage: T = ConjugateTree(2)
        sage: T.process([0, 1, 2])
        Traceback (most recent call last):
        ...
        ValueError: invalid word: letter 2 not in the alphabet {0, 1, ..., 1}
    """
    def __cinit__(self, *args, **kwds):
        r"""
        Set up empty buffers; ``__init__`` reads the arguments.

        TESTS::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: ConjugateTree(2).num_states()
            1
        """
        self.alphabet_size = 0
        self.dense = False
        self.max_letter = -1

        self.a_wbuf = array.array('i', [])
        self.wbuf = self.a_wbuf.data.as_ints
        self.wbuf_size = 0
        self.wbuf_capacity = 0
        self.a_wstart = array.array('i', [])
        self.wstart = self.a_wstart.data.as_ints
        self.a_wlen = array.array('i', [])
        self.wlen = self.a_wlen.data.as_ints
        self.nwords = 0
        self.words_capacity = 0

        self.nstates = 0
        self.capacity = 0
        self.a_dep = array.array('i', [])
        self.dep = self.a_dep.data.as_ints
        self.a_sl = array.array('i', [])
        self.sl = self.a_sl.data.as_ints
        self.a_parent = array.array('i', [])
        self.parent = self.a_parent.data.as_ints
        self.a_tword = array.array('i', [])
        self.tword = self.a_tword.data.as_ints
        self.a_tstart = array.array('i', [])
        self.tstart = self.a_tstart.data.as_ints
        self.a_tend = array.array('i', [])
        self.tend = self.a_tend.data.as_ints

        self.a_trans = array.array('i', [])
        self.trans = self.a_trans.data.as_ints
        self.a_fchild = array.array('i', [])
        self.fchild = self.a_fchild.data.as_ints
        self.a_nsib = array.array('i', [])
        self.nsib = self.a_nsib.data.as_ints

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
            ValueError: algorithm must be None, 'dense' or 'sparse'
        """
        cdef int n = alphabet
        cdef int r = reserve
        if n < 0:
            raise ValueError(f"alphabet (={alphabet}) must be a non-negative integer")
        if r < 0:
            raise ValueError(f"reserve (={reserve}) must be a non-negative integer")
        self.alphabet_size = n

        if algorithm is None:
            self.dense = 0 < n <= DENSE_MAX_ALPHABET
        elif algorithm == 'dense':
            if n == 0:
                raise ValueError("the 'dense' algorithm needs the alphabet")
            self.dense = True
        elif algorithm == 'sparse':
            self.dense = False
        else:
            raise ValueError("algorithm must be None, 'dense' or 'sparse'")

        self._reserve_nodes(r if r > 1 else 1)
        self._add_node()
        # NOTE: the root doubles as the transition of length one out of the
        # imaginary node -1 that _canonize starts from, hence its (-4, -3)
        self.dep[0] = 0
        self.sl[0] = -1
        self.parent[0] = -1
        self.tword[0] = 0
        self.tstart[0] = -4
        self.tend[0] = -3

    def __repr__(self):
        return "SuffixTree with {} states, {} leaves and {} implicit nodes".format(self.num_states(), len(self.leaves()), self.size())

    # ------------------------------------------------------------------
    # allocation
    # ------------------------------------------------------------------

    cdef int _reserve_nodes(self, int capacity) except -1:
        r"""
        Make room for ``capacity`` nodes and refresh every C view on the
        buffers, which ``array.resize`` invalidates.
        """
        if capacity <= self.capacity:
            return 0
        if self.dense and capacity > INT_MAX // self.alphabet_size:
            raise OverflowError("conjugate tree too large for a dense transition table")

        array.resize(self.a_dep, capacity)
        self.dep = self.a_dep.data.as_ints
        array.resize(self.a_sl, capacity)
        self.sl = self.a_sl.data.as_ints
        array.resize(self.a_parent, capacity)
        self.parent = self.a_parent.data.as_ints
        array.resize(self.a_tword, capacity)
        self.tword = self.a_tword.data.as_ints
        array.resize(self.a_tstart, capacity)
        self.tstart = self.a_tstart.data.as_ints
        array.resize(self.a_tend, capacity)
        self.tend = self.a_tend.data.as_ints

        if self.dense:
            array.resize(self.a_trans, capacity * self.alphabet_size)
            self.trans = self.a_trans.data.as_ints
        else:
            array.resize(self.a_fchild, capacity)
            self.fchild = self.a_fchild.data.as_ints
            array.resize(self.a_nsib, capacity)
            self.nsib = self.a_nsib.data.as_ints

        self.capacity = capacity
        return 0

    cdef int _reserve_words(self, int num_words, int num_letters) except -1:
        r"""
        Make room for ``num_words`` words and ``num_letters`` letters, and
        refresh the C views. Either bound may be ``0`` to leave it alone.
        """
        if num_words > self.words_capacity:
            array.resize(self.a_wstart, num_words)
            self.wstart = self.a_wstart.data.as_ints
            array.resize(self.a_wlen, num_words)
            self.wlen = self.a_wlen.data.as_ints
            self.words_capacity = num_words
        if num_letters > self.wbuf_capacity:
            array.resize(self.a_wbuf, num_letters)
            self.wbuf = self.a_wbuf.data.as_ints
            self.wbuf_capacity = num_letters
        return 0

    cdef int _add_node(self) except -1:
        r"""
        Append a node with nothing but its transitions initialized and return
        its index.
        """
        cdef int n = self.nstates
        cdef int c
        if n == self.capacity:
            self._reserve_nodes(2 * self.capacity if self.capacity else 8)
        self.nstates = n + 1
        # NOTE: -2 is the code for uninitialized, as in ConjugateTreeNaive
        self.dep[n] = -2
        self.sl[n] = -2
        self.parent[n] = -2
        self.tword[n] = -2
        self.tstart[n] = -2
        self.tend[n] = -2
        if self.dense:
            for c in range(self.alphabet_size):
                self.trans[n * self.alphabet_size + c] = -1
        else:
            self.fchild[n] = -1
            self.nsib[n] = -1
        return n

    cdef int _add_word(self, w) except -1:
        r"""
        Append the letters of ``w`` to the word buffer, checking them.
        """
        cdef int l = len(w)
        cdef int j, letter
        if l == 0:
            raise ValueError("empty word in input")
        if self.nwords == self.words_capacity:
            self._reserve_words(2 * self.words_capacity if self.words_capacity else 4, 0)
        if self.wbuf_size + l > self.wbuf_capacity:
            self._reserve_words(0, 2 * (self.wbuf_size + l))
        for j in range(l):
            letter = w[j]
            if letter < 0:
                raise ValueError("invalid word: must be made of non-negative integers")
            if self.alphabet_size and letter >= self.alphabet_size:
                raise ValueError(f"invalid word: letter {letter} not in the alphabet "
                                 f"{{0, 1, ..., {self.alphabet_size - 1}}}")
            self.wbuf[self.wbuf_size + j] = letter
            if letter > self.max_letter:
                self.max_letter = letter
        self.wstart[self.nwords] = self.wbuf_size
        self.wlen[self.nwords] = l
        self.wbuf_size += l
        self.nwords += 1
        return 0

    cdef void _pop_word(self) noexcept nogil:
        r"""
        Undo the last :meth:`_add_word`.
        """
        self.nwords -= 1
        self.wbuf_size = self.wstart[self.nwords]

    cdef void _truncate_word(self, int i, int size) noexcept nogil:
        r"""
        Keep only the first ``size`` letters of the ``i``-th word.
        """
        self.wlen[i] = size
        if i == self.nwords - 1:
            self.wbuf_size = self.wstart[i] + size

    # ------------------------------------------------------------------
    # words
    # ------------------------------------------------------------------

    cdef inline int _letter(self, int i, int k) noexcept nogil:
        cdef int l = self.wlen[i]
        k %= l
        if k < 0:
            k += l
        return self.wbuf[self.wstart[i] + k]

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
        return self.alphabet_size

    def algorithm(self):
        r"""
        Return how the children of a node are held, ``'dense'`` or
        ``'sparse'``.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: ConjugateTree().algorithm()
            'sparse'
            sage: ConjugateTree(8).algorithm()
            'dense'
            sage: ConjugateTree(256).algorithm()
            'sparse'
            sage: ConjugateTree(256, algorithm='dense').algorithm()
            'dense'
        """
        return 'dense' if self.dense else 'sparse'

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
        return self.nwords

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
        if j < 0 or j >= self.nwords:
            raise ValueError(f"i (={i}) must be the index of a word")
        cdef int start = self.wstart[j]
        return self.a_wbuf[start: start + self.wlen[j]]

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
        if j < 0 or j >= self.nwords:
            raise ValueError(f"i (={i}) must be the index of a word")
        return self.wlen[j]

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
        return [self.word(i) for i in range(self.nwords)]

    def letter(self, i, k):
        r"""
        Return the ``k``-th letter of the ``i``-th word, read cyclically.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0, 4, 2, 3])
            1
            sage: T.letter(0, 1)
            4
            sage: T.letter(0, 19)
            3
            sage: T.letter(0, -1)
            3
        """
        cdef int j = i
        if j < 0 or j >= self.nwords:
            raise ValueError(f"i (={i}) must be the index of a word")
        return self._letter(j, k)

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
        if j < 0 or j >= self.nwords:
            raise ValueError(f"i (={i}) must be the index of a word")
        cdef int a = k, b = p, t
        return [self._letter(j, t) for t in range(a, b)]

    # ------------------------------------------------------------------
    # transitions
    # ------------------------------------------------------------------

    cdef inline int _child(self, int s, int letter) noexcept nogil:
        cdef int t
        if self.dense:
            return self.trans[s * self.alphabet_size + letter]
        t = self.fchild[s]
        while t != -1:
            if self._letter(self.tword[t], self.tstart[t]) == letter:
                return t
            t = self.nsib[t]
        return -1

    cdef inline void _add_child(self, int s, int letter, int t) noexcept nogil:
        if self.dense:
            self.trans[s * self.alphabet_size + letter] = t
        else:
            self.nsib[t] = self.fchild[s]
            self.fchild[s] = t

    cdef void _replace_child(self, int s, int letter, int old, int new) noexcept nogil:
        r"""
        Put ``new`` where ``old`` sits among the children of ``s``.

        In the sparse representation ``old`` is found by its index and not by
        its letter, since the letter of a node is read off the label of the
        edge into it and a caller splitting that edge is about to move it.
        """
        cdef int prev, c
        if self.dense:
            self.trans[s * self.alphabet_size + letter] = new
            return
        prev = -1
        c = self.fchild[s]
        while c != old:
            prev = c
            c = self.nsib[c]
        self.nsib[new] = self.nsib[old]
        if prev == -1:
            self.fchild[s] = new
        else:
            self.nsib[prev] = new
        self.nsib[old] = -1

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
        if node < 0 or node >= self.nstates:
            raise ValueError(f"s (={s}) must be a node")
        cdef int c, t
        ans = {}
        if self.dense:
            for c in range(self.alphabet_size):
                t = self.trans[node * self.alphabet_size + c]
                if t != -1:
                    ans[c] = t
        else:
            # NOTE: sorted, so that the answer does not depend on the order
            # the siblings happen to sit in
            pairs = []
            t = self.fchild[node]
            while t != -1:
                pairs.append((self._letter(self.tword[t], self.tstart[t]), t))
                t = self.nsib[t]
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
        return self.nstates

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
        cdef int ans = 0
        cdef int s, k, p
        for s in range(self.nstates):
            k = self.tstart[s]
            p = self.tend[s]
            if p == -1:
                ans += 1
            else:
                ans += p - k
        return ans

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
        return [s for s in range(1, self.nstates) if self.tend[s] != -1]

    def leaves(self):
        r"""
        Return the leaves in this conjugate tree.

        Note that any further call to :meth:`process` might change the
        structure of the tree but not the word encoded by a given state. In
        particular, leaves remain leaves after an update.

        TESTS::

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
        return [s for s in range(1, self.nstates) if self.tend[s] == -1]

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
        """
        cdef int node = s
        if node < 0 or node >= self.nstates:
            raise ValueError
        cdef int i = self.tword[node]
        cdef int l = self.wlen[i]
        cdef int ans = (self.tstart[node] - self.dep[self.parent[node]]) % l
        if ans < 0:
            ans += l
        return (i, ans)

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
        if node < 0 or node >= self.nstates:
            raise ValueError
        if self.tend[node] != -1:
            raise ValueError(f"s(={s}) not a leaf")
        cdef int i = self.tword[node]
        cdef int k = self.tstart[node]
        cdef int letter = self._letter(i, k)
        cdef int ss = self.sl[self.parent[node]]

        if ss == -1:
            ss = 0
        else:
            ss = self._child(ss, letter)
        while self.tend[ss] != -1:
            k += self.tend[ss] - self.tstart[ss]
            letter = self._letter(i, k)
            ss = self._child(ss, letter)
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
        if node < 0 or node >= self.nstates:
            raise ValueError("s must be a node")
        # NOTE: no path[-1] here; this module is compiled with
        # wraparound=False, under which a negative index on a list is not
        # caught but read out of bounds
        path = []
        while node != 0:
            path.append(node)
            node = self.parent[node]
        ans = []
        cdef int i, k, p
        for node in reversed(path):
            i = self.tword[node]
            k = self.tstart[node]
            p = self.tend[node]
            if p == -1:
                p = self.wlen[i]
            ans.extend(self._slice(i, k, p))
        return ans

    def cyclically_sorted_leaves(self, angles):
        r"""
        Return the leaves sorted using a cyclic ordering of the alphabet given by ``angles``.

        The leaf of a conjugate comes before the leaf of another when the
        angle of its first letter is the smaller, ties being broken by the
        angle turned at each further letter.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0,1,1])
            1
            sage: T.process([0,1])
            1
            sage: T.cyclically_sorted_leaves([0, 1])
            [6, 1, 8, 4, 2]

        TESTS::

            sage: T.cyclically_sorted_leaves([])
            Traceback (most recent call last):
            ...
            ValueError: angles must be non-empty
            sage: T.cyclically_sorted_leaves([0, 2])
            Traceback (most recent call last):
            ...
            ValueError: angles must be a permutation of {0, 1}
            sage: T.cyclically_sorted_leaves([0])
            Traceback (most recent call last):
            ...
            ValueError: the letters of this tree do not fit in an alphabet of size 1
        """
        cdef int n = len(angles)
        if n == 0:
            raise ValueError("angles must be non-empty")
        cdef array.array a_ang = array.array('i', angles)
        cdef int *ang = a_ang.data.as_ints
        cdef int c
        for c in range(n):
            if ang[c] < 0 or ang[c] >= n:
                raise ValueError("angles must be a permutation of {%s}" % ", ".join(str(j) for j in range(n)))
        if self.max_letter >= n or (self.max_letter ^ 1) >= n:
            raise ValueError(f"the letters of this tree do not fit in an alphabet of size {n}")

        cdef array.array a_stack = array.clone(a_ang, self.nstates, False)
        cdef int *stack = a_stack.data.as_ints
        cdef array.array a_kids = array.clone(a_ang, n, False)
        cdef int *kids = a_kids.data.as_ints
        cdef array.array a_keys = array.clone(a_ang, n, False)
        cdef int *keys = a_keys.data.as_ints

        cdef int top = 0
        cdef int s, t, d, key, base
        leaves = []

        # the children of the root, ordered by the angle of their first letter
        d = 0
        if self.dense:
            for c in range(self.alphabet_size):
                t = self.trans[c]
                if t != -1:
                    kids[d] = t
                    keys[d] = ang[c]
                    d += 1
        else:
            t = self.fchild[0]
            while t != -1:
                kids[d] = t
                keys[d] = ang[self._letter(self.tword[t], self.tstart[t])]
                t = self.nsib[t]
                d += 1
        top = _push_sorted(stack, top, kids, keys, d)

        while top:
            top -= 1
            s = stack[top]
            if self.tend[s] == -1:
                leaves.append(s)
                continue
            # further down, the angle is measured from the reverse of the last
            # letter read
            base = ang[self._letter(self.tword[s], self.tend[s] - 1) ^ 1]
            d = 0
            if self.dense:
                for c in range(self.alphabet_size):
                    t = self.trans[s * self.alphabet_size + c]
                    if t != -1:
                        key = ang[c] - base
                        if key < 0:
                            key += n
                        kids[d] = t
                        keys[d] = key
                        d += 1
            else:
                t = self.fchild[s]
                while t != -1:
                    key = ang[self._letter(self.tword[t], self.tstart[t])] - base
                    if key < 0:
                        key += n
                    kids[d] = t
                    keys[d] = key
                    t = self.nsib[t]
                    d += 1
            top = _push_sorted(stack, top, kids, keys, d)

        return leaves

    def graph(self):
        r"""
        Return this conjugate tree as a directed graph, each edge labelled by
        the triple that describes the word it reads.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,1])
            1
            sage: T.graph()
            Digraph on 5 vertices
        """
        from sage.graphs.digraph import DiGraph
        G = DiGraph(self.num_states(), loops=False, multiedges=False)
        cdef int s
        for s in range(self.nstates):
            for t in self.transitions(s).values():
                G.add_edge(s, t, f"({self.tword[t]},{self.tstart[t]},{self.tend[t]})")
        return G

    def pprint(self):
        r"""
        Print the transitions of this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,1])
            1
            sage: T.pprint()
             0 --0(i=0, k=0)-->  1
             0 --1(array('i', [1]))-->  3
             3 --0(i=0, k=3)-->  4
             3 --1(i=0, k=2)-->  2
        """
        ans = []
        cdef int s, i, k, p
        for s in range(self.nstates):
            transitions = self.transitions(s)
            for letter in sorted(transitions):
                ss = transitions[letter]
                i = self.tword[ss]
                k = self.tstart[ss]
                p = self.tend[ss]
                if p != -1:
                    ans.append(f"{s:2} --{letter}({self.word(i)[k:p]})--> {ss:2}")
                else:
                    ans.append(f"{s:2} --{letter}(i={i}, k={k})--> {ss:2}")
        print("\n".join(ans))

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------

    cdef int _test_and_split(self, int s, int i, int k, int p, int letter) except -2:
        r"""
        Internal low-level function that checks whether upon reading one should
        create a branching in the tree.

        Given the canonical reference quadruple ``(s, i, k, p)`` this method tests
        whether adding ``letter`` creates a branching or whether the transition
        already exist. If it does not exist, ensure that the corresponding
        state is explicit.

        Return either ``-1`` if the transition exists or a non-negative integer ``s``
        corresponding to the node from which one needs to create a new transition.
        """
        cdef int t, ii, kk, index, lletter, ss, first
        if k < p:
            # implicit state
            # get the transition from s starting with word[i][k] and test
            # whether its (p - k)-th letter coincide with letter or not
            t = self._child(s, self._letter(i, k))
            ii = self.tword[t]
            kk = self.tstart[t]
            index = kk + p - k
            lletter = self._letter(ii, index)
            if letter == lletter:
                # the node already exists
                return -1
            # make the node explicit
            # the new node ss is the node made explicit
            # s ---> t becomes s --> ss --> t
            first = self._letter(ii, kk)
            ss = self._add_node()

            self.tword[ss] = ii
            self.tstart[ss] = kk
            self.tend[ss] = index
            self.parent[ss] = s
            self.dep[ss] = self.dep[s] + index - kk

            # NOTE: ss takes the place of t under s before the label of t is
            # shortened, since that label is where its first letter is read
            self._replace_child(s, first, t, ss)
            self.tstart[t] = index
            self.parent[t] = ss
            self._add_child(ss, lletter, t)

            return ss
        else:
            # explicit state
            if s == -1 or self._child(s, letter) != -1:
                # the node already exists
                return -1
            return s

    cdef void _canonize(self, int *sp, int i, int *kp, int p) noexcept nogil:
        r"""
        Canonize the quadruple ``(s, i, k, p)`` representing the (explicit or
        implicit) state obtained after reading ``word[i][k:p]`` from ``s``.

        Write the answer back into ``sp`` and ``kp``; ``i`` and ``p`` do not
        change.
        """
        cdef int s = sp[0]
        cdef int k = kp[0]
        cdef int ss, kk, pp
        if k >= p:
            # already explicit
            kp[0] = p
            return
        ss = 0 if s == -1 else self._child(s, self._letter(i, k))
        kk = self.tstart[ss]
        pp = self.tend[ss]
        while pp != -1 and pp - kk < p - k:
            k += pp - kk
            s = ss
            ss = self._child(s, self._letter(i, k))
            kk = self.tstart[ss]
            pp = self.tend[ss]
        if pp != -1 and pp - kk == p - k:
            # explicit
            sp[0] = ss
            kp[0] = p
        else:
            # implicit
            sp[0] = s
            kp[0] = k

    def canonize(self, s, i, k, p):
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
            sage: T.canonize(0, 0, 0, 0)
            (0, 0)
        """
        cdef int ss = s
        cdef int kk = k
        cdef int ii = i
        cdef int pp = p
        if ss < -1 or ss >= self.nstates:
            raise ValueError(f"s (={s}) must be a node")
        if ii < 0 or ii >= self.nwords:
            raise ValueError(f"i (={i}) must be the index of a word")
        self._canonize(&ss, ii, &kk, pp)
        return (ss, kk)

    cdef int _update(self, int *sp, int i, int *kp, int p) except -1:
        r"""
        Low-level internal function that updates by reading one letter.

        Here ``(s, i, k, p)`` should be the canonical reference pair of the
        active state from the previous state. Return the number of leaves
        that were created.
        """
        # (s, k, p): active state which is the first state along the boundary
        # path which is not an active leaf
        # r: closest branching from s (r is either s or its ancestor)
        cdef int s = sp[0]
        cdef int k = kp[0]
        cdef int letter = self._letter(i, p)
        cdef int old_r = 0
        cdef int created = 0
        cdef int r, rr
        r = self._test_and_split(s, i, k, p, letter)
        while r != -1:
            # create a leaf
            rr = self._add_node()
            created += 1
            self._add_child(r, letter, rr)
            self.parent[rr] = r
            self.tword[rr] = i
            self.tstart[rr] = p
            self.tend[rr] = -1
            if old_r != 0:
                self.sl[old_r] = r
            old_r = r
            # NOTE: s is never -1 here; _test_and_split returns -1 on the
            # imaginary node, which ends the loop
            s = self.sl[s]
            self._canonize(&s, i, &k, p)
            r = self._test_and_split(s, i, k, p, letter)

        if old_r != 0:
            self.sl[old_r] = s

        sp[0] = s
        kp[0] = k
        return created

    def process(self, w, check=True, hard_check=False):
        r"""
        Add the word ``w`` in this conjugate tree.

        The output value is an integer. Depending on its sign it encodes
        different information.

        - a positive ``exponent`` if the word ``w`` is not already present
          and ``exponent`` is the exponent of ``w`` (which is ``1`` if and only
          if ``w`` is primitive). In that case, the number of leaves in the tree
          increases by the period of ``w`` (which is its length divided by the
          exponent).

        - a non-negative ``-index`` if the word ``w``is already present and
          ``index`` is the index of the leaf corresponding to ``w`` in this
          conjugate tree

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
        cdef int i = self.nwords
        cdef int l = len(w)
        self._add_word(w)

        cdef int s = 0
        cdef int k = 0
        cdef int p = 0
        cdef int num_leaves = 0
        cdef int ss, ii, pp, exponent

        # To ensure that we find all conjugates we must create as many leaves
        # as the size of w (assuming it is primitive)
        while True:
            if p != k:
                ss = self._child(s, self._letter(i, k))
                ii = self.tword[ss]
                pp = self.tend[ss]
            else:
                ii = -1
                pp = -2
            num_leaves += self._update(&s, i, &k, p)
            if hard_check:
                self._check_structural()
            self._canonize(&s, i, &k, p + 1)

            # halt condition
            if num_leaves == l:
                # w is primitive
                break
            elif ii == i and p >= 2 * l:
                # w is non primitive
                break
            elif p >= l and num_leaves == 0 and ii != -1 and pp == -1 and l % self.wlen[ii] == 0:
                # w is conjugate to a power of the ii-th word
                break

            p += 1

        if num_leaves == 0:
            self._pop_word()
            if hard_check:
                self._check_bijection()
            return -ii
        else:
            if l % num_leaves:
                raise RuntimeError(f"the length (={l}) is not a multiple of the number of new leaves (={num_leaves})")
            exponent = l // num_leaves
            if exponent != 1:
                # NOTE: only store primitive words
                self._truncate_word(i, num_leaves)
            if hard_check:
                self._check_bijection()
            return exponent

    # ------------------------------------------------------------------
    # self-checks
    # ------------------------------------------------------------------

    def _check_structural(self):
        r"""
        Check the per-node invariants of this conjugate tree.

        Unlike :meth:`_check_bijection`, these invariants hold at every
        intermediate step of :meth:`process`, not only once it returns, since
        they say nothing about the leaves of the word currently being added.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,0,0,1])
            1
            sage: T._check_structural()
        """
        cdef int s, ss, i, k, p, n
        n = self.nstates

        for i in range(self.nwords):
            assert self.wlen[i] > 0, i

        assert self.parent[0] == -1
        assert self.tstart[0] == -4, self.tstart[0]
        assert self.tend[0] == -3, self.tend[0]

        for s in range(1, n):
            assert self.parent[s] >= 0, s
            assert self.tstart[s] >= 0, s
            assert self.tend[s] > -2, s
            assert self.tword[s] != -2, s

            if self.tend[s] != -1:
                # suffix link are only for internal nodes different from the root
                assert self.sl[s] != -2, s

        for s in range(n):
            transitions = self.transitions(s)
            for letter, ss in transitions.items():
                assert letter == self._letter(self.tword[ss], self.tstart[ss])
                assert self.parent[ss] == s

            # the leaves should correspond to the -1 states
            if s != 0:
                k = self.tstart[s]
                p = self.tend[s]
                assert p == -1 or p - k > 0, (s, k, p)

                assert (p == -1) == (not transitions)

                if p != -1:
                    # branching
                    ss = self.sl[s]
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

        Combines :meth:`_check_structural` and :meth:`_check_bijection`; only
        valid to call outside of a :meth:`process` call.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree
            sage: T = ConjugateTree()
            sage: T.process([0,1,0,0,1])
            1
            sage: T._check()
        """
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
        children = [self.transitions(s) for s in range(self.nstates)]

        # compute lexicographically sorted leaves
        leaves = []
        queue = [children[0][letter] for letter in sorted(children[0], reverse=reverse)]
        while queue:
            s = queue.pop()
            if self.tend[s] == -1:
                leaves.append(s)
            else:
                queue.extend(children[s][letter] for letter in sorted(children[s], reverse=reverse))

        pos = {}
        for i, s in enumerate(leaves):
            pos[s] = (xscale * self.dep[self.parent[s]] + 1, yscale * i)

        queue = set(range(self.nstates))
        queue.difference_update(leaves)
        while queue:
            treated = []
            for s in queue:
                assert children[s], "got a leaf!"
                if any(ss not in pos for ss in children[s].values()):
                    continue
                x = xscale * self.dep[s]
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
        if self.nwords == 1:
            colors = ["gainsboro"]
        if self.nwords <= 10:
            cmap = mpl.cm.tab10
        elif self.nwords <= 20:
            cmap = mpl.cm.tab20
        else:
            raise NotImplementedError
        if colors is None and cmap is not None:
            colors = [tuple(row[:3]) for row in cmap(range(self.nwords))]
        G = Graphics()
        for s in range(self.nstates):
            if self.tend[s] == -1:
                # leaf
                G += circle(pos[s], state_size, color=colors[self.tword[s]], fill=True, zorder=1)
            else:
                G += circle(pos[s], state_size, color="silver", fill=True, zorder=1)
                G += circle(pos[s], state_size, color="black", fill=False, zorder=2)
            G += text(str(s), pos[s], color="black", zorder=3)
            for ss in children[s].values():
                G += line2d([pos[s], pos[ss]], color="grey", zorder=0)
                mid = ((pos[s][0]+pos[ss][0])/2, (pos[s][1]+pos[ss][1])/2)
                if self.tend[ss] == -1:
                    label = str(self._letter(self.tword[ss], self.tstart[ss]))
                else:
                    label = ''.join(map(str, self._slice(self.tword[ss], self.tstart[ss], self.tend[ss])))
                G += text(label, mid, color="blue")
        G.axes(False)
        return G
