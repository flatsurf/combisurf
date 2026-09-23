r"""
Conjugate trees, reference implementation

This module holds :class:`ConjugateTreeNaive`, the pure Python conjugate
tree. It is the reference against which the fast Cython
:class:`~combisurf.conjugate_tree.ConjugateTree` of
:mod:`combisurf.conjugate_tree` is tested; everything else in the package
uses the fast one. See :mod:`combisurf.conjugate_tree` for what a conjugate
tree is.

The two classes answer the same questions with the same node numbering, so
either can be substituted for the other. This one stores the children of a
node in a dictionary and the per-node data in parallel Python lists, which
makes it easy to read and about twenty times slower.

EXAMPLES:

The main class from this module is :class:`ConjugateTreeNaive` which is initialized
with no argument::

    sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
    sage: T = ConjugateTreeNaive()

To populate a conjugate tree one uses the function :meth:`~ConjugateTreeNaive.process` that
takes as argument a word on non-negative integers (given as a list)::

    sage: T.process([0])
    1
    sage: T.process([0, 1, 0, 0, 1])
    1
    sage: T.process([1, 0, 1, 0])
    2
    sage: T.process([0, 1])
    -2

The output value of :meth:`~ConjugateTreeNaive.process` is either a pair ``(False,
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
"""

from combisurf.word import word_check, word_init


class ConjugateTreeNaive:
    r"""
    Tree structure to store all conjugates of a finite set of primitive words.

    The data structure works with words over non-negative integers.  The nodes
    are encoded with integers from 0 to the number of nodes minus one. The root
    always get the index ``0`` and created nodes gets the first available index
    (nodes are never deleted). In all algorithms, a node index is often denoted
    by a variable ``s``.

    Attributes

    - ``_words`` -- the list of (non-pairwise conjugate) primitive words that
      define this conjugate tree. This list might be updated when
      :meth:`process` is called

    - ``_depth`` -- list of depths of states

    - ``_transitions`` -- list of dictionaries that store children of each
      node. The keys are the first letter of the transition label and values are
      the target nodes.

    - ``_suffix_link`` -- pointer from internal states different from the root
      to their suffix obtained by removing the first letter

    - ``_ancestor`` -- list of ancestors

    - ``_transition_word``, ``_transition_start``, ``_transition_end`` -- lists
      that encode the information on a transition to a node ``s``. The
      associated variables are often denoted ``i``, ``k`` and ``p`` in the
      algorithms.
    """
    def __init__(self, alphabet=0, reserve=0, algorithm=None):
        r"""
        INPUT:

        - ``alphabet``, ``reserve``, ``algorithm`` -- ignored; they are
          accepted so that this class is a drop-in replacement for the
          Cython :class:`~combisurf.conjugate_tree.ConjugateTree`, which
          uses them to pick and to size its transition table.
        """
        # NOTE: -2 is a special code for uninitialized (see _add_node), when adding a node we reallocate accordingly
        self._words = []            # (primitive) words defining the tree
        self._depth = [0]           # internal state -> word length
        self._transitions = [{}]    # state -> (letter -> child)
        self._suffix_link = [-1]    # state -> state
        self._ancestor = [-1]       # state -> state
        self._transition_word = [0] # state -> word index
        self._transition_start = [-4]
        self._transition_end = [-3]

    def words(self):
        r"""
        Return the list of primitive words that define this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()
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
        return [w[:] for w in self._words]

    def num_words(self):
        r"""
        Return the number of primitive words that define this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
            sage: T = ConjugateTreeNaive()
            sage: T.process([0,1,0,1])
            2
            sage: T.num_words()
            1
        """
        return len(self._words)

    def word(self, i):
        r"""
        Return the ``i``-th primitive word of this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
            sage: T = ConjugateTreeNaive()
            sage: T.process([0,1,0,0,1])
            1
            sage: T.word(0)
            array('i', [0, 1, 0, 0, 1])
        """
        if i < 0 or i >= len(self._words):
            raise ValueError(f"i (={i}) must be the index of a word")
        return self._words[i][:]

    def word_length(self, i):
        r"""
        Return the length of the ``i``-th primitive word of this conjugate
        tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
            sage: T = ConjugateTreeNaive()
            sage: T.process([0,1,0,0,1])
            1
            sage: T.word_length(0)
            5
        """
        if i < 0 or i >= len(self._words):
            raise ValueError(f"i (={i}) must be the index of a word")
        return len(self._words[i])

    def letter(self, i, k):
        r"""
        Return the ``k``-th letter of the ``i``-th word, read cyclically.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
            sage: T = ConjugateTreeNaive()
            sage: T.process([0, 4, 2, 3])
            1
            sage: T.letter(0, 19)
            3
            sage: T.letter(0, -1)
            3
        """
        if i < 0 or i >= len(self._words):
            raise ValueError(f"i (={i}) must be the index of a word")
        return self._letter(i, k)

    def alphabet(self):
        r"""
        Return ``0``: this class is not told the alphabet.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
            sage: ConjugateTreeNaive(8).alphabet()
            0
        """
        return 0

    def algorithm(self):
        r"""
        Return ``'dict'``: this class holds the children of a node in a
        dictionary.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
            sage: ConjugateTreeNaive().algorithm()
            'dict'
        """
        return 'dict'

    def transitions(self, s):
        r"""
        Return the children of the node ``s`` as a dictionary mapping the
        first letter of a transition to its target, ordered by letter.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
            sage: T = ConjugateTreeNaive()
            sage: T.process([0,1,1])
            1
            sage: T.transitions(0)
            {0: 1, 1: 3}
            sage: T.transitions(1)
            {}
        """
        if s < 0 or s >= self.num_states():
            raise ValueError(f"s (={s}) must be a node")
        return dict(sorted(self._transitions[s].items()))

    def _pprint(self):
        r"""
        Print the transitions of this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
            sage: T = ConjugateTreeNaive()
            sage: T.process([0,1,1])
            1
            sage: T._pprint()
             0 --0(i=0, k=0)-->  1
             0 --1(array('i', [1]))-->  3
             3 --0(i=0, k=3)-->  4
             3 --1(i=0, k=2)-->  2
        """
        ans = []
        for s, transitions in enumerate(self._transitions):
            for letter in sorted(transitions):
                ss = transitions[letter]
                i = self._transition_word[ss]
                k = self._transition_start[ss]
                p = self._transition_end[ss]
                if p != -1:
                    ans.append(f"{s:2} --{letter}({self._words[i][k:p]})--> {ss:2}")
                else:
                    ans.append(f"{s:2} --{letter}(i={i}, k={k})--> {ss:2}")
        print("\n".join(ans))


    def __repr__(self):
        return "ConjugateTreeNaive with {} states, {} leaves and {} implicit nodes".format(self.num_states(), len(self.leaves()), self.size())

    def _leaf_shift(self, s):
        r"""
        Given a leaf with index ``s`` return the leaf corresponding to its shifted word.

        The function ``leaf_shift`` is a permutation of the leaves of this conjugate
        tree whose orbits represent conjugate words. There is no need for this function
        as each processing of a word provides a cycle of the created leaves (by
        increasing order).

        This function is kept only for testing purposes.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()

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
        if s < 0 or s >= len(self._ancestor):
            raise ValueError
        if self._transition_end[s] != -1:
            raise ValueError(f"s(={s}) not a leaf")
        assert self._suffix_link[s] == -2
        i = self._transition_word[s]
        k = self._transition_start[s]
        letter = self._letter(i, k)

        ss = self._ancestor[s]
        assert letter in self._transitions[ss] and self._transitions[ss][letter] == s

        ss = self._suffix_link[ss]
        if ss == -1:
            ss = 0
        else:
            assert letter in self._transitions[ss]
            ss = self._transitions[ss][letter]
        while self._transition_end[ss] != -1:
            k += self._transition_end[ss] - self._transition_start[ss]
            letter = self._letter(i, k)
            assert letter in self._transitions[ss]
            ss = self._transitions[ss][letter]
        return ss

    def leaf_as_conjugate(self, s):
        r"""
        Return the pair ``(i, k)`` such that the leaf ``s`` corresponds to the
        ``i``-th word shifted by ``k``.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()
            sage: T.process([0,0,0,1,0,2])
            1
            sage: T.process([1,2,1,2,1,1,2])
            1
            sage: T.process([0])
            1
            sage: [T.leaf_as_conjugate(s) for s in T.leaves()] == [(i, k) for i, w in enumerate(T.words()) for k in range(len(w))]
            True
        """
        if s < 0 or s >= len(self._ancestor):
            raise ValueError
        i = self._transition_word[s]
        k = self._transition_start[s]
        ss = self._ancestor[s]
        ans = (k - self._depth[ss]) % len(self._words[i])
        return (i, ans)

    def size(self):
        r"""
        Return the size of this conjugate tree.

        The size is the number of implicit states where each leaf accounts for
        1.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()
            sage: T.process([0,0,0,0,1])
            1
            sage: T.size()
            9
            sage: T.process([0,0,0,1])
            1
            sage: T.size()
            29

        The total size is the same if we input the two words in the opposite order::

            sage: T = ConjugateTreeNaive()
            sage: T.process([0,0,0,1])
            1
            sage: T.size()
            7
            sage: T.process([0,0,0,0,1])
            1
            sage: T.size()
            29
        """
        ans = 0
        for s in range(self.num_states()):
            i = self._transition_word[s]
            k = self._transition_start[s]
            p = self._transition_end[s]
            if p == -1:
                ans += 1
            else:
                ans += p - k
        return ans

    def num_states(self):
        r"""
        Return the number of states.
        """
        return len(self._transitions)

    def internal_states(self):
        r"""
        Return the internal states in this conjugate tree.

        Note that any further call to :meth:`process` might change the
        structure of the tree but not the word encoded by a given state. In
        particular, internal states remain internal states.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()
            sage: T.process([0,0,0,1])
            1
            sage: T.size()
            7

        """
        return [s for s in range(1, self.num_states()) if self._transitions[s]]

    def leaves(self):
        r"""
        Return the leaves in this conjugate tree.

        Note that any further call to :meth:`process` might change the
        structure of the tree but not the word encoded by a given state. In
        particular, leaves remain leaves after an update.

        TESTS::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()
            sage: T.process([0,1,1])
            1
            sage: T.process([0,1])
            1
            sage: len(T.leaves())
            5
        """
        return [s for s in range(1, self.num_states()) if not self._transitions[s]]

    def cyclically_sorted_leaves(self, order, pivot):
        r"""
        Return the leaves sorted by the order of the letters ``order``, the
        order below a node being turned by ``pivot``.

        See
        :meth:`~combisurf.conjugate_tree.ConjugateTree.cyclically_sorted_leaves`:
        the children of the root are visited by increasing ``order[c]`` and
        the children of an internal node whose label ends with the letter
        ``b`` by increasing ``(order[c] - pivot[b]) % n``.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()
            sage: T.process([0,1,1])
            1
            sage: T.process([0,1])
            1
            sage: T.cyclically_sorted_leaves([0, 1], [1, 0])
            [6, 1, 8, 4, 2]
        """
        n = len(order)
        leaves = []
        queue = [self._transitions[0][letter] for letter in sorted(self._transitions[0], key=lambda letter: order[letter], reverse=True)]
        while queue:
            s = queue.pop()
            if self._transition_end[s] == -1:
                leaves.append(s)
            else:
                i = self._transition_word[s]
                p = self._transition_end[s]
                base = pivot[self._letter(i, p - 1)]
                transitions = sorted(self._transitions[s], key=lambda letter: (order[letter] - base) % n, reverse=True)
                queue.extend(self._transitions[s][letter] for letter in transitions)

        return leaves

    def graph(self):
        r"""
        Return this conjugate tree as a directed graph, each edge labelled by
        the letters it reads, or by its first letter only if it ends at a
        leaf (such an edge reads an infinite periodic word).

        The letters are separated by commas once the tree has a letter
        larger than `9`, so that a label can be read back.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
            sage: T = ConjugateTreeNaive()
            sage: T.process([0,1,1])
            1
            sage: T.graph()
            Digraph on 5 vertices
            sage: sorted(T.graph().edges())
            [(0, 1, '0'), (0, 3, '1'), (3, 2, '1'), (3, 4, '0')]

        With a letter larger than `9`::

            sage: T = ConjugateTreeNaive()
            sage: T.process([12,1,2])
            1
            sage: T.process([12,1,1])
            1
            sage: sorted(T.graph().edges())
            [(0, 3, '2'),
             (0, 4, '12,1'),
             (0, 6, '1'),
             (4, 1, '2'),
             (4, 5, '1'),
             (6, 2, '2'),
             (6, 7, '1'),
             (6, 8, '12')]
        """
        from sage.graphs.digraph import DiGraph
        G = DiGraph(self.num_states(), loops=False, multiedges=False)
        for s in range(self.num_states()):
            for t in self._transitions[s].values():
                G.add_edge(s, t, self._edge_label(t))
        return G

    def _edge_label(self, t):
        r"""
        Return the label of the edge ending at the node ``t``: the letters it
        reads, or only its first letter if ``t`` is a leaf (such an edge
        reads an infinite periodic word). The letters are separated by
        commas once the tree has a letter larger than `9`.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive
            sage: T = ConjugateTreeNaive()
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
        if t <= 0 or t >= self.num_states():
            raise ValueError(f"t (={t}) must be a node")
        if self._transition_end[t] == -1:
            return str(self._letter(self._transition_word[t], self._transition_start[t]))
        sep = ',' if max((max(w) for w in self._words), default=-1) >= 10 else ''
        return sep.join(map(str, self._slice(self._transition_word[t], self._transition_start[t], self._transition_end[t])))

    def _add_node(self):
        r"""
        Internal low-level function that add a nodes and return its index.

        The function performs the necessary reallocation so that one can
        access ``self._transitions[i]``, etc where ``i`` is the index
        of the created node.
        """
        n = len(self._transitions)
        self._depth.append(-2)
        self._transitions.append({})
        self._suffix_link.append(-2)
        self._ancestor.append(-2)
        self._transition_start.append(-2)
        self._transition_end.append(-2)
        self._transition_word.append(-2)
        return n

    def _check_structural(self):
        r"""
        Check the per-node invariants of this conjugate tree.

        Unlike :meth:`_check_bijection`, these invariants hold at every
        intermediate step of :meth:`process`, not only once it returns, since
        they say nothing about the leaves of the word currently being added.
        """
        for w in self._words:
            assert word_check(w)
        n = len(self._transitions)
        assert len(self._suffix_link) == n
        assert len(self._ancestor) == n
        assert len(self._transition_word) == n
        assert len(self._transition_start) == n
        assert len(self._transition_end) == n

        assert self._ancestor[0] == -1
        assert self._transition_start[0] == -4, self._transition_start[0]
        assert self._transition_end[0] == -3, self._transition_end[0]

        for s in range(1, n):
            assert self._ancestor[s] >= 0, (s, self._ancestor)
            assert self._transition_start[s] >= 0, (s, self._transition_start)
            assert self._transition_end[s] > -2, (s, self._transition_end)
            assert self._transition_word[s] != -2, (s, self._transition_word)

            if s != 0 and self._transition_end[s] != -1:
                # suffix link are only for internal nodes different from the root
                assert self._suffix_link[s] != -2, (s, self._suffix_link)

        for s in range(n):
            for letter, ss in self._transitions[s].items():
                i = self._transition_word[ss]
                k = self._transition_start[ss]
                p = self._transition_end[ss]
                assert letter == self._letter(i, k)
                assert self._ancestor[ss] == s

            # the leaves should correspond to the -1 states
            if s != 0:
                k = self._transition_start[s]
                p = self._transition_end[s]
                assert p == -1 or p - k > 0, (s, k, p)

                assert (self._transition_end[s] == -1) == (not self._transitions[s])

                if self._transition_end[s] == -1:
                    # leaf
                    assert not self._transitions[s]
                else:
                    # branching
                    assert self._transitions[s]
                    ss = self._suffix_link[s]
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
        """
        assert [self.leaf_as_conjugate(s) for s in self.leaves()] == [(i, k) for i, w in enumerate(self.words()) for k in range(len(w))]

    def _check(self):
        r"""
        Check all invariants of this conjugate tree.

        Combines :meth:`_check_structural` and :meth:`_check_bijection`; only
        valid to call outside of a :meth:`process` call.
        """
        self._check_structural()
        self._check_bijection()

    def _test_and_split(self, s, i, k, p, letter):
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
        if k < p:
            # implicit state
            # get the transition from s starting with word[i][k] and test
            # whether its (p - k)-th letter coincide with letter or not
            t = self._transitions[s][self._letter(i, k)]
            ii = self._transition_word[t]
            # w = self._words[ii]
            kk = self._transition_start[t]
            assert kk >= 0
            index = kk + p - k
            assert index >= 0
            lletter = self._letter(ii, index)
            if letter == lletter:
                # the node already exists
                return -1
            else:
                # make the node explicit
                # the new node ss is the node made explicit
                # s ---> t becomes s --> ss --> t
                ss = self._add_node()

                self._transition_word[ss] = ii
                self._transition_start[ss] = kk
                self._transition_end[ss] = index
                self._transition_start[t] = index

                self._ancestor[t] = ss
                self._ancestor[ss] = s

                self._transitions[s][self._letter(ii, kk)] = ss
                self._transitions[ss][lletter] = t

                self._depth[ss] = self._depth[s] + index - kk

                return ss
        else:
            # explicit state
            if s == -1 or letter in self._transitions[s]:
                # the node already exists
                return -1
            else:
                return s

    def canonize(self, s, i, k, p):
        r"""
        Canonize the quadruple ``(s, i, k, p)`` representing
        the (explicit or implicit) state obtained after reading word[i][k:p]
        from s.

        Return a pair ``(s, k)`` (as ``i`` and ``p`` do not change).
        """
        assert s >= -1, s
        if k >= p:
            # already explicit
            return (s, p)
        else:
            ss = 0 if s == -1 else self._transitions[s][self._letter(i, k)]
            kk = self._transition_start[ss]
            pp = self._transition_end[ss]
            while pp != -1 and pp - kk < p - k:
                k += pp - kk
                s = ss
                ss = self._transitions[s][self._letter(i, k)]
                kk = self._transition_start[ss]
                pp = self._transition_end[ss]
            if pp != -1 and pp - kk == p - k:
                # explicit
                return (ss, p)
            else:
                # implicit
                return (s, k)

    def _update(self, s, i, k, p):
        r"""
        Low-level internal function that updates by reading one letter.

        INPUT:

        - ``s`` -- node
        - ``i`` -- index of a word
        - ``k``, ``p`` -- beginning and end of a slice in the ``i``-th word

        Here ``(s, i, k, p)`` should be the canonical reference pair of the
        active state from the previous state.
        """
        # (s, k, p): active state which is the first state along the boundary
        # path which is not an active leaf
        # r: closest branching from s (r is either s or its ancestor)
        letter = self._letter(i, p)
        old_r = 0
        created_leaves = []
        r = self._test_and_split(s, i, k, p, letter)
        while r != -1:
            assert k >= 0 and p >= 0, (k, p)
            assert r >= 0 and old_r >= 0, (r, old_r)
            # create a leaf
            rr = self._add_node()
            created_leaves.append(rr)
            self._transitions[r][letter] = rr
            self._ancestor[rr] = r
            self._transition_word[rr] = i
            self._transition_start[rr] = p
            self._transition_end[rr] = -1
            if old_r != 0:
                assert r != old_r
                self._suffix_link[old_r] = r
            old_r = r
            s, k = self.canonize(self._suffix_link[s], i, k, p)
            r = self._test_and_split(s, i, k, p, letter)

        if old_r != 0:
            assert old_r != s, (old_r, s)
            self._suffix_link[old_r] = s

        return s, k, created_leaves

    def _slice(self, i, k, p):
        r"""
        Return the slice from ``k`` to ``p`` of the ``i``-th word.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()
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
        w = self._words[i]
        return [w[j % len(w)] for j in range(k, p)]

    def _letter(self, i, k):
        r"""
        Return the ``k``-th letter of the ``i``-th word.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()
            sage: T.process([0, 4, 2, 3])
            1
            sage: T._letter(0, 1)
            4
            sage: T._letter(0, 19)
            3
        """
        return self._words[i][k % len(self._words[i])]

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
        """
        if not w:
            raise ValueError("empty word in input")
        if check:
            w = word_init(w)
        for letter in w:
            if letter < 0:
                raise ValueError("invalid word: must be made of non-negative integers")
        i = len(self._words)
        l = len(w)
        self._words.append(w)

        s = 0
        k = 0

        # To ensure that we find all conjugates we must create as many leaves
        # as the size rof w (assuming it is primitive)
        num_leaves = 0
        p = 0
        while True:
            if p != k:
                ss = self._transitions[s][self._letter(i, k)]
                ii = self._transition_word[ss]
                kk = self._transition_start[ss]
                pp = self._transition_end[ss]
            else:
                ii = -1
            s, k, created_leaves = self._update(s, i, k, p)
            num_leaves += len(created_leaves)
            if hard_check:
                self._check_structural()
            s, k = self.canonize(s, i, k, p + 1)

            # halt condition
            if num_leaves == l:
                # w is primitive
                break
            elif ii == i and p >= 2 * l:
                # w is non primitive
                break
            elif p >= l and num_leaves == 0 and ii != -1 and pp == -1 and l % len(self._words[ii]) == 0:
                # w is conjugate to a power of self._words[ii]
                break

            p += 1

        if num_leaves == 0:
            self._words.pop()
            if hard_check:
                self._check_bijection()
            return -ii
        else:
            assert len(w) % num_leaves == 0, (len(w), num_leaves)
            exponent = len(w) // num_leaves
            if exponent != 1:
                # NOTE: only store primitive words
                del self._words[-1][l//exponent:]
            if hard_check:
                self._check_bijection()
            return exponent

    def internal_state_word(self, s):
        r"""
        Return the word corresponding to the explicit state ``s``

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()
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
        if s < 0:
            raise ValueError("s must be a node")
        path = [s]
        while path[-1] != 0:
            path.append(self._ancestor[path[-1]])
        path.pop()
        ans = []
        for s in reversed(path):
            i = self._transition_word[s]
            k = self._transition_start[s]
            p = self._transition_end[s]
            if p == -1:
                p = len(self._words[i])
            ans.extend(self._slice(i, k, p))
        return ans

    def plot(self, state_size=.25, xscale=1, yscale=1, reverse=False):
        r"""
        Return a plot of this conjugate tree.

        EXAMPLES::

            sage: from combisurf.conjugate_tree_naive import ConjugateTreeNaive

            sage: T = ConjugateTreeNaive()
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
        # compute lexicographically sorted leaves
        leaves = []
        queue = [self._transitions[0][letter] for letter in sorted(self._transitions[0], reverse=reverse)]
        while queue:
            s = queue.pop()
            if self._transition_end[s] == -1:
                leaves.append(s)
            else:
                queue.extend(self._transitions[s][letter] for letter in sorted(self._transitions[s], reverse=reverse))

        pos = {}
        for i, s in enumerate(leaves):
            pos[s] = (xscale * self._depth[self._ancestor[s]] + 1, yscale * i)

        topological_order = []
        queue = set(range(len(self._transitions)))
        queue.difference_update(leaves)
        while queue:
            treated = []
            for s in queue:
                assert self._transitions[s], "got a leaf!"
                if any(ss not in pos for ss in self._transitions[s].values()):
                    continue
                x = xscale * self._depth[s]
                y = sum(pos[ss][1] for ss in self._transitions[s].values()) / len(self._transitions[s])
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
        if len(self._words) == 1:
            colors = ["gainsboro"]
        if len(self._words) <= 10:
            cmap = mpl.cm.tab10
        elif len(self._words) <= 20:
            cmap = mpl.cm.tab20
        else:
            raise NotImplementedError
        if colors is None and cmap is not None:
            colors = [tuple(row[:3]) for row in cmap(range(len(self._words)))]
        G = Graphics()
        for s in range(len(self._transitions)):
            if self._transition_end[s] == -1:
                # leaf
                G += circle(pos[s], state_size, color=colors[self._transition_word[s]], fill=True, zorder=1)
            else:
                G += circle(pos[s], state_size, color="silver", fill=True, zorder=1)
                G += circle(pos[s], state_size, color="black",fill=False, zorder=2)
            G += text(str(s), pos[s], color="black", zorder=3)
            for ss in self._transitions[s].values():
                G += line2d([pos[s], pos[ss]], color="grey", zorder=0)
                mid = ((pos[s][0]+pos[ss][0])/2, (pos[s][1]+pos[ss][1])/2)
                G += text(self._edge_label(ss), mid, color="blue")
        G.axes(False)
        return G
