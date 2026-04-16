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

The main class from this module is :class:`ConjugateTree` which is initialized
with no argument::

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
"""

# TODO: leaves should never be modified
VERBOSE = False
INDENT = 0
class ConjugateTree:
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
    def __init__(self):
        # NOTE: -2 is a special code for uninitialized (see _add_node), when adding a node we reallocate accordingly
        self._words = []            # (primitive) words defining the tree
        self._depth = [0]           # internal state -> word length
        self._transitions = [{}]    # state -> (letter -> child)
        self._suffix_link = [-1]    # state -> state
        self._ancestor = [-1]       # state -> state
        self._transition_word = [0] # state -> word index
        self._transition_start = [-4]
        self._transition_end = [-3]
        self._max_read = []

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
            [[0, 1], [0, 1, 0, 0, 1], [0, 1, 0, 1, 1]]
        """
        return [w[:] for w in self._words]

    def pprint(self):
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
        return "SuffixTree with {} states, {} leaves and {} implicit nodes".format(self.num_states(), len(self.leaves()), self.size())

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

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
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
        particular, leaves remain leaves.

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
        return [s for s in range(1, self.num_states()) if not self._transitions[s]]

    def graph(self):
        G = DiGraph(self.num_states(), loops=False, multiedges=False)
        for s in range(self.num_states()):
            for t in self._transitions[s].values():
                # w = self._words[self._transition_word[t]]
                i = self._transition_word[t]
                k = self._transition_start[t]
                p = self._transition_end[t]
                G.add_edge(s, t, f"({i},{k},{p})")
        return G

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

    def _check(self):
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
                    w0 = self.word(ss)
                    w1 = self.word(s)
                    assert w0 == w1[1:], (ss, w0, s, w1)

        assert [self.leaf_as_conjugate(s) for s in self.leaves()] == [(i, k) for i, w in enumerate(self.words()) for k in range(len(w))]

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
            if VERBOSE:
                print(" " * INDENT + f"[_test_and_split] _test_and_split(s={s}, i={i}, k={k}, p={p}, letter={letter})")
                print(" " * INDENT + f"[_test_and_split] implicit with target t={t}")
            ii = self._transition_word[t]
            # w = self._words[ii]
            kk = self._transition_start[t]
            assert kk >= 0
            index = kk + p - k
            assert index >= 0
            # TODO: remove check
            assert self._transition_end[t] == -1 or index < self._transition_end[t], (s, i, k, p, t, index)
            lletter = self._letter(ii, index)
            if VERBOSE:
                print(" " * INDENT + f"[_test_and_split] letter={letter} lletter={lletter}")
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

                if VERBOSE:
                    print(" " * INDENT + f"[_test_and_split] new node {ss} ({self.word(ss)}) from split between {s} ({self.word(s)}) and {t} ({self.word(t)}")

                return ss
        else:
            # explicit state
            if VERBOSE:
                print(" " * INDENT + f"[_test_and_split] _test_and_split(s={s}, i={i}, k={k}, p={p}, letter={letter}): explicit")
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
        if VERBOSE:
            print(" " * INDENT + f"[canonize] canonize(s={s}, i={i}, k={k}, p={p})")
        if k >= p:
            # already explicit
            # return (s, p)
            out = (s, p)
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
                # return (ss, p)
                if VERBOSE:
                    print(" " * INDENT + f"[canonize] 1: s={s} k={k} p={p} ss={ss} kk={kk} pp={pp}")
                out = (ss, p)
            else:
                # implicit
                # return (s, k)
                if VERBOSE:
                    print(" " * INDENT + f"[canonize] 2: s={s} k={k} p={p} ss={ss} kk={kk} pp={pp}")
                out = (s, k)

        if VERBOSE:
            print(" " * INDENT + f"[canonize] out={out}")

        # TODO: remove check
        s, k = out
        if k == p:
            return out
        elif k > p:
            raise RuntimeError
        else:
            assert s != -1
            t = self._transitions[s][self._letter(i, k)]
            kk = self._transition_start[t]
            pp = self._transition_end[t]
            assert pp == -1 or p - k < pp - kk, (k, p, kk, pp)
            return out

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
        global VERBOSE, INDENT

        letter = self._letter(i, p)
        old_r = 0
        created_leaves = []
        if VERBOSE:
            print(" " * INDENT + f"[update] process letter w[{p}]={letter}")
            print(" " * INDENT + f"[update] going through boundary path from active state s={s} k={k} p={p}")
            INDENT += 2
        r = self._test_and_split(s, i, k, p, letter)
        if VERBOSE:
            INDENT -= 2
        while r != -1:
            assert k >= 0 and p >= 0, (k, p)
            assert r >= 0 and old_r >= 0, (r, old_r)
            # create a leaf
            rr = self._add_node()
            created_leaves.append(rr)
            if VERBOSE:
                print(" " * INDENT + f"[update]   r={r} s={s} k={k}")
                print(" " * INDENT + f"[update]   new leaf rr={rr})")
            self._transitions[r][letter] = rr
            self._ancestor[rr] = r
            self._transition_word[rr] = i
            self._transition_start[rr] = p
            self._transition_end[rr] = -1
            if old_r != 0:
                if VERBOSE:
                    print(" " * INDENT + f"[update]   create suffix link {old_r} -> {r}")
                assert r != old_r
                self._suffix_link[old_r] = r
            old_r = r
            if VERBOSE:
                INDENT += 2
            s, k = self.canonize(self._suffix_link[s], i, k, p)
            r = self._test_and_split(s, i, k, p, letter)
            if VERBOSE:
                INDENT -= 2
        if VERBOSE:
            print(" " * INDENT + f"[update] end of loop: s={s} k={k} p={p} old_r={old_r}")

        if old_r != 0:
            if VERBOSE:
                print(" " * INDENT + f"[update] create suffix link {old_r} -> {s}")
            assert old_r != s, (old_r, s)
            self._suffix_link[old_r] = s

        return s, k, created_leaves

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
        w = self._words[i]
        return [w[j % len(w)] for j in range(k, p)]

    def _letter(self, i, k):
        r"""
        Return the ``k``-th letter of the ``i``-th word.

        EXAMPLES::

            sage: from combisurf.conjugate_tree import ConjugateTree

            sage: T = ConjugateTree()
            sage: T.process([0, 4, 2, 3])
            1
            sage: T._letter(0, 1)
            4
            sage: T._letter(0, 19)
            3
        """
        self._max_read[i] = max(self._max_read[i], k)
        return self._words[i][k % len(self._words[i])]

    def process(self, w, check=False):
        r"""
        Add the word ``w`` in this conjugate tree.

        The output value is an integer. Depending on its sign it encodes
        different informations.

        - a positive ``exponent`` if the word ``w`` is not already present
          and ``exponent`` is the exponent of ``w`` (which is ``1`` if and only
          if ``w`` is primitive). In that case, the number of leaves in the tree
          increases by the period of ``w`` (which is its length divided by the
          exponent).

        - a non-negative ``-index`` if the word ``w``is already present and
          ``index`` is the index of the leaf corresponding to ``w`` in this
          conjugate tree
        """
        w = list(map(int, w))
        for letter in w:
            if letter < 0:
                raise ValueError("invalid word: must be made of non-negative integers")
        global VERBOSE, INDENT
        i = len(self._words)
        l = len(w)
        self._words.append(w)
        self._max_read.append(-1)

        s = 0
        k = 0

        # To ensure that we find all conjugates we must create as many leaves
        # as the size rof w (assuming it is primitive)
        num_leaves = 0
        p = 0
        while True:
            if VERBOSE:
                print(f"[process] new loop with active state s={s} i={i} k={k} p={p}")
                INDENT += 2
            if p != k:
                ss = self._transitions[s][self._letter(i, k)]
                ii = self._transition_word[ss]
                kk = self._transition_start[ss]
                pp = self._transition_end[ss]
                if VERBOSE:
                    print(f"[process] endpoint of active state ss={ss} ii={ii} kk={kk} pp={pp}")
            else:
                ii = -1
            s, k, created_leaves = self._update(s, i, k, p)
            num_leaves += len(created_leaves)
            if VERBOSE:
                print(f"[process] {len(created_leaves)} new leaves, total={num_leaves}")
            if check:
                self._check()
            s, k = self.canonize(s, i, k, p + 1)
            if VERBOSE:
                INDENT -= 2

            # halt condition
            if num_leaves == l:
                # w is primitive
                if VERBOSE:
                    print(f"[process] enough leaves to determine conjugates at p={p} (i={i} ii={ii} len(w)={len(w)})")
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
            self._max_read.pop()
            return -ii
        else:
            assert len(w) % num_leaves == 0, (len(w), num_leaves)
            exponent = len(w) // num_leaves
            if exponent != 1:
                # NOTE: only store primitive words
                del self._words[-1][l//exponent:]
            return exponent

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
                if self._transition_end[ss] == -1:
                    label = str(self._letter(self._transition_word[ss], self._transition_start[ss]))
                else:
                    label = ''.join(map(str, self._slice(self._transition_word[ss], self._transition_start[ss], self._transition_end[ss])))
                G += text(label, mid, color="blue")
        G.axes(False)
        return G
