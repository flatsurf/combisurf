r"""
Word on non-negative integers and free group elements

The alphabet is always the non-negative integers. For words seen
as free group element, the letter `i ^ 1` is the inverse of `i`
(that is `2i` and `2i+1` are inverses of each other).
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

from combisurf.misc cimport str_to_int


def word_check(w):
    r"""
    Check that ``w`` is a valid word.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.word import word_check

        sage: word_check(array('i', []))
        True
        sage: word_check(array('d', []))
        False
        sage: word_check(array('i', [0]))
        True
        sage: word_check(array('i', [-1]))
        False
        sage: word_check("")
        False
    """
    if not isinstance(w, array.array):
        return False
    if w.typecode != 'i':
        return False
    if any(x < 0 for x in w):
        return False
    return True


def word_init(data=None):
    r"""
    Initialize a word from ``data``.

    EXAMPLES::

        sage: from combisurf.word import word_init
        sage: word_init([0, 1, 3, 2, 1])
        array('i', [0, 1, 3, 2, 1])
        sage: word_init("0,1,~2,~0")
        array('i', [0, 2, 5, 1])

    With no argument, returns the empty word::

        sage: word_init()
        array('i')
    """
    if data is None:
        return array.array('i')

    if isinstance(data, str):
        data = data.replace(' ', '')
        if data.startswith('(') or data.startswith('['):
            data = data[1:]
        if data.endswith(')') or data.endswith(')'):
            data = data[:len(data)-1]
        data = [str_to_int(x) for x in data.split(',')]
        data = [2 * x if x >= 0 else (2 * ~x + 1) for x in data]

    if isinstance(data, (array.array, tuple, list)):
        return array.array('i', data)
    else:
        raise TypeError("invalid argument")


def word_string(array.array w, edge_like=False, separator=', ', opening='[', closing=']'):
    r"""
    Return a string representing ``w``.

    INPUT:

    - ``edge_like`` -- (boolean, default ``False``) whether to print as
      half-edges and inverses

    - ``separator`` -- (str, default ``', '``) the string used to separate the elements
      in ``w``

    - ``opening`` -- (str, default ``'['``) the string used at the start

    - ``closing`` -- (str, default ``']'``) the string used at the end

    EXAMPLES::

        sage: from combisurf.word import word_init, word_string
        sage: w = word_init([0, 3, 1, 2])
        sage: word_string(w)
        '[0, 3, 1, 2]'
        sage: word_string(w, edge_like=True, separator=':', opening='', closing='')
        '0:~1:~0:1'
    """
    if edge_like:
        elt = lambda e: ('~%d' % (e // 2)) if e % 2 else '%d' % (e // 2)
    else:
        elt = str

    return opening + separator.join(map(elt, w)) + closing


def word_border_table(array.array u):
    r"""
    Return the border table of the word ``u``.

    A *border* of a word is a word which is both a proper prefix and a suffix
    of it. The entry ``b[i]`` of the border table is the length of the longest
    border of the prefix ``u[:i]``; the table has length ``len(u) + 1``.

    INPUT:

    - ``u`` -- a word

    EXAMPLES::

        sage: from combisurf.word import word_init, word_border_table

        sage: word_border_table(word_init([0, 0, 1, 0, 0, 1, 0]))
        array('i', [0, 0, 1, 0, 1, 2, 3, 4])

    The prefix ``[0, 0, 1, 0, 0]`` above has ``[0, 0]`` as longest border,
    whence the entry ``2`` at index ``5``. A word without repetition has a
    trivial border table::

        sage: word_border_table(word_init([0, 1, 2, 3]))
        array('i', [0, 0, 0, 0, 0])
        sage: word_border_table(word_init([0, 0, 0, 0]))
        array('i', [0, 0, 1, 2, 3])
        sage: word_border_table(word_init())
        array('i', [0])

    The smallest period of ``u`` is ``len(u) - b[len(u)]``::

        sage: u = word_init([0, 1, 2, 0, 1, 2, 0, 1])
        sage: b = word_border_table(u)
        sage: len(u) - b[len(u)]
        3

    .. SEEALSO::

        :func:`word_failure_table`
    """
    cdef Py_ssize_t m = len(u)
    cdef array.array b = array.clone(u, m + 1, False)
    cdef int * bb = b.data.as_ints

    bb[0] = 0
    if m == 0:
        return b
    bb[1] = 0

    cdef int * uu = u.data.as_ints
    cdef Py_ssize_t i
    cdef int k = 0
    for i in range(1, m):
        while k > 0 and uu[i] != uu[k]:
            k = bb[k]
        if uu[i] == uu[k]:
            k += 1
        bb[i + 1] = k

    return b


def word_failure_table(array.array u):
    r"""
    Return the failure table of the word ``u``, as consumed by
    :func:`word_find`.

    This is the strong variant of the border table of ``u``, in which a shift
    that would repeat a comparison already known to fail is skipped. It has
    length ``len(u)`` and is related to the border table ``b`` of
    :func:`word_border_table` by ``t[0] = -1`` and, for ``i >= 1``,

    - ``t[i] = b[i]`` if ``u[i] != u[b[i]]``,
    - ``t[i] = t[b[i]]`` otherwise.

    INPUT:

    - ``u`` -- a word

    EXAMPLES::

        sage: from combisurf.word import word_init, word_border_table, word_failure_table

        sage: u = word_init([0, 0, 1, 0, 0, 1, 0])
        sage: word_failure_table(u)
        array('i', [-1, -1, 1, -1, -1, 1, -1])

    It is indeed obtained from the border table by the rule above::

        sage: b = word_border_table(u)
        sage: t = [-1]
        sage: for i in range(1, len(u)):
        ....:     t.append(b[i] if u[i] != u[b[i]] else t[b[i]])
        sage: word_init(t) == word_failure_table(u)
        True

    ::

        sage: word_failure_table(word_init([0, 1, 2, 3]))
        array('i', [-1, 0, 0, 0])
        sage: word_failure_table(word_init([0, 0, 0, 0]))
        array('i', [-1, -1, -1, -1])
        sage: word_failure_table(word_init())
        array('i')

    .. SEEALSO::

        :func:`word_border_table`
    """
    cdef Py_ssize_t m = len(u)
    cdef array.array t = array.clone(u, m, False)
    if m == 0:
        return t

    cdef int * tt = t.data.as_ints
    cdef int * uu = u.data.as_ints

    tt[0] = -1
    cdef Py_ssize_t i
    cdef int cnd = 0
    for i in range(1, m):
        if uu[i] == uu[cnd]:
            tt[i] = tt[cnd]
        else:
            tt[i] = cnd
            while cnd >= 0 and uu[i] != uu[cnd]:
                cnd = tt[cnd]
        cnd += 1

    return t


def word_find(array.array u, array.array v, Py_ssize_t start=0, failure_table=None):
    r"""
    Return the lowest index in ``v`` at or after ``start`` where the word ``u``
    occurs, or ``-1`` if there is no such index.

    The conventions are the ones of :meth:`str.find`: the empty word occurs at
    ``start`` as soon as ``start <= len(v)`` and a negative ``start`` counts
    from the end of ``v``.

    INPUT:

    - ``u`` -- a word, the one searched for

    - ``v`` -- a word, the one searched in

    - ``start`` -- (integer, default ``0``) the index in ``v`` at which the
      search starts

    - ``failure_table`` -- (default: ``None``) the failure table of ``u``, as
      returned by :func:`word_failure_table`; when ``None`` it is computed.
      Pass it explicitly to search a single ``u`` in many ``v`` without
      recomputing it every time.

    EXAMPLES::

        sage: from combisurf.word import word_init, word_failure_table, word_find

        sage: word_find(word_init([1, 2]), word_init([0, 1, 2, 3]))
        1
        sage: word_find(word_init([2, 1]), word_init([0, 1, 2, 3]))
        -1

    The search can be restarted further to the right, which is how one
    enumerates all the occurrences::

        sage: u = word_init([0, 1])
        sage: v = word_init([0, 1, 0, 1, 0, 1])
        sage: word_find(u, v)
        0
        sage: word_find(u, v, 1)
        2
        sage: word_find(u, v, 3)
        4
        sage: word_find(u, v, 5)
        -1

    A negative ``start`` counts from the end of ``v``::

        sage: word_find(u, v, -3)
        4
        sage: word_find(u, v, -100)
        0

    As for :meth:`str.find`, the empty word occurs at ``start``::

        sage: word_find(word_init(), word_init([0, 1]), 1)
        1
        sage: word_find(word_init(), word_init([0, 1]), 2)
        2
        sage: word_find(word_init(), word_init([0, 1]), 3)
        -1

    When the same ``u`` is searched in many words, compute its failure table
    once and pass it along::

        sage: u = word_init([0, 1, 0])
        sage: t = word_failure_table(u)
        sage: [word_find(u, word_init(v), 0, t) for v in ([1, 0, 1, 0], [0, 1, 0], [1, 1])]
        [1, 0, -1]

    The usual variants are built on top of :func:`word_find` this way. All the
    occurrences, including the overlapping ones::

        sage: def occurrences(u, v):
        ....:     t = word_failure_table(u)
        ....:     res = []
        ....:     i = word_find(u, v, 0, t)
        ....:     while i != -1:
        ....:         res.append(i)
        ....:         i = word_find(u, v, i + 1, t)
        ....:     return res
        sage: occurrences(word_init([0, 1, 0]), word_init([0, 1, 0, 1, 0, 1, 0]))
        [0, 2, 4]

    and the number of them::

        sage: len(occurrences(word_init([0, 1, 0]), word_init([0, 1, 0, 1, 0, 1, 0])))
        3

    ALGORITHM:

    The Knuth-Morris-Pratt algorithm. Building the failure table of ``u`` takes
    time ``O(len(u))`` and the search then takes time ``O(len(v) - start)``, so
    that the whole computation is linear in ``len(u) + len(v)``. When the
    failure table is provided, the search alone is performed.

    That is the cost of a single search. Enumerating all the occurrences as
    above costs ``O(len(u))`` per occurrence on top of the scan, since each call
    restarts the automaton, and is not linear overall: on ``u = [0] * k`` inside
    ``v = [0] * n`` it takes time ``O(n k)``.

    .. TODO::

        Make the enumeration of the occurrences linear. What a call cannot
        recover from its arguments is the position ``k`` reached in ``u``, which
        a restart at ``i + 1`` sets back to ``0``; no choice of ``start`` gets
        around it, as resuming one period after a match is still ``i + 1`` on
        ``u = [0] * k``. The place for it is a ``word_occurrences`` that keeps
        ``k`` across the matches, resuming at the longest border of ``u``
        instead of at ``0``, rather than a further argument of
        :func:`word_find`, whose signature follows :meth:`str.find`.

        Note that this border is not in the failure table: that table has length
        ``len(u)``, because :func:`word_find` returns as soon as ``k`` reaches
        ``len(u)`` and never reads an entry there. It is
        ``word_border_table(u)[len(u)]``. Using the plain border for the resume
        and the failure table for the mismatches is correct, the strong variant
        only skipping shifts that would repeat a comparison known to fail, and
        there is no such comparison at a match.

    TESTS:

    The failure table, when provided, must match ``u``::

        sage: from combisurf.word import word_init, word_find

        sage: word_find(word_init([0, 1]), word_init([0, 1]), 0, word_init([-1]))
        Traceback (most recent call last):
        ...
        ValueError: failure_table must have the same length as u
        sage: word_find(word_init([0, 1]), word_init([0, 1]), 0, [-1, 0])
        Traceback (most recent call last):
        ...
        TypeError: Cannot convert list to array.array

    The agreement with :meth:`str.find` is checked exhaustively in
    ``test/test_word.py``.
    """
    cdef Py_ssize_t m = len(u)
    cdef Py_ssize_t n = len(v)

    if start < 0:
        start += n
        if start < 0:
            start = 0
    if start > n:
        return -1
    if m == 0:
        return start
    if m > n - start:
        return -1

    cdef array.array table
    if failure_table is None:
        table = word_failure_table(u)
    else:
        table = failure_table
        if len(table) != m:
            raise ValueError("failure_table must have the same length as u")

    cdef int * uu = u.data.as_ints
    cdef int * vv = v.data.as_ints
    cdef int * t = table.data.as_ints

    cdef Py_ssize_t j = start
    cdef int k = 0
    while j < n:
        if uu[k] == vv[j]:
            j += 1
            k += 1
            if k == m:
                return j - m
        else:
            k = t[k]
            if k < 0:
                k = 0
                j += 1

    return -1


def word_is_factor(array.array u, array.array v):
    r"""
    Return whether the word ``u`` is a factor of the word ``v``.

    A *factor* (also called a subword) is a contiguous subsequence. In
    particular the empty word is a factor of every word.

    This is ``word_find(u, v) != -1``. Use :func:`word_find` directly when the
    position of the occurrence is needed, when the search has to start further
    to the right, or when the same ``u`` is searched in many words and its
    failure table is worth computing only once.

    INPUT:

    - ``u`` -- a word, the one searched for

    - ``v`` -- a word, the one searched in

    EXAMPLES::

        sage: from combisurf.word import word_init, word_is_factor

        sage: u = word_init([1, 2])
        sage: word_is_factor(u, word_init([0, 1, 2, 3]))
        True
        sage: word_is_factor(u, word_init([0, 2, 1, 3]))
        False

    The empty word is a factor of every word and no non-empty word is a factor
    of the empty word::

        sage: word_is_factor(word_init(), word_init([0, 1]))
        True
        sage: word_is_factor(word_init(), word_init())
        True
        sage: word_is_factor(word_init([0]), word_init())
        False

    A word is a factor of itself::

        sage: w = word_init([0, 0, 1, 0, 0, 1, 0])
        sage: word_is_factor(w, w)
        True

    Repetitive words, on which the naive search is quadratic::

        sage: u = word_init([0] * 20 + [1])
        sage: word_is_factor(u, word_init([0] * 500))
        False
        sage: word_is_factor(u, word_init([0] * 500 + [1]))
        True

    To test whether ``u`` is a factor of ``v`` read cyclically, search in the
    concatenation of ``v`` with itself::

        sage: u = word_init([3, 0, 1])
        sage: v = word_init([0, 1, 2, 3])
        sage: word_is_factor(u, v)
        False
        sage: word_is_factor(u, v + v)
        True

    ALGORITHM:

    The Knuth-Morris-Pratt algorithm, see :func:`word_find`.

    TESTS:

    The agreement with the naive search is checked exhaustively in
    ``test/test_word.py``.
    """
    return word_find(u, v) != -1


def word_is_reduced(array.array w):
    r"""
    Return whether the free group word ``w`` is reduced.

    EXAMPLES::

        sage: from combisurf.word import word_init, word_is_reduced

        sage: word_is_reduced(word_init([0]))
        True
        sage: word_is_reduced(word_init([0, 1]))
        False
        sage: word_is_reduced(word_init([0, 2, 1]))
        True
    """
    if len(w) <= 1:
        return True

    cdef int i
    for i in range(len(w) - 1):
        if w[i] ^ 1 == w[i + 1]:
            return False

    return True


def word_is_cyclically_reduced(array.array w):
    r"""
    Return whether the free group word ``w`` is cyclically reduced.

    EXAMPLES::

        sage: from combisurf.word import word_init, word_is_cyclically_reduced

        sage: word_is_cyclically_reduced(word_init([0]))
        True
        sage: word_is_cyclically_reduced(word_init([0, 2, 3, 0]))
        False
        sage: word_is_cyclically_reduced(word_init([0, 2, 1]))
        False
    """
    if len(w) <= 1:
        return True

    cdef int i
    for i in range(len(w) - 1):
        if w.data.as_ints[i] ^ 1 == w.data.as_ints[i + 1]:
            return False

    return w.data.as_ints[0] ^ 1 != w.data.as_ints[len(w) - 1]


def word_reduce(array.array w):
    r"""
    EXAMPLES::

        sage: from combisurf.word import word_init, word_reduce

        sage: w = word_init([0, 2])
        sage: word_reduce(w)
        array('i', [0, 2])

        sage: w = word_init([0, 0, 2, 1, 1])
        sage: word_reduce(w)
        array('i', [0, 0, 2, 1, 1])

        sage: w = word_init([0, 0, 2, 3, 1, 1])
        sage: word_reduce(w)
        array('i')

        sage: w = word_init([0, 0, 1])
        sage: word_reduce(w)
        array('i', [0])

        sage: w = word_init([0, 2, 1, 0, 3, 1])
        sage: word_reduce(w)
        array('i')
    """
    if len(w) <= 1:
        return w
    cdef int i = 1
    ans = array.array('i', [w[0]])
    while i < len(w):
        if ans and w[i] ^ 1 == ans[len(ans) - 1]:
            ans.pop()
        else:
            ans.append(w[i])
        i += 1
    return ans


def word_cyclically_reduce(array.array w):
    r"""
    Return a cyclic reduction of ``w`` in the free group.

    EXAMPLES::

        sage: from combisurf.word import word_init, word_cyclically_reduce
        sage: w = word_init([0, 0, 2, 1, 1])
        sage: word_cyclically_reduce(w)
        array('i', [2])
        sage: w = word_init([0, 0, 1])
        sage: word_cyclically_reduce(w)
        array('i', [0])
        sage: w = word_init([0, 2, 0, 1, 0, 3, 1])
        sage: word_cyclically_reduce(w)
        array('i', [0])
        sage: word_cyclically_reduce(word_init())
        array('i')
    """
    if len(w) <= 1:
        return w
    ans = word_reduce(w)
    i = 0
    while i < len(ans) and ans[i] ^ 1 == ans[len(ans) - i - 1]:
        i += 1
    return ans[i:len(ans) - i]


def word_free_group_mul(array.array u, array.array v):
    r"""
    Return the multiplication of the words ``u`` and ``v`` in the free group.

    This is an enhanced version of concatenation where the output is reduced
    whenever both inputs ``u`` and ``v`` are.

    EXAMPLES::

        sage: from combisurf.word import word_init, word_free_group_mul
        sage: u = word_init([0, 2, 1])
        sage: v = word_init([0, 0])
        sage: word_free_group_mul(u, v)
        array('i', [0, 2, 0])
    """
    cdef int i = len(u) - 1
    cdef int j = 0
    while i >= 0 and j < len(v) and (u.data.as_ints[i] ^ 1 == v.data.as_ints[j]):
        i -= 1
        j += 1
    return u[:i+1] + v[j:]


def word_free_group_mul_inplace(array.array u, array.array v):
    r"""
    Update ``u`` by adding ``v`` to it and reducing.

    EXAMPLES::

        sage: from combisurf.word import word_init, word_free_group_mul_inplace
        sage: u = word_init([0, 2, 1])
        sage: v = word_init([0, 0])
        sage: word_free_group_mul_inplace(u, v)
        sage: u
        array('i', [0, 2, 0])
    """
    cdef int i = 0
    while i < len(u) and i < len(v) and u[len(u) - i - 1] ^ 1 == v[i]:
        i += 1
    del u[len(u) - i:]
    # TODO: useless copy!
    u.extend(v[i:])


def word_free_group_inverse(array.array w):
    r"""
    Return the inverse of ``w`` in the free group.

    EXAMPLES::

        sage: from combisurf.word import word_init, word_free_group_inverse
        sage: u = word_init([0, 2, 1, 3])
        sage: word_free_group_inverse(u)
        array('i', [2, 0, 3, 1])
        sage: word_free_group_inverse(word_init())
        array('i')
    """
    ans = array.clone(w, len(w), False)
    cdef int i
    for i in range(len(w)):
        ans[i] = w[len(w) - i - 1] ^ 1
    return ans
