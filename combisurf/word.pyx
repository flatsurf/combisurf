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


def word_init(data):
    r"""
    Initialize a word from ``data``.

    EXAMPLES::

        sage: from combisurf.word import word_init
        sage: word_init([0, 1, 3, 2, 1])
        array('i', [0, 1, 3, 2, 1])
        sage: word_init("0,1,~2,~0")
        array('i', [0, 2, 5, 1])
    """
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

        sage: from array import array
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


def fg_word_is_reduced(array.array w):
    r"""
    Return whether the free group word ``w`` is reduced.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.word import fg_word_is_reduced

        sage: fg_word_is_reduced(array('i', [0]))
        True
        sage: fg_word_is_reduced(array('i', [0, 1]))
        False
    """
    if len(w) <= 1:
        return True

    cdef int i
    for i in range(len(w) - 1):
        if w[i] ^ 1 == w[i + 1]:
            return False

    return True


def fg_word_is_cyclically_reduced(array.array w):
    r"""
    Return whether the free group word ``w`` is cyclically reduced.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.word import fg_word_is_cyclically_reduced
    """
    if len(w) <= 1:
        return True

    cdef int i
    for i in range(len(w) - 1):
        if w.data.as_ints[i] ^ 1 == w.data.as_ints[i + 1]:
            return False

    return w.data.as_ints[0] ^ 1 != w.data.as_ints[len(w) - 1]


def fg_word_reduce(array.array w):
    if len(w) <= 1:
        return w
    cdef int i = 1
    ans = array.array('i', [w[0]])
    while i < len(w):
        if ans and w[i] ^ 1 == ans[-1]:
            ans.pop()
        else:
            ans.append(w[i])
        i += 1
    return ans


def fg_word_mul(array.array u, array.array v):
    r"""
    Return the multiplication of the reduced words ``u`` and ``v``.
    """
    cdef int i = len(u) - 1
    cdef int j = 0
    while i >= 0 and j < len(v) and (u.data.as_ints[i] ^ 1 == v.data.as_ints[j]):
        i -= 1
        j += 1
    return u[:i+1] + v[j:]


def fg_word_inverse(array.array w):
    ans = array.clone(w, len(w), False)
    cdef int i
    for i in range(len(w)):
        ans[i] = w[len(w) - i - 1] ^ 1
    return ans


