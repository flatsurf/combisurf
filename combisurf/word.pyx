r"""
Word on non-negative integers and free group elements
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


def word_init(w):
    pass

def fg_word_init(w):
    pass


def fg_word_is_reduced(array.array w):
    r"""
    Return whether the free group word ``w`` is reduced.

    EXAMPLES::

        sage: word_is_reduced([0])
        True
        sage: word_is_reduced([0, 1])
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

        sage: from combisurf.free_group import word_is_cyclically_reduced
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
    while i >= 0 and j < len(v) and (u.data.as_ints[i] % 2 == v.data.as_ints[j] % 2 and u.data.as_ints[i] != v.data.as_ints[j]):
        i -= 1
        j += 1
    return u[:i+1] + v[j:]


def fg_word_inverse(array.array w):
    ans = array.clone(w, len(w), False)
    cdef int i
    for i in range(len(w)):
        ans[i] = w[len(w) - i - 1] ^ 1
    return ans


