r"""
Partial permutations on `\{0, 1, ..., n-1\}`.

TODO:

- In many situations, we need a bitarray of the size of
  the permutation (conjugation, composition, etc). But
  such array would better not be allocated each time the
  function is called.

- There are several functions that return a list of lists
  that could be encoded in a more compact form. Namely
  this only requires two attributes ``content`` and ``bounds``
  where each list is a slice content[bounds[i]:bounds[i+1]].
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
from math import log

import sage.all  # no-cython-lint
from sage.misc.prandom import shuffle, randint
from sage.arith.functions import lcm

from combisurf.misc cimport str_to_int


def argmin(l):
    r"""
    Return the position of the minimal element in the list ``l``.

    EXAMPLES::

        sage: from combisurf.permutation import argmin
        sage: argmin([3, 0, 1, 2])
        1
        sage: argmin([-1, 3, 5, -2, 50])
        3
    """
    if not l:
        raise ValueError("empty list")
    imin = 0
    jmin = l[0]
    for i, j in enumerate(l):
        if j < jmin:
            jmin = j
            imin = i
    return imin


def least_rotation(S):
    """
    Return the pair ``(index of smallest rotation, period)`` of the list ``S``.

    Implementation of Booth's algorithm.

    EXAMPLES::

        sage: from combisurf.permutation import least_rotation

        sage: least_rotation([1,0])
        (1, 2)
        sage: least_rotation([0,1,0])
        (2, 3)
        sage: least_rotation([0,1,1,0,1])
        (3, 5)
        sage: least_rotation([0,1,0,1,1])
        (0, 5)

    Adding some period to the above examples::

        sage: least_rotation([1,0]*4)
        (1, 2)
        sage: least_rotation([0,1,0]*5)
        (2, 3)
        sage: least_rotation([0,1,1,0,1]*3)
        (3, 5)
        sage: least_rotation([0,1,0,1,1]*2)
        (0, 5)
    """
    l = len(S)
    S = S + S
    f = [-1] * len(S)  # failure function
    k = 0              # least rotation of string found so far
    period = 0
    for j in range(1, len(S)):
        sj = S[j]
        i = f[j - k - 1]
        while i != -1 and sj != S[k + i + 1]:
            if sj < S[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if sj != S[k + i + 1]:  # if sj != S[k+i+1], then i == -1
            if sj < S[k]:  # better index
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    # NOTE: the loop below could probably be included in the above loop
    for period in range(l):
        if f[period+l] == l:
            return (k, period)
    return (k, l)

#####################################################################
# Initialization and conversion
#####################################################################

def perm_check(l, int n=-1, involution=None):
    r"""
    Checks that ``l`` is a partial permutation of `\{0, 1, ..., n-1\}`.

    INPUT:

    - ``n`` - integer (optional)

    EXAMPLES::

        sage: from combisurf.permutation import perm_check
        sage: from array import array

    Good permutations::

        sage: perm_check(array('i', [1, 0, 3, 2]), 4)
        True
        sage: perm_check(array('i', [-1]), 1)
        True
        sage: perm_check(array('i', [-1, 3, -1, 1]), 4)
        True

        sage: perm_check(array('i', [1,0,-1,-1,-1]), 2)
        True

    Bad permutations::

        sage: perm_check(array('i', [1, 0, 3, 2]), 3)
        False
        sage: perm_check(array('i', [2, 0]))
        False
        sage: perm_check(array('i', [1, 0, 1]))
        False

    With involution::

        sage: perm_check(array('i', [2,1,0]), involution=array('i', [2,1,0]))
        True
        sage: perm_check(array('i', [1,0,2]), involution=array('i', [2,1,0]))
        False
    """
    if not isinstance(l, array.array) or l.typecode != 'i':
        return False
    if n == -1:
        n = len(l)
    else:
        if len(l) < n:
            return False

    seen = [False] * n
    for i in range(n):
        if l[i] == -1:
            continue
        if l[i] < 0 or l[i] >= n or seen[l[i]]:
            return False
        seen[l[i]] = True

    if involution is not None:
        for i in range(n):
            if l[involution[i]] != involution[l[i]]:
                return False

    return True


def perm_trim(array.array p):
    r"""
    Clear trailing `-1` from the array ``p`` inplace.

    EXAMPLES::

        sage: from combisurf.permutation import perm_init, perm_trim
        sage: p = perm_init([3, 1, 0, 2, -1])
        sage: p
        array('i', [3, 1, 0, 2, -1])
        sage: perm_trim(p)
        sage: p
        array('i', [3, 1, 0, 2])
    """
    cdef int n = len(p) - 1
    while n >= 0 and p.data.as_ints[n] == -1:
        n -= 1
    array.resize(p, n + 1)


def perm_id(int n):
    r"""
    Return the identity permutation.

    EXAMPLES::

        sage: from combisurf.permutation import perm_id

        sage: perm_id(4)
        array('i', [0, 1, 2, 3])
    """
    return array.array('i', range(n))


def perm_init(data, int n=-1, edge_like=False, partial=False):
    """
    Return a permutation from the given data.

    If data is a list of integers, then they are considered to be
    the images of the permutation. If ``data`` is a list of list
    then each list in ``data`` is considered as a cycle. Finally,
    string input in cycle notation is allowed.

    EXAMPLES::

        sage: from combisurf.permutation import perm_init

    As a list of integers::

        sage: perm_init([3,2,1,4])
        array('i', [3, 2, 1, 4])
        sage: perm_init([3,1,-1,0])
        array('i', [3, 1, -1, 0])

    Cycle notation (not mentioned elements are considered to be fixed
    point)::

        sage: perm_init(([2,1],[3,4,0]))
        array('i', [3, 2, 1, 4, 0])
        sage: perm_init([[2,1],[3,4,0]])
        array('i', [3, 2, 1, 4, 0])

    As a string::

        sage: perm_init('(0,1)(3,2)')
        array('i', [1, 0, 3, 2])

        sage: perm_init("(2,4)")
        array('i', [0, 1, 4, 3, 2])
        sage: perm_init("(2,4)", partial=True)
        array('i', [-1, -1, 4, -1, 2])

    Edge-like::

        sage: perm_init("(0,~1,1)", edge_like=True)
        array('i', [3, 1, 0, 2])
        sage: perm_init("(0,~1,1)", edge_like=True, partial=True)
        array('i', [3, -1, 0, 2])

    Zerology::

        sage: perm_init([])
        array('i')
        sage: perm_init([[]])
        array('i')
        sage: perm_init('')
        array('i')
        sage: perm_init('()')
        array('i')
    """
    if isinstance(data, (array.array, tuple, list)):
        if not data:
            if n is not None:
                return array.array('i', range(n))
            else:
                return array.array('i', [])
        elif isinstance(data[0], (tuple, list)):
            return perm_from_cycles(data, n=n, edge_like=False, partial=partial)
        elif n == -1 or len(data) == n:
            return array.array('i', data)
        else:
            raise ValueError("invalid arguments (data={} n={})".format(data, n))

    if isinstance(data, str):
        c = str_to_cycles(data)
        return perm_from_cycles(c, n=n, edge_like=edge_like, partial=partial)


def perm_from_cycles(t, int n=-1, edge_like=False, partial=False):
    r"""
    Return a permutation on `[0, n-1]` from a list of cycles on `[0, n-1]`

    INPUT:

    - ``t`` - cycles

    - ``n`` - optional domain size

    - ``edge_like`` (optional boolean)

    - ``partial`` (optional boolean)

    EXAMPLES::

        sage: from combisurf.permutation import perm_from_cycles

        sage: perm_from_cycles([[1,3,5],[0,2,4],[6]])
        array('i', [2, 3, 4, 5, 0, 1, 6])

        sage: perm_from_cycles([])
        array('i')
        sage: perm_from_cycles([[],[]])
        array('i')

        sage: perm_from_cycles([[1,-2], [0,3]], n=8, edge_like=True)
        array('i', [6, 1, 3, 2, 4, 5, 0, 7])
        sage: perm_from_cycles([[1,-2], [0,3]], n=8, edge_like=True, partial=True)
        array('i', [6, -1, 3, 2, -1, -1, 0, -1])

        sage: perm_from_cycles([], n=4)
        array('i', [0, 1, 2, 3])
        sage: perm_from_cycles([], n=4, edge_like=True)
        array('i', [0, 1, 2, 3])
        sage: perm_from_cycles([], n=4, partial=True)
        array('i', [-1, -1, -1, -1])
    """
    if not any(tt for tt in t):
        if n != -1:
            if partial:
                return array.array('i', [-1] * n)
            else:
                return array.array('i', range(n))
        else:
            return array.array('i', [])

    if n == -1:
        n = max(map(max, t)) + 1
        if edge_like:
            n = max(n, -min(map(min, t)))
            n *= 2
    if edge_like and n % 2:
        raise ValueError("n={} must be even when edge_like=True".format(n))

    res = array.array('i', [-1] * n) if partial else array.array('i', range(n))

    for c in t:
        a = int(c[0])
        if edge_like:
            a = 2 * a if a >= 0 else (2 * ~a + 1)
        if not 0 <= a < n:
            raise ValueError("invalid input: entry {} is out of range".format(a))
        for j in range(1, len(c)):
            b = int(c[j])
            if edge_like:
                b = 2 * b if b >= 0 else (2 * ~b + 1)
            if not 0 <= b < n:
                raise ValueError("invalid input: entry {} is out of range".format(b))
            res[a] = b
            a = b
        b = int(c[0])
        if edge_like:
            b = 2 * b if b >= 0 else (2 * ~b + 1)
        if not 0 <= b < n:
            raise ValueError("invalid input: {} is out of range".format(b))
        res[a] = b

    return res


def str_to_cycles(s):
    """
    Return a list of cycles from a string.

    EXAMPLES::

        sage: from combisurf.permutation import str_to_cycles
        sage: str_to_cycles("(0,1)(3,2)")
        [[0, 1], [3, 2]]

    An integer ``i`` preceded by a ``"~"`` is interpreted as ``-i-1``::

        sage: str_to_cycles("(0,1,2)(~0,~1,~2)")
        [[0, 1, 2], [-1, -2, -3]]
    """
    if not isinstance(s, str):
        raise TypeError(f"s must be a string (got {type(s).__name__})")
    r = []
    if not s:
        return r
    if not s.startswith("(") or not s.endswith(")"):
        raise ValueError("invalid string to initialize a permutation")
    for c_str in s[1:len(s)-1].split(")("):
        if not c_str:
            continue
        r.append([str_to_int(c) for c in c_str.replace(' ', '').split(',')])
    return r


def str_to_cycles_and_data(s):
    r"""
    Return a list of cycles and data from a string.

    EXAMPLES::

        sage: from combisurf.permutation import str_to_cycles_and_data
        sage: str_to_cycles_and_data('(0:1,1:2)')
        ([[0, 1]], {0: 1, 1: 2})
        sage: str_to_cycles_and_data('(0:1,1:2)(3:0)')
        ([[0, 1], [3]], {0: 1, 1: 2, 3: 0})
        sage: str_to_cycles_and_data('(0,1:1)')
        ([[0, 1]], {0: 0, 1: 1})
        sage: str_to_cycles_and_data('(0:1,1)')
        ([[0, 1]], {0: 1, 1: 0})
        sage: str_to_cycles_and_data('(0:1,1:1)')
        ([[0, 1]], {0: 1, 1: 1})
    """
    r = []
    data = {}
    for c_str in s[1:-1].split(')('):
        if not c_str:
            continue
        cycle = []
        for c_dat in c_str.replace(' ', '').split(','):
            if not c_dat:
                raise ValueError('invalid input string')

            pos = c_dat.find(':')
            if pos == -1:
                j = 0
                c = c_dat
            elif c_dat.find(':', pos + 1) != -1:
                raise ValueError('invalid input string: two colons in a row in %s' % repr(c_dat))
            else:
                j = c_dat[pos + 1:]
                c = c_dat[:pos]

            if not c:
                raise ValueError('invalid input string: missing value in %s' % repr(c_dat))
            if c[0] == '~':
                c = c[1:]
                if not c:
                    raise ValueError('invalid input string: missing value in %s' % repr(c_dat))
                i = ~int(c)
            else:
                i = int(c)
            cycle.append(i)
            data[i] = int(j)
        r.append(cycle)
    return r, data


def perm_random(int n):
    r"""
    Return a random permutation.

    EXAMPLES::

        sage: from combisurf.permutation import perm_random, perm_check
        sage: perm_check(perm_random(13), 13)
        True
    """
    r = list(range(n))
    shuffle(r)
    return array.array('i', r)


def perm_random_centralizer(p):
    r"""
    Return a random permutation in the centralizer of ``p``.

    EXAMPLES::

        sage: from combisurf.permutation import *
        sage: p = perm_random(10)
        sage: q = perm_random_centralizer(p)
        sage: perm_compose(p, q) == perm_compose(q, p)
        True
    """
    if not p:
        return p

    cyc = perm_cycles(p)
    cyc.sort(key = lambda c: len(c))
    i = 0
    cdef array.array ans = array.array('i', [-1] * len(p))
    while i < len(cyc):
        j = i + 1
        s = len(cyc[i])
        while j < len(cyc) and len(cyc[j]) == s:
            j += 1

        # permutation of the cycles i, i+1, ..., j-1
        m = perm_random(j - i)

        for ii in range(i, j):
            jj = i + m[ii - i]
            shift = randint(0, s - 1)  # random shift
            for k in range(len(cyc[i])):
                ans[cyc[ii][k]] = cyc[jj][(k + shift) % s]

        # next round
        i = j

    return ans


def perm_random_conjugacy_class(c):
    r"""
    Return a random permutation with given conjugacy class ``c``.

    EXAMPLES::

        sage: from combisurf.permutation import perm_random_conjugacy_class, perm_cycle_type

        sage: p = perm_random_conjugacy_class([5, 2])
        sage: perm_cycle_type(p)
        [5, 2]

        sage: p = perm_random_conjugacy_class([7, 3, 3, 1])
        sage: perm_cycle_type(p)
        [7, 3, 3, 1]
    """
    n = sum(c)
    r = list(range(n))
    shuffle(r)
    p = array.array('i', [-1]*n)
    i = 0
    for k in c:
        # add a k-cycle following the list r
        for j in range(i, i+k-1):
            p[r[j]] = r[j+1]
        p[r[i+k-1]] = r[i]
        i += k
    return p

#####################################################################
# Serialization
#####################################################################

CHARS = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+-'
CHARS_INV = {j:i for i,j in enumerate(CHARS)}

def uint_base64_str(n, l=None):
    r"""
    EXAMPLES::

        sage: from combisurf.permutation import uint_base64_str

        sage: uint_base64_str(15)
        'f'
        sage: uint_base64_str(15, 3)
        '00f'
    """
    n = int(n)
    if n < 0:
        raise ValueError("negative input")
    s = ''
    while n:
        s = CHARS[n % 64] + s
        n //= 64
    if l is not None:
        if len(s) > l:
            raise ValueError
        else:
            s = CHARS[0] * (l - len(s)) + s
    return s


def uint_from_base64_str(s):
    r"""
    EXAMPLES::

        sage: from combisurf.permutation import uint_from_base64_str, uint_base64_str

        sage: uint_from_base64_str('mqb')
        91787
        sage: uint_base64_str(91787)
        'mqb'

        sage: uint_from_base64_str('00mqb')
        91787
    """
    n = 0
    d = 1
    for c in reversed(s):
        n += CHARS_INV[c] * d
        d *= 64
    return n


def perm_base64_str(p):
    r"""
    Make a canonical ASCII string out of ``p``.

    EXAMPLES::

        sage: from combisurf.permutation import perm_base64_str, perm_from_base64_str
        sage: from array import array

        sage: perm_base64_str([3,1,0,2])
        '3102'
        sage: s = perm_base64_str(range(2048))
        sage: s
        '00010203...v+v-'
        sage: perm_from_base64_str(s, 2048) == array('i', range(2048))
        True

        sage: perm_from_base64_str('vdh0keigmcjfpxtnrwsouyqba987654321zl', 36)
        array('i', [31, 13, 17, 0, 20, 14, 18, 16, 22, 12, 19, 15, 25, 33, 29, 23, 27, 32, 28, 24, 30, 34, 26, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 35, 21])

    Partial permutation::

        sage: perm_base64_str([-1, 3, -1, 1])
        '~3~1'
    """
    n = len(p)
    l = int(log(n, 64)) + 1 # number of digits used for each entry
    return ''.join('~' * l if i == -1 else uint_base64_str(i, l) for i in p)


def perm_from_base64_str(s, n):
    r"""
    EXAMPLES::

        sage: from combisurf.permutation import perm_from_base64_str, perm_base64_str
        sage: from array import array

        sage: p = array('i', [3,0,2,1])
        sage: s = perm_base64_str(p)
        sage: perm_from_base64_str(s, 4) == p
        True

        sage: r = list(range(3000))
        sage: shuffle(r)
        sage: p = array('i', r)
        sage: perm_from_base64_str(perm_base64_str(p), 3000) == p
        True

    Partial permutation::

        sage: perm_from_base64_str('~3~1', 4)
        array('i', [-1, 3, -1, 1])
    """
    l = int(log(n, 64)) + 1 # number of digits used for each entry
    if len(s) != n * l:
        raise ValueError('wrong length')
    return array.array('i', (-1 if s[i:i+l] == '~' * l else uint_from_base64_str(s[i:i+l]) for i in range(0, len(s), l)))

#####################################################################
# Boolean properties
#####################################################################

def perm_is_one(array.array p, int n=-1):
    r"""
    Return whether ``p`` is the identity permutation.
    """
    if n == -1:
        n = len(p)
    cdef int i
    for i in range(n):
        if p.data.as_ints[i] != i:
            return False
    return True

#####################################################################
# Cycles and action
#####################################################################

def perm_dense_cycles(array.array p, int n=-1):
    r"""
    Return an array of length ``n`` labelling each of the first ``n`` points of
    ``p`` by the index of the cycle it belongs to.

    The cycles are numbered by consecutive integers starting from zero,
    following the order in which they are met while scanning the points
    ``0, 1, ..., n - 1``. Inactive points, encoded by ``-1`` in ``p``, get the
    label ``-1``.

    INPUT:

    - ``p`` -- a permutation

    - ``n`` -- (default: ``-1``) only use the first ``n`` points of ``p``; if
      ``-1`` use them all. It must lie between ``0`` and ``len(p)``, and ``p``
      must map ``[0, n)`` to itself, since a cycle leaving that range has no
      label to be given; a ``ValueError`` is raised otherwise.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_dense_cycles

        sage: p = array('i', [1,3,4,0,5,7,6,2])
        sage: perm_dense_cycles(p)
        array('i', [0, 0, 1, 0, 1, 1, 2, 1])

    The same array could also be constructed from :func:`perm_cycles` as
    follows::

        sage: from combisurf.permutation import perm_cycles
        sage: perm_cycles(p)
        [[0, 1, 3], [2, 4, 5, 7], [6]]
        sage: ans = array('i', [-1] * 8)
        sage: for i, c in enumerate(perm_cycles(p)):
        ....:     for j in c:
        ....:         ans[j] = i
        sage: ans
        array('i', [0, 0, 1, 0, 1, 1, 2, 1])

    The labels are consecutive, whatever the position of the cycles::

        sage: perm_dense_cycles(array('i', [1,0,3,2,5,4]))
        array('i', [0, 0, 1, 1, 2, 2])

    Inactive points, encoded by ``-1``, get the label ``-1`` and are not
    counted::

        sage: perm_dense_cycles(array('i', [1,0,-1,4,3]))
        array('i', [0, 0, -1, 1, 1])
        sage: perm_dense_cycles(array('i', [2,-1,0]))
        array('i', [0, -1, 0])

    TESTS:

    With ``n`` only the first ``n`` points are scanned, and they must be
    stable under ``p``::

        sage: perm_dense_cycles(array('i', [1,0,3,2]), 2)
        array('i', [0, 0])
        sage: perm_dense_cycles(array('i', [1,0,3,2]), 0)
        array('i')
        sage: perm_dense_cycles(array('i', [3,1,2,0]), 2)
        Traceback (most recent call last):
        ...
        ValueError: p does not map [0, 2) to itself
        sage: perm_dense_cycles(array('i', [1,0,3,2]), 9)
        Traceback (most recent call last):
        ...
        ValueError: n (=9) must be between 0 and len(p) (=4)

    .. SEEALSO::

        :func:`perm_cycles`
    """
    if n == -1:
        n = len(p)
    elif n < 0 or n > len(p):
        raise ValueError(f"n (={n}) must be between 0 and len(p) (={len(p)})")

    cdef array.array res = array.array('i', [-1] * n)
    cdef int * pp = p.data.as_ints
    cdef int * rr = res.data.as_ints
    cdef int i, j, k = 0
    for i in range(n):
        if pp[i] == -1 or rr[i] != -1:
            continue
        j = i
        while rr[j] == -1:
            rr[j] = k
            j = pp[j]
            if j < 0 or j >= n:
                raise ValueError(f"p does not map [0, {n}) to itself")
        k += 1
    return res


def perm_dense_cycle_positions(array.array p, int n=-1):
    r"""
    Return an integral array of length ``n`` whose element at index ``i`` is the
    position of ``i`` in its cycle.

    A cycle is read from the smallest point it contains, which is the one at
    position zero, in the order given by ``p``. Inactive points, encoded by
    ``-1`` in ``p``, get the value ``-1``.

    INPUT:

    - ``p`` -- a permutation

    - ``n`` -- (default: ``-1``) only use the first ``n`` points of ``p``; if
      ``-1`` use them all. It must lie between ``0`` and ``len(p)``, and ``p``
      must map ``[0, n)`` to itself, since the position of a point outside that
      range is not recorded; a ``ValueError`` is raised otherwise.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_dense_cycle_positions

        sage: p = array('i', [1,3,4,0,5,7,6,2])
        sage: perm_dense_cycle_positions(p)
        array('i', [0, 1, 0, 2, 1, 2, 0, 3])

    The same array could also be constructed from :func:`perm_cycles` as
    follows::

        sage: from combisurf.permutation import perm_cycles
        sage: ans = array('i', [-1] * 8)
        sage: for c in perm_cycles(p):
        ....:     for pos, j in enumerate(c):
        ....:         ans[j] = pos
        sage: ans
        array('i', [0, 1, 0, 2, 1, 2, 0, 3])

    Inactive points, encoded by ``-1``, gets value ``-1``::

        sage: perm_dense_cycle_positions(array('i', [1,0,-1,4,3]))
        array('i', [0, 1, -1, 0, 1])

        sage: perm_dense_cycle_positions(array('i', [2,-1,0]))
        array('i', [0, -1, 1])

    TESTS:

    With ``n`` only the first ``n`` points are scanned, and they must be
    stable under ``p``::

        sage: perm_dense_cycle_positions(array('i', [1,0,3,2]), 2)
        array('i', [0, 1])
        sage: perm_dense_cycle_positions(array('i', [1,0,3,2]), 0)
        array('i')
        sage: perm_dense_cycle_positions(array('i', [3,1,2,0]), 2)
        Traceback (most recent call last):
        ...
        ValueError: p does not map [0, 2) to itself
        sage: perm_dense_cycle_positions(array('i', [1,0,3,2]), 9)
        Traceback (most recent call last):
        ...
        ValueError: n (=9) must be between 0 and len(p) (=4)

    .. SEEALSO::

        :func:`perm_dense_cycles`
    """
    if n == -1:
        n = len(p)
    elif n < 0 or n > len(p):
        raise ValueError(f"n (={n}) must be between 0 and len(p) (={len(p)})")

    cdef array.array res = array.array('i', [-1] * n)
    cdef int * pp = p.data.as_ints
    cdef int * rr = res.data.as_ints
    cdef int i, j, k
    for i in range(n):
        if pp[i] == -1 or rr[i] != -1:
            continue
        j = i
        k = 0
        while rr[j] == -1:
            rr[j] = k
            j = pp[j]
            if j < 0 or j >= n:
                raise ValueError(f"p does not map [0, {n}) to itself")
            k += 1
    return res


def perm_cycles(array.array p, singletons=True, int n=-1):
    r"""
    Return the cycle decomposition of ``p`` as a list of lists.

    INPUT:

    - ``p`` -- the permutation

    - ``singletons`` -- bool (default: ``True``) - return or not the singletons

    - ``n`` -- (optional) only use the first ``n`` elements of the permutation ``p``

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_cycles

        sage: perm_cycles(array('i', [0,2,1]))
        [[0], [1, 2]]
        sage: perm_cycles(array('i', [0,2,1]), False)
        [[1, 2]]

        sage: perm_cycles(array('i', [2,-1,0]))
        [[0, 2]]

        sage: perm_cycles(array('i', [2,0,1,-1,-1]), n=3)
        [[0, 2, 1]]
    """
    if n == -1:
        n = len(p)
    elif n < 0 or n > len(p):
        raise ValueError

    cdef array.array seen = array.clone(p, n, True)
    cdef list res = []
    cdef list cycle
    cdef int i, j

    for i in range(n):
        if seen.data.as_ints[i] or p.data.as_ints[i] == -1:
            continue
        if p.data.as_ints[i] == i and not singletons:
            continue
        cycle = []
        j = i
        while not seen.data.as_ints[j]:
            seen[j] = True
            cycle.append(j)
            j = p.data.as_ints[j]
            assert 0 <= j < n, (j, n)
        res.append(cycle)

    return res


def perm_cycles_lengths(array.array p, int n=-1):
    r"""
    Return the array of orbit sizes.

    INPUT:

    - ``p`` -- the permutation

    - ``n`` -- (optional) only use the first ``n`` elements of the permutation ``p``

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_cycles_lengths

        sage: perm_cycles_lengths(array('i', [0,2,1]))
        array('i', [1, 2, 2])
        sage: perm_cycles_lengths(array('i', [2,-1,0]))
        array('i', [2, -1, 2])
        sage: perm_cycles_lengths(array('i', [2,0,1,-1,-1]), n=3)
        array('i', [3, 3, 3])
    """
    if n == -1:
        n = len(p)
    elif n < 0 or n > len(p):
        raise ValueError

    cdef array.array lengths = array.array('i', [-1] * n)
    cdef int i, j, m

    for i in range(n):
        if lengths.data.as_ints[i] != -1 or p.data.as_ints[i] == -1:
            continue
        j = i
        m = 0
        while lengths.data.as_ints[j] == -1:
            lengths.data.as_ints[j] = 0
            j = p.data.as_ints[j]
            m += 1
        while lengths.data.as_ints[j] == 0:
            lengths.data.as_ints[j] = m
            j = p[j]

    return lengths


def perm_num_cycles(array.array p, n=-1):
    r"""
    Return the number of cycles of the permutation ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_num_cycles

        sage: perm_num_cycles(array('i', [1,2,3,0]))
        1
        sage: perm_num_cycles(array('i', [0,2,3,1]))
        2
        sage: perm_num_cycles(array('i', [3,2,1,0]))
        2
        sage: perm_num_cycles(array('i', [3,1,2,0]))
        3
        sage: perm_num_cycles(array('i', [0,1,2,3]))
        4
    """
    if n == -1:
        n = len(p)
    cdef array.array seen = array.clone(p, n, True)
    ans = 0
    cdef int i, j
    for i in range(n):
        if seen.data.as_ints[i] or p.data.as_ints[i] == -1:
            continue
        ans += 1
        j = i
        while not seen.data.as_ints[j]:
            seen.data.as_ints[j] = 1
            j = p.data.as_ints[j]
    return ans


def perm_cycle_type(array.array p, int n=-1):
    r"""
    Return the lengths of the cycles of the permutation ``p`` in size of
    decreasing order.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_cycle_type
        sage: perm_cycle_type(array('i', [1,2,3,0]))
        [4]
        sage: perm_cycle_type(array('i', [0,2,3,1]))
        [3, 1]
        sage: perm_cycle_type(array('i', [3,2,1,0]))
        [2, 2]
        sage: perm_cycle_type(array('i', [3,1,2,0]))
        [2, 1, 1]
        sage: perm_cycle_type(array('i', [0,1,2,3]))
        [1, 1, 1, 1]
    """
    if n == -1:
        n = len(p)
    cdef array.array seen = array.clone(p, n, True)
    cdef list c = []
    cdef int i, j, k
    for i in range(n):
        if seen.data.as_ints[i] or p.data.as_ints[i] == -1:
            continue
        k = 0
        j = i
        while not seen.data.as_ints[j]:
            seen.data.as_ints[j] = 1
            k += 1
            j = p.data.as_ints[j]
        c.append(k)
    c.sort(reverse=True)
    return c


def perm_order(array.array p, int n=-1):
    r"""
    Return the multiplicative order of the permutation ``p``.

    EXAMPLES::

        sage: from combisurf.permutation import perm_init, perm_order
        sage: p = perm_init("(1,3)(2,4,6)(5)")
        sage: perm_order(p)
        6
    """
    return lcm(perm_cycle_type(p))


def perm_is_involution(array.array p, int n=-1):
    r"""
    Return whether ``p`` is a non-trivial involution.

    EXAMPLES::

        sage: from combisurf.permutation import perm_init, perm_is_involution
        sage: perm_is_involution(perm_init('(0)(1)(2)'))
        False
        sage: perm_is_involution(perm_init('(0)(1,2)'))
        True
        sage: perm_is_involution(perm_init('(0,2)(1)'))
        True
        sage: perm_is_involution(perm_init('(0,2,1)'))
        False

        sage: perm_is_involution(perm_init('(0,2)', partial=True))
        True
        sage: perm_is_involution(perm_init('(0,2,4)', partial=True))
        False
    """
    if n == -1:
        n = len(p)
    cdef int i, j, is_identity = 1
    for i in range(n):
        if p.data.as_ints[i] == -1:
            continue
        j = p.data.as_ints[i]
        is_identity &= (i == j)
        if p.data.as_ints[j] != i:
            return False
    return not is_identity


def perm_cycles_to_string(list cycles, edge_like=False, data=None):
    r"""
    Return a string representing a list of cycles.

    INPUT:

    - ``cycles`` -- list of cycles

    - ``edge_like`` -- optional boolean (default ``False``)

    - ``data`` -- optional data
    """
    if edge_like:
        if data:
            elt = lambda e: ('~%d' % (e // 2) if e % 2 else '%d' % (e // 2)) + ((':%d' % data[e]) if data[e] else '')
        else:
            elt = lambda e: ('~%d' % (e // 2)) if e % 2 else '%d' % (e // 2)
    elif data:
        elt = lambda e: ('%d' % e) + ((':%d' % data[e]) if data[e] else '')
    else:
        elt = str

    return ''.join(map(lambda x: '(' + ','.join(map(elt, x)) + ')', cycles))


def perm_cycle_string(array.array p, singletons=True, n=-1, edge_like=False):
    r"""
    Return a string representing the cycle decomposition of `p`

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_cycle_string

        sage: perm_cycle_string(array('i', [0, 2, 1]))
        '(0)(1,2)'
        sage: perm_cycle_string(array('i', [0, 2, 1]), False)
        '(1,2)'
        sage: perm_cycle_string(array('i', [0, 1, 2]), False)
        '()'
    """
    c = perm_cycles(p, singletons, n)
    return '()' if not c else perm_cycles_to_string(c, edge_like)


def perm_orbit(array.array p, int i):
    r"""
    Return the forward orbit of ``i`` under the permutation ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_orbit

        sage: perm_orbit(array('i', [0,3,1,2]), 2)
        [2, 1, 3]
    """
    if i < 0 or i >= len(p):
        raise ValueError("permutation index out of range")
    cdef int j
    cdef list res = [i]
    j = p.data.as_ints[i]
    while j != -1 and j != i:
        res.append(j)
        j = p.data.as_ints[j]
    return res


def perm_are_in_same_orbit(array.array p, int i, int j):
    r"""
    Return whether ``i`` and ``j`` belong to the same orbit of the permutation ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_are_in_same_orbit


        sage: p = array('i', [5, 3, 0, 4, 6, 2, 1])
        sage: perm_are_in_same_orbit(p, 0, 0)
        True
        sage: perm_are_in_same_orbit(p, 0, 1)
        False
    """
    if i < 0 or j < 0 or i >= len(p) or j >= len(p):
        raise ValueError("permutation indices out of range")

    if i == j:
        return True
    k = p.data.as_ints[i]
    while k !=i and k != j:
        k = p.data.as_ints[k]
    return k == j


def perm_orbit_size(array.array p, int i):
    r"""
    Return the size of the orbit of ``i`` under the permutation ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_orbit_size

        sage: perm_orbit_size(array('i', [0,3,1,2]), 2)
        3
        sage: perm_orbit_size(array('i', [0,3,1,2]), 0)
        1
    """
    if i < 0 or i >= len(p):
        raise ValueError
    cdef int j, s
    s = 1
    j = p.data.as_ints[i]
    while j != i:
        s += 1
        j = p.data.as_ints[j]
    return s


def perm_preimage(array.array p, int i):
    r"""
    Return the preimage of ``i`` under ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_init, perm_preimage

        sage: p = perm_init("(0,3,1,5)(2,4)")
        sage: perm_preimage(p, 3)
        0
        sage: perm_preimage(p, 2)
        4
    """
    cdef int j = i
    while p.data.as_ints[j] != i:
        j = p.data.as_ints[j]
    return j


def perm_on_array(array.array dest, array.array src, array.array p, int n=-1):
    r"""
    Inplace action of permutation on list-like objects.

    INPUT:

    - ``dest``, ``src`` -- destination and source (could be the same array in which
      case the operation is done inplace)

    - ``p`` -- permutation

    - ``n`` -- (optional) size of permutation

    - ``swap`` -- (optional) a swap function

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import *

        sage: l = array('i', [0, 1, 2, 3, 4])
        sage: p = array('i', [4, 2, 3, 0, 1])
        sage: perm_on_array(l, l, p)
        sage: l
        array('i', [3, 4, 1, 2, 0])

        sage: src = array('i', [0, 1, 2, 3, 4])
        sage: dest = array('i', [-1] * 5)
        sage: perm_on_array(dest, src, p)
        sage: dest
        array('i', [3, 4, 1, 2, 0])
        sage: src
        array('i', [0, 1, 2, 3, 4])
    """
    cdef array.array seen
    cdef int i, j, tmp
    if n == -1:
        n = len(p)
    if len(dest) < n:
        raise ValueError("invalid dest")
    if len(src) < n:
        raise ValueError("invalid src")
    if dest is src:
        seen = array.clone(p, n, True)
        for i in range(n):
            if seen.data.as_ints[i]:
                continue
            seen.data.as_ints[i] = 1
            j = p.data.as_ints[i]
            if j == -1:
                continue
            while not seen.data.as_ints[j]:
                tmp = src.data.as_ints[i]
                src.data.as_ints[i] = src.data.as_ints[j]
                src.data.as_ints[j] = tmp
                seen.data.as_ints[j] = 1
                j = p.data.as_ints[j]
                if j == -1:
                    raise ValueError("invalid permutation p")
    else:
        for i in range(n):
            j = p.data.as_ints[i]
            if not 0 <= j < n:
                raise ValueError("partial permutation p={}".format(p))
            dest.data.as_ints[j] = src.data.as_ints[i]


def perm_on_list(l, array.array p, int n=-1, swap=None):
    r"""
    Inplace permutation action on list like objects.

    EXAMPLES:

    Permutation action on matrix rows::

        sage: from combisurf.permutation import perm_init, perm_compose, perm_on_list

        sage: m1 = matrix(ZZ, 5, 5, 1)
        sage: m2 = matrix(ZZ, 5, 5, 1)
        sage: m = matrix(ZZ, 5, 5, 1)
        sage: p1 = perm_init([4,1,3,2,0])
        sage: p2 = perm_init([1,0,3,4,2])
        sage: perm_on_list(m1, p1, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: perm_on_list(m2, p2, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: perm_on_list(m, perm_compose(p1, p2), swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: m == m2 * m1
        True

        sage: m1 = matrix(ZZ, 5, range(25))
        sage: m2 = matrix(ZZ, 5, range(137, 87, -2))
        sage: m = m1 * m2
        sage: p = perm_init([4, 1, 3, 2, 0])
        sage: perm_on_list(m1, p, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: perm_on_list(m1, p, swap=sage.matrix.matrix0.Matrix.swap_columns)
        sage: perm_on_list(m2, p, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: perm_on_list(m2, p, swap=sage.matrix.matrix0.Matrix.swap_columns)
        sage: perm_on_list(m, p, swap=sage.matrix.matrix0.Matrix.swap_rows)
        sage: perm_on_list(m, p, swap=sage.matrix.matrix0.Matrix.swap_columns)
        sage: m == m1 * m2
        True
    """
    cdef int i, j
    if n == -1:
        n = len(p)
    cdef array.array seen = array.clone(p, n, True)
    for i in range(n):
        if seen.data.as_ints[i]:
            continue
        seen.data.as_ints[i] = 1
        j = p.data.as_ints[i]
        if j == -1:
            continue
        while not seen.data.as_ints[j]:
            if swap:
                swap(l, i, j)
            else:
                l[i], l[j] = l[j], l[i]
            seen.data.as_ints[j] = 1
            j = p.data.as_ints[j]
            if j == -1:
                raise ValueError("invalid permutation p")


def perm_on_edge_array(array.array dest, array.array src, array.array p, int n=-1):
    r"""
    Inplace action of permutation on list-like objects.

    INPUT:

    - ``dest``, ``src`` -- destination and source (could be the same array in which
      case the operation is done inplace)

    - ``p`` -- permutation

    - ``n`` -- (optional) size of permutation

    - ``swap`` -- (optional) a swap function

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import *

        sage: p = array('i', [6, 7, 3, 2, 9, 8, 0, 1, 5, 4])
        sage: l = array('i', [0, 1, 2, 3, 4])
        sage: perm_on_edge_array(l, l, p)
        sage: l
        array('i', [3, 1, 4, 0, 2])

        sage: src = array('i', [0, 1, 2, 3, 4])
        sage: dest = array('i', [-1] * 5)
        sage: perm_on_edge_array(dest, src, p)
        sage: dest
        array('i', [3, 1, 4, 0, 2])
        sage: src
        array('i', [0, 1, 2, 3, 4])
    """
    cdef array.array seen
    cdef int i, j, tmp
    if n == -1:
        n = len(p)
    if n % 2:
        raise ValueError("n must be even")
    if len(dest) < n // 2:
        raise ValueError("invalid dest")
    if len(src) < n // 2:
        raise ValueError("invalid src")
    if len(p) < n:
        raise ValueError("invalid p")
    if dest is src:
        seen = array.clone(p, n // 2, True)
        for i in range(n // 2):
            if seen.data.as_ints[i]:
                continue
            seen.data.as_ints[i] = 1
            j = p.data.as_ints[2 * i] // 2
            assert 0 <= j < n // 2
            if p.data.as_ints[2 * i + 1] != -1 and p.data.as_ints[2 * i] != p.data.as_ints[2 * i + 1] ^ 1:
                raise ValueError("inconsistent relabelling: p[2*{}] = {}, p[2*{}+1] = {}".format(i, p.data.as_ints[2 * i],
                                                                                               i, p.data.as_ints[2 * i + 1]))
            while not seen.data.as_ints[j]:
                tmp = src.data.as_ints[i]
                src.data.as_ints[i] = src.data.as_ints[j]
                src.data.as_ints[j] = tmp
                seen.data.as_ints[j] = 1
                j = p.data.as_ints[2 * j] // 2
                assert 0 <= j < n // 2
    else:
        for i in range(n // 2):
            j = p.data.as_ints[2 * i] // 2
            assert 0 <= j < n
            dest.data.as_ints[j] = src.data.as_ints[i]


def perm_on_edge_list(array.array p, a, int n=-1, swap=None):
    r"""
    Inplace action of permutation on list-like objects.

    INPUT:

    - ``p`` - permutation of even length such that {p(2i), p(2i+1)} = {2j, 2j+1}

    - ``a`` - list, array of length half the length of the permutation

    - ``n`` - (optional) size of permutation

    - ``swap`` - (optional) a swap function

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import *

        sage: l = [0,1,2]
        sage: p = array('i', [3,2,0,1,5,4])
        sage: perm_on_edge_list(p, l)
        sage: l
        [1, 0, 2]
    """
    if n == -1:
        n = len(p)
    if n % 2:
        raise ValueError("n must be even")
    if len(a) < n // 2:
        raise ValueError("invalid argument a")
    if len(p) < n:
        raise ValueError("invalid argument p")
    cdef array.array seen = array.clone(p, n // 2, True)
    cdef int i, j
    for i in range(n // 2):
        if seen.data.as_ints[i]:
            continue
        seen.data.as_ints[i] = 1
        j = p.data.as_ints[2 * i] // 2
        assert 0 <= j < n // 2
        assert p.data.as_ints[2 * i + 1] == -1 or p.data.as_ints[2 * i] == p.data.as_ints[2 * i + 1] ^ 1
        while not seen.data.as_ints[j]:
            if swap:
                swap(a, i, j)
            else:
                tmp = a[i]
                a[i] = a[j]
                a[j] = tmp
            seen.data.as_ints[j] = 1
            assert p.data.as_ints[2 * j + 1] == -1 or p.data.as_ints[2 * j] == p.data.as_ints[2 * j + 1] ^ 1
            j = p.data.as_ints[2 * j] // 2
            assert 0 <= j < n // 2


#####################################################################
# Group operations
#####################################################################

def perm_invert(array.array p, int n=-1):
    r"""
    Return the inverse of the permutation ``l``.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_invert

        sage: perm_invert(array('i', [0, 3, 1, 2]))
        array('i', [0, 2, 3, 1])

        sage: perm_invert(array('i', [2, -1, 5, 0, -1, 3]))
        array('i', [3, -1, 0, 5, -1, 2])
    """
    if n == -1:
        n = len(p)
    cdef array.array res = array.clone(p, n, False)
    cdef int i
    for i in range(n):
        if p.data.as_ints[i] == -1:
            res.data.as_ints[i] = -1
        else:
            res.data.as_ints[p.data.as_ints[i]] = i
    return res


def perm_compose(array.array p1, array.array p2, int n=-1):
    r"""
    Return the product ``p1 p2``.

    In the notation ``p1 p2`` we use the right action, in other words
    ``p1`` is applied first.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_compose

        sage: perm_compose(array('i', [0,2,1]), array('i', [0,2,1]))
        array('i', [0, 1, 2])
        sage: perm_compose(array('i', [-1,2,3,1]), array('i', [-1,2,1,3]))
        array('i', [-1, 1, 3, 2])

        sage: perm_compose(array('i', [1,0,2,-1,-1]), array('i', [2,1,0,-1]), 3)
        array('i', [1, 2, 0])
    """
    if n == -1:
        n = len(p1)
    cdef array.array r = array.clone(p1, n, False)
    cdef int i
    for i in range(n):
        if p1.data.as_ints[i] == -1:
            r.data.as_ints[i] = -1
        else:
            r.data.as_ints[i] = p2.data.as_ints[p1.data.as_ints[i]]
    return r


perm_compose_00 = perm_compose


# TODO: do something less stupid
# (do it for each cycle independently, detecting if needed the period)
def perm_pow(array.array p, int k, int n=-1):
    r"""
    Return the power of the permutation ``p``.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_pow

        sage: perm_pow(array('i', [3, 0, 1, 2]), 2)
        array('i', [2, 3, 0, 1])
        sage: perm_pow(array('i', [3, 0, 1, 2]), -1)
        array('i', [1, 2, 3, 0])
    """
    if n == -1:
        n = len(p)
    if k == 0:
        return perm_id(n)

    cdef array.array q
    if k < 0:
        p = perm_invert(p, n)
        k = -k

    q = array.copy(p)
    k -= 1
    while k:
        q = perm_compose(q, p)
        k -= 1
    return q


def perm_compose_10(array.array p1, array.array p2, int n=-1):
    r"""
    Return the product ``p1^(-1) p2``

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_compose, perm_invert, perm_compose_10

        sage: p1 = array('i', [0,5,2,1,3,4])
        sage: p2 = array('i', [3,1,5,4,2,0])
        sage: perm_compose_10(p1, p2) == perm_compose(perm_invert(p1), p2)
        True
        sage: shuffle(p1)
        sage: shuffle(p2)
        sage: perm_compose_10(p1, p2) == perm_compose(perm_invert(p1), p2)
        True
    """
    if n == -1:
        n = len(p1)
    cdef int i
    cdef array.array r = array.clone(p1, n, False)
    for i in range(n):
        r.data.as_ints[p1.data.as_ints[i]] = p2.data.as_ints[i]
    return r


def perm_compose_01(array.array p1, array.array p2, int n=-1):
    r"""
    Return the product ``p1 p2^(-1)``

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_compose, perm_invert, perm_compose_01

        sage: p1 = array('i', [0,5,2,1,3,4])
        sage: p2 = array('i', [3,1,5,4,2,0])
        sage: perm_compose_01(p1, p2) == perm_compose(p1, perm_invert(p2)) # not tested
        True
        sage: shuffle(p1)
        sage: shuffle(p2)
        sage: perm_compose_01(p1, p2) == perm_compose(p1, perm_invert(p2)) # not tested
        True
    """
    raise NotImplementedError


def perm_compose_11(array.array p1, array.array p2, int n=-1):
    r"""
    Return the product `p1^(-1) p2^(-1)`.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_compose, perm_invert, perm_compose_11

        sage: p1 = array('i', [0,5,2,1,3,4])
        sage: p2 = array('i', [3,1,5,4,2,0])
        sage: perm_compose_11(p1, p2) == perm_compose(perm_invert(p1), perm_invert(p2))
        True
        sage: shuffle(p1)
        sage: shuffle(p2)
        sage: perm_compose_11(p1, p2) == perm_compose(perm_invert(p1), perm_invert(p2))
        True

    TESTS::

        sage: from combisurf.permutation import perm_invert, perm_compose
        sage: from itertools import permutations

        sage: for p1 in permutations(range(4)):
        ....:     for p2 in permutations(range(4)):
        ....:         q1 = perm_compose_11(array('i', p1), array('i', p2))
        ....:         q2 = perm_compose(perm_invert(array('i', p1)), perm_invert(array('i', p2)))
        ....:         assert q1 == q2, (p1, p2, q1, q2)
    """
    if n == -1:
        n = len(p1)
    cdef array.array r = array.clone(p1, n, False)
    for i in range(n):
        r.data.as_ints[p1.data.as_ints[p2.data.as_ints[i]]] = i
    return r


def perm_conjugate(array.array p1, array.array p2, int n=-1):
    r"""
    Conjugate ``p1`` by ``p2``.

    Let call ``res`` the output of this function. If ``p1`` was
    mapping ``a`` to ``b`` then ``res`` will map ``p2[a]``
    to ``p2[b]``.

    EXAMPLES::

        sage: from array import array
        sage: from combisurf.permutation import perm_conjugate, perm_random

        sage: p1 = perm_random(23)
        sage: p2 = perm_random(23)
        sage: res = perm_conjugate(p1, p2)
        sage: all(res[p2[i]] == p2[p1[i]] for i in range(23))
        True

        sage: perm_conjugate(array('i', [8, 7, 2, 3, 4, 5, 6, 1, 0]), array('i', [0, 2, 4, 5, 6, 7, 3, 1]))
        Traceback (most recent call last):
        ...
        ValueError: invalid input
    """
    if n == -1:
        n = len(p1)
    if len(p2) < n:
        raise ValueError("invalid input")

    cdef array.array res = array.clone(p1, n, False)
    cdef int i
    for i in range(n):
        res[i] = -1
    for i in range(n):
        if p1.data.as_ints[i] == -1 or p2.data.as_ints[i] == -1:
            continue
        if not (0 <= p1.data.as_ints[i] < n) or not (0 <= p2.data.as_ints[i] < n):
            raise ValueError("invalid input")
        else:
            res.data.as_ints[p2.data.as_ints[i]] = p2.data.as_ints[p1.data.as_ints[i]]
    return res


def perm_conjugate_transposition_inplace(array.array p, int i, int j, int pi_inv=-1, int pj_inv=-1, int n=-1):
    if n == -1:
        n = len(p)
    if len(p) < n:
        raise ValueError("invalid input")
    if i < 0 or i >= n:
        raise ValueError(f"i (={i}) out of range")
    if j < 0 or j >= n:
        raise ValueError(f"j (={j}) out of range")

    cdef int pi = p.data.as_ints[i]
    cdef int pj = p.data.as_ints[j]

    if pi_inv == -1:
        pi_inv = perm_preimage(p, i)
    elif pi_inv < 0 or pi_inv >= n:
        raise ValueError(f"pi_inv (={pi_inv}) out of range")
    elif p.data.as_ints[pi_inv] != i:
        raise ValueError(f"invalid input: pi_inv={pi_inv} mapped to {p.data.as_ints[pi_inv]} and not i={i}")

    if pj_inv == -1:
        pj_inv = perm_preimage(p, j)
    elif pj_inv < 0 or pj_inv >= n:
        raise ValueError(f"pj_inv (={pj_inv}) out of range")
    elif p.data.as_ints[pj_inv] != j:
        raise ValueError(f"invalid input: pj_inv={pj_inv} mapped to {p.data.as_ints[pj_inv]} and not j={j}")

    if p.data.as_ints[i] == j and p.data.as_ints[j] == i:
        #  ij)
        pass
    elif p.data.as_ints[i] == j:
        # (.ij.)
        p.data.as_ints[pi_inv] = j
        p.data.as_ints[j] = i
        p.data.as_ints[i] = pj
    elif p.data.as_ints[j] == i:
        # (.ji.)
        p.data.as_ints[pj_inv] = i
        p.data.as_ints[i] = j
        p.data.as_ints[j] = pi
    else:
        if pi == i:
            # (i)
            p.data.as_ints[j] = j
        else:
            # (.i.)
            assert pi_inv != j
            assert pi != j
            p.data.as_ints[pi_inv] = j
            p.data.as_ints[j] = pi

        if pj == j:
            # (j)
            p.data.as_ints[i] = i
        else:
            # (.j.)
            assert pj_inv != i
            assert pj != i
            p.data.as_ints[pj_inv] = i
            p.data.as_ints[i] = pj

#####################################################################
# Transitivity test
#####################################################################

def perms_transitive_components(p, int n=-1):
    r"""
    Return the list of transitive components of the subgroup generated by the
    permutations ``p``.

    INPUT:

    - ``p`` -- a list of permutations given as lists

    EXAMPLES::

        sage: from combisurf.permutation import perms_transitive_components

        sage: perms_transitive_components([[1,0,2,3],[0,1,3,2]])
        [(0, 1), (2, 3)]

        sage: perms_transitive_components([[2,3,0,1]])
        [(0, 2), (1, 3)]

        sage: perms_transitive_components([[3,1,2,0], [0,3,2,1], [0,1,3,2]])
        [(0, 1, 2, 3)]
    """
    if n == -1:
        n = len(p[0])
    seen = [-1] * n
    cc_num = 0
    for i in range(n):
        if seen[i] != -1:
            continue

        todo = [i]
        seen[i] = cc_num
        while todo:
            j = todo.pop()
            for pp in p:
                k = pp[j]
                if seen[k] == -1:
                    todo.append(k)
                    seen[k] = cc_num
        cc_num += 1

    return [tuple(i for i in range(n) if seen[i] == j) for j in range(cc_num)]


def perms_are_transitive(p, int n=-1):
    """
    Test whether the group generated by the permutations in ``p`` is transitive.

    The function assumes that the list of partial permutations act on
    the same domain (ie the -1 occur at the same positions).

    INPUT:

    - ``p`` - a list of permutations of `[0, n-1]`

    EXAMPLES::

        sage: from combisurf.permutation import perms_are_transitive
        sage: perms_are_transitive([[0,1,2],[0,2,1]])
        False
        sage: perms_are_transitive([[0,1,2],[1,2,0]])
        True

        sage: p0 = [0,2,3,1,7,5,6,4]
        sage: p1 = [7,1,2,3,4,5,6,0]
        sage: p2 = [6,1,2,3,4,5,0,7]
        sage: p3 = [1,0,2,3,4,5,6,7]
        sage: p4 = [0,1,4,5,2,3,6,7]
        sage: perms_are_transitive([p0,p1,p2,p3,p4])
        True
        sage: perms_are_transitive([p0,p1,p2,p3])
        False

    Example with partial permutations::

        sage: p0 = [2, -1, 0, 3, 4]
        sage: p1 = [0, -1, 3, 2, 4]
        sage: p2 = [0, -1, 4, 3, 2]
        sage: perms_are_transitive([p0, p1, p2])
        True
        sage: perms_are_transitive([p0, p1])
        False
    """
    if not p:
        raise ValueError("empty list")

    p0 = p[0]
    if n == -1:
        n = len(p0)

    # compute the connected component of the first active integer
    cc0 = [True if j == -1 else False for j in p0]
    if all(cc0):
        return True

    todo = [next(i for i, j in enumerate(p0) if j != -1)]
    cc0[todo[0]] = True
    while todo:
        j = todo.pop()
        for pp in p:
            k = pp[j]
            if cc0[k] is False:
                todo.append(k)
                cc0[k] = True

    return all(cc0)


def perm_edge_orbits(array.array p, int ne=-1):
    r"""
    Return the list of orbits of the permutation group generated by ``p``.

    The function assumes that the list of partial permutations act on
    the same domain (ie the -1 occur at the same positions).

    EXAMPLES::

        sage: from combisurf.permutation import perm_init, perm_edge_orbits

        sage: p = perm_init([2, 6, 4, -1, 0, -1, 8, -1, 1, -1])
        sage: perm_edge_orbits(p)
        [[0, 1, 2, 3, 4]]
    """
    if ne == -1:
        if len(p) % 2:
            raise ValueError("invalid permutation")
        ne = len(p) // 2

    if len(p) < 2 * ne:
        raise ValueError("inconsistent data")

    # compute the connected component of 0
    viewed = [False] * ne
    cc_list = []
    cdef int i
    for i in range(ne):
        if viewed[i] or p.data.as_ints[2 * i] == -1:
            continue
        cc = [i]
        todo = [i]
        viewed[i] = True
        while todo:
            i = todo.pop()
            j = p.data.as_ints[2 * i] // 2
            if not viewed[j]:
                viewed[j] = True
                cc.append(j)
                todo.append(j)
            j = p[2 * i + 1]
            if j != -1:
                j //= 2
                if not viewed[j]:
                    viewed[j] = True
                    cc.append(j)
                    todo.append(j)
        cc.sort()
        cc_list.append(cc)

    return cc_list
