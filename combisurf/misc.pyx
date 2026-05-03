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

cpdef Py_hash_t array_hash(int[:] a):
    cdef Py_hash_t acc

    cdef Py_hash_t _PyTuple_HASH_XXPRIME_1 = <Py_hash_t>11400714785074694791ULL
    cdef Py_hash_t _PyTuple_HASH_XXPRIME_2 = <Py_hash_t>14029467366897019727ULL
    cdef Py_hash_t _PyTuple_HASH_XXPRIME_5 = <Py_hash_t>2870177450012600261ULL

    cdef Py_ssize_t lane
    cdef Py_ssize_t l = len(a)
    cdef Py_ssize_t i

    acc = _PyTuple_HASH_XXPRIME_5;
    for i in range(l):
        lane = a[i]
        acc += lane * _PyTuple_HASH_XXPRIME_2
        acc = ((acc << 31) | (acc >> 33))
        acc *= _PyTuple_HASH_XXPRIME_1

    acc += l ^ (_PyTuple_HASH_XXPRIME_5 ^ 3527539UL)

    if acc == <Py_hash_t> - 1:
        acc = 1546275796

    return acc


cpdef str_to_int(c):
    r"""
    Return a Python integer corresponding to the string ``c`` possibly starting with ``"~"``

    EXAMPLES::

        sage: from combisurf.misc import str_to_int
        sage: str_to_int("3")
        3
        sage: str_to_int("~2")
        -3
    """
    if not isinstance(c, str):
        raise TypeError(f"c must be a string (got {type(c).__name__})")
    if not c:
        raise ValueError("empty string")
    if c[0] == "~":
        c1 = c[1:]
        if not c1:
            raise ValueError(f"invalid string c (={c}) to initialize a half-edge")
        return ~int(c1)
    elif not c:
        raise ValueError(f"invalid string c (={c}) to initialize a half-edge")
    return int(c)
