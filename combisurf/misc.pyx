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


