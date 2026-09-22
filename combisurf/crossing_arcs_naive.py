r"""
Counting crossing chords of a circle, reference implementation

This module holds :func:`crossing_arcs_sweep` and
:func:`crossing_arcs_sweep_sorted`, the pure Python versions of
:func:`combisurf.crossing_arcs.crossing_arcs_sweep` and
:func:`combisurf.crossing_arcs.crossing_arcs_sweep_sorted`, and
:func:`crossing_arcs_double_sum`, which computes the same number by an
`O(n^2)` double sum over the endpoints. They are the reference against which
the Cython sweep is tested in ``test/test_geometric_intersection.py``;
everything else in the package uses the Cython sweep, which is faster than
both at every ``n``. See :mod:`combisurf.crossing_arcs` for the quantity
computed.

EXAMPLES::

    sage: from combisurf.crossing_arcs_naive import crossing_arcs_sweep, crossing_arcs_double_sum
    sage: n = 4
    sage: arcs = {2 * n + 0: [1, 0], 3 * n + 1: [0, 1]}
    sage: crossing_arcs_sweep(n, arcs), crossing_arcs_double_sum(n, arcs)
    (1, 1)
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

from combisurf.partial_sums_naive import PartialSumsNaive


def crossing_arcs_double_sum(n, arcs, symmetric=False):
    r"""
    Return the weighted number of pairs of crossing arcs in ``arcs``, by a
    double sum over the endpoints.

    The arcs are chords between the positions ``0``, ..., ``n - 1`` of a
    circle. Two of them, `A = (i_0, j_0)` and `B = (i_1, j_1)`, cross when
    their endpoints strictly interleave, that is `i_0 < i_1 < j_0 < j_1`; arcs
    that share an endpoint do not cross. Each arc carries a `u`-weight and a
    `v`-weight and the returned value is

    .. MATH::

        S = \sum_{A, B} u(A) v(B) + v(A) u(B)

    over the ordered pairs `(A, B)` of crossing arcs with `i_0 < i_1`.

    INPUT:

    - ``n`` -- positive integer, the number of positions on the circle

    - ``arcs`` -- dictionary; the arc from ``first`` to ``last`` (with
      ``first < last``) is the key ``last * n + first`` and its value is the
      pair ``[u, v]`` of its weights; each arc appears at most once

    - ``symmetric`` -- boolean (default: ``False``); whether the `u`-weight
      and the `v`-weight are known to be equal on every arc; not used here,
      accepted for the signature of :func:`crossing_arcs_sweep`

    This costs `O(n^2)` whatever the number of arcs. See
    :func:`combisurf.crossing_arcs.crossing_arcs_sweep` for the same quantity
    in `O(|\text{arcs}| \log(n))`.

    EXAMPLES::

        sage: from combisurf.crossing_arcs_naive import crossing_arcs_double_sum
        sage: n = 4
        sage: crossing_arcs_double_sum(n, {2 * n + 0: [1, 0], 3 * n + 1: [0, 1]})
        1
        sage: crossing_arcs_double_sum(n, {2 * n + 0: [1, 1], 3 * n + 1: [1, 1]}, True)
        2

    Arcs that share an endpoint do not cross::

        sage: crossing_arcs_double_sum(n, {2 * n + 0: [1, 1], 3 * n + 2: [1, 1]})
        0
    """
    Nu = [[0] * n for _ in range(n)]
    Nv = [[0] * n for _ in range(n)]
    for key, (u, v) in arcs.items():
        last, first = divmod(key, n)
        Nu[first][last] = u
        Nv[first][last] = v

    # NOTE: below is a O(n^2) time version of the two following O(n^4) time sums
    #     sum(Nu[i0][j0] * Nv[i1][j1]
    #         for i0 in range(n)
    #         for j0 in range(i0 + 1, n)
    #         for i1 in range(i0 + 1, j0)
    #         for j1 in range(j0 + 1, n))
    #
    #     sum(Nu[i1][j1] * Nv[i0][j0]
    #         for i0 in range(n)
    #         for j0 in range(i0 + 1, n)
    #         for i1 in range(i0 + 1, j0)
    #         for j1 in range(j0 + 1, n))
    #
    # We optimize the computation of the first sum by transforming Nu and
    # Nv to contain partial sums in respectively i0 and j1 respectively
    # (O(n^2) time).  Then we do a double sum in i1, j0 (O(n^2) time). We
    # reverse the role of Nu and Nv to handle the second sum.
    Nu1 = [l[:] for l in Nu]
    for j in range(n):
        for i in range(j - 1):
            Nu1[i + 1][j] += Nu1[i][j]
    Nv1 = [l[:] for l in Nv]
    for i in range(n):
        for j in range(n - 1, i + 1, -1):
            Nv1[i][j - 1] += Nv1[i][j]

    Nv2 = [l[:] for l in Nv]
    for j in range(n):
        for i in range(j - 1):
            Nv2[i + 1][j] += Nv2[i][j]
    Nu2 = [l[:] for l in Nu]
    for i in range(n):
        for j in range(n - 1, i + 1, -1):
            Nu2[i][j - 1] += Nu2[i][j]

    return sum(Nu1[i1 - 1][j0] * Nv1[i1][j0 + 1] + Nv2[i1 - 1][j0] * Nu2[i1][j0 + 1]
               for i1 in range(1, n - 2) for j0 in range(i1 + 1, n - 1))


def crossing_arcs_sweep(n, arcs, symmetric=False):
    r"""
    Return the weighted number of pairs of crossing arcs in ``arcs``, by a
    sweep over the arcs.

    The input and the output are those of :func:`crossing_arcs_double_sum`.

    The arcs `B = (i_1, j_1)` are visited by increasing right endpoint
    `j_1`. When `B` is visited, the arcs already inserted in the partial sums
    are exactly the arcs `A = (i_0, j_0)` with `j_0 < j_1`, and among those
    `A` crosses `B` if and only if `i_0 < i_1` and not `j_0 \le i_1` (the
    second condition implies the first since `i_0 < j_0`). So an arc `A` is
    inserted as `+u(A)` at position `i_0 + 1` and `-u(A)` at position `j_0`,
    and the sum of the positions `0, \ldots, i_1` is the total `u`-weight of
    the arcs crossing `B` from the left. The same goes for `v`. The arcs with
    the same right endpoint `j_1` are all visited before any of them is
    inserted, since they do not cross each other.

    This is the algorithm of
    :func:`combisurf.crossing_arcs.crossing_arcs_sweep`, written with the
    plain array of :class:`~combisurf.partial_sums_naive.PartialSumsNaive`
    rather than a Fenwick tree, so it costs `O(|\text{arcs}| n)`.

    When ``symmetric`` is ``True`` the `v`-weights are not read, and a single
    partial sums structure is used.

    EXAMPLES::

        sage: from combisurf.crossing_arcs_naive import crossing_arcs_sweep
        sage: n = 4
        sage: crossing_arcs_sweep(n, {2 * n + 0: [1, 0], 3 * n + 1: [0, 1]})
        1
        sage: crossing_arcs_sweep(n, {2 * n + 0: [1, 1], 3 * n + 1: [1, 1]}, True)
        2
        sage: crossing_arcs_sweep(n, {2 * n + 0: [1, 1], 3 * n + 2: [1, 1]})
        0

    It agrees with the double sum::

        sage: from combisurf.crossing_arcs_naive import crossing_arcs_double_sum
        sage: n = 9
        sage: arcs = {}
        sage: for _ in range(20):
        ....:     first, last = sorted(sample(range(n), 2))
        ....:     arcs[last * n + first] = [randint(0, 3), randint(0, 3)]
        sage: crossing_arcs_sweep(n, arcs) == crossing_arcs_double_sum(n, arcs)
        True
        sage: for weights in arcs.values():
        ....:     weights[1] = weights[0]
        sage: crossing_arcs_sweep(n, arcs, True) == crossing_arcs_double_sum(n, arcs)
        True
    """
    keys = sorted(arcs)
    num = len(keys)
    Fu = PartialSumsNaive(n)
    if symmetric:
        Fv = Fu
    else:
        Fv = PartialSumsNaive(n)
    S = 0
    k = 0
    while k < num:
        last = keys[k] // n
        k0 = k
        while k < num and keys[k] // n == last:
            first = keys[k] - last * n
            u, v = arcs[keys[k]]
            if symmetric:
                S += 2 * u * Fu.partial_sum(0, first + 1)
            else:
                S += v * Fu.partial_sum(0, first + 1) + u * Fv.partial_sum(0, first + 1)
            k += 1
        for t in range(k0, k):
            first = keys[t] - last * n
            u, v = arcs[keys[t]]
            Fu.update(first + 1, u)
            Fu.update(last, -u)
            if not symmetric:
                Fv.update(first + 1, v)
                Fv.update(last, -v)
    return S


def crossing_arcs_sweep_sorted(n, ukeys, uweights, vkeys=None, vweights=None):
    r"""
    Return the weighted number of pairs of crossing arcs, the arcs being given
    as two sorted lists, one for the `u`-weights and one for the `v`-weights.

    This is :func:`crossing_arcs_sweep` on the arcs ``arcs`` with
    ``arcs[ukeys[k]][0] = uweights[k]``, ``arcs[vkeys[k]][1] = vweights[k]``
    and the weights missing from these lists equal to zero. When ``vkeys``
    and ``vweights`` are ``None`` the `v`-weights are the `u`-weights, which
    is the ``symmetric`` case of :func:`crossing_arcs_sweep`. The keys are
    those of :func:`crossing_arcs_double_sum`, and each list of keys is
    strictly increasing.

    This is the algorithm of
    :func:`combisurf.crossing_arcs.crossing_arcs_sweep_sorted`: the arcs are
    visited group by group, a group being the arcs with a given right
    endpoint, by merging the two lists. The `u`-arcs of a group are queried
    against the `v`-weights inserted so far and the `v`-arcs against the
    `u`-weights, then all of them are inserted. It is written with the plain
    array of :class:`~combisurf.partial_sums_naive.PartialSumsNaive`, so it
    costs `O((|ukeys| + |vkeys|) n)`.

    EXAMPLES::

        sage: from combisurf.crossing_arcs_naive import crossing_arcs_sweep_sorted
        sage: n = 4
        sage: crossing_arcs_sweep_sorted(n, [2 * n + 0], [1], [3 * n + 1], [1])
        1
        sage: crossing_arcs_sweep_sorted(n, [2 * n + 0, 3 * n + 1], [1, 1])
        2
        sage: crossing_arcs_sweep_sorted(n, [2 * n + 0], [1], [3 * n + 2], [1])
        0

    It agrees with the double sum::

        sage: from combisurf.crossing_arcs_naive import crossing_arcs_double_sum
        sage: n = 9
        sage: ukeys = sorted(set(b * n + a for a, b in (sorted(sample(range(n), 2)) for _ in range(15))))
        sage: vkeys = sorted(set(b * n + a for a, b in (sorted(sample(range(n), 2)) for _ in range(15))))
        sage: uweights = [randint(1, 3) for _ in ukeys]
        sage: vweights = [randint(1, 3) for _ in vkeys]
        sage: arcs = {k: [0, 0] for k in ukeys + vkeys}
        sage: for k, u in zip(ukeys, uweights):
        ....:     arcs[k][0] = u
        sage: for k, v in zip(vkeys, vweights):
        ....:     arcs[k][1] = v
        sage: S = crossing_arcs_sweep_sorted(n, ukeys, uweights, vkeys, vweights)
        sage: S == crossing_arcs_double_sum(n, arcs)
        True
        sage: S = crossing_arcs_sweep_sorted(n, ukeys, uweights)
        sage: S == crossing_arcs_double_sum(n, {k: [u, u] for k, u in zip(ukeys, uweights)})
        True
    """
    Fu = PartialSumsNaive(n)
    if vkeys is None:
        # the v-arcs are the u-arcs, and the two partial sums are the same
        S = 0
        k = 0
        while k < len(ukeys):
            last = ukeys[k] // n
            k0 = k
            while k < len(ukeys) and ukeys[k] // n == last:
                S += 2 * uweights[k] * Fu.partial_sum(0, ukeys[k] - last * n + 1)
                k += 1
            for t in range(k0, k):
                Fu.update(ukeys[t] - last * n + 1, uweights[t])
                Fu.update(last, -uweights[t])
        return S

    Fv = PartialSumsNaive(n)
    lu = len(ukeys)
    lv = len(vkeys)
    S = 0
    iu = iv = 0
    while iu < lu or iv < lv:
        if iv == lv or (iu < lu and ukeys[iu] < vkeys[iv]):
            last = ukeys[iu] // n
        else:
            last = vkeys[iv] // n
        ju = iu
        while ju < lu and ukeys[ju] // n == last:
            S += uweights[ju] * Fv.partial_sum(0, ukeys[ju] - last * n + 1)
            ju += 1
        jv = iv
        while jv < lv and vkeys[jv] // n == last:
            S += vweights[jv] * Fu.partial_sum(0, vkeys[jv] - last * n + 1)
            jv += 1
        for t in range(iu, ju):
            Fu.update(ukeys[t] - last * n + 1, uweights[t])
            Fu.update(last, -uweights[t])
        for t in range(iv, jv):
            Fv.update(vkeys[t] - last * n + 1, vweights[t])
            Fv.update(last, -vweights[t])
        iu = ju
        iv = jv
    return S
