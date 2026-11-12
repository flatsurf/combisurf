r"""
Reduced maps.

A reduced map is a map with a single vertex and a single face.
"""

import numbers

from sage.all import catalan_number, binomial, falling_factorial, randrange, choice

from combisurf.oriented_map import OrientedMap


def uniform_reduced_map(g):
    r"""
    Return a uniform oriented map of genus ``g`` with a single vertex and a single face.

    EXAMPLES::

        sage: from extra.reduced_map import uniform_reduced_map

        sage: for g in range(2, 10):
        ....:     m = uniform_reduced_map(3)
        ....:     assert m.num_vertices() == 1
        ....:     assert m.num_faces() == 1
        ....:     assert m.genus() == 3
        ....:     assert m.is_connected()
    """
    if not isinstance(g, numbers.Integral):
        raise TypeError("g must be a non-negative integer")
    if g < 0:
        raise ValueError("g must be a non-negative integer")
    if g == 0:
        return OrientedMap()
    _, trace = split_and_join(g)
    vp = vertex_permutation_from_trace(trace)
    return OrientedMap(vp)


def split_and_join_with_failure(g):
    r"""
    Run the split-and-join Markov chain on partitions starting from [4g]
    and fails if a [0] appears.

    OUTPUT: a pair ``(validity, trace)`` where ``validity`` is a boolean that
    is ``True`` if it is a complete trace and ``trace`` is the sequence of
    split and join instructions.
    """
    partition = [4*g]
    trace = []
    for n in range(4 * g, 2, -2):
        # choose at random a split or a join with respect to the atom sizes
        # n*(n-1) ordered pairs
        assert sum(partition) == n
        i = randrange(n)
        j = randrange(n - 1)
        if j >= i:
            j += 1
        else:
            i, j = j, i
        assert 0 <= i < j < n
        tot = 0
        pi = pj = -1
        for k, part in enumerate(partition):
            if tot <= i < tot + part:
                pi = k
            if tot <= j < tot + part:
                pj = k
            tot += part
        assert pi != -1 and pj != -1
        if pi == pj:
            # split
            part = partition[pi]
            s = randrange(1, part)
            trace.append(("s", s, part - s))
            if s == 1 or s == part - 1:
                # failure: split trivial piece
                return False, trace
            partition[pi] = s - 1
            partition.append(part - s - 1)
        else:
            # join
            trace.append(("j", partition[pi], partition[pj]))
            if partition[pi] == partition[pj] == 1:
                # failure: join [1] and [1]
                return False, trace
            partition[pi] += partition.pop(pj) - 2

    assert partition == [1, 1]
    trace.append(("j", 1, 1))
    return True, trace


def split_and_join(g):
    r"""
    Run the split-and-join Markov chain on partitions starting from [4g]
    until it finds a valid trace.

    OUTPUT: a pair ``(number_of_tries, trace)`` where ``number_of_tries`` is the number
    of run of the Markov chain that had to be performed before finding a valid trace.
    """
    tries = 1
    is_valid, trace = split_and_join_with_failure(g)
    while not is_valid:
        is_valid, trace = split_and_join_with_failure(g)
        tries += 1
    return tries, trace


def vertex_permutation_from_trace(trace):
    r"""
    Given a sequence of split and join of partitions, build an associated involution without fixed
    points giving a unicellular map.
    """
    # the edge permutations is (0 1)(2 3) ... (4g-2 4g-1)
    # we build the vertex permutation such that multiplied on the right
    # by the edge permutation it gives a 4g-cycle
    # join: (i A) (j B) (i j) = (i A j B)
    # split: (i A j B) (i j) = (i A) (j B)
    assert len(trace) % 2 == 0
    g = len(trace) // 2
    cycles = [[0], [1]]
    assert trace.pop() == ("j", 1, 1)
    n = 1
    while trace:
        op, s0, s1 = trace.pop()
        if op == "s":
            # invert the split [s0+s1] -> [s0-1,s1-1]
            cycle_pos0 = [i for i, c in enumerate(cycles) if len(c) == s0 - 1]
            assert cycle_pos0, (op, s0, s1, cycles)
            i0 = choice(cycle_pos0)
            cycle0 = cycles.pop(i0)

            cycle_pos1 = [i for i, c in enumerate(cycles) if len(c) == s1 - 1]
            assert cycle_pos1, (op, s0, s1, cycles)
            i1 = choice(cycle_pos1)
            cycle1 = cycles.pop(i1)

            p0 = randrange(len(cycle0))
            p1 = randrange(len(cycle1))
            cycles.append([2*n] + cycle0[p0:] + cycle0[:p0] + [2*n+1] + cycle1[p1:] + cycle1[:p1])

        elif op == "j":
            # invert the join [s0], [s1] -> [s0+s1-2]
            cycle_pos = [i for i, c in enumerate(cycles) if len(c) == s0 + s1 - 2]
            assert cycle_pos, (op, s0, s1, cycles)
            i = choice(cycle_pos)
            cycle = cycles.pop(i)
            p = randrange(len(cycle))
            cycle = cycle[p:] + cycle[:p] 
            cycles.append([2*n] + cycle[:s0-1])
            cycles.append([2*n+1] + cycle[s0-1:])

        else:
            raise ValueError("invalid trace")

        n += 1

    return cycles
