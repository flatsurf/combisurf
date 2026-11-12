r"""
Canonical labelling for 3-constellations.
"""

import collections

# try to do partition refinement: as soon as the leftmost atom is a singleton we
# should stop and return
def initial_partition(constellation):
    # make a partition based on cycle type
    # + whether each (x, p[i][x]) belongs to identical or different p[j] orbits
    dense_cycles = []
    for p in constellation:
        c = perm_cycles(p)
        dense_cycles.append(perm_dense_cycles(p))

    # could already look at intersection of p[0]-orbits, p[1]-orbits, p[2]-orbits, etc



# make something compatible with
# (sortable) half-edge/vertex/edge/face data
# data can be seen as half-edge partitions that can probably be refined by
# iterating the group elements
def refine(partition, perm):
    r"""
    EXAMPLES::

        sage: from combisurf.canonical_labelling import refine
        sage: refine([[0,1,2,3],[4,5,6]], [4,1,2,3,0,5,6])
        [[1, 2, 3], [0], [4], [5, 6]]
        sage: refine([[0,1,2,3],[4,5,6]], [0,1,2,3,4,5,6])
        [[0, 1, 2, 3], [4, 5, 6]]
    """
    n = len(perm)
    partition_new = collections.defaultdict(list)
    partition_index = [-1] * n

    while True:
        print(f"partition={partition}")
        assert sum(len(part) for part in partition) == n
        for i, part in enumerate(partition):
            for j in part:
                partition_index[j] = i

        for i in range(n):
            partition_new[partition_index[i], partition_index[perm[i]]].append(i)

        if len(partition) == len(partition_new):
            return partition
        print(f"partition_new={partition_new}")
        partition = [partition_new[x] for x in sorted(partition_new)]
        partition_new.clear()


def good_starts(p, q):
    r"""
    - ``c`` -- constellation
    """
    n = len(vp)

    vertices = perm_cyles(vp)
    by_degrees = collections.defaultdict(list)
    for v in vertices:
        by_degrees[len(v)].extend(v)

    P = [by_degrees[d] for d in sorted(by_degrees)]
    for i, part in P:
        for j in part:
            Pindex[j] = i

