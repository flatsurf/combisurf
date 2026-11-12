

def markov_chain(g):
    r"""

    The proba to terminate is 1/(2g+1)::

        sage: M = markov_chain(8)
        sage: M
        Digraph on 2782 vertices
        sage: sum(prod(M.edge_label(path[i], path[i+1]) for i in range(len(path) - 1)) for path in M.all_paths((4,), (0,)))
        1/3
        sage: sum(prod(M.edge_label(path[i], path[i+1]) for i in range(len(path) - 1)) for path in M.all_paths((8,), (0,)))
        1/5
        sage: sum(prod(M.edge_label(path[i], path[i+1]) for i in range(len(path) - 1)) for path in M.all_paths((12,), (0,)))
        1/7
        sage: sum(prod(M.edge_label(path[i], path[i+1]) for i in range(len(path) - 1)) for path in M.all_paths((16,), (0,)))
        1/9

    Less costly computation::

        sage: layer = set([(32,)])
        sage: weight = {(32,): 1}
        sage: while layer:
        ....:     new_layer = set().union(u for v in layer for u in M.neighbors_out(v))
        ....:     for v in new_layer:
        ....:         weight[v] = sum(weight[u] * M.edge_label(u, v) for u in M.vertices_in(v))
        ....:     layer = new_layer
    """
    M = DiGraph(loops=False, multiedges=False)
    todo = set([(4*g,)])

    def add_transition(u, v, p):
        if M.has_edge(u, v):
            M.set_edge_label(u, v, M.edge_label(u, v) + p)
        else:
            if not M.has_vertex(v) and v and 0 not in v:
                todo.add(v)
            M.add_edge(u, v, p)

    while todo:
        p = todo.pop()
        s = sum(p)
        num_pairs = ZZ(s * (s - 1) / 2)
        total = 0

        # split
        for i, part in enumerate(p):
            pp = list(p)
            pp.pop(i)
            for j in range(part - 1):
                new_p = pp[:]
                new_p.append(j)
                new_p.append(part - j - 2)
                new_p.sort()
                count = part / ZZ(2)
                total += count
                add_transition(p, tuple(new_p), count / num_pairs)

        # join
        for i in range(len(p)):
            part_i = p[i]
            for j in range(i + 1, len(p)):
                part_j = p[j]
                new_p = list(p)
                new_p.pop(j)
                new_p.pop(i)
                new_p.append(part_i + part_j - 2)
                new_p.sort()
                count = part_i * part_j
                total += count
                add_transition(p, tuple(new_p), count / num_pairs)

        assert total == num_pairs

    return M


def epsilon_prefactor(g, n):
    if not isinstance(g, numbers.Integral):
        raise TypeError
    g = int(g)
    if g < 0:
        raise ValueError
    return sum(falling_factorial(n + 1, 2*g - 2 + len(p)) / prod(factorial(mult) * (2 * part + 1)**mult for part, mult in p.to_exp_dict().items()) for p in Partitions(g)) / 2**(2*g)
