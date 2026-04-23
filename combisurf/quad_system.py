r"""
Homotopy test on maps

This module introduced three classes: :class:`QuadSystems`, :class:`Geodesic` and :class:`Walk`.
"""
# ****************************************************************************
#  This file is part of combisurf
#
#       Copyright (C) 2026 Oscar Fontaine
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


from combisurf import OrientedMap
from collections import deque



def labels(G):
    r"""
    Assign a label (an integer) to each vertex and each face of G.
    Compute for each half-edge the vertex and the face it belongs. 
    """
    vertices = G.vertices()
    faces = G.faces()
    hedge_to_vertex = [None for _ in G.half_edges()]
    hedge_to_face = [None for _ in G.half_edges()]
    for i in range(len(vertices)):
        for e in vertices[i]:
            hedge_to_vertex[e] = i
    for i in range(len(faces)):
        for e in faces[i]:
            hedge_to_face[e] = i
    return hedge_to_vertex, hedge_to_face


def tree_co_tree(G):
    r"""
    Compute a tree/co-tree decomposition of G.
    Return a list res of length the number of edges of G.
    res[e] is 0 if the edge e is in the tree, 1 if e is in the co-tree and 2 otherwise.
    """
    hedge_to_vertex, hedge_to_face = labels(G)
    vp = G.vertex_permutation(copy=False)
    T_found = [False for _ in G.vertices()]
    T_found[0] = True
    res = [2 for _ in G.edge_indices()]
    e = 0
    path = [e]

    if not T_found[hedge_to_vertex[1]]:
        T_found[hedge_to_vertex[1]]=True
        res[0]=0
        path.append(1)
        e=vp[1]
    else:
        e=vp[0]
    
    while len(path) > 0:
        e1 = G._ep(e)
        if e == path[-1]:
            path.pop()
            e = vp[e1]
        elif not T_found[hedge_to_vertex[e1]]:
            T_found[hedge_to_vertex[e1]] = True
            res[e // 2] = 0
            path.append(e1)
            e = vp[e1]
        else:
            e = vp[e]

    fp = G.face_permutation(copy=False)
    C_found = [False for _ in G.faces()]
    C_found[0] = True
    e = fp[0]
    path = [0]
    
    while len(path) > 0:
        e1 = G._ep(e)
        if e == path[-1]:
            path.pop()
            e = fp[e1]
        elif res[e//2] == 0:
            e = fp[e]
        elif not C_found[hedge_to_face[e1]]:
            C_found[hedge_to_face[e1]] = True
            res[e//2] = 1
            path.append(e1)
            e = fp[e1]
        else:
            e = fp[e]
        
    return res


def tree_contraction(G, treecotree):
    r"""
    Compute the OrientedMap F obtained from G by:
        - contracting all the edges in the tree of treecotree
        - removing all the edges in the co-tree of treecotree
    Compute a list correspondence of length the number of half-edges of G such that correspondence[e] is:
        - None if e is an edge of the tree of treecotree
        - d such that fp[d]=e in G after contracting the edge of the tree if e is an edge of the cotree
        - the corresponding edge of F otherwise
    Compute a list recor such that recor[f] for f an half-edge in F is the corresponding half_edge in G
    Compute a list rank such that rank[e] is:
        - None if e is an edge of the tree of treecotree
        - the index of e around the vertex if e is an edge of the cotree
        - the total number of edge of the cotree between e and the next edge in F

    INPUT:

    treecotree should be a list of lenght the number of edges in G where G[e] is:
        - 0 if e is in the tree
        - 1 if e is in the co-tree
        - 2 otherwise
    """
    fp = G.face_permutation(copy=False)
    correspondence = [None for _ in G.half_edges()]
    rank = [None for _ in G.half_edges()]
    nfp = [None for _ in range(4*G.genus())]
    recor = [None for _ in range(4*G.genus())]
    i = 0
    while treecotree[i] != 2:
        i += 1
    last=0
    correspondence[2 * i] = 0
    correspondence[2 * i + 1] = 1
    recor[0] = 2 * i
    recor[1] = 2 * i + 1
    e = fp[2 * i]
    free_index = 2
    r = 0
    
    while e != 2 * i:
        e1 = G._ep(e)
        if treecotree[e // 2] == 1:
            correspondence[e] = last
            rank[e] = r
            r += 1
            e = fp[e1]
            
        elif treecotree[e//2] == 0:
            e = fp[e]
        else:
            if correspondence[e] == None:
                correspondence[e] = free_index
                correspondence[e1] = free_index+1
                recor[free_index] = e
                recor[free_index + 1] = e1
                free_index += 2
            rank[e] = r
            r = 0
            nfp[last] = correspondence[e]
            last = correspondence[e]
            e = fp[e]
    rank[2 * i] = r
    nfp[last] = 0
    F = OrientedMap(fp=nfp)
    return F, correspondence, recor, rank



def turn_remove(s):
    r"""
    Remove a turn at the end of s if there is one.
    """
    if s:
        (t, num) = s.pop()
        if num > 1:
            s.append((t, num-1))

def turn_remove_left(s):
    r"""
    Remove a turn at the beginning of s if there is one.
    """
    if s:
        (t, num) = s.popleft()
        if num > 1:
            s.appendleft((t, num-1))

def turn_add(s, turn, num):
    r"""
    Add num turns at the end of s.
    """
    if num<=0:
        return
    if len(s) == 0:
        s.append((turn,num))
    else:
        if s[-1][0] == turn:
            (a,b) = s.pop()
            s.append((a, b + num))
        else:
            s.append((turn, num))

def turn_add_left(s, turn, num):
    r"""
    Add num turns at the beginning of s.
    """
    if num<=0:
        return
    if len(s) == 0:
        s.append((turn,num))
    else:
        if s[0][0] == turn:
            (a,b) = s.popleft()
            s.appendleft((a, b + num))
        else:
            s.appendleft((turn, num))
    

def turn_modif(s, mod, d):
    if len(s) == 0:
        return
    (t, n) = s.pop()
    r = (t + mod) % d
    if n > 1:
        s.append((t, n-1))
        s.append((r, 1))
    else:
        turn_add(s, r, 1)

def turn_modif_left(s, mod, d):
    if len(s) == 0:
        return
    (t, n) = s.popleft()
    r = (t + mod) % d
    if n > 1:
        s.appendleft((t, n-1))
        s.appendleft((r, 1))
    else:
        turn_add_left(s, r, 1)


def bracket_removal(Q, geo, s, positive, length, d):
    r"""
    Remove a bracket of length ``length`` at the end of geo with turn sequence s.
    """
    fp=Q.face_permutation(copy=False)
    if positive:
        l = deque([])
        for i in range(length + 1, 0, -1):
            edge = geo[-i]
            edge1 = Q._ep(edge)
            l.append(fp[fp[edge1]])
        for i in range(length + 2):
            geo.pop()
        geo = geo.extend(l)
        s.pop()
        if length != 0:
            s.pop()
        turn_modif(s, -1, d)
        turn_add(s, d - 2, length)
    else:
        l=deque([])
        for i in range(length + 1, 0, -1):
            edge = fp[fp[geo[-i]]]
            edge1 = Q._ep(edge)
            l.append(edge1)
        for _ in range(length + 2):
            geo.pop()
        geo = geo.extend(l)
        s.pop()
        if length != 0:
            s.pop()
        turn_modif(s, 1, d)
        turn_add(s, 2, length)


def bracket_removal_left(Q, geo, s, positive, length, d):
    r"""
    Remove a bracket of length ``length`` at the beginning of geo with turn sequence s.
    """
    fp=Q.face_permutation(copy=False)
    if positive:
        l = deque([])
        for i in range(length + 1):
            edge = geo[i]
            edge1 = Q._ep(edge)
            l.appendleft(fp[fp[edge1]])
        for i in range(length + 2):
            geo.popleft()
        geo = geo.extendleft(l)
        s.popleft()
        if length != 0:
            s.popleft()
        turn_modif_left(s, -1, d)
        turn_add_left(s, d - 2, length)
    else:
        l=deque([])
        for i in range(length + 1):
            edge = fp[fp[geo[i]]]
            edge1 = Q._ep(edge)
            l.appendleft(edge1)
        for _ in range(length + 2):
            geo.popleft()
        geo = geo.extendleft(l)
        s.popleft()
        if length != 0:
            s.popleft()
        turn_modif_left(s, 1, d)
        turn_add_left(s, 2, length)


def test_KMP(u, v):
    r"""
    Test if u is a subword of v in O(|u|+|v|).
    """
    if len(u) == 0:
        return True
    elif len(v) == 0:
        return False
    cnd = 0
    T = [-1]
    for i in range(1, len(u)):
        if u[i] == u[cnd]:
            T.append(u[cnd])
        else:
            T.append(cnd)
            while cnd>=0 and u[i] == u[cnd]:
                cnd = T[cnd]
            cnd += 1
    j = 0
    k = 0
    res = False
    while j < len(v) and not res:
        if u[k] == v[j]:
            j += 1
            k += 1
            if k == len(u):
                res = True
        else:
            k = T[k]
            if k == -1:
                k += 1
                j += 1
    return res





class QuadSystem:

    #TODO : Documentation

    def __init__(self, G, treecotree=None, check=True):
        r"""
        Methods:
            _origin_map: the original OrientedMap
            _genus: the genus of the underlying surface
            _quad: the quad system
            _proj: the projection from _origin_map half-edges to path of length 0 or 2 of _quad
            _turn: a list of edges that gives a coefficient to each half-edge around each vertex corresponding at the turn
        """

        if check:
            G._check()
            G._assert_connected()
            if G.has_folded_edge():
                raise NotImplementedError
            if G.genus() == 0:
                raise ValueError("Cannot look for homotopy in genus 0")
            elif G.genus() == 1:
                raise ValueError("Quad systems and homotopy are not effective in genus 1")

        if G.is_mutable():
            G = G.copy(mutable=False)

        self._origin_map = G
        self._genus = G.genus()

        if treecotree == None:
            treecotree = tree_co_tree(G)
        
        F,cor, _, _ = tree_contraction(G,treecotree)
        vp = F.vertex_permutation(copy=False)
        fp = F.face_permutation(copy=False)
        qvp = [None for _ in range(8 * self._genus)]
        remember = [None for _ in range(4 * self._genus)]
        remember2 = [None for _ in range(4 * self._genus)]
        
        e = fp[0]
        last_edge = 0
        next_edge = 1
        remember[0] = vp[0]
        remember2[0] = 0
        while e != 0:
            qvp[2 * last_edge] = 2 * next_edge
            remember[next_edge] = vp[e]
            remember2[e] = next_edge
            e = fp[e]
            last_edge = next_edge
            next_edge += 1
        qvp[-2] = 0
        for ne in range(4*self._genus):
            qvp[2 * ne + 1] = 2 * remember2[remember[ne]] + 1
        self._quad = OrientedMap(vp=qvp)
        
        qfp = self._quad.face_permutation(copy=False)
        cor2 = []
        for e in G.half_edges():
            if treecotree[e//2] == 0:
                cor2.append([])
            elif treecotree[e//2] == 1:
                start = fp[cor[e]]
                e1 = G._ep(e)
                end = fp[cor[e1]]
                das = 2 * remember2[start] + 1
                dae = 2 * remember2[end]
                cor2.append([das,dae])
            else:
                edge = cor[e]
                da = 2 * remember2[edge] + 1
                cor2.append([da, qvp[da - 1]])
        self._proj = cor2

        self._turn = list(self._quad.half_edges())
        label = 0
        for v in self._quad.vertices():
            for e in v:
                self._turn[e] = label
                label += 1


    def _check(self):
        self._origin_map._check()
        if self._genus < 2:
            raise ValueError("Too low genus for a quad system")
        self._quad._check()

    def __eq__(self, other):
        return (self._origin_map == other._origin_map) and (self._quad == other._quad) and (self._proj == other._proj)

    def __repr__(self, *args, **kwds):
        return self._quad.__repr__(*args, **kwds)

    def turn(self, h1, h2):
        l1 = self._turn[h1]
        l2 = self._turn[h2]
        if l1 // (4 * self._genus) != l2 // (4 * self._genus):
            raise ValueError("The two half_edges does not belong to the same vertex")
        nl1 = l1 % (4 * self._genus)
        nl2 = l2 % (4 * self._genus)
        if nl1 <= nl2:
            return nl2 - nl1
        else:
            return 4 * self._genus + nl2 - nl1


    def rotate_list(self):
        res = []
        d = 4 * self._genus
        for h in self._quad.half_edges():
            res.append([])
            current = h
            for _ in range(d):
                res[h].append(current)
                current = self._quad._vp[current]
        return res
            



class Geodesic:

    # TODO : Documentation + Change name ?

    def __init__(self, Q, geo=None, turn=None, check=False):
        r"""
        Methods:
            _quadsystem: the underlying quad system
            _geodesic: the canonical geodesic representative in the quad system (as a deque !)
            _turn_sequence: the turn sequence associated to _geodesic

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~3,~5,~1,~6,~2,~4,~7)", "(0,~7,6,~1)(~0,7,~4,3)(1,~5,4,~2)(2,~6,5,~3)")
            sage: p = Geodesic(Q)
            sage: TestSuite(p).run()
        """
        self._quadsystem = Q
        if geo == None or not geo:
            self._geodesic = deque([])
            self._turn_sequence = deque([])
        elif turn == None:
            if check:
                raise NotImplementedError
            self._geodesic = deque(geo)
            self._turn_sequence = deque([])
            h0 = geo[0]
            for index in range(1, len(geo)):
                h1 = geo[index]
                turn_add(self._turn_sequence, Q.turn(Q._quad._ep(h0), h1), 1)
                h0 = h1
        else:
            if check:
                raise NotImplementedError
            self._geodesic = deque(geo)
            self._turn_sequence = deque(turn)
                

    def __eq__(self, other):
        return (self._quadsystem == other._quadsystem) and (self._geodesic == other._geodesic)

    def __repr__(self, *args, **kwds):
        return f"Geodesic \"{self._geodesic}\" with turns \"{self._turn_sequence}\""

    def __len__(self):
        r"""
        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~3,~5,~1,~6,~2,~4,~7)", "(0,~7,6,~1)(~0,7,~4,3)(1,~5,4,~2)(2,~6,5,~3)")
            sage: p = Geodesic(Q)
            sage: len(p)
            0
            sage: q = Geodesic(Q, [2, 5, 8, 1])
            sage: len(q)
            4
        """
        return len(self._geodesic)
    
    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~3,~5,~1,~6,~2,~4,~7)", "(0,~7,6,~1)(~0,7,~4,3)(1,~5,4,~2)(2,~6,5,~3)")
            sage: p = Geodesic(Q)
            sage: list(iter(p))
            []
            sage: q = Geodesic(Q, [2, 5, 8, 1])
            sage: list(iter(q))
            [2, 5, 8, 1]
        """
        yield from self._geodesic

    def add_edge(self, e):
        r"""
        Add an edge at the end of self and perform simplification to maintain it geodesical.

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~3,~5,~1,~6,~2,~4,~7)", "(0,~7,6,~1)(~0,7,~4,3)(1,~5,4,~2)(2,~6,5,~3)")
            sage: p = Geodesic(Q)
            sage: p.add_edge(12)
            sage: p.add_edge(5)
            sage: p.add_edge(8)
            sage: p.add_edge(1)
            sage: p
            Geodesic "deque([12, 5, 8, 1])" with turns "deque([(1, 1), (2, 2)])"
            sage: p.add_edge(2)
            sage: p
            Geodesic "deque([10, 1, 12])" with turns "deque([(6, 2)])"
            sage: p.add_edge(13)
            sage: p
            Geodesic "deque([10, 1])" with turns "deque([(6, 1)])"
            sage: p.add_edge(14)
            sage: p.add_edge(9)
            sage: p
            Geodesic "deque([10, 7])" with turns "deque([(7, 1)])"
            sage: p.add_edge(8)
            sage: p.add_edge(15)
            sage: p
            Geodesic "deque([10, 1])" with turns "deque([(6, 1)])"
            sage: p.add_edge(14)
            sage: p.add_edge(5)
            sage: p.add_edge(2)
            sage: p
            Geodesic "deque([10, 7, 10])" with turns "deque([(7, 1), (2, 1)])"
        """

        Q = self._quadsystem._quad
        fp = Q.face_permutation(copy=False)
        d = 4 * self._quadsystem._genus
        if len(self._geodesic) == 0:
            self._geodesic.append(e)
        else:
            pre = self._geodesic[-1]
            newturn = self._quadsystem.turn(Q._ep(pre), e)
            
            if newturn == 0: #spur
                self._geodesic.pop()
                turn_remove(self._turn_sequence)
            elif len(self._turn_sequence) >= 1 and newturn == 1 and self._turn_sequence[-1][0] == 1: # positive bracket of length one
                bracket_removal(Q, self._geodesic, self._turn_sequence, True, 0, d)
            elif len(self._turn_sequence) >= 2 and newturn == 1 and self._turn_sequence[-1][0] == 2 and self._turn_sequence[-2][0] == 1:
                # positive bracket of length more than one
                bracket_removal(Q, self._geodesic, self._turn_sequence, True, self._turn_sequence[-1][1], d)
            elif len(self._turn_sequence) >= 1 and newturn == d - 1 and self._turn_sequence[-1][0] == d - 1: # negative bracket of length one
                bracket_removal(Q, self._geodesic, self._turn_sequence, False, 0, d)
            elif len(self._turn_sequence) >= 2 and newturn == d - 1 and self._turn_sequence[-1][0] == d - 2 and self._turn_sequence[-2][0] == d - 1:
                # negative bracket of length more than one
                bracket_removal(Q, self._geodesic, self._turn_sequence, False, self._turn_sequence[-1][1], d)
            else:
                self._geodesic.append(e)
                turn_add(self._turn_sequence, newturn, 1)


    def add_edge_left(self, e):
        r"""
        Add an edge at the beginning of self and perform simplification to maintain it geodesical.

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~3,~5,~1,~6,~2,~4,~7)", "(0,~7,6,~1)(~0,7,~4,3)(1,~5,4,~2)(2,~6,5,~3)")
            sage: p = Geodesic(Q)
            sage: p.add_edge_left(2)
            sage: p.add_edge_left(1)
            sage: p.add_edge_left(8)
            sage: p.add_edge_left(5)
            Geodesic "deque([5, 8, 1, 2])" with turns "deque([(2, 2), (1, 1)])"
            sage: p.add_edge_left(12)
            sage: p
            Geodesic "deque([10, 1, 12])" with turns "deque([(6, 2)])"
            sage: p.add_edge_left(11)
            sage: p
            Geodesic "deque([1, 12])" with turns "deque([(6, 1)])"
            sage: p.add_edge_left(6)
            sage: p.add_edge_left(9)
            sage: p
            Geodesic "deque([15, 12])" with turns "deque([(7, 1)])"
            sage: p.add_edge_left(8)
            sage: p.add_edge_left(7)
            sage: p
            Geodesic "deque([1, 12])" with turns "deque([(6, 1)])"
            
            sage: p.add_edge_left(6)
            sage: p.add_edge_left(11)
            sage: p.add_edge_left(2)
            sage: p
            Geodesic "deque([4, 15, 12])" with turns "deque([(2, 1), (7, 1)])"
        """
        
        Q = self._quadsystem._quad
        fp = Q.face_permutation(copy=False)
        d = 4 * self._quadsystem._genus
        if len(self._geodesic) == 0:
            self._geodesic.appendleft(e)
        else:
            pre = self._geodesic[0]
            newturn = self._quadsystem.turn(Q._ep(e), pre)
            
            if newturn == 0: #spur
                self._geodesic.popleft()
                turn_remove_left(self._turn_sequence)
            elif len(self._turn_sequence) >= 1 and newturn == 1 and self._turn_sequence[0][0] == 1: # positive bracket of length one
                bracket_removal_left(Q, self._geodesic, self._turn_sequence, True, 0, d)
            elif len(self._turn_sequence) >= 2 and newturn == 1 and self._turn_sequence[0][0] == 2 and self._turn_sequence[1][0] == 1:
                # positive bracket of length more than one
                bracket_removal_left(Q, self._geodesic, self._turn_sequence, True, self._turn_sequence[0][1], d)
            elif len(self._turn_sequence) >= 1 and newturn == d - 1 and self._turn_sequence[0][0] == d - 1: # negative bracket of length one
                bracket_removal_left(Q, self._geodesic, self._turn_sequence, False, 0, d)
            elif len(self._turn_sequence) >= 2 and newturn == d - 1 and self._turn_sequence[0][0] == d - 2 and self._turn_sequence[1][0] == d - 1:
                # negative bracket of length more than one
                bracket_removal_left(Q, self._geodesic, self._turn_sequence, False, self._turn_sequence[0][1], d)
            else:
                self._geodesic.appendleft(e)
                turn_add_left(self._turn_sequence, newturn, 1)


    def origin_simplification(self):
        r"""
        Perform the geodesic simplification at the based point

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~3,~5,~1,~6,~2,~4,~7)", "(0,~7,6,~1)(~0,7,~4,3)(1,~5,4,~2)(2,~6,5,~3)")
            sage: p = Geodesic(Q, [5, 12, 3, 4])
            sage: p.origin_simplification()
            sage: p
            Geodesic "deque([12, 3])" with turns "deque([(7, 1)])"
            
            sage: p = Geodesic(Q, [2, 5, 8, 1])
            sage: p.origin_simplification()
            sage: p
            Geodesic "deque([10, 1])" with turns "deque([(6, 1)])"
            
            sage: p = Geodesic(Q, [0, 9, 4, 3])
            sage: p.origin_simplification()
            sage: p
            Geodesic "deque([2, 11])" with turns "deque([(2, 1)])"

            sage: p = Geodesic(Q, [6, 13, 4, 7, 12, 5])
            sage: p.origin_simplification()
            sage: p
            Geodesic "deque([13, 4, 7, 10])" with turns "deque([(4, 2), (2, 1)])"

            sage: p = Geodesic(Q, [10, 1, 12, 3, 14, 7, 2, 13])
            sage: p.origin_simplification()
            sage: p
            Geodesic "deque([5, 8, 1, 14, 7, 2])" with turns "deque([(2, 2), (7, 1), (2, 1), (6, 1)])"

            sage: p = Geodesic(Q, [1, 12, 3, 14, 7, 2, 13, 10])
            sage: p.origin_simplification()
            sage: p
            Geodesic "deque([14, 7, 2, 5, 8, 1])" with turns "deque([(2, 1), (6, 1), (2, 3)])"

            sage: p = Geodesic(Q, [0, 11, 14, 1, 4, 15])
            sage: p.origin_simplification()
            sage: p
            Geodesic "deque([4, 11])" with turns "deque([(5, 1)])"
        """

        Q = self._quadsystem._quad
        fp = Q.face_permutation(copy=False)
        d = 4 * self._quadsystem._genus
        geo = self._geodesic
        s = self._turn_sequence
        if len(geo) == 0 or len(geo) % 2 == 1:
            return
        first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])
        
        while first_turn == 0 and geo: #spur
            geo.pop()
            geo.popleft()
            turn_remove(s)
            turn_remove_left(s)
            first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])

        if first_turn == 1 and len(s) == 1 and s[0][0] == 2: # whole geodesic bracket
            geo.pop()
            geo.popleft()
            l = deque([])
            for e in geo:
                e1 = Q._ep(e)
                l.append(fp[fp[e1]])
            geo.clear()
            geo.extend(l)
            s.append((d - 2, s[0][1] - 2))
            s.popleft()
            
        elif first_turn == d - 1 and len(s) == 1 and s[0][0] == d - 2: # whole geodesic bracket
            geo.pop()
            geo.popleft()
            l = deque([])
            for e in geo:
                e1 = Q._ep(e)
                l.append(fp[fp[e1]])
            geo.clear()
            geo.extend(l)
            s.append((2, s[0][1] - 2))
            s.popleft()

        if first_turn == 1 and len(s) >= 1 and s[-1][0] == 1: # bracket ending at the origin
            geo.popleft()
            turn_remove_left(s)
            bracket_removal(Q, geo, s, True, 0, d)
            if geo:
                first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])
                while first_turn == 0 and geo:
                    geo.pop()
                    geo.popleft()
                    turn_remove(s)
                    turn_remove_left(s)

        elif first_turn == 1 and len(s) >= 2 and s[-1][0] == 2 and s[-2][0] == 1: # bracket ending at the origin
            geo.popleft()
            turn_remove_left(s)
            bracket_removal(Q, geo, s, True, s[-1][1], d)
            if geo:
                first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])
                while first_turn == 0 and geo:
                    geo.pop()
                    geo.popleft()
                    turn_remove(s)
                    turn_remove_left(s)

        elif first_turn == d - 1 and len(s) >= 1 and s[-1][0] == d - 1: # bracket ending at the origin
            geo.popleft()
            turn_remove_left(s)
            bracket_removal(Q, geo, s, False, 0, d)
            if geo:
                first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])
                while first_turn == 0 and geo:
                    geo.pop()
                    geo.popleft()
                    turn_remove(s)
                    turn_remove_left(s)
                    
        elif first_turn == d - 1 and len(s) >= 2 and s[-1][0] == d - 2 and s[-2][0] == d - 1: # bracket ending at the origin
            geo.popleft()
            turn_remove_left(s)
            bracket_removal(Q, geo, s, False, s[-1][1], d)
            if geo:
                first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])
                while first_turn == 0 and geo:
                    geo.pop()
                    geo.popleft()
                    turn_remove(s)
                    turn_remove_left(s)

        if first_turn == 1 and len(s) >= 1 and s[0][0] == 1: # bracket starting at the origin
            geo.pop()
            turn_remove(s)
            bracket_removal_left(Q, geo, s, True, 0, d)
            if geo:
                first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])
                while first_turn == 0 and geo:
                    geo.pop()
                    geo.popleft()
                    turn_remove(s)
                    turn_remove_left(s)

        elif first_turn == 1 and len(s) >= 2 and s[0][0] == 2 and s[1][0] == 1: # bracket starting at the origin
            geo.pop()
            turn_remove(s)
            bracket_removal_left(Q, geo, s, True, s[0][1], d)
            if geo:
                first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])
                while first_turn == 0 and geo:
                    geo.pop()
                    geo.popleft()
                    turn_remove(s)
                    turn_remove_left(s)

        elif first_turn == d - 1 and len(s) >= 1 and s[0][0] == d - 1: # bracket starting at the origin
            geo.pop()
            turn_remove(s)
            bracket_removal_left(Q, geo, s, False, 0, d)
            if geo:
                first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])
                while first_turn == 0 and geo:
                    geo.pop()
                    geo.popleft()
                    turn_remove(s)
                    turn_remove_left(s)
                    
        elif first_turn == d - 1 and len(s) >= 2 and s[0][0] == d - 2 and s[1][0] == d - 1: # bracket starting at the origin
            geo.pop()
            turn_remove(s)
            bracket_removal_left(Q, geo, s, False, s[0][1], d)
            if geo:
                first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])
                while first_turn == 0 and geo:
                    geo.pop()
                    geo.popleft()
                    turn_remove(s)
                    turn_remove_left(s)

        if first_turn == 2 and len(s) >= 3: # bracket containing the origin
            while first_turn == 2:
                geo.append(geo.popleft())
                (t, n) = s.popleft()
                if n != 1:
                    s.appendleft((t, n - 1))
                (t2, n2) = s.pop()
                if t2==2:
                    s.append((2, n2 + 1))
                else:
                    s.append((t2, n2))
                    s.append((2, 1))
                first_turn = t
            if first_turn == 1 and s[-2][0] == 1:
                bracket_removal(Q, geo, s, True, s[-1][1], d)
                geo.popleft()
                turn_remove_left(s)

        elif first_turn == d - 2 and len(s) >= 3: # bracket containing the origin
            while first_turn == d - 2:
                geo.append(geo.popleft())
                (t, n) = s.popleft()
                if n != 1:
                    s.appendleft((t, n - 1))
                (t2, n2) = s.pop()
                if t2==2:
                    s.append((2, n2 + 1))
                else:
                    s.append((t2, n2))
                    s.append((2, 1))
                first_turn = t
            if first_turn == d - 1 and s[-2][0] == d - 1:
                bracket_removal(Q, geo, s, False, s[-1][1], d)
                geo.popleft()
                turn_remove_left(s)

    def canonical(self):
        r"""
        Perform the right push to the cyclic geodesic path to make it canonical.

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~3,~5,~1,~6,~2,~4,~7)", "(0,~7,6,~1)(~0,7,~4,3)(1,~5,4,~2)(2,~6,5,~3)")
            
            sage: p = Geodesic(Q, [6, 3])
            sage: p.canonical()
            sage: p
            Geodesic "deque([13, 8])" with turns "deque([(6, 1)])"
            
            sage: p = Geodesic(Q, [0, 11, 14, 7, 10, 9, 4, 15])
            sage: p.canonical()
            sage: p
            Geodesic "deque([11, 4, 3, 14, 5, 2, 9, 2])" with turns "deque([(5, 1), (6, 3), (7, 1), (3, 1), (5, 1)])"

            sage: p = Geodesic(Q, [0, 11, 14, 7, 10, 15])
            sage: p.canonical()
            sage: p
            Geodesic "deque([9, 4, 3, 14, 5, 2])" with turns "deque([(6, 4), (7, 1)])"
        """

        Q = self._quadsystem._quad
        fp = Q.face_permutation(copy=False)
        vp = Q.vertex_permutation(copy=False)
        d = 4 * self._quadsystem._genus
        geo = self._geodesic
        s = self._turn_sequence
        
        if len(geo) % 2 == 1:
            raise ValueError("This is not a cyclic geodesic path")
        elif len(geo) == 0:
            return

        self.origin_simplification()
        first_turn = self._quadsystem.turn(Q._ep(geo[-1]), geo[0])
        
        if first_turn == 2 and len(s) == 1 and s[0][0] == 2: # whole geodesic of turn 2
            l = deque([])
            for e in geo:
                e1 = Q._ep(e)
                l.append(fp[fp[e1]])
            geo.clear()
            geo.extend(l)
            m = s[0][1]
            s.clear()
            s.append((d - 2, m))
            
        else:
            i = 0
            finish = False
            while i <= len(geo):
                if first_turn == 1: # left L
                    if s[0][0] == 2 and s[-1][0] == 2:
                        (t1, n1) = s.pop()
                        (t2, n2) = s.popleft()
                    elif s[0][0] == 2:
                        n1 = 0
                        (t2, n2) = s.popleft()
                    elif s[-1][0] == 2:
                        n2 = 0
                        (t1, n1) = s.pop()
                    else:
                        n1 = 0
                        n2 = 0
                        
                    if n1 != 0 and n2 != 0:
                        first_turn = d - 3
                    elif n1 != 0 or n2 != 0:
                        first_turn = d - 2
                    else:
                        first_turn = d - 1
                        
                    if len(s) == 0: # whole path bracket. Should have been cleared at origin_simplification
                        raise ValueError("The path was not geodesic")
                    elif len(s) == 1 and s[0][1] == 1: # whole geodesic is a L
                        (t3, n3) = s.pop()
                        finish = True
                        if t3 == 3:
                            turn_add(s, d - 2, n2)
                            turn_add(s, d - 1, 1)
                            turn_add(s, d - 2, n1)
                            l = deque([])
                            e = geo[0]
                            for j in range(n2):
                                e = geo[1 + j]
                                e1 = Q._ep(e)
                                l.append(fp[fp[e1]])
                            x = fp[vp[e1]]
                            l.append(x)
                            l.append(fp[x])
                            for j in range(n1 + 1, 1, -1):
                                e = geo[-j]
                                e1 = Q._ep(e)
                                l.append(fp[fp[e1]])
                            geo.clear()
                            geo.extend(l)
                        else:
                            turn_add(s, d - 2, n2 - 1)
                            if n2 != 0:
                                turn_add(s, d - 1, 1)
                            turn_add(s, t3 - 2, 1)
                            if n1 != 0:
                                turn_add(s, d - 1, 1)
                            turn_add(s, d - 2, n1 - 1)
                            l = deque([])
                            e = geo[1]
                            e1 = Q._ep(e)
                            for j in range(n2):
                                e = geo[1 + j]
                                e1 = Q._ep(e)
                                l.append(fp[fp[e1]])
                            x = vp[e1]
                            l.append(Q._ep(x))
                            y = geo[n2 + 1]
                            l.append(Q.previous_at_vertex(y))
                            for j in range(n1 + 1, 1, -1):
                                e = geo[-j]
                                e1 = Q._ep(e)
                                l.append(fp[fp[e1]])
                            geo.clear()
                            geo.extend(l)
                            
                    else: # general L
                        turn_modif(s, -1, d)
                        if n1 != 0:
                            turn_add(s, d - 1, 1)
                            turn_add(s, d - 2, n1 - 1)
                        turn_modif_left(s, -1, d)
                        if n2 != 0:
                            turn_add_left(s, d - 1, 1)
                            turn_add_left(s, d - 2, n2 - 1)
                        l = deque()
                        e = geo.popleft()
                        e1 = Q._ep(e)
                        for _ in range(n2):
                            e = geo.popleft()
                            e1 = Q._ep(e)
                            l.appendleft(fp[fp[e1]])
                        geo.appendleft(Q.previous_in_face(e1))
                        geo.extendleft(l)
                        l = deque()
                        e = geo.pop()
                        e1 = Q._ep(e)
                        for _ in range(n1):
                            e = geo.pop()
                            e1 = Q._ep(e)
                            l.appendleft(fp[fp[e1]])
                        l.appendleft(fp[e1])
                        geo.extend(l)
                i+=1
                geo.rotate()
                turn_add_left(s, first_turn, 1)
                (t1, n1) = s.pop()
                first_turn = t1
                if n1 != 1:
                    s.append((t1, n1 - 1))




class Walk:

    # TODO : Documentation + Change name ?

    def __init__(self, Q, walk):
        r"""
        Methods:
            _quadsystem: the underlying quad system
            _walk: the walk in the original OrientedMap (as a list)
            _geodesic: the canonical geodesic representative of walk
        """

        self._quadsystem = Q
        self._walk = walk
        self._geodesic = Geodesic(Q)
        for e in self._walk:
            for f in Q._proj[e]:
                self._geodesic.add_edge(f)
        self._geodesic.canonical()

    def __eq__(self, other):
        return (self._quadsystem == other._quadsystem) and (self._walk == other._walk)

    def is_homotopic(self, other):
        r"""
        Return whether self and other are freely homotopic.
        
        EXAMPLES::
        
            sage: from combisurf import OrientedMap, QuadSystem, Geodesic, Walk
            sage: m = OrientedMap(vp=[[0, 2, 4, 6],[7, 8, 5], [9, 10, 12, 11], [3, 15, 1, 13, 14]])
            sage: Q = QuadSystem(m)
            sage: w1 = Walk(Q, [])
            sage: w2 = Walk(Q, [4, 8, 11, 9, 7])
            sage: w3 = Walk(Q, [6, 5])
            sage: w1.is_homotopic(w2)
            True
            sage: w1.is_homotopic(w3)
            False
            sage: w4 = Walk(Q, [6, 8, 11, 9, 7])
            sage: w5 = Walk(Q, [2, 15])
            sage: w3.is_homotopic(w4)
            True
            sage: w3.is_homotopic(w5)
            False
            sage: w6 = Walk(Q, [11])
            sage: w3.is_homotopic(w6)
            True

        """
        if self._quadsystem != other._quadsystem:
            raise ValueError("The quadsystems are different")
        if len(self._geodesic) != len(other._geodesic):
            return False

        c = self._geodesic._geodesic.copy()
        c.extend(c)
        return test_KMP(other._geodesic._geodesic, c)


    def simplicity(self, oriented_map, edge_to_vertex=None, quadsystem=None):

        r"""
        Return whether self lift into a simple walk in the universal covering.
        """
        
        if not self._walk:
            return True
        elif len(self._walk)==2 and oriented_map._ep(self._walk[0]) == self._walk[1]:
            return False
        if quadsystem is None:
            quadsystem = QuadSystem(oriented_map)
        if edge_to_vertex is None:
            edge_to_vertex = oriented_map.half_edges_to_vertices()

        star = StarShapedSpace(quadsystem, edge_to_vertex[self._walk[0]])
        seen_value = {} 
        previous_vertex = 0
        first_values = (0, edge_to_vertex[self._walk[0]])
        seen_value[first_values] = True

        for i in range(len(self._walk)):
            edge = self._walk[i]
            for elt in quadsystem._proj[edge]:
                previous_vertex = star.insert_edge(previous_vertex, elt)
            map_vertex = edge_to_vertex[oriented_map._ep(edge)]
            if seen_value.get((previous_vertex, map_vertex)) is None:
                seen_value[(previous_vertex, map_vertex)] = True
                #print(seen_value)
            elif i == len(self._walk) - 1 and (previous_vertex, map_vertex) == first_values:
                return True
            else:
                #print(seen_value)
                return False
        return True




class LazyGeodesic:

    # TODO : Documentation. Class of geodesic that contains only the first and the last edge and the sequence of turn. Should not be used to test homotopy ! Only contractibility.
    
    def __init__(self, Q, geo=None, turn=None, check=False):
        r"""
        Methods:
            _quadsystem: the underlying quad system
            _first: the first edge
            _last: the last edge
            _turn_sequence: the turn sequence associated to _geodesic

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, LazyGeodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~3,~5,~1,~6,~2,~4,~7)", "(0,~7,6,~1)(~0,7,~4,3)(1,~5,4,~2)(2,~6,5,~3)")
            sage: p = LazyGeodesic(Q)
            sage: TestSuite(p).run()
        """
        self._quadsystem = Q
        if geo == None or not geo:
            self._first = None
            self._last = None
            self._turn_sequence = deque([])
        elif turn == None:
            if check:
                raise NotImplementedError
            self._first = geo[0]
            self._last = geo[-1]
            self._turn_sequence = deque([])
            h0 = geo[0]
            for index in range(1, len(geo)):
                h1 = geo[index]
                turn_add(self._turn_sequence, Q.turn(Q._quad._ep(h0), h1), 1)
                h0 = h1
        else:
            if check:
                raise NotImplementedError
            self._first = geo[0]
            self._first = geo[-1]
            self._turn_sequence = deque(turn)
                

    def __eq__(self, other):
        return (self._quadsystem == other._quadsystem) and (self._first == other._first) and (self._last == other._last) and (self._turn_sequence == other._turn_sequence)

    def is_empty(self):
        return self._first == None
    
    def __repr__(self, *args, **kwds):
        return f"LazyGeodesic starting with \"{self._first}\", ending with \"{self._last}\" and with turns \"{self._turn_sequence}\""

    def length(self):
        r"""
        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~3,~5,~1,~6,~2,~4,~7)", "(0,~7,6,~1)(~0,7,~4,3)(1,~5,4,~2)(2,~6,5,~3)")
            sage: p = LazyGeodesic(Q)
            sage: p.length()
            0
            sage: q = LazyGeodesic(Q, [2, 5, 8, 1])
            sage: q.length()
            4
        """
        if self._first == None:
            return 0
        n = 1
        for elt in self._turn_sequence:
            n += elt[1]
        return n


    def add_edge(self, e):
        r"""
        Add an edge at the end of self and perform simplification to maintain it geodesical.

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, LazyGeodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6],[7, 8, 5], [9, 10, 12, 11], [3, 15, 1, 13, 14]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~5,~6,~7,~4,~1,~2,~3)", "(0,~3,2,~1)(~0,7,~6,5)(1,~4,3,~2)(4,~7,6,~5)")
            sage: p = LazyGeodesic(Q)
            sage: p.add_edge(12)
            sage: p.add_edge(5)
            sage: p.add_edge(6)
            sage: p.add_edge(11)
            sage: p
            LazyGeodesic starting with "12", ending with "11" and with turns "deque([(4, 1), (1, 1), (2, 1)])"
            sage: p.add_edge(12)
            sage: p
            LazyGeodesic starting with "12", ending with "14" and with turns "deque([(3, 1), (6, 1)])"
            sage: p.add_edge(15)
            sage: p
            LazyGeodesic starting with "12", ending with "3" and with turns "deque([(3, 1)])"
        """

        Q = self._quadsystem._quad
        vp = Q.vertex_permutation(copy=False)
        fp = Q.face_permutation(copy=False)
        d = 4 * self._quadsystem._genus
        if self._first == None:
            self._first = e
            self._last = e
        else:
            pre = self._last
            newturn = self._quadsystem.turn(Q._ep(pre), e)
            
            if newturn == 0 and not self._turn_sequence: #spur
                self._first = self._last = None
            elif newturn == 0: #spur
                cur = self._last 
                (t, n) = self._turn_sequence.pop()
                if n > 1:
                    self._turn_sequence.append((t, n-1))
                for i in range(d - t):
                    cur = vp[cur]
                self._last = Q._ep(cur)
            elif len(self._turn_sequence) >= 1 and newturn == 1 and self._turn_sequence[-1][0] == 1: # positive bracket of length one
                self._turn_sequence.pop()
                turn_modif(self._turn_sequence, -1, d)
                self._last = Q.previous_in_face(Q._ep(e))
            elif len(self._turn_sequence) >= 2 and newturn == 1 and self._turn_sequence[-1][0] == 2 and self._turn_sequence[-2][0] == 1:
                # positive bracket of length more than one
                _, n = self._turn_sequence.pop()
                self._turn_sequence.pop()
                turn_modif(self._turn_sequence, -1, d)
                turn_add(self._turn_sequence, d - 2, n)
                self._last = Q.previous_in_face(Q._ep(e))
            elif len(self._turn_sequence) >= 1 and newturn == d - 1 and self._turn_sequence[-1][0] == d - 1: # negative bracket of length one
                self._turn_sequence.pop()
                turn_modif(self._turn_sequence, 1, d)
                self._last = Q._ep(fp[e])
            elif len(self._turn_sequence) >= 2 and newturn == d - 1 and self._turn_sequence[-1][0] == d - 2 and self._turn_sequence[-2][0] == d - 1:
                # negative bracket of length more than one
                _, n = self._turn_sequence.pop()
                self._turn_sequence.pop()
                turn_modif(self._turn_sequence, 1, d)
                turn_add(self._turn_sequence, 2, n)
                self._last = Q._ep(fp[e])
            else:
                self._last = e
                turn_add(self._turn_sequence, newturn, 1)


    def add_edge_left(self, e):
        r"""
        Add an edge at the beginning of self and perform simplification to maintain it geodesical.

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, LazyGeodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6],[7, 8, 5], [9, 10, 12, 11], [3, 15, 1, 13, 14]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,1,2,3,4,5,6,7)(~0,~5,~6,~7,~4,~1,~2,~3)", "(0,~3,2,~1)(~0,7,~6,5)(1,~4,3,~2)(4,~7,6,~5)")
            sage: p = LazyGeodesic(Q)
            sage: p.add_edge_left(12)
            sage: p.add_edge_left(5)
            sage: p.add_edge_left(6)
            sage: p.add_edge_left(11)
            sage: p
            LazyGeodesic starting with "11", ending with "12" and with turns "deque([(6, 1), (7, 1), (4, 1)])"
            sage: p.add_edge_left(12)
            sage: p
            LazyGeodesic starting with "14", ending with "12" and with turns "deque([(2, 1), (5, 1)])"
            sage: p.add_edge_left(15)
            sage: p
            LazyGeodesic starting with "3", ending with "12" and with turns "deque([(5, 1)])"
        """
        
        Q = self._quadsystem._quad
        vp = Q.vertex_permutation(copy=False)
        fp = Q.face_permutation(copy=False)
        d = 4 * self._quadsystem._genus
        if self._first == None:
            self._first = e
            self._last = e
        else:
            pre = self._first
            newturn = self._quadsystem.turn(Q._ep(e), pre)
            
            if newturn == 0 and not self._turn_sequence: #spur
                self._first = self._last = None
            elif newturn == 0: #spur
                cur = Q._ep(self._first)
                (t, n) = self._turn_sequence.popleft()
                if n > 1:
                    self._turn_sequence.appendleft()
                for i in range(t):
                    cur = vp[cur]
                self._first = cur
            elif len(self._turn_sequence) >= 1 and newturn == 1 and self._turn_sequence[0][0] == 1: # positive bracket of length one
                self._turn_sequence.popleft()
                turn_modif_left(self._turn_sequence, -1, d)
                self._first = Q.previous_at_vertex(e)
            elif len(self._turn_sequence) >= 2 and newturn == 1 and self._turn_sequence[0][0] == 2 and self._turn_sequence[1][0] == 1:
                # positive bracket of length more than one
                _, n = self._turn_sequence.popleft()
                self._turn_sequence.popleft()
                turn_modif_left(self._turn_sequence, -1, d)
                turn_add_left(self._turn_sequence, d - 2, n)
                self._first = Q.previous_at_vertex(e)
            elif len(self._turn_sequence) >= 1 and newturn == d - 1 and self._turn_sequence[0][0] == d - 1: # negative bracket of length one
                self._turn_sequence.popleft()
                turn_modif_left(self._turn_sequence, 1, d)
                self._first = vp[e]
            elif len(self._turn_sequence) >= 2 and newturn == d - 1 and self._turn_sequence[0][0] == d - 2 and self._turn_sequence[1][0] == d - 1:
                # negative bracket of length more than one
                _, n = self._turn_sequence.popleft()
                self._turn_sequence.popleft()
                turn_modif_left(self._turn_sequence, 1, d)
                turn_add_left(self._turn_sequence, 2, n)
                self._first = vp[e]
            else:
                self._first = e
                turn_add_left(self._turn_sequence, newturn, 1)


class StarShapedSpace:

    # TODO

    def __init__(self, Q, root):

        self._vertices = [root] #??
        self._quadsystem = Q # the underlying quadsystem
        self._inedges = [{}] # for each vertex a dictionnary containing the entering edges
        self._outedges = [[]] # for each vertex the list of (at most 2) outedges
        self._rotate_list = Q.rotate_list() 

    
    def add_vertex(self, vertex):
        
        self._vertices.append(vertex)
        self._inedges.append({})
        self._outedges.append([])
        return len(self._vertices) - 1

    
    def add_edge(self, vertex1, vertex2, edge):
        r"""
        Add an edge from vertex1 to vertex2 that project to edge
        """
        Q = self._quadsystem
        if not self.contains_edge(vertex1, edge):
            self._outedges[vertex1].append((edge, vertex2))
        if not self.contains_edge(vertex2, Q._quad._ep(edge)):
            self._inedges[vertex2][Q._quad._ep(edge)] = vertex1


    def contains_edge(self, vertex, edge):
        
        result = False
        for elt in self._outedges[vertex]:
            if elt[0] == edge:
                result = True
        return result or (not self._inedges[vertex].get(edge) is None)
            

    def opposite_vertex(self, vertex, edge):
        for elt in self._outedges[vertex]:
            if elt[0] == edge:
                return elt[1]
        for elt in self._inedges[vertex]:
            if elt == edge:
                return self._inedges[vertex][elt]
        raise ValueError("There is no such edge")

        
    def rotate(self, vertex, edge, turn):
        if self.contains_edge(vertex, edge):
            return self._rotate_list[edge][turn]
        else:
            raise ValueError("There is no such edge")

    def turn(self, vertex, edge):
        turn = None
        out_edge = None
        Q = self._quadsystem
        d = 4 * Q._genus
        for elt in self._outedges[vertex]:
            newturn = Q.turn(edge, elt[0])
            if newturn > d//2:
                newturn = newturn - d
            if turn == None or (abs(turn) > abs(newturn)):
                turn = newturn
                out_edge = elt[0]
        return turn, out_edge


    def insert_edge(self, vertex, edge):  
        if self.contains_edge(vertex, edge):
            return self.opposite_vertex(vertex, edge)
        turn, out_edge = self.turn(vertex, edge)
        opposite_edge = self._quadsystem._quad._ep(edge)
        new_vertex = self.add_vertex(opposite_edge%2)
        if turn != None and abs(turn) == 1:
            old_vertex = self.opposite_vertex(vertex, out_edge)
            new_edge = self._rotate_list[out_edge][-turn]
            self.insert_edge(old_vertex, new_edge)
            square_vertex = self.opposite_vertex(old_vertex, new_edge)
            square_edge = self.rotate(square_vertex, self._quadsystem._quad._ep(new_edge), -turn)
            self.add_edge(new_vertex, square_vertex, self._quadsystem._quad._ep(square_edge))
        self.add_edge(vertex, new_vertex, edge)
        return new_vertex
