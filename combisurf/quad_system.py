r"""
Homotopy test on maps

This module introduced two classes: :class:`QuadSystem` and :class:`Geodesic`.
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
    r"""
    Modify the last turn of s by mod.
    """
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
    r"""
    Modify the first turn of s by mod.
    """
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
    Test if ``u`` is a subword of ``v`` in ``O(|u|+|v|)``.
    """
    if len(u) == 0:
        return True
    elif len(v) == 0:
        return False
    cnd = 0
    T = [-1]
    for i in range(1, len(u)):
        if u[i] == u[cnd]:
            T.append(T[cnd])
        else:
            T.append(cnd)
            while cnd>=0 and u[i] != u[cnd]:
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

    r"""
    A quadsystem associated to an oriented map.

    A quadsystem is a reduced form of a map with two vertices and such that all faces have degree 4. The quadsystem is given as the corresponding map and the projection from the original map to the quadsystem.

    """

    def __init__(self, m, forest=None, coforest=None, mutable=False, check=True):
        r"""
        INPUT:

        - ``m`` -- the original oriented map

        - ``treecotree`` -- ``None`` or the tree-cotree decomposition used to construct the quadsystem

        - ``check`` (boolean, default ``True``) -- whether to perform consistency checks of
          the data
        """

        if check:
            m._check()
            m._assert_connected()
            if m.has_folded_edge():
                raise NotImplementedError
            if m.genus() == 0:
                raise ValueError("Cannot look for homotopy in genus 0")
            elif m.genus() == 1:
                raise ValueError("Quad systems and homotopy are not effective in genus 1")
            # TODO: we should check that forest and coforest are both forest and coforest if provided.
        
        if m.is_mutable():
            m = m.copy(mutable=False)
        
        self._origin_map = m
        self._genus = m.genus()

        quad, proj = m.quad_system(forest=forest, coforest=coforest, relabel=True, mapping=True, mutable=mutable, check=check)

        self._quad = quad
        self._proj = proj

        self._turn = list(self._quad.half_edges())
        label = 0
        for v in self._quad.vertices(): # we label each edge with a integer such that the difference between label[e] and label[e'] is the number of turns between e and e'.
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
        r"""
        Compute the number of turns from ``h1`` to ``h2`` around a vertex. Raise an error if ``h1`` and ``h2``do not belong to the same vertex. 
        """
        l1 = self._turn[h1]
        l2 = self._turn[h2]
        if l1 // (4 * self._genus) != l2 // (4 * self._genus):
            raise ValueError("The two half_edges do not belong to the same vertex")
        nl1 = l1 % (4 * self._genus)
        nl2 = l2 % (4 * self._genus)
        if nl1 <= nl2:
            return nl2 - nl1
        else:
            return 4 * self._genus + nl2 - nl1


    def rotate_list(self):
        r"""
        Compute the rotate_list of the quadsystem, i.e. the list ``rotate`` such that ``rotate[e][t]`` is the half-edge obtained in the quadsystem by turning ``t`` time around a vertex starting from the half-edge ``e``.
        """
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
    r"""
    A geodesic in a quadsystem

    A geodesic is a shortest walk in its free homotopy class. A geodesic is seen as a walk encoded by a deque with a given turn sequence.
    """

    def __init__(self, Q, geo=None, turn=None, check=False):
        r"""
        INPUT:
            
            - ``Q`` -- the underlying quadsystem
            - ``geo`` -- ``None`` or the walk corresponding to the geodesic
            - ``turn`` -- ``None`` or the turn sequence corresponding to ``geo``
            - ``check`` (boolean, default ``False``) -- whether to perform consistency checks of
          the data

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,4,5,3,1,2,6,7)(~0,~3,~6,~4,~1,~7,~5,~2)", "(0,~2,1,~4)(~0,7,~1,3)(2,~5,4,~6)(~3,5,~7,6)")
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
        Return the length of the geodesic.
        
        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,4,5,3,1,2,6,7)(~0,~3,~6,~4,~1,~7,~5,~2)", "(0,~2,1,~4)(~0,7,~1,3)(2,~5,4,~6)(~3,5,~7,6)")
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
            OrientedMap("(0,4,5,3,1,2,6,7)(~0,~3,~6,~4,~1,~7,~5,~2)", "(0,~2,1,~4)(~0,7,~1,3)(2,~5,4,~6)(~3,5,~7,6)")
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
            OrientedMap("(0,4,5,3,1,2,6,7)(~0,~3,~6,~4,~1,~7,~5,~2)", "(0,~2,1,~4)(~0,7,~1,3)(2,~5,4,~6)(~3,5,~7,6)")
            sage: p = Geodesic(Q)
            sage: p.add_edge(12)
            sage: p.add_edge(9)
            sage: p.add_edge(6)
            sage: p.add_edge(9)
            sage: p
            Geodesic "deque([12, 9, 6, 9])" with turns "deque([(1, 1), (2, 2)])"
            sage: p.add_edge(10)
            sage: p
            Geodesic "deque([4, 15, 4])" with turns "deque([(6, 2)])"
            sage: p.add_edge(5)
            sage: p
            Geodesic "deque([4, 15])" with turns "deque([(6, 1)])"
            sage: p.add_edge(0)
            sage: p.add_edge(7)
            sage: p
            Geodesic "deque([4, 3])" with turns "deque([(5, 1)])"
            sage: p.add_edge(6)
            sage: p.add_edge(1)
            sage: p
            Geodesic "deque([4, 15])" with turns "deque([(6, 1)])"
            sage: p.add_edge(12)
            sage: p.add_edge(1)
            sage: p.add_edge(14)
            sage: p
            Geodesic "deque([4, 11, 2])" with turns "deque([(7, 1), (2, 1)])"
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
            OrientedMap("(0,4,5,3,1,2,6,7)(~0,~3,~6,~4,~1,~7,~5,~2)", "(0,~2,1,~4)(~0,7,~1,3)(2,~5,4,~6)(~3,5,~7,6)")
            sage: p = Geodesic(Q)
            sage: p.add_edge_left(10)
            sage: p.add_edge_left(9)
            sage: p.add_edge_left(6)
            sage: p.add_edge_left(9)
            sage: p
            Geodesic "deque([9, 6, 9, 10])" with turns "deque([(2, 2), (1, 1)])"
            sage: p.add_edge_left(12)
            sage: p
            Geodesic "deque([4, 15, 4])" with turns "deque([(6, 2)])"
            sage: p.add_edge_left(5)
            sage: p
            Geodesic "deque([15, 4])" with turns "deque([(6, 1)])"
            sage: p.add_edge_left(2)
            sage: p.add_edge_left(7)
            sage: p
            Geodesic "deque([1, 4])" with turns "deque([(5, 1)])"
            sage: p.add_edge_left(6)
            sage: p.add_edge_left(3)
            sage: p
            Geodesic "deque([15, 4])" with turns "deque([(6, 1)])"
            sage: p.add_edge_left(10)
            sage: p.add_edge_left(3)
            sage: p.add_edge_left(14)
            sage: p
            Geodesic "deque([0, 13, 4])" with turns "deque([(2, 1), (7, 1)])"
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
            OrientedMap("(0,4,5,3,1,2,6,7)(~0,~3,~6,~4,~1,~7,~5,~2)", "(0,~2,1,~4)(~0,7,~1,3)(2,~5,4,~6)(~3,5,~7,6)")

            Remove a spur at origin.
            sage: p = Geodesic(Q, [5, 12, 1, 4])
            sage: p.origin_simplification()
            sage: p
            Geodesic "deque([12, 1])" with turns "deque([(6, 1)])"

            TODO: Make new examples and tests.
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

        elif first_turn == 1 and len(s) >= 1 and s[0][0] == 1: # bracket starting at the origin
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

        elif first_turn == 2 and (len(s) >= 2 or (len(s) == 1 and s[0][0] != 2)): # bracket containing the origin
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
            if first_turn == 1:
                self.origin_simplification() # the bracket now starts at the origin so we can recursively call origin_simplification

        elif first_turn == d - 2 and (len(s) >= 2 or (len(s) == 1 and s[0][0] != d - 2)): # bracket containing the origin
            while first_turn == d - 2:
                geo.append(geo.popleft())
                (t, n) = s.popleft()
                if n != 1:
                    s.appendleft((t, n - 1))
                (t2, n2) = s.pop()
                if t2==d-2:
                    s.append((d-2, n2+1))
                else:
                    s.append((t2, n2))
                    s.append((d-2, 1))
                first_turn = t
            if first_turn == d - 1:
                self.origin_simplification() # the bracket now starts at the origin so we can recursively call origin_simplification

    def canonical(self):
        r"""
        Perform the right push to the cyclic geodesic path to make it canonical.

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, Geodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,4,5,3,1,2,6,7)(~0,~3,~6,~4,~1,~7,~5,~2)", "(0,~2,1,~4)(~0,7,~1,3)(2,~5,4,~6)(~3,5,~7,6)")
            
            sage: p = Geodesic(Q, [7, 4])
            sage: p.canonical()
            sage: p
            Geodesic "deque([14, 9])" with turns "deque([(6, 1)])"
            
            TODO: Make new examples and tests.
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
                            e1 = Q._ep(e)
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




class LazyGeodesic:
    r"""
    A geodesic in a quadsystem

    A lazy geodesic is a geodesic but only the first and the last edge of the geodesic is maintained together with the turn sequence. It cannot be used to test homotopy but it can be used for contractibility.
    """
    
    def __init__(self, Q, geo=None, turn=None, check=False):
        r"""
        INPUT:
            
            - ``Q`` -- the underlying quadsystem
            - ``geo`` -- ``None`` or the walk corresponding to the geodesic
            - ``turn`` -- ``None`` or the turn sequence corresponding to ``geo``
            - ``check`` (boolean, default ``False``) -- whether to perform consistency checks of
          the data

        EXAMPLES::

            sage: from combisurf import OrientedMap, QuadSystem, LazyGeodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,4,5,3,1,2,6,7)(~0,~3,~6,~4,~1,~7,~5,~2)", "(0,~2,1,~4)(~0,7,~1,3)(2,~5,4,~6)(~3,5,~7,6)")
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
            self._last = geo[-1]
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

            sage: from combisurf import OrientedMap, QuadSystem, LazyGeodesic
            sage: m = OrientedMap(vp=[[0, 2, 4, 6], [5, 8, 10, 12], [3, 11, 13, 7, 1, 9]])
            sage: Q = QuadSystem(m)
            sage: Q
            OrientedMap("(0,4,5,3,1,2,6,7)(~0,~3,~6,~4,~1,~7,~5,~2)", "(0,~2,1,~4)(~0,7,~1,3)(2,~5,4,~6)(~3,5,~7,6)")
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
            OrientedMap("(0,3,4,5,1,6,7,2)(~0,~6,~4,~5,~1,~3,~7,~2)", "(0,~2,7,~3)(~0,2,~7,6)(1,~5,4,~6)(~1,5,~4,3)")
            
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
            OrientedMap("(0,3,4,5,1,6,7,2)(~0,~6,~4,~5,~1,~3,~7,~2)", "(0,~2,7,~3)(~0,2,~7,6)(1,~5,4,~6)(~1,5,~4,3)")
            
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
    r"""
    A star shaped space in the universal covering of a quadsystem.

    A star shaped space is a pointed convex set in the universal covering of a quadsystem. It thus contains all geodesics from a point in it to the root. Thus all edges are oriented towards the root.

    """

    def __init__(self, Q, root):
        r"""
        INPUT:

            - ``Q`` -- the quadsystem
            - ``root`` -- the vertex in ``Q`` whose lift is the root of the star shaped space

        """

        self._vertices = [root] # the projection of each vertex of the StarShapedSpace into the quadsystem
        self._quadsystem = Q # the underlying quadsystem
        self._inedges = [{}] # for each vertex a dictionary containing the entering edges
        self._outedges = [[]] # for each vertex the list of (at most 2) outedges
        self._rotate_list = Q.rotate_list() 

    
    def add_vertex(self, vertex):
        r"""
        Add a vertex to the star shaped space with no edges.
        """
        
        self._vertices.append(vertex)
        self._inedges.append({})
        self._outedges.append([])
        return len(self._vertices) - 1

    
    def add_edge(self, vertex1, vertex2, edge):
        r"""
        Add an edge from vertex1 to vertex2 that project to ``edge``
        """
        Q = self._quadsystem
        if not self.contains_edge(vertex1, edge):
            self._outedges[vertex1].append((edge, vertex2))
        if not self.contains_edge(vertex2, Q._quad._ep(edge)):
            self._inedges[vertex2][Q._quad._ep(edge)] = vertex1


    def contains_edge(self, vertex, edge):
        r"""
        Test whether the star shaped space contains the edge ``edge`` at vertex ``vertex``.
        """
        
        result = False
        for elt in self._outedges[vertex]:
            if elt[0] == edge:
                result = True
        return result or (not self._inedges[vertex].get(edge) is None)
            

    def opposite_vertex(self, vertex, edge):
        r"""
        Return the vertex opposite to ``edge`` starting at ``vertex``.
        """
        for elt in self._outedges[vertex]:
            if elt[0] == edge:
                return elt[1]
        for elt in self._inedges[vertex]:
            if elt == edge:
                return self._inedges[vertex][elt]
        raise ValueError("There is no such edge")

        
    def rotate(self, vertex, edge, turn):
        r"""
        Rotate around ``vertex`` starting from ``edge`` and starting from ``edge``
        """
        if self.contains_edge(vertex, edge):
            return self._rotate_list[edge][turn]
        else:
            raise ValueError("There is no such edge")

    def turn(self, vertex, edge):
        r"""
        Return the smallest turn from ``edge`` to an outedge of ``vertex``
        """
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
        r"""
        Insert an edge in the star shaped space and add vertices and edges to maintain it star shaped.
        """
        
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
