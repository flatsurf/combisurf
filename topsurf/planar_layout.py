"""
Layout of planar maps

We use a form a force-directed graph drawing taking care of the
underlying planar embedding.

Coordinates are handled as numpy array of floating points.

"""

import math
import numpy as np

from sage.rings.real_double import RDF

from topsurf.permutation import perm_orbit

def orientation(x, y, i, j, k, tol=1e-9):
    r"""
    Return the position of (x[k],y[k]) wrt to the line through (x[i],y[i]) and (x[j],y[j])

    EXAMPLES::

        sage: from topsurf.planar_layout import orientation

    Third point on the left::

        sage: assert orientation([0,1,0],[0,0,1],0,1,2) == 1
        sage: assert orientation([0,1,2],[0,0,1],0,1,2) == 1
        sage: assert orientation([0,1,-2],[0,0,1],0,1,2) == 1

    Third point on the right::

        sage: assert orientation([0,1,0],[0,0,-1],0,1,2) == -1
        sage: assert orientation([0,1,2],[0,0,-1],0,1,2) == -1
        sage: assert orientation([0,1,-2],[0,0,-1],0,1,2) == -1
    """
    val = (x[j] - x[i]) * (y[k] - y[j]) - (y[j] - y[i]) * (x[k] - x[j])
    if abs(val) < tol:
        return 0
    elif val > 0:
        return 1
    else:
       return -1

def is_on_segment(x, y, i, j, k, tol=1e-9):
    r"""
    Return whether (x[k], y[k]) is on the segment [(x[i],y[i]), (x[j],y[j])] assuming they are colinear.
    """
    # TODO: make a difference between
    # - big intersection
    # - small or almost intersection (possibly due to numerical noise)
    # - no intersection

        # small or almost intersection
    xm = min(x[i], x[j])
    xM = max(x[i], x[j])
    ym = min(y[i], y[j])
    yM = max(y[i], y[j])

    # no intersection
    if (x[k] < xm - tol or x[k] > xM + tol or
        y[k] < ym - tol or y[k] > yM + tol):
        return 0

    # big intersection
    if (xm + tol <= x[k] <= xM - tol or ym + tol <= y[k] <= yM - tol):
        return 2

    # small or almost intersection
    return 1


def are_segments_intersecting(x, y, i0, i1, j0, j1, tol=1e-9):
    r"""
    EXAMPLES::

        sage: from topsurf.planar_layout import are_segments_intersecting
        sage: are_segments_intersecting([-1,1,0,0],[0,0,-1,1],0,1,2,3)
        3

        sage: are_segments_intersecting([-1,1,0,0],[-2,-2,-1,1],0,1,2,3)
        0
        sage: are_segments_intersecting([-1,1,0,0],[2,2,-1,1],0,1,2,3)
        0
        sage: are_segments_intersecting([-1,1,-2,-2],[0,0,-1,1],0,1,2,3)
        0
        sage: are_segments_intersecting([-1,1,2,2],[0,0,-1,1],0,1,2,3)
        0

        sage: are_segments_intersecting([0,1,0,-1], [0,0,0,0], 0,1,2,3)
        1

        sage: are_segments_intersecting([0,2,1,2], [0,0,0,2], 0,1,2,3)
        2
        sage: are_segments_intersecting([0,2,1,-2], [0,0,0,0], 0,1,2,3)
        2
        sage: are_segments_intersecting([0,2,-2,1], [0,0,0,0], 0,1,2,3)
        2
    """
    # degenerate segments
    if -tol <= x[i0] - x[i1] <= tol and -tol <= y[i0] - y[i1] <= tol:
        return -1
    if -tol <= x[j0] - x[j1] <= tol and -tol <= y[j0] - y[j1] <= tol:
        return -1

    o1 = orientation(x, y, i0, i1, j0)
    o2 = orientation(x, y, i0, i1, j1)
    o3 = orientation(x, y, j0, j1, i0)
    o4 = orientation(x, y, j0, j1, i1)

    # stable intersection
    if o1 * o2 == -1 and o3 * o4 == -1:
        return 3

    # TODO: we could have one small and one big intersection that should result in one big!!!
    # degenerate intersection
    if not o1:
        a = is_on_segment(x, y, i0, i1, j0)
        if a:
            return a

    if not o2:
        a = is_on_segment(x, y, i0, i1, j1)
        if a:
            return a
    if not o3:
        a = is_on_segment(x, y, j0, j1, i0)
        if a:
            return a
    if not o4:
        a = is_on_segment(x, y, j0, j1, i1)
        if a:
            return a

    return 0


def angle(x, y, i, j, k, tol=1e-9):
    r"""
    Return the angle between (z[i], z[j]) and (z[i], z[k]) as a multiple of 2pi (ie a number between 0 and 1)

    EXAMPLES::

        sage: from topsurf.planar_layout import angle
        sage: angle([0, 1, 1], [0, 0, 1], 0, 1, 2)
        0.125
        sage: angle([0, 1, 0], [0, 0, 1], 0, 1, 2)
        0.25
        sage: angle([0, 1, -1], [0, 0, 1], 0, 1, 2)
        0.375
        sage: angle([0, 1, -1], [0, 0, 0], 0, 1, 2)
        0.5
        sage: angle([0, 1, -1], [0, 0, -1], 0, 1, 2)
        0.625
        sage: angle([0, 1, 0], [0, 0, -1], 0, 1, 2)
        0.75
        sage: angle([0, 1, 1], [0, 0, -1], 0, 1, 2)
        0.875
    """
    x1 = x[j] - x[i]
    y1 = y[j] - y[i]
    x2 = x[k] - x[i]
    y2 = y[k] - y[i]
    if ((-tol <= x1 <= tol and -tol <= y1 <= tol) or
        (-tol <= x2 <= tol and -tol <= y2 <= tol)):
        raise ValueError("degenerate vector")

    dot = x1 * x2 + y1 * y2
    det = x1 * y2 - y1 * x2
    angle = math.atan2(det, dot) / math.pi / 2
    return angle if angle >= 0 else angle + 1


class PlanarLayout:
    r"""
    Layout of planar maps.

    Vertices are indexed from 0 to n-1 and are split into
    - nodes (the node of the original graph)
    - subdivisions (all of which have degree 2)

    EXAMPLES::

        sage: from topsurf import OrientedMap
        sage: from topsurf.planar_layout import PlanarLayout
        sage: m = OrientedMap(vp="(0,~2)(~0,3,1)(~1,~5,2)(~3,4)(~4,5)")
        sage: pl = PlanarLayout(m)
        sage: pl.plot()    # default embedding from sage
        Graphics object consisting of ... graphics primitives
        sage: pl.refine()  # perform some refinement steps
        sage: pl.plot()    # nicer plot
        Graphics object consisting of ... graphics primitives
    """
    def __init__(self, m, root=None, check=True):
        if root is None:
            root = len(m._vp) - 2

        G, embedding, root_edge, edge_vertices = m.graph(root=root, subdivide=2)
        pos = G.layout_planar(on_embedding=embedding, external_face=root_edge)
        self._root = root
        self._cm = m
        self._half_edge_to_vertex = [-1] * (2 * len(m._vp))
        for i, v in enumerate(m.vertices()):
            for h in v:
                self._half_edge_to_vertex[h] = i
        self._graph = G
        self._edge_vertices = edge_vertices
        self._pos = np.zeros((2, len(pos)), dtype=float)
        self._num_nodes = m.num_vertices()
        self._num_subdivisions = len(pos) - self._num_nodes
        self._embedding = embedding
        # NOTE: sage embeddings is specified by clockwise ordering. We make it counter-clockiwse
        # below
        for neighbors in self._embedding.values():
            neighbors.reverse()
        for i, (xi, yi) in pos.items():
            self._pos[0][i] = xi
            self._pos[1][i] = yi

        if check and not self.is_valid():
            raise RuntimeError

        # force parameters
        self._repulsion_power = 1.0
        self._spring_length = 1.0
        self._spring_coupling = 2.0
        self._center_pull_coupling = 3.0
        self._stretch_force_power = 1.0

    def plot_force(self, force, *args, **kwds):
        from sage.plot.arrow import arrow2d
        from sage.plot.graphics import Graphics
        ans = Graphics()
        for v in range(self._graph.num_verts()):
            p = self._pos[:, v]
            f = force[:, v]
            ans += arrow2d(p, p+f, *args, **kwds)
        return ans

    def vertices(self):
        return self._graph.vertices()

    def vertices_on_half_edge(self, h, include_start=True, include_end=False):
        r"""
        EXAMPLES::

            sage: from topsurf import OrientedMap
            sage: from topsurf import OrientedMap
            sage: from topsurf.planar_layout import PlanarLayout
            sage: m = OrientedMap(vp="(0,~2)(~0,3,1)(~1,~5,2)(~3,4)(~4,5)")
            sage: pl = PlanarLayout(m)
            sage: pl.vertices_on_half_edge(0)
            [0, 5, 6]
            sage: pl.vertices_on_half_edge(1)
            [1, 6, 5]
        """
        vertices = self._edge_vertices[h // 2]
        if h % 2:
            vertices = vertices[::-1]
        if include_start and include_end:
            return vertices
        elif include_start:
            return vertices[:-1]
        elif include_end:
            return vertices[1:]
        else:
            return vertices[1:-1]

    def internal_faces(self):
        r"""
        Return the internal faces of the map as a list of underlying SageMath graph vertices.
        """
        ans = []
        for f in self._cm.faces():
            if self._root in f:
                continue
            ans.append([v for h in f for v in self.vertices_on_half_edge(h, include_start=True, include_end=False)])
        return ans

    def edges(self):
        ans = []
        for e, verts in enumerate(self._edge_vertices):
            for i in range(len(verts) - 1):
                ans.append(verts[i:i + 2])
        return ans

    def refine(self, num_steps=20):
        step_size = 1.0
        pos = self._pos
        if not self.is_valid():
            raise RuntimeError
        for step in range(num_steps):
            print(f"step={step}")
            spring = self.spring_force()
            stretch = self.stretch_force()
            repulsion = self.repulsion_force()

            spring /= np.linalg.norm(spring)
            stretch /= np.linalg.norm(stretch)
            repulsion /= np.linalg.norm(repulsion)
            f = spring + stretch + repulsion

            while not self.is_valid(self._pos + f):
                f /= 2

            self._pos += f

    def is_valid(self, pos=None, tol=1e-9, verbose=False):
        r"""
        Test whether the embedding does make sense, that is

        - no edge intersection

        .. TODO::

            This function does not check that the embedding is the one given by
            the underlying oriented map.
        """
        if pos is None:
            pos = self._pos

        for u in range(self._pos.shape[1]):
            for v in range(u):
                if (-tol <= self._pos[0,v] - self._pos[0,u] <= tol and
                    -tol <= self._pos[1,v] - self._pos[1,u] <= tol):
                    if verbose:
                        print(f"vertices u={u} and v={v} are too close")
                    return False

        # check vertex ordering
        for u in self.vertices(internal=True, subdivisions=False):
            edges = self._embedding[u]
            a = sum(angle(self._pos[0], self._pos[1], u, edges[i], edges[(i+1) % len(edges)]) for i in range(len(edges)))
            a_round = round(a)
            if abs(a_round - a) > tol:
                raise RuntimeError
            if a_round != 1.0:
                if verbose:
                    print(f"wrong embedding at u={u}, got total angle={a}")
                return False

        # check face embedding
        for face in self.internal_faces():
            for i1 in range(len(face)):
                u1 = face[i1]
                v1 = face[(i1 + 1) % len(face)]
                for i2 in range(i1):
                    u2 = face[i2]
                    v2 = face[(i2 + 1) % len(face)]
                    ans = are_segments_intersecting(pos[0], pos[1], u1, v1, u2, v2)
                    if ans == -1:
                        raise RuntimeError
                    elif ans >= 2:
                        if verbose:
                            print(f"big intersection between ({u1} ({pos[:,u1]}), {v1} ({pos[:,v1]})) and ({u2} ({pos[:,u2]}) ,{v2} ({pos[:,v2]}))")
                        return False
                    if u1 == u2 or u1 == v2 or v1 == u2 or v1 == v2:
                        if ans != 1:
                            if verbose:
                                print(f"adjacent edges ({u1} ({pos[:,u1]}), {v1} ({pos[:,v1]})) and ({u2} ({pos[:,u2]}) ,{v2} ({pos[:,v2]}))")
                            return False
                    else:
                        if ans != 0:
                            if verbose:
                                print(f"disjoint edges ({u1} ({pos[:,u1]}), {v1} ({pos[:,v1]})) and ({u2} ({pos[:,u2]}) ,{v2} ({pos[:,v2]}))")
                            return False
        return True

    def plot(self):
        r"""
        Make a plot of the planar version of the graph.
        """
        from sage.plot.graphics import Graphics
        from sage.plot.point import point2d
        from sage.plot.line import line2d
        from sage.plot.text import text
        G = Graphics()
        G += point2d([(self._pos[0][i], self._pos[1][i]) for i in range(self._num_nodes)],
                     color="red",
                     pointsize=30, zorder=10)
        G += point2d([(self._pos[0][self._num_nodes+i], self._pos[1][self._num_nodes+i]) for i in range(self._num_subdivisions)],
                     color="blue",
                     pointsize=10, zorder=10)
        for (u, v, label) in self._graph.edges():
            start = (self._pos[0][u], self._pos[1][u])
            end = (self._pos[0][v], self._pos[1][v])
            G += line2d([start, end], color="black", zorder=0)

        for h in range(len(self._cm._vp)):
            u, v = self.vertices_on_half_edge(h)[:2]
            G += text(f"{h}", (self._pos[:, u] + self._pos[:, v]) / 2, color="black")

        return G

    # NOTE: forces
    # each force is a 2d vector for each vertex of the map
    def repulsion_force(self):
        r"""
        Return the repulsion force as a numpy array.

        In each face we consider a repulsion force which for each pair of edges in the face

        The repulsion force is a sum 
        """
        # NOTE: in force.js this corresponds to repulsionForce which then calls
        # repulsionForceFace, then repulsionForceEdgeEdge, then repulsionForceNodeLink
        # It involves the springLength and springCoupling parameters
        # The central repulsionForceNodeLink(n,l,strength,calcForce) is
        #
        #     var scale = repulsionPower * energy / distdiff;
        #     var lton = [ n.pos.minus(l[0].pos).normalize(), n.pos.minus(l[1].pos).normalize() ];
        #     var l0tol1 = l[1].pos.minus(l[0].pos).normalize();
        #     n.force.addVec( lton[0].plus(lton[1]).mult(scale) );
        #     l[0].force.subVec( lton[0].minus(l0tol1).mult(scale) );
        #     l[1].force.subVec( lton[1].plus(l0tol1).mult(scale) );
        #
        force = np.zeros((2, self._graph.num_verts()), dtype=float)

        for f in self.internal_faces():
            for i in range(len(f)):
                ii = (i + 1) % len(f)
                u0 = f[i]
                p0 = self._pos[:, u0]
                u1 = f[ii]
                p1 = self._pos[:, u1]
                assert u0 != u1
                assert (p0 != p1).any()
                u0_to_u1 = p1 - p0
                u0_to_u1 /= np.linalg.norm(u0_to_u1)

                for j in range(len(f)):
                    if j == i or j == ii:
                        continue
                    v = f[j]
                    if v == u0 or v == u1:
                        # NOTE: the same face could be adjacent several times to the same vertex
                        continue

                    q = self._pos[:, v]
                    assert (q != p0).any()
                    assert (q != p1).any()

                    u0_to_v = q - p0
                    u0_to_v /= np.linalg.norm(u0_to_v)
                    u1_to_v = q - p1
                    u1_to_v /= np.linalg.norm(u1_to_v)

                    distdiff = np.linalg.norm(p0 - q) + np.linalg.norm(p1 - q) - np.linalg.norm(p1 - p0)
                    scale = self._repulsion_power / distdiff

                    force[:, v] += scale * (u0_to_v + u1_to_v)

                    # NOTE: we ignore the u0,u1 force here!
                    # force[u0] -= scale * 
                    # force[u1] -=

        return force

    def spring_force(self):
        r"""
        Return the spring force as a numpy array.

        The spring force tries to make each edge have unit length.

        .. TODO::

            it seems reasonable to have spring length depending on the edge
            position. On the boundary, edges could be expected to have length
            1. But in the k-th layer they should probably be much smaller
        """
        # NOTE: in force.js this is computed independently for each subdivided edge
        # NOTE: in force.js this corresponds to springForce which then calls
        # springforceEdge
        force = np.zeros((2, self._graph.num_verts()), dtype=float)
        for u, v in self.edges():
            u_to_v = self._pos[:,v] - self._pos[:,u]
            length = np.linalg.norm(u_to_v) - self._spring_length
            fv = self._spring_coupling * length * u_to_v
            force[:, u] += u_to_v
            force[:, v] -= u_to_v
        return force

    def center_pull_force(self):
        # NOTE: in force.js this corresponds to centerPullForce which then calls
        # centerPullForceVertex
        force = np.zeros((2, self._graph.num_verts()), dtype=float)
        for v in self.vertices():
            radial = self._center_pull_coupling * self._pos[:, v]
            force[:, v] -= radial
        return force

    def stretch_force(self):
        r"""
        For each corner, try to move the two ends so that the angle fits to the target.

        For vertices only adjacent to the external face: try to make all angles
        equal For vertices alternating between internal faces and external
        face: try to make the sum of internal angles and the external angle
        proportional
        """
        # NOTE: in force.js this corresponds to stretchForce
        force = np.zeros((2, self._graph.num_verts()), dtype=float)

        bdry_half_edges = [0] * len(self._cm._vp)
        for h in perm_orbit(self._cm._fp, self._root):
            bdry_half_edges[h] = 1

        # compute target angles as multiple of 2 pi
        # TODO: should be stored once and for all
        target_angles = [.5] * len(self._cm._vp)

        # vertices of the original map
        for v in self._cm.vertices():
            num_external = sum(bdry_half_edges[h] for h in v)
            num_internal = len(v) - num_external
            if num_external == 0 or num_internal == 0:
                for h in v:
                    target_angles[h] = 1. / len(v)
            else:
                # pick half of the angle for external and half for internal
                for h in v:
                    if bdry_half_edges[h]:
                        target_angles[h] = .5 / num_external
                    else:
                        target_angles[h] = .5 / num_internal

            assert abs(sum(target_angles[h] for h in v) - 1.) < 0.000001

        for h in range(len(self._cm._vp)):
            if target_angles[h] == 1.:
                # ignore degree one vertex
                continue

            # make corners close to their target
            hh = self._cm.next_at_vertex(h)
            h_verts = self.vertices_on_half_edge(h, include_start=True, include_end=True)
            hh_verts = self.vertices_on_half_edge(hh, include_start=True, include_end=True)

            assert h_verts[0] == hh_verts[0]

            v = h_verts[0]
            v_next = h_verts[1]
            v_prev = hh_verts[1]

            left = self._pos[:, v_prev] - self._pos[:, v]
            left /= np.linalg.norm(left)
            right = self._pos[:, v_next] - self._pos[:, v]
            right /= np.linalg.norm(right)
            angle = np.angle((left[0] + 1j * left[1]) / (right[0] + 1j * right[1])) / 2. / np.pi
            if angle < 0:
                angle += 1
            assert 0 < angle < 1

            left_perp = np.array([left[1], -left[0]])
            right_perp = np.array([-right[1], right[0]])
            force[:, v_prev] += -(target_angles[h] - angle) * left_perp
            force[:, v_next] += -(target_angles[h] - angle) * right_perp

            # make subdivision corner close to pi
            for i in range(1, len(h_verts) - 1):
                v_prev = h_verts[i - 1]
                v = h_verts[i]
                v_next = h_verts[i + 1]
                left = self._pos[:, v_prev] - self._pos[:, v]
                left /= np.linalg.norm(left)
                right = self._pos[:, v_next] - self._pos[:, v]
                right /= np.linalg.norm(right)
                angle = np.angle((left[0] + 1j * left[1]) / (right[0] + 1j * right[1])) / 2. / np.pi
                if angle < 0:
                    angle += 1
                assert 0 < angle < 1
                energy = (angle - .5) / angle / (1 - angle)


                left_perp = np.array([left[1], -left[0]])
                right_perp = np.array([-right[1], right[0]])
                force[:, v_prev] += -(target_angles[h] - angle) * left_perp
                force[:, v_next] += -(target_angles[h] - angle) * right_perp

        return force

