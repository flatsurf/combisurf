import pytest


def test_geometric_intersection():
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection

    torus = OrientedMap(fp="(0,1,~0,~1)")
    octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")

    torus_gi = GeometricIntersection(torus)
    octagon_gi = GeometricIntersection(octagon)

    for (cmap, u, v, expected_intersection) in [
        (torus_gi, [0], [0,2], 1),
        (torus_gi, [0,2], [0,2,0], 1),
        (torus_gi, [0,2,0], [0,2,0,0,2], 1),
        (torus_gi, [0,3,1,2], [0,2,0,0,2], 0),
        (torus_gi, [0,0,2,2], [1], 2),
        (torus_gi, [0,2,0,0,3], [0], 0),
        (torus_gi, [0,2,0,0,3], [2], 3),
        (torus_gi, [0], None, 0),
        (torus_gi, [0,2,0,0,2], None, 0),
        (torus_gi, [0,0,2,2], None, 1),
        (torus_gi, [0,2,2,0], None, 1),
        (torus_gi, [0,2,0,2,0,0], None, 1),
        (torus_gi, [0,2,0,0,2,0,0,2,0,2], None, 1),
        (torus_gi, [0,2,0,0,2,0,2,0], None, 0),
        (torus_gi, [0,2,0,3], None, 1),
        (torus_gi, [1,3,1,2], None, 1),
        (torus_gi, [0,0,2,2,2], None, 2),
        (torus_gi, [0,0,2,2,2,2], None, 3),
        (torus_gi, [0,0,3,3,0,3], None, 2),
        (octagon_gi, [0,6], [6,0,2], 0)]:
        computed_intersection = cmap.primitive_curve_geometric_intersection(u, v)
        assert computed_intersection == expected_intersection, (cmap, u, v, expected_intersection, computed_intersection)
