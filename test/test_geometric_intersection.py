import pytest


def test_geometric_intersection():
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection

    torus = OrientedMap(fp="(0,1,~0,~1)")
    octagon = OrientedMap(fp="(0,1,2,3,~0,~1,~2,~3)")

    torus_gi = GeometricIntersection(torus)
    octagon_gi = GeometricIntersection(octagon)

    for (cmap, u, v, expected_intersection) in [
        (torus_gi, [0], [2], 1),
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
        (torus_gi, [0,2,2,0,2,2,2,2], None, 1),
        (torus_gi, [0,2,0,0,2,0,0,2,0,2], None, 1),
        (torus_gi, [0,2,0,0,2,0,2,0], None, 0),
        (torus_gi, [0,2,0,3], None, 1),
        (torus_gi, [1,3,1,2], None, 1),
        (torus_gi, [0,0,2,2,2], None, 2),
        (torus_gi, [0,0,2,2,2,2], None, 3),
        (torus_gi, [0,0,3,3,0,3], None, 2),
        (octagon_gi, [0,6], [6,0,2], 0)]:
        if v is None:
            computed_intersection = cmap.geometric_intersection([u], None)
        else:
            computed_intersection = cmap.geometric_intersection([u], [v])
        assert computed_intersection == expected_intersection, (cmap, u, v, expected_intersection, computed_intersection)


def test_torus_mcg():
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection

    torus = OrientedMap(fp="(0,1,~0,~1)")
    gi = GeometricIntersection(torus)

    f0 = [[2], [3], [0], [1]]
    f1 = [[1], [0], [2], [3]]
    f2 = [[0, 2], [3, 1], [2], [3]]
    f3 = [[0], [1], [2, 0], [1, 3]]
    f = [f0, f1, f2, f3]
    def apply_mcg(f, w):
        ww = []
        for i in w:
            ww.extend(f[i])
        return ww

    for w in [[0], [0, 0, 2, 2], [0, 0, 0, 2, 2], [0, 0, 2, 2, 0, 0, 2, 2], [0, 3, 1, 2], [0, 2, 0, 0, 3]]:
        intersection = gi.geometric_intersection([w])
        for s in [[0], [1], [2], [3], [0, 1], [0, 2], [1, 2], [2, 2, 3, 0], [3, 2, 1]]:
            ww = w[:]
            for i in s:
                ww = apply_mcg(f[i], ww)
            assert gi.geometric_intersection([ww]) == intersection


def test_geometric_intersection_multilinearity():
    import itertools
    from combisurf import OrientedMap
    from combisurf.geometric_intersection import GeometricIntersection

    torus = OrientedMap(fp="(0,1,~0,~1)")
    gi = GeometricIntersection(torus)

    # non-primitivity self-intersections
    for u in [[0, 0, 2], [0, 2, 1, 3], [0, 0, 2, 2], [0, 0, 2, 2, 1, 1, 3]]:
        intersection = gi.geometric_intersection([u])
        assert gi.geometric_intersection([u * 2]) == 4 * intersection + 1
        assert gi.geometric_intersection([u, u]) == 4 * intersection
        assert gi.geometric_intersection([u * 2, u]) == 9 * intersection + 1
        assert gi.geometric_intersection([u * 3]) == 9 * intersection + 2

    wlist = [[0], [1], [0, 1], [0, 0, 1], [0, 1, 1], [0, 0, 1, 1]]
    Q = [[gi.geometric_intersection([w0], [w1]) for w1 in wlist] for w0 in wlist]
    for ucoeffs in [[3,5,1,0,2,4], [1,3,0,1,2,1], [5,1,2,3,0,3]]:
        ulist = [w * mult for w, mult in zip(wlist, ucoeffs)]
        for vcoeffs in [[0,1,0,2,0,3], [1,1,1,1,0,2], [2,6,4,1,5,3]]:
            vlist = [w * mult for w, mult in zip(wlist, vcoeffs)]
            ans0 = sum(ucoeffs[i] * vcoeffs[j] * Q[i][j] for i in range(6) for j in range(6))
            ans1 = gi.geometric_intersection(ulist, vlist)
            assert ans0 == ans1, gi
