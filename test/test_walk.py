import pytest


def test_constructor():
    from combisurf import OrientedMap
    from combisurf import Walk

    m = OrientedMap(vp = [2, 12, 4, 5, 6, 11, 0, 8, 10, 13, 1, 18, 14, 20, 7, 16, 15, 3, 17, 9, 21, 19], mutable=False)

    Walk(m, [])
    Walk(m, [14, 16, 18])
    Walk(m, [0, 7, 6, 10])
    Walk(m, [18, 20, 21])
    Walk(m, [13, 7, 2, 18])

    with pytest.raises(ValueError, match="""edge 22 does not belong to the oriented map OrientedMap("(0,1,2,3)(~0,6,7,~3,4,5)(~1,~2,~5,9,~8)(~4,~6,10,~10,~9)(~7,8)", "(0,5,~2,1,~8,~7,6,~4,~3,2,~1)(~0,3,7,8,9,~10,~6)(4,~9,~5)(10)")."""):
        Walk(m, [15, 22])
    with pytest.raises(ValueError, match="edge 3 does not follow edge 0."):
        Walk(m, [0, 3])
    with pytest.raises(ValueError, match="edge 20 does not follow edge 10."):
        Walk(m, [11, 14, 15, 10, 20, 13, 7])



def test_is_homotopy():
    from combisurf import OrientedMap
    from combisurf import Walk

    m = OrientedMap(vp = [2, 12, 4, 5, 6, 11, 0, 8, 10, 13, 1, 18, 14, 20, 7, 16, 15, 3, 17, 9, 21, 19], mutable=False)

    w1 = Walk(m, [])
    w2 = Walk(m, [20])
    w3 = Walk(m, [1, 6, 14, 16, 18, 21, 13])
    w4 = Walk(m, [12, 13])
    
    w5 = Walk(m, [7, 2, 17, 15])
    w6 = Walk(m, [7, 0])
    
    w7 = Walk(m, [8, 19, 17, 15])
    w8 = Walk(m, [10, 17, 15])
    w9 = Walk(m, [8, 21, 19, 17, 15])

    assert w1.is_homotopic(w2)
    assert w1.is_homotopic(w3)
    assert w1.is_homotopic(w4)
    assert w3.is_homotopic(w4)

    assert not w4.is_homotopic(w6)
    assert w6.is_homotopic(w6)
    assert not w5.is_homotopic(w6)

    assert w7.is_homotopic(w8)
    assert w7.is_homotopic(w9)
    assert w9.is_homotopic(w8)
    assert w8.is_homotopic(w9)
    assert w8.is_homotopic(w8)


def test_simplicity():
    from combisurf import OrientedMap
    from combisurf import Walk

    m = OrientedMap(vp = [2, 12, 4, 5, 6, 11, 0, 8, 10, 13, 1, 18, 14, 20, 7, 16, 15, 3, 17, 9, 21, 19], mutable=False)

    w1 = Walk(m, [])
    w2 = Walk(m, [7, 6])
    w3 = Walk(m, [16])
    w4 = Walk(m, [18, 13, 14, 16, 11, 7, 2])
    w5 = Walk(m, [8, 21, 19, 17, 15])
    w6 = Walk(m, [2, 17, 15, 12, 9, 7])
    w7 = Walk(m, [2, 17, 15, 12, 9, 7, 4, 3, 0, 10, 5])
    w8 = Walk(m, [2, 17, 15, 12, 9, 7, 4, 3, 0, 10, 5, 0])
    w9 = Walk(m, [15, 7, 2, 17, 15, 12, 9, 7, 4, 3, 0, 10, 5, 0, 12, 19, 17])

    assert w1.simplicity()
    assert not w2.simplicity()
    assert w3.simplicity()
    assert w4.simplicity()
    assert not w5.simplicity()
    assert w6.simplicity()
    assert w7.simplicity()
    assert not w8.simplicity()
    assert not w9.simplicity()