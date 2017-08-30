#!/usr/bin/env python3


def ncart(am):
    return ((am+1)*(am+2))//2


def iterate_gaussian(lmn):
    am = lmn[0] + lmn[1] + lmn[2]
    if lmn[2] >= am:
        return None
    if lmn[2] < (am - lmn[0]):
        return (lmn[0], lmn[1]-1, lmn[2]+1)
    else:
        return (lmn[0]-1, am-lmn[0]+1, 0)


def all_cartesian_components(am):
    all_cart = []
    lmn = (am, 0, 0)

    while lmn:
        all_cart.append(lmn)
        lmn = iterate_gaussian(lmn)

    return all_cart
