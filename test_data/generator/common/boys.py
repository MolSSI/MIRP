#!/usr/bin/env python3

from mpmath import mp


def calculate_boys(m, t):
    zero_mp = mp.mpf(0)
    one_mp = mp.mpf(1)
    two_mp = mp.mpf(2)
    m_mp = mp.mpf(m)
    t_mp = mp.mpf(t)

    if t_mp == zero_mp:
        F_mp = one_mp/(two_mp * m_mp + one_mp)
    else:
        M = m_mp + (one_mp / two_mp)
        F_mp = mp.gammainc(M, zero_mp, t_mp) * one_mp/(two_mp * mp.power(t_mp, M))
    return (m, t_mp, F_mp)


def calculate_boys_list(m_list, t_list):
    ret = []
    for t in t_list:
        for m in m_list:
            ret.append(calculate_boys(m, t))
    return ret


def write_file(f, boys_data):
    pass
