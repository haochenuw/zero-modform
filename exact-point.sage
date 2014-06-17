### find the exact z ###

from itertools import *


def mindiff(lst):
    return min([abs(a-b) for a,b in combinations(lst,2)])


def taulist(Fj,N,prec):
    G = Gamma0(N)
    cosreps = G.coset_representatives()
    C = ComplexField(prec)
    jlist = Fj.roots()
    jroots = [a for a,b in jlist]
    taulist = []
    eps = mindiff(jroots)/20
    R.<x> = C[]
    for j,m in jlist:
        y = (256*(1-x)^2 - j*x^2).roots(multiplicities = False)[0]
        lam = (x^2-x-y).roots(multiplicities = False)[0]
        tau = hyp(tau,prec)
        found_it = False
        for g in cosreps:
            if found_it: break
            else:
                for a,b in jlist:
                    if abs(gp.ellj(C(N*g.acton(tau)))- ) < eps:
                        taulist.append(g.acton(tau),m)
                        found_it = True
                        break
    if len(taulist) != len(jlist):
        raise RuntimeError('failed to compute , find %s roots, expected %s'%(len(taulist),len(jlist)))
    return taulist