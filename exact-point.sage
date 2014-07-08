### find the exact z ###
### dependence: depend on exceptional-points.sage to compute the exceptional points.
### implicitly, depend on zero-modform.sage to provide the zero polynomial Fj.
from itertools import *
from mpmath import *
import numpy as np
load('exceptional-points.sage')

def twodargmin(lst1,lst2,used2):
    """
    (i,j) such that abs(lst[i]-lst2[j]) is the minimum among
    (i,j) such that j is not in used2
    """

    v = [min([abs(a-lst2[j]) for j in range(len(lst2)) if j not in used2]) for a in lst1]
    i1 = np.argmin(v)
    #print v,i1
    x = lst1[i1]
    t = [j for j in range(len(lst2)) if j not in used2 and abs(x-lst2[j]) == min(v)]
    assert len(t) > 0
    return (i1,t[0])


def mindiff(lst):
    if len(lst) < 2:
        # void
        return 1
    return min([abs(a-b) for a,b in combinations(lst,2)])

def hyp(z,prec):
    mp.prec = prec
    C = ComplexField(prec)
    tau = hyp2f1(1/2,1/2,1,1-z)/hyp2f1(1/2,1/2,1,z)
    return C(I)*(C(tau.real)+C(tau.imag)*I)


def taus(jlist,N,prec):
    """
    return taulist such that j( taulist[i] ) = jlist[i]
    """
    t = cputime()
    G = Gamma0(N)
    cosreps = list(G.coset_reps())
    C = ComplexField(prec)
    R.<x> = C[]
    #jlist = R(Fj).roots(multiplicities = False)
    taulist = []
    eps = 1
    verbose('eps = %s '%eps)
    t1 = cputime()
    for j in jlist:
        y = R(256*(1-x)^3 - j*x^2).roots(multiplicities = False)[0]
        lam = R(x-x^2-y).roots(multiplicities = False)[0]
        tau = hyp(lam,prec)
        taulist.append(tau)
    verbose('2nd step done, found raw taus, took %s seconds'%cputime(t1))
    zlist = []
    count = 0
    num = len(taulist)
    return taulist
    """
    # use field smaller precision since gp.ellj is slow
    C1 = CDF
    for tau,m in taulist:
        count = count +1
        verbose('adjusting the %s th point (total # of points = %s)'%(count,num))
        found = False
        t2 = cputime()
        for g in cosreps:
            if found: break
            else:
                for a in jlist:
                    #verbose('trying the combination with g = %s, a = %s'%(g,ComplexField(10)(a)))
                    if abs(C1(gp.ellj(C1(N*g.acton(tau)))) - C1(a)) < eps:
                        zlist.append((C(g.acton(tau)),m))
                        found = True
                        break
        if not found:
            zlist.append((tau,'failed'))
        verbose('done(successful = %s )adjusting for the %s th point, took %s seconds'%(found, count, cputime(t2)))
    #if len(zlist) != R(Fj).degree():
    #    raise RuntimeError('failed to compute , find %s points, expected %s'%(len(zlist),R(Fj).degree()))
    verbose('The whole computation took %s seconds'%cputime(t))
    return zlist
    """


#to-do: factor the above code, since it's not correct


def from_tau_to_z(taulist,jlist,N,eps,prec):
    verbose('eps = %s'%eps)
    C = ComplexField(prec)
    pari.set_real_precision(prec)
    G = Gamma0(N)
    cosreps = list(G.coset_reps())
    n = len(taulist)
    verbose('n = %s'%n)
    used = []
    zlist = []
    for i in range(n):
        tau = taulist[i]
        verbose('the %s th tau is equal to %s'%(i+1,CDF(tau)))
        cands = []
        count = 0
        t = cputime()
        for g in cosreps:
            cands.append(C(pari(C(N*g.acton(tau))).ellj()))
            #verbose('computed a j-invariant')
        # cands = [C(pari(C(N*g.acton(tau))).ellj()) for g in cosreps] # apparently slow
        verbose('the list of j(Ngtau) for the %s th tau is computed with %s seconds'%(i+1, cputime(t)))
        j,k = twodargmin(cands,jlist,used) # find the indexes that minimizes the difference between gtau and j
        verbose('the index of the argmin: %s, %s'%(j,k))
        mindist = abs(cands[j]- jlist[k])
        verbose('the minimal difference between the two list element is %s'%mindist)
        if mindist < eps:
            used.append(k)
            g = cosreps[j]
            zlist.append(C(g.acton(tau)))
            verbose('the indexes of used j: %s'%str(used))
        else:
            raise RuntimeError('not able to identify the %s th point '%(i+1))
    return zlist


def from_tau_to_z_onepoint(tau,jlist,N,eps,used,prec):
    C = ComplexField(prec)
    verbose('tau is equal to %s'%CDF(tau))
    cands = [C(gp.ellj(C(N*g.acton(tau)))) for g in cosreps]
    verbose('the list of j(Ngtau) for the %s th tau is computed '%(i+1))
    j,k = twodargmin(cands,jlist,used) # find the indexes that minimizes the difference between gtau and j
    verbose('the index of the argmin: %s, %s'%(j,k))
    mindist = abs(cands[j]- jlist[k])
    verbose('the minimal difference between the two list element is %s'%CDF(mindist))
    if mindist < eps:
        used.append(k)
        g = cosreps[j]
        zlist.append(g.acton(tau))
        verbose('the indexes of used j: %s'%str(used))
    else:
        raise RuntimeError('not able to identify the %s th point '%(i+1))
    return zlist


def evaluate_sym(lst):
    """
    given a list of complex numbers with length n,
    evaluate all basic symmetric polynomials its entries.
    """
    n = len(lst)
    Sym = SymmetricFunctions(QQ)
    e = Sym.elementary()
    return [e[i].expand(n)(lst) for i in range(1,n+1)]

def adjust(lst):
    """
    multiply by -1 every other entry
    """
    return [lst[i]*(-1)**(i+1) for i in range(0,len(lst))]



def exact_zeros(Fj,N,prec,eps):
    """
    return an approximation in complex upper half plane representatives
    of a well-defined set of points on Y0(N) , where Fj is the polynomial
    satisfied by their j-invariant. See readme.md for details.

    prec -- precision of computation
    eps -- allowed error to detect exceptional points.
    """
    try:
        S.<x> = QQ[]
        Fj = S(Fj)
    except: pass
    if not Fj.is_irreducible():
        raise NotImplementedError('Fj must be irreducible.')
    if not is_prime(N):
        raise NotImplementedError('N must be a prime.')
    C = ComplexField(prec)
    R.<x> = C[]
    jlist = R(Fj).roots(multiplicities = False)
    exceptional_jlist = [j for j,f in exc_j(N,prec)]
    j0 = jlist[0]
    i = closest(j0,exceptional_jlist)
    if abs(exceptional_jlist[i][0] - j0) < eps:
        exceptional = True
        form = exceptional_jlist[i][0]
        a,b,c = form[0],form[1],form[2]
        d = gcd(gcd(a,b),c)
        D = form.discriminant() // d^2
        quadforms =  BinaryQF_reduced_representatives(D,primitive_only = True) # find all reps primitive forms of that discriminant
        zs = []
        for f in quadforms:
            pass
            # do something: mainly just adjust
    else:
        exceptional = False
        taulist = taus(Fj,N,prec,jlist)
        zs = from_tau_to_z(taulist,jlist,N,prec,eps)
    return zs





