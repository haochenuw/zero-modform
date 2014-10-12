### find the exact z ###
### dependence: depend on exceptional-points.sage to compute the exceptional points.
### implicitly, depend on zero-modform.sage to provide the zero polynomial Fj.


### added 7/12: in fact if the polynomial in question is CM then we don't need to mess with it.
# so it suffices to get a list of CM polynomials that are possible and forget about the exceptional points.
# what kind of CM is possible? Is there an upper bound on discriminant?

# or we could just...? formula for the discriminant of the whole thing? Probably not, but I could bound the degree, so...?






from itertools import *
from mpmath import mp, hyp2f1
import numpy as np
#load('exceptional-points.sage')

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


def taus(jlist,prec):
    """
    return taulist such that j( taulist[i] ) = jlist[i]
    """
    t = cputime()
    C = ComplexField(prec)
    R.<x> = C[]
    #jlist = R(Fj).roots(multiplicities = False)
    taulist = []
    t1 = cputime()
    for j in jlist:
        y = R(256*(1-x)^3 - j*x^2).roots(multiplicities = False)[0]
        lam = R(x-x^2-y).roots(multiplicities = False)[0]
        tau = hyp(lam,prec)
        taulist.append(tau)
    verbose('2nd step done, found raw taus, took %s seconds'%cputime(t1))
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




def from_tau_to_z(taulist,jlist,N,eps,prec):
    """
    given a list of taus, for each tau find g such that j(g(tau)) and j(Ngtau) are both in jlist
    by 'in jlist' I mean epsilon-close to an element in jlist.


    Output:

    zlist -- a list of z's : complex zeros.
    pairs -- a pairing data of integers from 0 to len(jlist)-1.
    """
    verbose('eps = %s'%eps)
    C = ComplexField(prec)
    pari.set_real_precision(prec)
    G = Gamma0(N)
    cosreps = list(G.coset_reps())
    n = len(taulist)
    verbose('n = %s'%n)
    used = []
    result = []
    pairs = []



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
            pairs.append((i,k))
            used.append(k)
            g = cosreps[j]
            result.append(C(g.acton(tau)))
            verbose('the indexes of used j: %s'%str(used))
        else:
            raise RuntimeError('not able to identify the %s th point '%(i+1))
    return result, pairs


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



def exact_points(F,N,prec,eps=1):
    """
    return an approximation in complex upper half plane representatives
    of a well-defined set of points on Y0(N) , where Fj is the polynomial
    satisfied by their j-invariant. See readme.md for details.

    prec -- precision of computation
    eps -- allowed error to detect exceptional points.
    """
    try:
        S.<x> = QQ[]
        F = S(F)
    except:
        pass
    if not F.is_irreducible():
        raise NotImplementedError('F must be irreducible.')
    C = ComplexField(prec)
    R.<x> = C[]
    jlist = R(F).roots(multiplicities = False)
    verbose('j invariants computed')
    taulist = taus(jlist,prec)
    verbose('complex candidate taus computed')
    zs = from_tau_to_z(taulist,jlist,N,eps,prec)
    return zs


from sage.schemes.elliptic_curves.heegner import * # need ringclassfield

### checking if a polynomial is CM, not tested yet.
def is_cm(f):
    """
    check if f is a CM polynomial, i.e. if f = H_O(x) for some...
    INPUT:
        f -- an irreducible polynomial in QQ[x]
    """
    R.<x> = QQ[]
    try:
        f = R(f)
    except:
        raise ValueError('f must be an intege polynomial')
    df = f.degree()
    if not f.is_irreducible():
        raise NotImplementedError
    D = ZZ(f.disc())
    v = sorted(D.prime_factors())
    verbose('prime factors of discriminant: %s'%str(v))
    possibleList = []
    for n in product_combinations(v): # since the ramified primes are all factors of D
        K.<a> = QuadraticField(-n)
        dK = K.discriminant()
        verbose('dK = %s'%dK)
        cond = 1
        while True:
            Kc = RingClassField(dK,cond)
            deg = Kc.degree_over_K()
            if deg < df:
                cond += 1
            elif deg == df:
                possibleList.append((dK,cond))
                break
            else: # deg > df, the discriminant dK is impossible here.
                break
    if len(possibleList) == 0:
        return False
    else:
        for dK,c in possibleList:
            disc = dK*cond^2
            h_disc = hilbert_class_polynomial(disc)
            if h_disc == f:
                return True, disc
        return False


import itertools as it

def product_combinations(lst):
    n = len(lst)
    pick = list(it.product(range(2), repeat=n))
    result = []
    for p in pick:
        p = zip(lst,list(p))
        result.append(prod([a**e for a, e in p]))
    return result

def gamma0n_good_rep(z,N):
    """
    z -- a complex number representing a point [z] on X_0(N)

    OUTPUT:
        (z',sign) such that z' \in H, sign = +-1, and [z'] = z
        or if sign = -1, [z'] = [wN(z)]
    """
    C = z.parent()
    prec = C.prec()
    R = RealField(prec)
    half = R(1/2)
    const = C(1/sqrt(N))
    sign = 1
    N  = C(N)
    while True:
        if abs(z.real_part()) < half:
            if abs(z) > const:
                return (z,sign)
            else:
                z = -1/(N*z)
                sign = -sign
        else:
            re = z.real_part()
            z = z - C(floor(re))
            if abs(z.real_part()) > half:
                z = z - 1
    return (z,sign)

def x_poly_prec(x,deg,prec,name='t'):
    #x_reduced_prec = x.n(prec=floor(prec*0.8))
    poly_factorization = ZZ[name](x.algdep(deg)).factor()
    verbose('polynomial factorization is %s'%poly_factorization)
    # Test that full precision x is a root of one of the factors
    # to precision greater than self._prec*0.9
    found_correct_factor=False
    RF = RealField(prec)
    n = len(poly_factorization)-1
    while n>=0 and found_correct_factor==False:
        poly = poly_factorization[n][0]
        newprec = floor(-log(abs(poly(x)))/log(RF(2)))
        verbose('real precision is: %s'%newprec)
        if newprec >= floor(prec*0.9):
            result = poly
            found_correct_factor=True
            n=0
        n -= 1
    # Raise error if we didn't find a polynomial factor for which
    # x is a root
    if found_correct_factor==False:
        s = "Defining polynomial does not stabilize. Try using greater precision."
        raise FloatingPointError(s)
    return result


def map_to_curve(E,zlist):
    """
    return the complex point on the curve E such that
    P = \sum_{z in zlist} phi(z)
    """
    phi = E.modular_parametrization()
    Pcomplex = sum([phi.map_to_complex_numbers(z) for z in zlist])
    return E.elliptic_exponential(Pcomplex)




