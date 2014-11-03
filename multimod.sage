"""
computes a relationship F(r,u) = 0 with integer coefficients
when r, u are modular functions on X0(N) with known valence
and known q-expansion to arbitrary precision.
"""

from sage.matrix.matrix_integer_dense import _lift_crt


def Bbound(N,n, C = 2, k = 20):
    """
    return a naive upper bound for the q^n coefficient of f^N
    where f = sum a_m q^m with |a_m| <= C m^k
    """
    return RealField(100)((n/N)**(N*k)*(C**N)*binomial(n+N-1,N-1))


def ell_points(N):
    """
    return a tuple (x,y) of number of
    elliptic points of order 2 and 3 on Gamma0(N)
    """
    v = N.prime_divisors()
    if Mod(N,4) == 0:
        x = 0
    else:
        x = prod([(1 + kronecker(-1,p)) for p in v if p != 2])
    if Mod(N,9) == 0:
        y = 0
    else:
        y = prod([(1 + kronecker(-3,p)) for p in v])
    return (x,y)


def ABC(k,d,N):
    """
    the exponents on on delta, E4 and E6 on the denominator, used in the wt_factor function
    """
    etwo, ethree = ell_points(N)
    return(ZZ(k*d/12-ethree*k/3-etwo*k/4),ethree*k, etwo*k/2)

def wt_factor(k,d,N,prec):
    """
    see Ono's weierstrass point paper for details and an example.
    """
    E4 = eisenstein_series_qexp(4, prec)
    E6 = eisenstein_series_qexp(6, prec)
    delta = delta_qexp(prec)
    A,B,C = ABC(k,d,N)
    verbose('A,B,C = %s,%s,%s'%(A,B,C))
    return delta**A*E4**B*E6**C

def normalize(g):
    """
    setting the first nonzero coefficient to be 1
    """
    return g/g.padded_list()[g.valuation()]

def get_primes(P,M):
    """
    return a list of consecutive primes, starting from P
    such that the product is bigger than M
    """
    result = [P]
    prod = P
    for p in Primes():
        if prod > M:
            return result
        if p > P:
            prod *= p
            result.append(p)





### Computing a (very crude) upper bound on the absolute value
# of the coefficients of the zeropoly, to be used in the multimodular algorithm.

def cn(n,d): # the maximal n-th coefficient the norm can take
    """
    here d = the index of Gamma0(N) in SL_2(Z)
    """
    return RR((2**d)*binomial(n+d-1,d-1)*(n/d*N)**(k*d/2))



def jqbound(m,n):
    """
    upperbound on the q^n coeff of (jq)^m
    """
    return RealField(100)(1200*e**(4*pi*sqrt(m*n)))

def coef_bound(g,G):
    """
    bound on the coefficient of F_f
    where we have W = F_f(j), W has valuation -(2g-2)
    such that W*q^{2g-2} = sum a_n q^n, with a_0 = 1
    j = q^-1 + ...
    """
    result = [1]
    prec_small = 2*g-2+20 # for safety, in fact using 2g-2 was enough.
    E4 = eisenstein_series_qexp(4, prec_small)
    E6 = eisenstein_series_qexp(6, prec_small)
    delta = delta_qexp(prec_small)
    jlow = delta/normalize(E4**3)
    q = jlow.parent().gen()
    jq = QQ[[q]](jlow*q)  # = j*q, a power series, not a Laurent series
    deg = 2*g-2 # the degree of the polynomial
    a = (jq)**deg+G
    #print type(a)
    while deg > 0:
        n = 2*g-2-deg+1
        a = a + (a.padded_list()[n])*jq**(deg-1)
        result.append(a.padded_list()[n])
        deg =  deg-1
    return max([RR(a) for a in result])

@cached_method
def h8qbounds(M,N):
    """
    upperbound on the q^n coeff of (h8q)^m
    where 1 <= m <= M and 1 <= n <= N
    where h8 = EtaProduct(8,{8:8,2:4,4:-12})
    """
    result = [[0 for x in xrange(N)] for x in xrange(M)]
    prec = N + 10
    h8 = EtaProduct(8,{8:8,2:4,4:-12}).qexp(prec)
    R = h8.parent()
    h8q = R(h8*q)
    currentPower = h8q
    for i in range(1,M+1):
        result[i] == [abs(a) for a in currentPower.padded_list()[:N]]
        currentPower *= h8q
    return result



def coef_bound_new(g,lstG,useh8 = False):
    """
    the upper bound on the coefficients
    """
    result = [1]
    deg = 2*g-2
    for i in range(1,deg+1):
        ai = RR(lstG[i])
        lstnew = [jqbound(2*g-2-j,i-j) for j in range(i)]
        assert len(result) == len(lstnew)
        adj = RR(sum([a*b for a,b in zip(lstnew,result)]))
        result.append(ai+adj)
        #print 'aterm = ', aterm
        #print 'adj =', adj
        print 'current large = ', result[-1]
    assert len(result) == deg+1
    M = result[-1]
    if useh8:
        M = M**{1/12}*2**(4*g-3)
    return RealField(100)(M)

def find_B_bound(N,k,d,prec,useh8 = False,useh16 = False):
    """
    d = [SL_2(Z): Gamma_0(N)]
    k = the weight of the modular form (usually 2)
    """

    g = Gamma0(N).genus()
    wtf = wt_factor(k,d,N,prec)

    val = wtf.valuation(); valNf = val - (2*g-2);

    R = wtf.parent(); q = R.gen()

    wtf_trunc = R(wtf.padded_list()[val:val+2*g-2+5]).add_bigoh(val+2*g-2+5); wtfprime = R(1/wtf_trunc);

    v = R([abs(a) for a in wtfprime.padded_list()[:2*g-2+1]]);

    Nfbound = R([0]*valNf + [cn(ZZ(n),d) for n in range(valNf,valNf+2*g-2+10)])
    bigG = (Nfbound*v); bigG = bigG/q**bigG.valuation()
    # bigG is the power series that majorizes F_f(q)
    coefB = coef_bound(g,bigG)
    # coefB is the upper bound on the coefficient of F when we write F_f = F(j).
    if useh8:
        deg = 2*g-2
        C = (66.2)^deg * (2^(deg/12) -1);
        coefB = RealField(100)(C*(coefB)^(1.0/12))
    if useh16:
        deg = 2*g-2
        C = (66.2)^deg * (2^(deg/12) -1);
        coefB = RealField(100)(C*(coefB)^(1.0/12))
        coefB = RealField(100)(sqrt(2)*2*coefB^(0.5))
    return coefB




def degree(l,dega,degb):
    """
    obtain the highest degree of x and y appearing in the polynoial l represented as a list.
    """
    assert len(l) == (dega+1)*(degb+1)
    maxa,maxb =0,0
    for i in range(len(l)):
        if l[i] != 0:
            maxa = max(maxa,l[i]//(degb+1))
            maxb = max(maxb,l[i]%(degb+1))




# for chow-heegner purposes, or for any occasion where the zeros of a
# function is not known to be integral under the other function.

import os

def relation_zz(r,u,degr,degu,padding,description):
    matrix = []
    num_terms = (degu+1)*(degr+1)
    prec = num_terms + padding
    prec_small = prec - padding//2
    t = cputime()
    verbose('prec_small = %s'%prec_small)
    R = r.parent()
    q = R.gen()
    pows_of_r, pows_of_u = [R(1)],[R(1)]
    while len(pows_of_r) < degr+1:
        s = pows_of_r[-1]
        pows_of_r.append(s*r)
    while len(pows_of_u) < degu+1:
        m = pows_of_u[-1]
        pows_of_u.append(m*u)

    verbose('powers of r and u made')

    verbose('making the matrix ...')

    for a in range(degr+1):
        for b in range(degu+1):
            matrix.append(R(pows_of_r[a]*pows_of_u[b]).add_bigoh(prec_small).padded_list()[:prec_small])
    Mp = Matrix(QQ,num_terms,prec_small,matrix)
    save(matrix, os.path.join(os.environ['HOME'],'poly-relation','results','matrix-%s'%description))

    verbose('matrix made, took %s seconds'%cputime(t))

    verbose('computing the kernel...')
    Ker = Mp.kernel()
    verbose('got here')

    K = list(Ker.basis_matrix())

    verbose('dimension of kernel =  %s'%len(K))
    if len(K) == 0: # kernel must be one-dimensional.
        raise ValueError('got trivial kernel')


    elif len(K) > 1:
        verbose('the dimension of the kernel is: %s'%len(K))
        if not os.path.exists('debug'):
            os.makedirs('debug')
        save(K,'debug/kernel-%s'%description)
        #save(K, os.path.join(os.environ['HOME'],'poly-relation','debug','kernel-%s'%description))
        raise ValueError('kernel is greater than one dimensional, please debug')

    k = K[0]
    kmat = Matrix(k[0].parent(),1,len(k),list(k))

    if not os.path.exists('results'):
        os.makedirs('results')

    save(kmat,'results/kmat-%s'%description)
    save(list(kmat), os.path.join(os.environ['HOME'],'poly-relation','results','relation-%s'%description))
    verbose('the whole computation took %s seconds'%cputime(t))

    kk = list(kmat)[0]
    R.<r,u> = PolynomialRing(QQ,2)
    F = R(0)
    for a in range(degr+1):
        for b in range(degu+1):
                F += kk[a*(degu+1)+b]*r**a*u**b
    save(F,'results/F-%s'%description)
    return F


# one prime at a time




def zeropoly_modp(p,degr,degu,r,u,padding,firstTime = True):

    try:
        Mp = load(os.path.join(os.environ['HOME'],'poly-relation','results','Mp-%s-%s'%(E.label(),p)))
        verbose('already computed.')
        return Mp
    except:
        pass

    num_terms = (degu+1)*(degr+1)
    prec = num_terms + padding
    prec_small = prec - padding//2
    verbose('prec_small = %s'%prec_small)
    t = cputime()
    result = []
    verbose('computing modulo the prime %s'%p)

    # convert everything modulo p and computing r
    # p = 0 means we are working over QQ


    q = var('q')
    Rp = GF(p)[[q]]
    q = Rp.gen()
    """
    if not zzfirst:
        RpL = Rp.laurent_series_ring()
        deltap = Rp(delta)
        E4p = normalize(Rp(E4))
        u0p = deltap/(E4p**3)
        if useh8:
            up = Rp(EtaProduct(8,{8:4,1:-4,2:2,4:-2}).q_expansion(prec))
        else:
            up = u0p
        du0p = Rp(u0p.derivative())
        fp = Rp(f)
        gp = deltap(q = q^2)/deltap
        gnew =  Rp((512*gp-1)/(256*gp+1)) #g = (g0-1/512)/(g0+1/256)
        rp = Rp(-fp*u0p*gnew/(q*du0p))
        if useh4:
            hp = RpL(EtaProduct(4,{1:8,4:-8}).q_expansion(prec) + 32)
            rp = RpL(rp)*hp
    """
    RpL = Rp.laurent_series_ring()
    up = RpL(u)
    rp = RpL(r)
    pows_of_r, pows_of_u = [RpL(1)],[RpL(1)]



    verbose('rp and up series made')

    # now make powers of rp and up.
    while len(pows_of_r) < degr+1:
        s = pows_of_r[-1]
        pows_of_r.append(s*rp)
    while len(pows_of_u) < degu+1:
        m = pows_of_u[-1]
        pows_of_u.append(m*up)

    verbose('powers of rp and up made')

    verbose('making the matrix ...')

    for a in range(degr+1):
        for b in range(degu+1):
            result.append(Rp(pows_of_r[a]*pows_of_u[b]).add_bigoh(prec_small).padded_list()[:prec_small])
    Mp = Matrix(GF(p),num_terms,prec_small,result)
    verbose('matrix made for the prime p = %s'%p)

    verbose('computing the kernel...')
    t1 = cputime()
    K = list(Mp.kernel().basis_matrix())
    if len(K) == 0: # kernel must be one-dimensional.
        raise ValueError('got trivial kernel')


    elif len(K) > 1:
        verbose('the dimension of the kernel is: %s'%len(K))
        save(K, os.path.join(os.environ['HOME'],'poly-relation','debug','kernel-%s-%s'%(E.label(),p)))
        raise ValueError('kernel is greater than one dimensional for the prime p = %s, please debug'%p)


    #k = K[0][:degu+1]
    k = K[0]

    kmat = Matrix(k[0].parent(),1,len(k),list(k))
    save(kmat, os.path.join(os.environ['HOME'],'poly-relation','results','Mp-%s-%s'%(E.label(),p)))
    verbose('the whole computation for p = %stook %s seconds'%(p,cputime(t)))
    return kmat

modulus = 11


@parallel(ncpus = 10)
def _multimodular_zeropoly(plist,degr,degu,r,u,padding,mod):
    multiMat = []
    for p in plist:
        if Mod(p,modulus) == mod:
            kmat = zeropoly_modp(p,degr,degu,r,u,padding)
            multiMat.append(kmat)
            kmat = None
    return multiMat

def multimodular_zeropoly(plist,degr,degu,r,u,padding):
    multiMat = []
    inputs = [(plist,degr,degu,r,u,padding,a) for a in range(modulus) if gcd(a,modulus) == 1]
    for output in _multimodular_zeropoly(inputs):
        multiMat += output[1]
    return multiMat

def get_xinv(E,prec):
    e = list(gp.elltaniyama(E,prec))[0];
    x = var('x')
    e = ZZ[[x]]([gp.polcoeff(e,a) for a in range(-2,prec-1)]).add_bigoh(prec+1)
    x = e.parent().gen()
    return x^2/e



