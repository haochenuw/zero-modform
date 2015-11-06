from sage.matrix.matrix_integer_dense import _lift_crt


def zero_poly(f,level = None,weight = None, algorithm = 'multimodular'):
    """
    Input
        f -- a modular form, usually the newform attached to some elliptic curve.
        level -- the level of f (Must be prime).
        weight -- the weight of f.
        algorithm:
            'multimodular' -- default. Use the multimodular method to compute norm.
            'integer' -- compute everything over ZZ

    Output:
        The polynomial satified by the j-invariants of zeros f(z)(dz)^(weight)/2.
    """
    if level is None:
        level = f.level()
    if weight is None:
        weight = f.weight()

    if not is_prime(level):
        raise ValueError('level ( = %s) must be prime.'%level)

    p = level
    k = weight

    t = cputime()

    # compute various precisions in the computations of series involved
    v = precs(p,k)

    verbose('Computing %s Fourier coefficients of f...'%(v[2] +10))
    try:
        f = f.qexp(v[2]+10)
    except:
        pass

    verbose('Done computing the expansion of f, took %s seconds'%(cputime(t)))
    verbose('precision of f = %s'%f.prec())

    # computing the norm of f
    if algorithm == 'multimodular':
        Nf = Norm(f,p,k)
    else:
        Nf = Norm_int(f,p,k)

    verbose('done computing the norm of f, took %s seconds'%cputime(t))

    verbose('computing Fq := F_{E,j}(j(q))...')
    Fq = normalize(Nf)/normalize(weight_factor(p,k))

    # Knowing that the Laurent series Fq is a polynomial in the j-invariant,
    # we recover the polynomial F_{E,j}.
    L = Fq.truncate_laurentseries(1).coefficients()
    alist = []

    prec_low = v[0]


    # computing the q-expansion of the j-invariant.
    E4 = eisenstein_series_qexp(4, prec_low+1)
    delta = delta_qexp(prec_low+1)
    j_inv = normalize(E4**3)/delta
    jinvs = [j_inv**0,j_inv]
    while len(jinvs) < len(L):
        a =  jinvs[-1]
        jinvs.append(a*j_inv)

    verbose('j-invariants computed')

    # recognize the coefficients. i.e. find the polynomial F s.t. F(j(q)) = Fq
    for i in range(prec_low):
        d = prec_low-1-i
        fd = get_principal_part(jinvs[d],prec_low)
        ad = L[i]/fd[i]
        alist.append(ad)
        L = [L[j]-ad*fd[j] for j in range(prec_low)]

    F = QQ[x](alist[::-1])
    verbose('zero polynomial computed')
    return F


def critjpoly_modl(E,l = None):
    """
    compute the
    critical j-polynomial modulo l of an elliptic curve E.
    conductor must be prime p. and l = 1 modulo p.
    """
    p = E.conductor()
    weight = 2
    v = precs(p,weight)

    # first time the preparation
    t = cputime()
    f = E.modular_form().qexp(v[2]+10)
    print 'terms of expansions used is %s'%f.prec()

    # also get the estimate on the number of primes used.

    bigN = coef_bound(f,p,weight)
    llist = mod1_primes(p, 2*bigN)
    nprimes = len(llist)

    print 'the number of primes used = %s'%nprimes

    if l is None:
        l = llist[0]
    print 'the prime l = %s'%l
    phip = cyclotomic_polynomial(p)

    print 'prep time = %s'%cputime(t)

    t = cputime()
    # time the norm computation mod one l.
    normfmodl = norm_mod_l(f,p,phip,l)
    print 'norm time = %s'%cputime(t)


    # print 'time estimate', RR(t1) + RR(t2)*RR(nprimes)

    t = cputime()
    Rmodl.<q> = GF(l)[[]]
    Normf = Rmodl(f.truncate(f.prec()//p))*Rmodl(normfmodl.list())
    print 'part 2.5, get Norm: %s'%cputime(t)
    Fq = normalize(Normf)/Rmodl(normalize(weight_factor(p,weight)))
    # get the principal part + constant term.
    L = Fq.shift(-Fq.valuation()).power_series().padded_list(-Fq.valuation()+1)
    print 'part 3, get weight factors and L: %s'%cputime(t)

    # powers of j-invariants modulo l
    t = cputime()
    prec_low = v[0]
    if len(L) != prec_low:
        raise ValueError('number of coefficients is incorrect. Please debug.')
    RmodlL = Rmodl.laurent_series_ring()
    E4 = Rmodl(eisenstein_series_qexp(4, prec_low+1))
    delta = Rmodl(delta_qexp(prec_low+1))
    j_inv = RmodlL(normalize(E4**3)/delta)
    jinvs = [j_inv**0,j_inv]
    while len(jinvs) < len(L):
        a =  jinvs[-1]
        jinvs.append(a*j_inv)
    print 'part 4, get j invariants: %s'%cputime(t)

    t = cputime()
    alist = []
    for i in range(prec_low):
        d = prec_low-1-i
        fd = get_principal_part(jinvs[d],prec_low,modp = True)
        ad = L[i]/fd[i]
        alist.append(ad)
        L = [L[j]-ad*fd[j] for j in range(prec_low)]
    print 'part 5, recognization: %s'%cputime(t)

    return GF(l)[x](alist[::-1])
    # print 'computation for prime l took %s seconds'%cputime(t)

def zero_poly_modp(f,p,k):
    '''
    This function computes the zero poly F_{E,j} modulo p,
    where p is the level.
    (Not directly used for F_{E,j} computation, but can provide a sanity check.)
    '''
    v = precs(p,k)
    try:
        f = f.qexp(v[1]+1) # added one since the way qexp works.
    except:
        pass
    fmodp = f.change_ring(GF(p))
    Nfmodp = normalize(fmodp^2)
    wf = weight_factor(p,k)
    #verbose('wf has valuation %s'%(wf.valuation()))
    wfmodp = wf.change_ring(GF(p))
    #verbose('wfmop = %s'%str(wfmodp))
    Ffmodp = Nfmodp/wfmodp
    # In this case we don't need to compute the norm, since we have Nf = f^2(mod p)
    val = -Ffmodp.valuation()
    q = Ffmodp.parent().gen()
    shifted = GF(p)[[q]](Ffmodp*(q**val))
    L = shifted.padded_list()
    L = L[:val+1]

    #L = Ffmodp.truncate_laurentseries(1).coefficients()
    alist = []
    prec_low = precs(p,k)[0]
    assert len(L) == prec_low
    # computing the j-invariant
    E4 = eisenstein_series_qexp(4, prec_low+1)
    delta = delta_qexp(prec_low+1)
    j_invmodp = (normalize(E4**3)/delta).change_ring(GF(p))
    jinvs = [j_invmodp**0,j_invmodp]
    while len(jinvs) < prec_low:
        a =  jinvs[-1]
        jinvs.append(a*j_invmodp)
    assert len(jinvs) == prec_low
    for i in range(prec_low):
        d = prec_low-1-i
        fd = get_principal_part(jinvs[d],prec_low,modp = True)
        ad = L[i]/fd[i]
        alist.append(ad)
        L = [L[j]-ad*fd[j] for j in range(prec_low)]
    T.<x> = GF(p)[]
    F = T(alist[::-1])
    verbose('zero polynomial computed and saved.')
    return F


def normalize(g):
    """
    return the scaled power series so that its leading coefficient is one.
    """
    return g/g.padded_list()[g.valuation()]



def weight_factor(p,k,prec = None):
    if prec is None:
        prec = precs(p,k)[0]
    A,B,C = weight_index(p,k)
    E4 = eisenstein_series_qexp(4, prec)
    E6 = eisenstein_series_qexp(6, prec)
    delta = delta_qexp(prec+1)
    verbose('adjustment factors: %s, %s, %s'%(A,B,C))
    return delta**A*E4**B*E6**C


def precs(p,weight):
    """
    returns the precisions needed for computing the zero
    polynomial of a modular form f.
    return a tuple (k(g-1)+2,kg+1,p(kg+1)). Ordered from
    small to large, and they are the precisions needed
    for (Delta,E4,E6), (f), (f(q^{1/p})), respectively.
    We are going to add 1 for safety.
    """
    k = weight
    g = Gamma0(p).genus()
    return (k*(g-1)+ 1 , k*g+1+1 , p*(k*g+1 + 1))

def coef_bound(f,p,k):
    """
    Used for the multimodular algorithm.
    bound the largest possible coefficient of
    a product of p modular forms of level p and weight k. with degree prec, where if f = \sum a_nq^n
    then each term in the product is f = \sum a_n b_n q^n
    with |b_n| = 1 for all n.
    """
    prec = f.prec()
    verbose('the precision of the power series = %s'%prec)
    verbose('weight = %s'%k)
    v = f.padded_list(prec)
    C = max([abs(v[n])/RR(n^(k/2)) for n in range(1,prec)])
    verbose('the constant C used is %s'%RR(C))
    return RR(binomial(prec + p-1, p-1)*(RR(C)**p)*(prec//p)^(p*k/2))


def coef_bound_weightfree(f,p):
    """
    given any cusp form f of level divisible by p,
    compute a bound on the abs value of coefficients of norm(f),
    defined as
        norm(f) = \prod_{i mod p} f((z+i)/p).
    p doesn't have to be a prime.
    This bound does not depend on weight of f. It may or may not give a worse bound than coef_bound.
    Also, this bound is more general: it does not require f to be an eigenform.
    """
    M = f.prec()
    v = list(f.polynomial())
    # we compute a constant d(f) such that |a_n(f)| \leq (d(f))^n.
    df = RR(max([abs(RR(v[n]))**(1.0/n) for n in range(1,len(v))]))
    verbose('df = %s'%df)
    return RR(binomial(M + p-1, p-1)*(df**M))




def mod1_primes(p,N):
    """
    p -- a prime
    N -- a (large) integer
    return a list of primes l such that l = 1 (mod p)
    and the product of all primes l in the list is larger than N.
    This is used for multimodular algorithms.
    """
    a = 1
    result = []
    prod = 1
    while prod < N:
        a = a+p
        if is_prime(a):
            result.append(a)
            prod *= a
    return result

def get_principal_part(f,n,modp = False):
    """
    returns the principal part + constant part of the
    coefficients of a Laurent series f
    if the length l < n, put (n-l) zeros in front so that the length of output is
    always n.
    """
    if modp:
        q = f.parent().gen()
        GF = f.base_ring()
        val = -f.valuation()
        L = GF[[q]](f*q**val).padded_list()
        s = L[:val+1]
    else:
        s = f.truncate_laurentseries(1).coefficients()
    return  [0]*(n-len(s)) + s



def norm_mod_l(f,p,phip,l):
    """
    Take a modular form f of level p with integer coefficents,
    this function computes \prod_i (f(zeta_p^i q^1/p))
    modulo a prime l that is congruent to 1 modulo p.

    Input: f -- a modular form with integer coefficients
    """

    # get a p-th root of unity in the finite field F_l.
    phipl = phip.change_ring(GF(l))
    v = phipl.roots(multiplicities = False)
    v.sort()
    z = v[0]

    prec_big = f.prec()
    prec_med = (prec_big-10)//p
    print("Number of multiplications to perform on powerseries mod %s with precision %s : %s"%(l, prec_big,p))
    L = f.padded_list() # raised everything to pth power q = q'^p

    # construct the finite field,
    Fl = z.parent()
    Rl.<x> = Fl[[]]
    F = Rl(1)
    Lbar = [Fl(a) for a in L]

    t = cputime()
    for k in range(p):
        tmp = Rl([z**(i*k)*Lbar[i] for i in range(prec_big)])
        F = F * tmp
        F = F.truncate(prec_big)
    print("The norm computation for the prime %s is performed within %s seconds"%(l, cputime(t)))
    L = F.padded_list()[0::p]
    return Matrix(GF(l),1,prec_med,L[:prec_med])


@parallel(ncpus = 32)
def _norm(f,p,phip,llist,a):
    Matlist = []
    for l in llist:
        verbose('l = %s'%l)
        if Mod(l,64) == a:
            t = cputime()
            Matlist.append(norm_mod_l(f,p,phip,l))
    return Matlist

def norm(f,p,k, proof = False, ncpus = 64): # to-do: add parallel = false
    """
    Given a modular form of weight k and level p.
    Use multimodular algorithm (computing mod primes, and use CRT to lift back)
    to compute the norm of a modular form f of level p.
    ASSUMPTION: f has integer coefficients.
    """
    verbose('computing the norm of f, with level %s and weight %s...'%(p,k))
    Matlist = []
    count = 0
    prec_big = f.prec()
    prec_med = (prec_big)//p

    bigN1 = floor(coef_bound_weightfree(f,p))
    bigN2 = floor(coef_bound(f,p,k))
    bigN = min(bigN1,bigN2)
    if proof is False:
        bigN = RR((bigN)**(0.75))
    phip = cyclotomic_polynomial(p)
    verbose('multimodular bound is %s'%bigN.n())
    llist = mod1_primes(p,2*bigN) #multiply by 2 since we are inside the interval [-bigN, bigN]
    verbose('we are using %s primes'%len(llist))
    verbose('the primes are %s'%str(llist))

    inputs = [(f,p,phip,llist,a) for a in range(ncpus*2) if gcd(a,ncpus*2) == 1]
    count = 0
    for output in _norm(inputs):
        count += 1
        Matlist += output[1]
        print 'output[1] = %s'%output[1]
        print 'count = %s'%count
    print 'Number of vectors mod p obtained = %s'%len(Matlist)

    M = _lift_crt(Matrix(ZZ, 1, prec_med), Matlist)
    R.<q> = QQ[[]]
    return R(M.list()).add_bigoh(prec_med)

def Norm(f,p,k):
    nf = norm(f,p,k)
    prec_med = nf.prec()
    R.<q> = QQ[[]]
    return R(nf*f.truncate(prec_med)).add_bigoh(prec_med)


# algorithm for general square free level N #

def base_prec(N,k):
    """
    the precision B needed for Norm_N(f), i.e., we need
    to compute N_N(f)(q) = \cdots + q^(B) + \cdots
    """
    g = Gamma0(N).genus()
    return (k*(g-1)+1, k*(g-1) + 1 + 2**sigma(N,0) + 10)


def weight_index(N,k):
    """
    the exponent of D, E4, E6 on the denominator
    """
    v = N.prime_divisors()
    B = prod([(1 + kronecker(-3,p)) for p in v])*k
    C = prod([(1 + kronecker(-1,p)) for p in v if p != 2])*k//2
    A = (k*Gamma0(N).index()-4*B-6*C)/12
    verbose('A,B,C = %s,%s,%s'%(A,B,C))
    try:
        return (ZZ(A),ZZ(B),ZZ(C))
    except:
        return None


def Norm_comp(f,N,weight,algorithm = 'multimodular'):
    """
    Given a modular form of weight k and level N(N = square free),
    using multimodular algorithm(computing mod primes, and use CRT to lift back)
    to compute the norm of f, with accuracy upto q^prec.
    ASSUMPTION: f has integer coefficients
    """
    v = N.prime_divisors()
    tmp = f
    for p in v:
        if algorithm == 'multimodular':
            tmp = Norm(tmp,p,weight)
        else:
            tmp = Norm_int(tmp,p,weight)
        weight = weight*(p+1)
    return tmp



def zero_poly_comp(f,level = None, weight = None, algorithm = 'multimodular'):
    """
    zero-polynomial for composite square free level

    Input:
        f -- a modular form.
        level -- the level of f.
        weight -- the weight of f.
        algorithm -- 'multimodular' is default.

    Output:
        the polynomial with roots the j-invariants of zeros of f(z)(dz)^(weight/2)

    EXAMPLE::
    sage:E = EllipticCurve('26a1')
    sage:f = E.modular_form()
    sage:F = zero_poly_comp(f,26,2,'critical'); F.factor()
    (x-1728)^2
    """
    if level is None:
        level = f.level()
        N = level
    if weight is None:
        weight = f.weight()
        k = weight

    if not level.is_squarefree():
        raise ValueError('N must be square free')

    prec_low,prec_med = base_prec(N,k)
    prec_high = sigma(N,1)*prec_med
    verbose('medium and high precs = %s, %s'%(prec_med,prec_high))
    try:
        f = f.qexp(prec_high+10)
    except:
        pass

    # computing norm
    Nf = Norm_comp(f,N,k,algorithm = algorithm)
    verbose('valuation of Norm(f) = %s'%Nf.valuation())
    verbose('done computing the norm')


    Fq = normalize(Nf)/normalize(weight_factor(N,k))

    L = Fq.truncate_laurentseries(1).coefficients()
    alist = []
    if len(L) != prec_low:
        raise ValueError('Length of principal part is off( = %s). Please debug.'%len(L))
    # computing the j-invariant

    n = len(L)

    E4 = eisenstein_series_qexp(4, n+1)
    delta = delta_qexp(n+1)
    j_inv = normalize(E4**3)/delta
    jinvs = [j_inv**0,j_inv]
    while len(jinvs) < n:
        a =  jinvs[-1]
        jinvs.append(a*j_inv)
    verbose('j-invariants computed')
    for i in range(n):
        d = prec_low-1-i
        fd = get_principal_part(jinvs[d],n)
        ad = L[i]/fd[i]
        alist.append(ad)
        L = [L[j]-ad*fd[j] for j in range(n)]
    F = QQ[x](alist[::-1])
    verbose('zero polynomial computed.')
    return F



def norm_int(f,m,use_cc = False):
    return Norm_int(f,m,use_cc, multiply_by_f = False)


def Norm_int(f,m,use_cc = False,multiply_by_f = True):
    """
    Computes the "norm" of a modular form.
    when m = p is prime, this is
    \prod f|_A for A in SL_2(ZZ)/Gamma0(p)
    recommened prec = (g^3-g+1)*(p+1), to get enough precision for later use.
    Input:
        f -- a power series
        m -- a positive integer
        use_cc -- if true, uses complex floating point numbers

    """
    verbose('Computing a %s-norm'%m)
    if use_cc:
        C = ComplexField(f.prec()//3)
        zetam = C.zeta(m)
    else:
        C.<zetam> = CyclotomicField(m)

    # verbose('parent of series = %s'%f.base_ring())
    try:
        q = f.parent().gen()
        f = QQ[[q]](f)
    except:
        raise ValueError("f must have rational coefficients")
    f = f.change_ring(C)
    R = f.parent()
    q = R.gen()
    prec = f.prec()
    verbose('prec of series = %s'%prec)

    verbose("Number of multiplications to perform on power series with %s terms is %s"%(prec,m))
    L = f.padded_list(prec) # raised everything to pth power q = q'^p
    F = R(1)
    for k in range(m):
        t = cputime()

        #tmp = R(sum([z**(ZZ(Mod(i*k,p)))*(q**i)*L[i] for i in range(prec)]))
        #verbose("Creating poly took %s seconds"%cputime(t))
        #t = cputime()

        tmp = R([zetam**(ZZ(Mod(i*k,m)))*L[i] for i in range(prec)]).add_bigoh(prec)
        #verbose("Creating poly tmp way took %s seconds"%cputime(t))

        F = F * tmp # k = 0,1,...,p-1
        verbose("The %s th multiplication is performed within %s seconds"%(k+1, cputime(t)))
    # convert back to q
    verbose('prec of F = %s'%F.prec())
    newprec = ZZ(F.prec()//m)
    newF = F.padded_list()[0::m]
    result = R(newF).add_bigoh(newprec)
    if multiply_by_f:
        result*f
    return result
