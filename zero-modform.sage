from sage.matrix.matrix_integer_dense import _lift_crt


def zero_poly(f,p,k,description,deg = 1,check =False,proof = True):
    """
    Input: f -- the power series expansion of a modular form
           k -- the weight of f
    Output: the polynomial satified by the j-invariants of f(z)(dz)^k/2
    """
    v = precs(p,k)
    try:
        f = f.qexp(v[2]+10) # added one since the way qexp works.
    except:
        pass
    Nf = Norm(f,p,k,deg,check,proof)
    verbose('done computing the norm')
    Fq = normalize(Nf)/normalize(weight_factor(p,k))
    L = Fq.truncate_laurentseries(1).coefficients()
    verbose("L = %s"%str(L))
    alist = []
    prec_low = precs(p,k)[0]
    verbose('prec_low =  %s'%prec_low)
    assert len(L) == prec_low
    # computing the j-invariant
    E4 = eisenstein_series_qexp(4, prec_low+1)
    delta = delta_qexp(prec_low+1)
    j_inv = normalize(E4**3)/delta
    jinvs = [j_inv**0,j_inv]
    while len(jinvs) < prec_low:
        a =  jinvs[-1]
        jinvs.append(a*j_inv)
    verbose('j-invariants computed')
    assert len(jinvs) == prec_low
    for i in range(prec_low):
        d = prec_low-1-i
        fd = get_principal_part(jinvs[d],prec_low)
        ad = L[i]/fd[i]
        alist.append(ad)
        L = [L[j]-ad*fd[j] for j in range(prec_low)]
    #print 'L = ', L
    T.<x> = QQ[]
    F = T(alist[::-1])
    save(F, os.path.join(os.environ['HOME'],'critical-point','zero-modform','F%s-%s-%s'%(p,k,description)))
    verbose('zero polynomial computed and saved.')
    return F


# the following is a simple version that computes the zero poly modulo p,
# where p is the level.

def zero_poly_modp(f,p,k,description):
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
    save(F, os.path.join(os.environ['HOME'],'critical-point','zero-modform','F%s-%s-%s'%(p,k,description)))
    verbose('zero polynomial computed and saved.')
    return F


def normalize(g):
    return g/g.padded_list()[g.valuation()]




def weight_factor(p,k):
    prec_low = precs(p,k)[0]
    A,B,C = weight_index(p,k)
    E4 = eisenstein_series_qexp(4, prec_low)
    E6 = eisenstein_series_qexp(6, prec_low)
    delta = delta_qexp(prec_low+1)
    verbose('adjustment factors: %s, %s, %s'%(A,B,C))
    return delta**A*E4**B*E6**C


def precs(p,k):
    """
    returns the precisions needed for computing the zero
    polynomial of a modular form f.
    return a tuple (k(g-1)+2,kg+1,p(kg+1)). Ordered from
    small to large, and they are the precisions needed
    for (Delta,E4,E6), (f), (f(q^{1/p})), respectively.
    We are going to add 1 for safety
    """
    g = Gamma0(p).genus()
    return (k*(g-1)+ 1 , k*g+1+1 , p*(k*g+1 + 1))

def coef_bound(f,p,k,deg,new = True):
    """
    bound the largest possible coefficient of
    a product of p modular forms with degree prec, where if f = \sum a_nq^n
    then each term in the product is f = \sum a_n b_n q^n
    with |b_n| = 1 for all n.
    """
    prec = f.prec()
    v = f.padded_list()
    if new:
        C = max([abs(v[n])/n^(k/2) for n in range(1,prec)])
    else: C = 2
    verbose('the constant C used is %s'%C)
    return RR(binomial(prec + p-1, p-1)*((C*deg*prec//p)^(k*p/2)))




def mod1_primes(p,N):
    """
    return a list of primes l such that l = 1 (mod p)
    and the product of all l in the list is > N
    used for CRT lift
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
    if the length l <n, put n-l zeros in there.
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
    z -- a root of \Phi_p(mod l) in GF(l)
    """
    phipl = phip.change_ring(GF(l))
    v = phipl.roots(multiplicities = False)
    v.sort()
    z = v[0]

    prec_big = f.prec()
    prec_med = (prec_big-pad)//p
    verbose("Number of multiplications to perform on powerseries with precision %s : %s"%(prec_big,p))
    L = f.padded_list() # raised everything to pth power q = q'^p
    Fl = z.parent()
    Rl.<x> = Fl[[]]
    F = Rl(1)
    Lbar = [Fl(a) for a in L]
    l = Fl.characteristic()
    verbose('Computing modulo the prime %s'%l)
    t = cputime()
    for k in range(p):
        tmp = Rl([z**(i*k)*Lbar[i] for i in range(prec_big)])
        F = F * tmp
        F = F.truncate(prec_big)
    verbose("The multiplication for the prime %s is performed within %s seconds"%(l, cputime(t)))
    #verbose("the old and new of coefficients for prime %s: %s, %s"%(l,len(F.padded_list()),len(F.padded_list()[0::p])))
    L = F.padded_list()[0::p]
    return Matrix(GF(l),1,prec_med,L[:prec_med])


@parallel(ncpus = 12)
def _norm(f,p,phip,llist,a):
    Matlist = []
    for l in llist:
        if Mod(l,21) == a:
            Matlist.append(norm_mod_l(f,p,phip,l))
    return Matlist

def norm(f,p,k,deg,check,proof):
    """
    Given a modular form of weight k and level p.
    sing multimodular algorithm(computing mod primes, and use CRT to lift back)
    to compute the norm of a modular form f of level p
    ASSUMPTION: f has integer coefficients
    """

    Matlist = []
    count = 0
    prec_big = f.prec()
    prec_med = (prec_big-pad)//p
    bigN = coef_bound(f,p,k,deg,new)
    if check:
        bigN = floor(bigN^(1.2))
    if not proof:
        verbose('proof = False')
        bigN = RealField(200)(sqrt(bigN)) # take its square root. Then the result could be double checked by mod p (?)
    phip = cyclotomic_polynomial(p)
    verbose('the bound on the size of coefficents of the norm is %s'%bigN)
    llist = mod1_primes(p,2*bigN) #multiply by 2 since we are inside the interval [-bigN, bigN]
    verbose('we are using %s primes'%len(llist))
    verbose('the primes are %s'%str(llist))

    inputs = [(f,p,phip,llist,a) for a in range(21) if gcd(a,21) == 1]
    for output in _norm(inputs):
        Matlist += output[1]

    # save the matrix list to a file
    save(Matlist, os.path.join(os.environ['HOME'],'critical-point','debug','modlcoeffs%s'%p))
    M = _lift_crt(Matrix(ZZ, 1, prec_med), Matlist)
    R.<q> = QQ[[]]
    return R(M.list()).add_bigoh(prec_med)

def Norm(f,p,k,deg,check = False,proof = True):
    R.<q> = QQ[[]]
    nf = norm(f,p,k,deg,check,proof)
    prec_med = nf.prec()
    return R(nf*f.truncate(prec_med)).add_bigoh(prec_med)

# algorithm for general square free level N #

def base_prec(N,k):
    """
    the precision B needed for N_N(f), i.e., we need
    to compute N_N(f)(q) = \cdots + q^(B) + \cdots
    """
    g = Gamma0(N).genus()
    return (k*(g-1)+1, k*(g-1) + 1 + sigma(N,0)*k//2)


def weight_index(N,k):
    """
    the exponent of D, E4, E6 on the denominator
    """
    v = N.prime_divisors()
    B = prod([(1 + kronecker(-3,p)) for p in v])*k
    C = prod([(1 + kronecker(-1,p)) for p in v if p != 2])*k//2
    A = (k*(sigma(N,1))-4*B-6*C)/12
    verbose('A,B,C = %s,%s,%s'%(A,B,C))
    try:
        return (ZZ(A),ZZ(B),ZZ(C))
    except:
        return None


def Norm_comp(f,N,weight,deg,check =False):
    """
    Given a modular form of weight k and level N(N = square free),
    using multimodular algorithm(computing mod primes, and use CRT to lift back)
    to compute the norm of f, with accuracy upto q^prec.
    ASSUMPTION: f has integer coefficients
    """
    v = N.prime_divisors()
    tmp = f
    for p in v:
        tmp = Norm(tmp,p,weight,deg,check)
        weight = weight*(p+1)
    return tmp


pad = 10 # the padding to make sure we have computed enough

# to-do: return to Norm to deal with the issue of precision information


def zero_poly_comp(f,N,k,description,deg =1,check = False):
    """
    zero-polynomial for composite square free level

    EXAMPLE::
    sage:E = EllipticCurve('26a1')
    sage:f = E.modular_form()
    sage:F = zero_poly_comp(f,26,2,'critical'); F.factor()
    (x-1728)^2
    """
    if not N.is_squarefree():
        raise NotImplementedError('N must be square free')
    prec_low,prec_med = base_prec(N,k)
    prec_high = sigma(N,1)*prec_med
    verbose('prec_high = %s'%prec_high)
    try:
        f = f.qexp(prec_high+pad+1)
    except:
        pass
    Nf = Norm_comp(f,N,k,deg,check)
    num_terms = prec_low
    verbose('valuation of Nf = %s'%Nf.valuation())
    verbose('done computing the norm')
    verbose('debug: type of N %s, k %s'%(type(N),type(k)))
    verbose('Nf = %s'%Nf)
    Fq = normalize(Nf)/normalize(weight_factor(N,k))
    L = Fq.truncate_laurentseries(1).coefficients()
    verbose("length of L is %s"%len(L))
    alist = []
    verbose('number of terms =  %s'%num_terms)
    assert len(L) == num_terms
    # computing the j-invariant
    E4 = eisenstein_series_qexp(4, prec_low+1)
    delta = delta_qexp(prec_low+1)
    j_inv = normalize(E4**3)/delta
    jinvs = [j_inv**0,j_inv]
    while len(jinvs) < num_terms:
        a =  jinvs[-1]
        jinvs.append(a*j_inv)
    verbose('j-invariants computed')
    for i in range(num_terms):
        d = prec_low-1-i
        fd = get_principal_part(jinvs[d],num_terms)
        ad = L[i]/fd[i]
        alist.append(ad)
        L = [L[j]-ad*fd[j] for j in range(num_terms)]
    print 'L = ', L
    T.<x> = QQ[]
    F = T(alist[::-1])
    save(F, os.path.join(os.environ['HOME'],'critical-point','zero-modform','F%s-%s-%s'%(N,k,description)))
    verbose('zero polynomial computed and saved.')
    return F

