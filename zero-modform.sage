from sage.matrix.matrix_integer_dense import _lift_crt



def zero_poly(f,p,k,description):
    """
    Input: f -- the power series expansion of a modular form
           k -- the weight of f
    Output: the polynomial satified by the j-invariants of f(z)(dz)^k/2
    """
    v = precs(p,k)
    try:
        f = f.qexp(v[2])
    except:
        pass
    Nf = norm_multimodular(f,p,k)
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
    print 'L = ', L
    T.<x> = QQ[]
    F = T(alist[::-1])
    save(F, os.path.join(os.environ['HOME'],'critical-point','zero-modform','F%s-%s-%s'%(p,k,description)))
    return F

def normalize(g):
    return g/g.padded_list()[g.valuation()]

def weight_index(p,k):
    """
    return the tuple (A,B,C) in the writeup in misc math project/zero-polynomial
    of modform. Here p = level, k = weight.
    """
    B = (1 + kronecker(-3,p))*k
    C = (1 + kronecker(-1,p))*k/2
    A = (k*(p+1)-4*B-6*C)/12
    try:
        return (ZZ(A),ZZ(B),ZZ(C))
    except:
        return None

# this weight_index function is one of the main functions
# that should be revised if we want to generalize this to square free level.



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
    return (k*(g-1)+ 1 , k*g+1 + 1, p*(k*g+1 + 1))

def coef_bound(p,k):
    """
    bound the largest possible coefficient of
    a product of p modular forms with degree prec, where if f = \sum a_nq^n
    then each term in the product is f = \sum a_n b_n q^n
    with |b_n| = 1 for all n.
    """
    prec= precs(p,k)[1] # the medium precision, needed for the norm of f
    return RR(binomial(prec + p-1, p-1)*(prec^(k*p/2)))


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

def get_principal_part(f,n):
    """
    returns the principal part + constant part of the
    coefficients of a Laurent series f
    if the length l <n, put n-l zeros in there.
    """
    s = f.truncate_laurentseries(1).coefficients()
    return  [0]*(n-len(s)) + s



def norm_mod_l(f,p,z):
    """
    Take a modular form f of level  p with integer coefficents,
    this function computes \prod_i (f(zeta_p^i q^1/p))
    modulo a prime l that is congruent to 1 modulo p.

    Input: f -- a modular form with integer coefficients
    z -- a root of \Phi_p(mod l) in GF(l)
    """

    prec_big = f.prec()
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
    new_prec = ZZ(prec_big//p)
    verbose("The multiplication for the prime %s is performed within %s seconds"%(l, cputime(t)))
    return F.padded_list()[0::p]


def norm_multimodular(f,p,k):
    """
    Given a modular form of weight k and level p.
    sing multimodular algorithm(computing mod primes, and use CRT to lift back)
    to compute the norm of a modular form f of level p
    ASSUMPTION: f has integer coefficients
    """
    bigN = coef_bound(p,k)
    v = precs(p,k)
    prec_big = v[2]
    phip = cyclotomic_polynomial(p)
    verbose('the bound on the size of coefficents of the norm is %s'%bigN)
    llist = mod1_primes(p,2*bigN) #multiply by 2 since we are inside the interval [-bigN, bigN]
    verbose('we are using %s primes'%len(llist))
    Matlist = []
    prec_med = v[1]
    count = 0
    for l in llist:
        phipl = phip.change_ring(GF(l))
        v = phipl.roots(multiplicities = False)
        v.sort()
        z = v[0]
        Nfbarl = norm_mod_l(f,p,z)
        Matlist.append(Matrix(GF(l),1,prec_med,Nfbarl[:prec_med]))
        count += 1
        verbose('computation done for the %s-th prime'%(count))
    M = _lift_crt(Matrix(ZZ, 1, prec_med), Matlist)
    R.<q> = QQ[[]]
    return R(M.list())*f.truncate(prec_med)









