def normalize(g):
    """
    setting the first nonzero coefficient to be 1
    """
    return g/g.padded_list()[g.valuation()]

def q_exp_eta(etaElement,prec):
    R,q = PowerSeriesRing(QQ, 'q').objgen()
    pr = R(1)
    f = etaElement
    v = f._keys
    eta = qexp_eta(R, prec)
    for d in v:
        if f.r(d) != 0:
            pr *= (eta.truncate(prec//d + 1)(q = q**d).add_bigoh(prec))**f.r(d)
    return pr*q**(f._sumDR / ZZ(24))*( R(1).add_bigoh(prec))


def sub(u,N):
    prec = u.prec()
    q = u.parent().gen()
    return u.truncate(prec//N + 1)(q = q^N).add_bigoh(prec)

def r_series_twisted(E,prec):
    #first compute f|w_N
    f = E.modular_form()
    eigwN = f.atkin_lehner_eigenvalue()
    fq = eigwN*f.qexp(prec)


    # then compute u/du where u = 1/j
    def normalize(g):
        """
        setting the first nonzero coefficient to be 1
        """
        return g/g.padded_list()[g.valuation()]


    E4 = eisenstein_series_qexp(4, prec)
    E6 = eisenstein_series_qexp(6, prec)
    delta = delta_qexp(prec)

    #save(f, os.path.join(os.environ['HOME'],'poly-relation','results','modform-%s'%E.label()))

    N = E.conductor()
    R = fq.parent(); q = R.gen()
    u0 = R(delta/normalize(E4)**3)  # u0 = 1/j
    u0wN  = sub(u0,N)
    verbose('1')
    du0 = R(u0.derivative())
    du0wN = N*q^N*sub(du0,N)

    du0wN = du0wN/N # don't want the extra N factor
    # then compute g
    g = sub(delta,N)/sub(delta,N//2)*(2**12)
    g = R(g)
    gnew =  normalize(R((g-512)/(g+256))) # level = 2 , div(gnew) = [i] - [rho]
    verbose('2')
    r = R(-fq*u0wN*gnew/(du0wN))
    RL = R.laurent_series_ring()
    verbose('3')


    # at last compute h

    h4 = EtaProduct(4,{1:-8,4:8})
    h4inv = EtaProduct(4,{1:8,4:-8})
    h4wN = sub(RL(q_exp_eta(h4inv,prec)),N//4)/(4**4)
    verbose('4')
    h = RL(h4wN*32 + 1) # div(h) = wN([(i+1)/2]) - 0
    h *= 8 # to make it integral
    r =  (RL(r)*(h))
    return r

def make_pows(r,degr):
    """
    input:
    r -- a power series
    a -- the power
    output:
    a list  [1,r,..,r^a]
    """
    prec = r.prec()
    R =  r.parent()
    pows_of_r = [R(1)]
    while len(pows_of_r) < degr+1 :
        s = pows_of_r[-1]
        pows_of_r.append(s*r)
    return pows_of_r



def poly_relation(r,u,degr,degu,description):
    rows = degr + 1
    cols = degu + 1
    M=[]
    for row in xrange(rows): M += [[0]*cols]
    M[degr][0] = -1
    M[0][degu] = 1
    _,c = find_ldgtm(r)
    if c== -1:
        r = -r
    rs = make_pows(r,degr)
    us = make_pows(u,degu)


    remainder = (-rs[degr] + us[degu]).truncate(1)

    while remainder != 0:
        # suppose remainder starts with cq^-d
        d,c  = find_ldgtm(remainder)
        if d > 0:
            raise ValueError('got positive valuation, please debug')
        elif d == 0:
            M[0][0] == -c
            break
        else:
            # solve x,y such that 228x + 143y = d
            a,b  = solve_dioph(degu,degr,abs(d))
            verbose('a,b,c,d = %s,%s,%s,%s'%(a,b,c,d))
            remainder -= (c*rs[a]*us[b]).truncate(1)
            M[a][b] =  -c


    save(M, 'results/polyrel-%s'%description)
    return M


def _polyrel(r,u,degr,degu,remainder, Max):
    for i in xrange(Max):
        if remainder ==0 :
            break
        else:
        # suppose remainder starts with cq^-d
            d,c  = find_ldgtm(remainder)
            if d >= 0:
                raise ValueError('got positive valuation, please debug')
            # solve x,y such that 228x + 143y = d
            a,b  = solve_dioph(degu,degr,abs(d))
            verbose('a,b,c,d = %s,%s,%s,%s'%(a,b,-c,d))
            remainder -= (c*r[a]*u[b]).truncate(1)
    return remainder


def solve_dioph(a,b,c):
    """
    knowing that ax+by = c has a unique pair of nonnegative integer solution, find a pair of solution.
    """
    assert gcd(a,b) == 1
    assert c <= a*b
    for y in range(0,a+1):
        if a.divides(c-b*y):
            return ((c-b*y)//a, y)
    return None

def find_ldgtm(f):
    """
    f is a laurent series
    """
    d = f.valuation()
    R,q = f.parent(), f.parent().gen()
    Rp = R.power_series_ring()
    c = Rp(q**(-d)*f).padded_list()[0]
    return (d,c)



def atkin_lehner_eta(etaElement):
    """
    returns the *normalized* result of atkin-lehner involution acting on
    this eta product.
    using eta(-1/z) = sqrt(-iz)eta(z)
    """
    R,q = PowerSeriesRing(QQ, 'q').objgen()
    pr = R(1)
    f = etaElement
    N = f.level()
    v = {}
    for d in f._keys:
        v[N/d] = f.r(d)
    return EtaProduct(N,v)




def yang_product(N,positive = True):
    """
    Input: N — the level
    Output: an eta product u of level N whose divisor is
    D - m[oo], where D >= 0 has Supp(D) = Allcusps - oo.
    i.e., every cusp not equal to oo is a zero of u, and the only pole of u is at infinity. and we want to find such u with minimum degree m.
    """

    p = MixedIntegerLinearProgram(maximization = False)
    b = p.new_variable(integer = True)
    A = divisor_matrix_eta(N)
    number_of_etas = len(EtaGroup(N).basis())
    number_of_cusps = len(Gamma0(N).cusps())

    assert A.dimensions() == (number_of_etas,number_of_cusps)

    if positive:
        minorder = 1
    else:
        minorder = 0
    # we made it so that the first column always correspond to the cusp oo.
    for i in range(number_of_cusps):

        c = A.column(i)
        if i > 0:
            p.add_constraint(sum([b[j]*c[j] for j in range(number_of_etas)]),min = minorder) # zero at all other cusps
        else:
            p.add_constraint(sum([b[j]*c[j] for j in range(number_of_etas)]),max = -1) # Infinity is a pole.
            p.set_objective(-sum([b[j]*c[j] for j in range(number_of_etas)])) # we want to minimize the degree

        p.set_min(b[i],None) # so that b[i] can take negative values


    p.solve(objective_only = False)
    x = p.get_values(b)


    v = EtaGroup(N).basis()
    return prod([v[i]**x[i] for i in range(len(v)) if x[i] >0])/prod([v[i]**(-x[i]) for i in range(len(v)) if x[i] <0])


def divisor_matrix_eta(N):
    """
    return two things: a matrix of divisors of eta generators, and the index of infinity.
    """
    M = []
    for a in EtaGroup(N).basis():
        M.append([a.order_at_cusp(b) for b in AllCusps(N)])
    return Matrix(ZZ,M)









