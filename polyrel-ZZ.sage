
class rFunction():
    """
    The r4 series and related modular functions on X_0(N)
    when 4 \mid N
    """
    def __init__(self,E,n):
        """
        n -- a divisor of the level N.
        """
        if n != 4 and n!=9:
            raise ValueError('n must be 4 or 9.')
        N = E.conductor()
        if Mod(N,n) != 0:
            raise ValueError('N( = %s) must be divisible by %s'%n)
        self.E = E
        self.N = N
        self.n = n

    def poles(self):
        """
        returns a dictionary with keys the cusps and values the order of poles.
        """
        N = self.N
        result = {}
        v = cusp_diagram(N)

        n = self.n
        if n == 4:
            for c in v.keys():
                if c.is_gamma0_equiv(Cusp(oo),4):
                    result[c] = v[c]-1
                else:
                    result[c] = 0
        elif n == 9:
            for c in v.keys():
                if c.is_gamma0_equiv(Cusp(0),9):
                    result[c] = v[c]//9 -1
                elif  c.is_gamma0_equiv(Cusp(oo),9):
                    result[c] = v[c] - 1
                else:
                    result[c] = 0
        else:
            raise NotImplementedError
        return result

    def zeros(self):
        """
        the extra zeros outside the zeros of omega.
        """
        N = self.N
        n = self.n
        result = {}

        if n == 4:
            for c in Gamma0(N).cusps():
                if not c.is_gamma0_equiv(Cusp(oo),4):
                    result[c] = 1
                else:
                    result[c] = 0
        elif n == 9:
            for c in Gamma0(N).cusps():
                if c.is_gamma0_equiv(Cusp(1/3),9) or c.is_gamma0_equiv(Cusp(2/3),9):
                    result[c] = 1
                else:
                    result[c] = 0
        else:
            raise NotImplementedError
        return result

    def degree(self):
        """
        degree of this modular function. This is independent of dprime.
        """
        N = self.N
        deg = sum([b for a,b in self.poles().items()])
        if deg - 2*Gamma0(N).genus() + 2 - sum([b for a,b in self.zeros().items()]) != 0:
            raise ValueError("zeros and poles do not match")
        return deg


    def find_etas(self,dprime):
        """
        finds a pair (h1,h2) of eta products on X0(N) of weigt 0, such that
        (1) h1 has only pole at infinity
        (2) h1*r4 has only pole at oo.
        (2) h2 has only pole at oo, and only zero P_d. (the sum of all cusps of denom d)

        also return the divisor D sucht that div(r1h) = div(w) + D - m[\infty]
        """
        v = self.poles()
        N = self.N
        w = change_keys(v,N)
        h1 = yang_product(N,goal ='majorize',Dmin = w)



        degh1 =  ZZ(h1.degree())

        denom  = N// dprime

        h2 = yang_product(N,goal = 'concentrate',denom = denom)



        degh2 = ZZ(h2.degree())
        verbose('degh1 = %s, degh2 = %s'%(degh1,degh2))
        verbose('divisor of h2 = %s'%h2.divisor())

        if gcd(degh1,degh2) == 1:
            verbose('degrees = %s,%s'%(h1.degree(),h2.degree()))
            return h1, h2
        else:
            verbose('extra work being done...')

            # need some extra work]
            w1 = change_keys(cusp_diagram(N),N)

            g = Gamma0(N).genus()
            c = adjust_to_coprime(degh1,degh2,g+1)
            while True:
                verbose('c = %s'%c)
                if gcd(degh1 + c,degh2) == 1:
                    h3 = yang_product(N,goal = 'prescribe_degree',deg = c)
                    if h3 is not None:
                        h1new = h1*h3
                        verbose('degrees = %s,%s'%(h1new.degree(),h2.degree()))
                        break
                    else:
                        verbose('no eta function of such degree')
                c += 1
        verbose('divisor of h1 = %s'%h1new.divisor())
        return h1new, h2

    def q_exp(self,prec):
        """
        Div(r4) = Z_\omega + (cusps ~0 or 1/2 mod Gamma0(4)) - (pi_4^*(oo) - (cusps ~oo mod Gamma0(4)))


        Only works for 4 | N.
        If power_series = False, return r4 such that Div0r4 = Z_w + cusps~0 + cusps ~1/2
        If True, return Div0r4 = Zw + cusps~1/2 + cusps ~oo

        """
        E = self.E
        N = self.N
        n = self.n

        f = E.modular_form()
        fq = f.qexp(prec)

        E4 = eisenstein_series_qexp(4, prec)
        E6 = eisenstein_series_qexp(6, prec)
        delta = delta_qexp(prec)

        verbose('E4,E6,delta are computed.')


        R = fq.parent(); q = R.gen()
        u0 = R(delta/normalize(E4)**3)  # u0 = 1/j
        du0 = R(u0.derivative())

        # different recipies for different n:
        if n == 4:
            g2 = delta/sub(delta,2)
            gnew =  normalize(((g2-512)/(g2+256)).power_series()) # level = 2 , div(gnew) = [i] - [rho]

            r = R(fq*u0*gnew/(du0*q))
            RL = R.laurent_series_ring()

            h4 = RL(q_exp_eta(EtaProduct(4,{1:8,4:-8}),prec) + 32)
            r =  RL(r)*(h4)
        elif n == 9:
            r =R(fq*u0/(du0*q))
            RL = R.laurent_series_ring()


            m3 =  q_exp_eta(EtaGroup(3).basis()[0],prec)
            g3 =  RL((m3**2 - 486*m3 - 19683)/(m3+243))
            verbose('level 3 adjustment computed')
            m9 = q_exp_eta(EtaGroup(9).basis()[0],prec)
            g9 = RL(1/(m9+9))
            verbose('level 9 adjustment computed')
            r = RL(r)*g3*g9

            # we add this h9 so that the divisor of r9
            # is div(omega) + cusps - pi9^*([0]+[infty]),so it is non-vanishing at oo.

            h9 = q_exp_eta(EtaGroup(9).basis()[1],prec)

            r = r*RL(h9)
            verbose('r computed')
        else:
            raise NotImplementedError

        return r

    def get_expansions(self,dprime, padding =30):
        """
        return the expansions of rh1 and h2
        also return the extra "order" to be subtracted from our list
        """
        h1, h2 = self.find_etas(dprime)

        degh1, degh2 = ZZ(h1.degree()), ZZ(h2.degree())
        verbose('degh1 = %s, degh2 = %s'%(degh1,degh2))

        assert gcd(degh1,degh2) == 1

        prec = (degh1+1)*(degh2+1)+padding

        verbose('prec = %s'%prec)

        rq = self.q_exp(prec)

        R = rq.parent()

        h1q,h2q = R(q_exp_eta(h1,prec)), R(q_exp_eta(h2,prec))

        verbose('expansions of h1 and h2 obtained.')
        rfinal = rq*h1q
        hfinal = h2q

        verbose('final expansions obtained')

        import os
        if not os.path.exists('results'):
            os.makedirs('results')
        save(rfinal,'results/rfinal-%s'%degh1)
        save(hfinal,'results/hfinal-%s'%degh2)

        N = self.N
        c = CuspFamily(N,dprime,'1')
        orderFromr = self.zeros()[Cusp(1/(N//dprime))] - self.poles()[Cusp(1/(N//dprime))]
        orderFromh1 = h1.order_at_cusp(c)


        denom = N//dprime

        verbose('order from r = %s'%orderFromr)
        verbose('order from h1 = %s'%orderFromh1)
        order = orderFromr + orderFromh1


        return rfinal,hfinal,degh1,degh2, order

    @cached_method
    def order(self,dprime = 4):
        """
        return the order of vanishing of omega at cusps of denom dprime.
        """
        r,h,degr,degh,order = self.get_expansions(dprime)
        verbose('expansions obtained. Now computing a relation...')
        P = yang_relation(r,h,degh,degr,'test')
        x, y = P.parent().gens()
        f0 = P(x = 0)
        verbose("f0 = %s"%f0)
        N = self.N
        denom = N//dprime
        sizeOfOrbit = euler_phi(gcd(denom,dprime))
        totalDeg = f0.gcd(y**f0.degree()).degree()
        if Mod(totalDeg, sizeOfOrbit) != 0:
            raise ValueError('wrong')
        return (totalDeg//sizeOfOrbit ) -order



def cusp_diagram(N):
    G = Gamma0(N)
    v = list(G.cusps())
    result = dict((a,0) for a in v)
    for g in G.coset_reps():
        try:
            c = Cusp(g[0][0]/g[1][0])
        except: # division by zero
            c = Cusp(Infinity)
        for a in v:
            if c.is_gamma0_equiv(a,N):
                result[a] +=1
    return result


def non_unitary_divisors(N):
    return [d for d in N.divisors() if gcd(d,N//d) > 1]

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


def normalize(g):
    """
    setting the first nonzero coefficient to be 1
    """
    return g/g.padded_list()[g.valuation()]


def r4_series(E,prec,power_series = False):
    """
    when there's no elliptic point
    Div(r4) = Z_\omega + (cusps ~0 or 1/2 mod Gamma0(4)) - (pi_4^*(oo) - (cusps ~oo mod Gamma0(4)))


    Only works for 4 | N.
    If power_series = False, return r4 such that Div0r4 = Z_w + cusps~0 + cusps ~1/2
    If True, return Div0r4 = Zw + cusps~1/2 + cusps ~oo

    """
    N = E.conductor()
    assert Mod(N,4) == 0
    f = E.modular_form()
    fq = f.qexp(prec)

    E4 = eisenstein_series_qexp(4, prec)
    E6 = eisenstein_series_qexp(6, prec)
    delta = delta_qexp(prec)

    verbose('E4,E6,delta are computed.')


    N = E.conductor()
    R = fq.parent(); q = R.gen()
    u0 = R(delta/normalize(E4)**3)  # u0 = 1/j
    du0 = R(u0.derivative())

    verbose('1')

    g = delta/sub(delta,2)
    gnew =  normalize(((g-512)/(g+256)).power_series()) # level = 2 , div(gnew) = [i] - [rho]

    verbose('2')
    r = R(fq*u0*gnew/(du0*q))
    RL = R.laurent_series_ring()
    if not power_series:
        h4 = RL(q_exp_eta(EtaProduct(4,{1:8,4:-8}),prec) + 32)
        r =  (RL(r)*(h4))
    else:
        h4 = R(1 + 32*q_exp_eta(EtaProduct(4,{1:-8,4:8}),prec))
        r =  (RL(r)*(h4)).power_series()
    return r



def invert_eta(h):
    """
    h  = EtaGroupElement
    return the inverse of h
    EXAMPLES::

        sage: f = EtaGroup(100).basis()[0]; g = invert_eta(f)
        sage: g*f
        Eta product of level 100 : 1
    """
    N = h.level()
    v = N.divisors()
    dic = {}
    for d in v:
        dic[d] = -h.r(d)
    return EtaProduct(N,dic)

def degr4(N):
    """
    the degree of the function r4.
    """
    G = Gamma0(N)
    d = 0
    for c in G.cusps():
        if c.is_gamma0_equiv(Cusp(0),4) or c.is_gamma0_equiv(Cusp(1/2),4):
            d += 1
    return 2*G.genus()-2 + d


def r4_poles(N):
    result = {}
    v = cusp_diagram(N)
    for c in v.keys():
        if c.is_gamma0_equiv(Cusp(oo),4):
            result[c] = v[c]-1
        else:
            result[c] = 0
    return result


def r4_zeros(N):
    result = {}
    for c in Gamma0(N).cusps():
        if not c.is_gamma0_equiv(Cusp(oo),4):
            result[c] = 1
    return result

def r_series_twisted(E,prec):
    #first compute f|w_N
    f = E.modular_form()
    eigwN = f.atkin_lehner_eigenvalue()
    fq = eigwN*f.qexp(prec)



    # then compute u/du where u = 1/j



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



import os

def yang_relation(r,u,degr,degu,description):
    verbose('degr = %s'%degr)
    verbose('degu = %s'%degu)
    rows = degr + 1
    cols = degu + 1
    M=[]
    for row in xrange(rows): M += [[0]*cols]
    M[degr][0] = -1
    M[0][degu] = 1
    _,c = find_ldgtm(r)
    _,d = find_ldgtm(u)
    assert abs(c) == 1
    assert d == 1
    if c== -1:
        r = -r
    rs = make_pows(r,degr)
    us = make_pows(u,degu)


    remainder = (-rs[degr] + us[degu]).truncate(1)

    R.<x,y> = PolynomialRing(ZZ,2)
    F = -x**degr+y^degu
    while remainder != 0:
        # suppose remainder starts with cq^-d
        d,c  = find_ldgtm(remainder)
        if d > 0:
            raise ValueError('got positive valuation, please debug')
        elif d == 0:
            M[0][0] == -c
            F += -c
            break
        else:
            # solve x,y such that degu x + degr y = d
            a,b  = solve_dioph(degu,degr,abs(d))
            verbose('a,b,c,d = %s,%s,%s,%s'%(a,b,c,d))
            remainder -= (c*rs[a]*us[b]).truncate(1)
            M[a][b] =  -c
            F += -c*x**a*y**b


    if not os.path.exists('results'):
        os.makedirs('results')
    save(M, 'results/polyrel-matrix-%s'%description)
    save(F,'results/polyrel-%s'%description)


    return F


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


def r0(E,prec,twisted = False):
    """
    returns fj(j-1728)/j'
    """

    f = E.modular_form().qexp(prec)
    N = E.conductor()
    R = f.parent(); q = R.gen()
    E4 = eisenstein_series_qexp(4, prec)
    E6 = eisenstein_series_qexp(6, prec)
    delta = delta_qexp(prec)
    u = R(delta/normalize(E4)**3)  # u0 = 1/j
    du = R(u.derivative())
    power =1
    if twisted:
        u = sub(u,N).add_bigoh(prec)
        du = sub(du,N).add_bigoh(prec)
        power = N
        R = R.laurent_series_ring()
    r = R(f*(1-1728*u)/du)
    return R(r/q**power)


def valr0(E):
    N = E.conductor()
    G = Gamma0(N)
    return ZZ(G.index() - len(G.cusps()))


def find_h(N):
    """
    find an eta product satisfying (1) only pole is oo,
    (2) order at 0 is > 1/2*deg and coprime to deg.
    """
    h0 = yang_product(N)
    h1 = yang_product(N, goal ='at_zero')
    degh0, degh1 = ZZ(h0.degree()), ZZ(h1.degree())
    ord0h0,ord0h1 = h0.order_at_cusp(CuspFamily(N,N)), h1.order_at_cusp(CuspFamily(N,N)) # cusp 0
    k = 0
    m = 0
    print h0.divisor(), h1.divisor()
    print degh0, degh1,ord0h0,ord0h1
    assert gcd(ord0h0, ord0h1) == 1 or gcd(ord0h0,degh0) == 1
    while True:
        ord0h = ord0h0 + m*ord0h1
        degh = degh0 + m*ord0h1
        if ord0h > 1/ZZ(2)*degh and gcd(ord0h,degh) == 1:
            return h0*h1**m
        else:
            m+=1


def change_keys(dictionary,N):
    """
    change keys from cusp to divisors
    """
    v = dictionary
    result = {}
    for c in v.keys():
        if c != Infinity:
            d = c.denominator()
            if result.has_key(d):
                result[d] = max(v[c],result[d])
            else:
                result[d] = v[c]
        else:
            result[N] = v[c]
    return result


def formalsum_to_dict(formalsum):
    v = list(formalsum)
    result = {}
    for t in v:
        a = t[0]
        c = t[1]
        result[c] = a
    return result



def yang_product(N,goal = 'small_deg',Dmin = None,denom = None, deg = None,exclude = None):
    """
    Input:
        N — the level
        goal -- the type of eta product we want
            'concentrate': m([Pd]-[oo]), d = denom.
            'majorize': an eta product u of level N whose divisor is D - m[oo], where D majorizes Dmin.
            'small_deg': want a eta product of divisor
            'prescribe_degree': want the degree to be deg
            'majorize-prescribe_degree': do both.
            'concentrate_exclude': an extension of concentrate: provide a list of denominators where the eta product can have zero.

        Dmin -- a dictionary representing the divisor we want to majorize, only needed in 'majorize' mode.
    Output:
        The corresponding eta product

    EXAMPLES::
        sage: h = yang_product(8); h.divisor()
        sage: -(Inf) + (0)

    concentrate example::
        sage: h = yang_product(8,goal = 'concentrate',denom = 2); h.divisor()
        sage: (c_{4}) - (Inf)

        sage: h = yang_product(32,goal = 'concentrate',denom = 4); h.divisor()
        sage: -2*(Inf) + (c_{8,2}) + (c_{8,1})

    concentrate_exclude::
        sage: yang_product(162,goal = 'concentrate_exclude',denom = 18, exclude = []).divisor()
        sage: (c_{9,6}) + (c_{27,1}) + (c_{27,2}) + (c_{81}) + (c_{9,4}) + (c_{9,1}) - 9*(Inf) + (c_{9,5}) + (c_{9,2}) + (c_{9,3})

    majorize-prescribe_degree::
        sage: yang_product(162,goal = 'majorize-prescribe_degree',Dmin = {1:1,27:2},deg = 24).divisor()
        sage: 3*(c_{54,1}) + 3*(c_{54,2}) + 3*(0) + 3*(c_{6,1}) + 3*(c_{3,2}) + 3*(c_{2}) + 3*(c_{3,1}) + 3*(c_{6,2}) - 24*(Inf)
    """

    p = MixedIntegerLinearProgram(maximization = False)
    b = p.new_variable(integer = True)
    A = divisor_matrix_eta(N)
    number_of_etas = len(EtaGroup(N).basis())
    number_of_cusps = len(Gamma0(N).cusps())

    assert A.dimensions() == (number_of_etas,number_of_cusps)

    verbose('goal = %s'%goal)
    verbose('Dmin = %s'%Dmin)
    # we made it so that the first column always correspond to the cusp oo.
    for i in range(number_of_cusps):
        c = A.column(i)
        if i > 0:
            p.add_constraint(sum([b[j]*c[j] for j in range(number_of_etas)]),min = 0)

            if goal == 'concentrate':
                f = AllCusps(N)[i]
                d =  f.level()//f.width()
                if d != denom:
                    p.add_constraint(sum([b[j]*c[j] for j in range(number_of_etas)]),max = 0)

            elif goal == 'concentrate_exclude':
                if exclude is None:
                    exclude = non_unitary_divisors(N)
                f = AllCusps(N)[i]
                d =  f.level()//f.width()
                if d == denom:
                    p.add_constraint(sum([b[j]*c[j] for j in range(number_of_etas)]),min = 1)
                elif d in exclude:
                    p.add_constraint(sum([b[j]*c[j] for j in range(number_of_etas)]),max = 0)
            elif goal == 'majorize' or goal == 'majorize-prescribe_degree':
                f = AllCusps(N)[i]
                d =  f.level()//f.width() # the denominator of this cusp. Note that this width is not the usual definition of width
                try:
                    dmin = Dmin[d]
                except:
                    dmin = 0
                p.add_constraint(sum([b[j]*c[j] for j in range(number_of_etas)]),min = dmin) # zero at all other cusps
            else:
                if not goal == 'small_deg' and not goal == 'prescribe_degree':
                    raise ValueError('invalid goal %s'%goal)
        else: # the cusp oo
            p.add_constraint(sum([b[j]*c[j] for j in range(number_of_etas)]),max = -1) # Infinity is a pole.
            p.set_objective(-sum([b[j]*c[j] for j in range(number_of_etas)])) # we want to minimize the degree
            if goal == 'prescribe_degree' or goal == 'majorize-prescribe_degree':
                p.add_constraint(sum([b[j]*c[j] for j in range(number_of_etas)]),max = -deg)
                p.add_constraint(sum([b[j]*c[j] for j in range(number_of_etas)]),min = -deg)
        p.set_min(b[i],None) # so that b[i] can take negative values

    try:
        p.solve(objective_only = False)
        x = p.get_values(b)
        verbose('x = %s'%x)
    except:
        verbose('no solution')
        return None
    v = EtaGroup(N).basis()
    one = EtaProduct(N,{})
    d = [invert_eta(v[i])**(-x[i]) for i in range(len(v)) if x[i] <0] +[one]
    n = [v[i]**x[i] for i in range(len(v)) if x[i] >0]+ [one]

    return prod(n)*prod(d)



def divisor_matrix_eta(N):
    """
    return two things: a matrix of divisors of eta generators, and the index of infinity.
    """
    M = []
    for a in EtaGroup(N).basis():
        M.append([a.order_at_cusp(b) for b in AllCusps(N)])
    return Matrix(ZZ,M)



def find_etas(N):
    """
    finds a pair (hr,h) of eta products, such that
    (1) hr has only pole at infinity
    (2) hr*r4 has only pole at oo.
    (2) h has only pole at oo, and the degree of h is the smallest among such

    """

    assert Mod(N,4) == 0
    v = r4_poles(N)
    w = change_keys(v,N)
    hr = yang_product(N,goal ='majorize',Dmin = w)
    valh1 =  ZZ(hr.degree())
    h = yang_product(N)
    valh2 = ZZ(h.degree())
    print valh1,valh2

    return hr, h


def get_expansions(E,hr,h,padding = 30):
    """
    given a pair of functions relatively small degree
    such that zeros of r contains zeros of omega + some cusps
    and h is a small degree eta product.
    """

    valr, valh = ZZ(hr.degree()), ZZ(h.degree())
    verbose('valr = %s, varh = %s'%(valr,valh))
    prec = (valr+1)*(valh+1)+padding
    verbose('prec = %s'%prec)
    r = r4_series(E,prec)
    hrq = q_exp_eta(hr,prec)
    hq = q_exp_eta(h,prec)


    rfinal = r*hrq
    hfinal = hq
    save(rfinal,'r-%s'%valr)
    save(hfinal,'h-%s'%valh)
    return r*hrq, hq



# have (a,b) > 1. Want to find congruence conditions on c >= M such that (a+c,b) = 1.
# using the fact that if c = 1 mod b then (c,b) = 1.
def adjust_to_coprime(a,b,M):
    assert b>0
    for i in range(M,M+b):
        if gcd(a+i,b) == 1:
            return i
    return None


