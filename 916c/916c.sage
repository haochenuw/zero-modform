

# E = 916c and today it's going to reveal itself to me.


# r = load('r.sobj')
# u = load('u.sobj')

# we assume div_oo(r) = 228 oo and div_oo(u) = 143 oo
# and the leading terms are +1
# and r,u  are two lists storing powers of r and powers of u.






def q_exp(etaElement,prec):
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
    h4wN = sub(RL(q_exp(h4inv,prec)),N//4)/(4**4)
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



def poly_relation(r,u,degr,degu):
    rows = degr + 1
    cols = degu + 1
    M=[]
    for row in xrange(rows): M += [[0]*cols]
    M[degr][0] = -1
    M[0][degu] = 1
    remainder = -r[degr] + u[degu]
    #print type(r[degr])
    #print type(u[degu])
    #print 'remainder = ', remainder
    #print type(remainder)
    while remainder != 0:
        # suppose remainder starts with cq^-d
        d,c  = find_ldgtm(remainder)
        assert d < 0
        # solve x,y such that 228x + 143y = d
        a,b  = solve_dioph(degu,degr,abs(d))
        remainder -= c*r[a]*u[b]
        M[a][b] =  -c
        print 'a,b,c,d = ', a,b,c,d


    save(M, os.path.join(os.environ['HOME'],'poly-relation','results','coefMatrix-916c'))
    return M

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


#def faithful_digits(a,b):
    # how many digits shallwe compute for r^a u^b
    # to make sure the principal part is correct?
    # consider this: r has polar part 228 and u has 143
    # so r^a*u^b has a q^(-C) where C = 228a+143b to it. and we also know C <= 228*143 =M
    # then we only need the normalized part of r and u to have precision q^(C+1) is enough.
    # so we need precision M+1 of both r and u.
    # then forget about the precision of r^a and u^b: they take care of themselves.


degu = 228
degr = 143
padding = 500
prec = padding+(degu+1)*(degr+1)

print('prec = %s'%prec)

try:
    powsOfr = load(os.path.join(os.environ['HOME'],'poly-relation','results','powsOfr-916c'))
    powsOfu = load(os.path.join(os.environ['HOME'],'poly-relation','results','powsOfu-916c'))
except:
    r = r_series_twisted(EllipticCurve('916c1'),prec)
    assert r.valuation() == -228
    save(r, os.path.join(os.environ['HOME'],'poly-relation','results','rnew-916c'))


    v = EtaGroup(916).basis()
    u = v[3]*v[2]*v[4]/v[1];
    assert u.degree() == 143
    print u.divisor()
    u = q_exp(u,prec)
    assert u.valuation() == -143
    save(u, os.path.join(os.environ['HOME'],'poly-relation','results','unew-916c'))

    print('r = %s'%r.truncate(1))
    print('u = %s'%u.truncate(1))

    powsOfr = make_pows(r,degr)
    powsOfu = make_pows(u,degu)


print('computing powers done')

save(powsOfr, os.path.join(os.environ['HOME'],'poly-relation','results','powsOfr-916c'))
save(powsOfu, os.path.join(os.environ['HOME'],'poly-relation','results','powsOfu-916c'))

M = poly_relation(powsOfr,powsOfu,degr,degu)
