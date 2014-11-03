def mod_poly(M1,N):
    x = var('x')
    y = var('y')
    result = 0
    nrows = len(M1)
    ncols = len(M1[0])
    for a in range(nrows):
        for b in range(min(ncols,N)): # if u^N appears then we don't care.
            if M1[a][b] != 0:
                result += M1[a][b]*x**a*y**b
    R = PolynomialRing(ZZ,2,'xy',order = 'lex')
    return R(result)


def lift(M,r0,u):
    n = r0.prec()
    u = u.truncate(n+1).add_bigoh(n+1)
    #print 'r0 = ', r0
    f = mod_poly(M,n+1) # trimming all parts where degu >= n+1
    #print 'f = ',f
    K = r0.padded_list()[0].parent()
    T.<b> = K[]
    R = r0.parent()
    RT = R.change_ring(T)
    q = RT.gen()
    r = RT(r0.padded_list()+[b])
    feval = f(x=r,y=u)
    val = feval.valuation()
    assert val >= n # did not mess up the old solution
    fb = feval.padded_list()[val] # a linear function of b
    b0 = fb.roots(multiplicities = False)[0]
    return R(r0.padded_list()+[b0]).add_bigoh(n+1)



def get_series(M,u,r0,N):
    """
    Input -- the first term r0 of the power sereis r.
    """
    r = r0
    for n in range(1,N):
        r = lift(M,r,u)
    return r


def get_lead_terms(M):
    R.<x> = PolynomialRing(ZZ)
    f = R(mod_poly(M,1))
    K.<a> = f.splitting_field()
    fK = f.change_ring(K)
    v = fK.roots(multiplicities = False)
    verbose('minpoly of a is %s'%a.minpoly())
    return v



def get_r_qexps(M,u,maximum = 10):
    N = min(maximum,u.prec())
    """
    Given a polyrelation M between r and u
    and a sufficiently precise q-expansion for u.
    Obtain all expansions of f at big cusps N \mid d^2.
    """
    result = []
    v = get_lead_terms(M)
    print v
    for r0 in v:
        K = r0.parent()
        q = var('q')
        r = K[[q]](r0).add_bigoh(1)
        print r
        result.append(get_series(M,u,r,N))
    return result