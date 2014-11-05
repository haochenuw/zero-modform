def get_leading_terms(E,F,d,hq):
    """
    E -- the elliptic curve
    d -- the denominator of the cusp in question
    F -- the polynomial relation F(r,h) = 0
    hq -- the second modular function h.
    """
    N = E.conductor()
    assert Mod(d^2,N) == 0
    dprime = ZZ(N//d)

    K = CyclotomicField(dprime,'zeta_%s'%dprime)
    verbose('K = %s'%K)
    T.<b> = K[]
    R.<q> = T[[]]
    rq = R([b])
    Sn = F(r = rq, h = hq) # means modulo q^{n+1}
    vn = Sn.valuation()
    # I think vn must be 0, since F(b,1) = 0 can't have infinitely many solution.
    func = Sn.padded_list()[Sn.valuation()]
    return [K[[q]]([a]).add_bigoh(1) for a in func.roots(multiplicities = False)]


def lifts(F,rq,hq):
    """
    suppose we know rq (mod q^prec)
    want to lift to rq (mod q^prec+1)
    """
    n = rq.prec()
    verbose('the series to be lifted is rq = %s'%rq)
    assert n < hq.prec()
    verbose('lifiting from mod q^%s to mod q^%s'%(n,n+1))
    R = rq.parent()
    K = R.base_ring()
    T.<b> = K[]
    RT = R.change_ring(T)
    rqnew= RT(rq.padded_list()+[b])
    verbose('rqnew = %s'%rqnew)
    Feval = F(r=rqnew,h=hq)
    val = Feval.valuation()
    verbose('Feval = %s'%Feval)
    # did not mess up the old solution
    fb = Feval.padded_list()[val] # a polynomial in one variable b.
    verbose('fb = %s'%fb)
    blist = fb.roots(multiplicities = False)
    lifts = [R(rq.padded_list()+[b0]).add_bigoh(n+1) for b0 in blist]
    verbose('lifting done.')
    verbose('the lifts of rq are %s'%lifts)
    return lifts