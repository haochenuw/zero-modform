def get_leading_terms(F,K):
    """
    Input:

    F -- the polynomial relation F(r,h) = 0
    K -- a number field. We are only looking for the solutions in K[[q]].
    Output:

    the solutions to F(r,hq = 0) as power series. Removing all the conjugate duplicates.
    """
    T.<b> = K[]
    R.<q> = T[[]]
    rq = R([b])
    r, u = F.parent().gens()

    Sn = F(r = rq,u = 0) # means modulo q^{n+1}
    vn = Sn.valuation()
    # I think vn must be 0, since F(b,1) = 0 can't have infinitely many solution.
    func = Sn.padded_list()[Sn.valuation()]
    v = func.roots(multiplicities = False)
    v1 = remove_conjugates(K,v)
    return [K[[q]]([a]).add_bigoh(1) for a in v1]

def remove_conjugates(K,v):
    """
    K -- a Galois number field
    v -- a list of elements in K.

    Output:
    v' -- a list of v, where we keep one representative of each conjugacy class.

    Examples::
        sage: K.<zeta3> = CyclotomicField(3)
        sage: v = [zeta3, -zeta3, -zeta3 - 1, zeta3 + 1]
        sage: remove_conjugates(K,v)
        [-zeta3, zeta3]
    """
    vs = Set([])
    for a in v:
        if not a in K:
            a = K(a)
        if Set(a.galois_conjugates(K)).intersection(vs) == Set([]):
            vs = vs.union(Set([a]))
    return list(vs)




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
    r, u = F.parent().gens()
    Feval = F(r=rqnew,u=hq)
    val = Feval.valuation()
    # did not mess up the old solution
    fb = Feval.padded_list()[val] # a polynomial in one variable b.
    blist = fb.roots(multiplicities = False)
    lifts = [R(rq.padded_list()+[b0]).add_bigoh(n+1) for b0 in blist]
    verbose('lifting done.')
    verbose('the lifts of rq are %s'%lifts)
    return lifts