# various helper functions


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


def newton_polygon(v):
    """
    Input: v — a list of points on ZZ^2.
    Output: the newton-polygon of v: given as a sublist [p1,..,ps]of v, such that the newton polygon is [p1p2] \cup []… \cup [ps-1ps]
    """
    if len(v) <= 1:
        return v
    result = []
    from numpy import argmin
    start = v[0]
    print 'start = %s'%str(start)
    result.append(start)
    slopes = [(pt[1]-start[1])/(pt[0]-start[0]) for pt in v[1:]]
    verbose('slopes =%s'%slopes)
    t = argmin(slopes)
    return result + newton_polygon(v[t+1:])

def newton_polygon_helper(M,degr,degu):
    """
    get all the points for the newton polygon algorithm. Input is a polynomial in 2 variables.
    The valuation here is ord_u. 
    """
    v = []
    for a in range(degr+1):
        for b in range(degu+1):
            if M[a*(degu+1)+b] != 0:
                v.append((a,b))
                break
    a,b = v[-1]
    w = [(a-c,d-b) for c,d in v]
    w.reverse()
    return w