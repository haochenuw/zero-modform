def canonical_rep_cusp(k,N):
    """
    require N | (gcd(k,N))^2. i.e. 1/k is a width one cusp.
    Output a pair (d,l), where d is the denominator, (l,N) = 1
    and 1/k = 1/ld
    """
    d = gcd(k,N)
    l = k//d

    dprime = N//d
    l = ZZ(Mod(l,dprime))
    while gcd(l,N) > 1:
        l += dprime
    return (d,l)


def cusps_of_denom_d(d,N):
    """
    same assumption that d \mid N and N | d^2
    """
    result = []
    dprime = N//d
    for i in range(1,dprime+1):
        if gcd(i,dprime) == 1:
            result.append(canonical_rep_cusp(i*d,N))
    return result

def add_cusp(tup1,tup2,N):
    a,b = tup1
    c,d = tup2
    return canonical_rep_cusp(a*b + c*d,N)


def add_orbits(d1,d2,N):
    result = []
    v1,v2 = cusps_of_denom_d(d1,N),cusps_of_denom_d(d2,N)
    for a in v1:
        for b in v2:
            result.append(add_cusp(a,b,N))
    return result

def do_they_commute(p,N,d):
    dprime = ZZ(N/d)
    cusps = [a for a in Gamma0(N).cusps() if a.denominator() == d]
    ms = [alpha(cusp) for cusp in cusps]
    #print ms
    v = hecke_matrices(p,N)
    lst1,lst2 = [],[]
    for a in v:
        for m in ms:
            lst1.append(m*a)
            lst2.append(a*m)
    return same_mod_gamma0n(lst1,lst2,N)

def hecke_matrices(p,N):
    """

    """
    results = [matrix(QQ,[[1,j],[0,p]]) for j in range(p)]
    if Mod(N,p) == 0:
        return results
    else:
        return results + [matrix(QQ,[[p,0],[0,1]])]

def alpha(cusp):
    """
    Input: a cusp c for Gamma_0(N)
    Output: a matrix alpha \in SL_2(\bZ) s.t. alpha(oo) = c

    """
    x = cusp
    if x == Cusp(Infinity):
        return matrix.identity(2)
    else:
        c,d  = x.numerator(), x.denominator()
        _,u,v = xgcd(c,d)
        return matrix([[c,-v],[d,u]])


def same_mod_gamma0n(lst1,lst2,N):
    same = True
    matches = [False for _ in range(len(lst1))]
    for k in range(len(lst1)):
        a = lst1[k]
        #print 'a = ', a
        for b in lst2:
            x = b*(~a)
            if Mod(x[1][0],N) == 0 and all([x[i][j] in ZZ for i in range(2) for j in range(2)]):
                #print 'matching = ', b
                matches[k] = True
    return all(matches)