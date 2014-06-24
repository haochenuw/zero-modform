#####################################################################################
# for the detailed definition of exceptional points, see the writeup in Math.Works
# basically it's  the preimage of singular points of the plane model Z_0(p) of the modular curve Y_0(p)
# under the rational map Y_0(p) \to {\Phi_p(X,Y) = 0} = Z_0(p).
#####################################################################################


# reference: http://ac.els-cdn.com/S1071579701903193/1-s2.0-S1071579701903193-main.pdf?_tid=3d4e38b4-fbbf-11e3-b320-00000aacb362&acdnat=1403628633_ec0267dd80a8d17066269839dca8817b


# local reference: singularities-of-modular-curve.pdf

load('atkin-lehner.sage')


def exc_points(p):
    '''
    INPUT:
    p -- a prime number
    OUTPUT:
    the full list of exceptional points on Y_0(p)
    '''
    result = []
    for t in range(1, 2*p): # loop over the trace of \alpha
        if t != p:
            D = t^2- 4*p^2
            quadforms =  BinaryQF_reduced_representatives(D,primitive_only = False)
            for f in quadforms:
                result = result + find_exc_points(f,t,p)
    return result

def find_exc_points(f,t,p):
    D = f.discriminant()
    R.<x> = ZZ[]
    K.<tau> = NumberField(R(f.polynomial()(y=1)))
    a = f[0]; b = f[1] # tau = [-b+sqrt(D)]/2a
    sqrtD = 2*a*tau + b
    A = mat(a,b,t,D) # the matrix of alpha acting on [1,tau], where alpha = (t + sqrtD) /2
    assert A.det() == p^2
    assert A.trace() == t

    A = A.change_ring(ZZ)
    D,P,Q = A.smith_form()
    assert d[0][0] == 1
    assert d[1][1] == p^2

    glist = [] # the generators of cyclic p-subgroup of C/Z+Z*tau
    glist.append((P[1][0]+P[1][1]*tau)/p)
    Qinv = ~Q
    f1 = Qinv[0][0] + Qinv[0][1]*tau
    glist.append(f1/p)
    for g in glist:
        result.append(adjust(tau,g,p))
    return result


def adjust(tau,g,p):
    c,d = vector([p*g.matrix()[0][0], p*g.matrix()[0][1]])
    a,b = xgcd(c,d)[1],xgcd(c,d)[2]
    assert a*d - b*c == 1
    return (a+b*tau)/(c+d*tau)






