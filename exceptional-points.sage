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
    result = []
    D = f.discriminant()
    R.<x> = QQ[]
    a, b = f[0],f[1]
    K.<tau> = NumberField(R(f.polynomial()(y=1)/a))
    # tau = [-b+sqrt(D)]/2a
    sqrtD = 2*a*tau + b
    A = mat(a,b,t,D) # the matrix of alpha acting on [1,tau], where alpha = (t + sqrtD) /2
    verbose('A = %s'%str(A))
    assert A.det() == p^2
    assert A.trace() == t

    A = A.change_ring(ZZ)
    D,P,Q = A.smith_form()
    verbose('det(D) = %s'%D.det())
    assert D[0][0] == 1

    glist = [] # the generators of cyclic p-subgroup of C/Z+Z*tau
    glist.append((P[1][0]+P[1][1]*tau)/p)
    Qinv = ~Q
    f1 = Qinv[0][0] + Qinv[0][1]*tau
    glist.append(f1/p)
    for g in glist:
        result.append(adjust(tau,g,p))
    return result

def mat(a,b,t,D):
    """
    return the matrix of (t+sqrtD)/2 acting on column vector [1, (-b+sqrtD)/2a]
    """
    return matrix([[(t+b)/2,a],[(D-b^2)/(4*a),(t-b)/2]])



def adjust(tau,g,p):
    c,d = vector([p*g.matrix()[0][0], p*g.matrix()[0][1]])
    a,b,c,d = lift_to_sl2z(c,d,0)
    assert a*d - b*c == 1
    return (a+b*tau)/(c+d*tau)

#p = 37
#set_verbose(1)
#z = exc_points(p)






