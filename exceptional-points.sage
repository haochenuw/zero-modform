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
    for t in corrected_interval(p):
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

# ******* auxiliary functions, to avoid overcount ********
def adjust_i(p):
    assert Mod(p,4) == 1
    K.<i> = QuadraticField(-1)
    a = K.elements_of_norm(p)[0]
    return sorted([abs((a^2).trace()), abs((i*a^2).trace())])

def adjust_rho(p):
    assert Mod(p,3) == 1
    K.<w> = NumberField(x^2+x+1)
    a = K.elements_of_norm(p)[0]
    return sorted([abs((a^2).trace()), abs((w*a^2).trace()),abs((w^2*a^2).trace())])

def corrected_interval(p):
    v = list(range(1,2*p))
    v.remove(p)
    if Mod(p,4) == 1:
        s1 = adjust_i(p)
        verbose('adjustment for D = -4 : %s'%str(s1))
        v.remove(s1[1])

    if Mod(p,3) == 1:
        s2 = adjust_rho(p)
        verbose('adjustment for D = -3 : %s'%str(s2))
        v.remove(s2[1])
        v.remove(s2[2])
    return v


#p = 37
#set_verbose(1)
#z = exc_points(p)






