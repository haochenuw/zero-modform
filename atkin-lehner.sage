from sage.structure.sage_object import SageObject

##################################################################
#             Hao Chen (chenh123@uw.edu)                         #
##################################################################


#################################
# part 0: chow-heenger index    #
#################################


def chow_heegner_index(lab1,lab2):
    """
    returns an integer m = 2**n
    dividing the index of P_{E,F}.
    When P_{E,F} is torsion, sometimes
    prints a message, sometimes not.
    """
    E = EllipticCurve(lab1)
    F = EllipticCurve(lab2)
    phi = E.modular_parametrization()
    psi = F.modular_parametrization()
    N,N1 = E.conductor(),F.conductor()
    if N1 != N:
        raise NotImplementedError
    qlist0 = []
    for d in N.divisors():
        if d != 1 and gcd(d,N/d) == 1:
            qlist0.append(d)
    #print 'list of possible q\'s:', str(qlist0)
    wbarlist = []
    greedypoints = set([])
    for q in qlist0:
        wq = atkin_lehner(q,N)
        _,signE,PE = wq.induced_map(E)
        #print 'induced map of wq on E: q = ' , q, signE,PE
        wqbar,signF,PF = wq.induced_map(F)
        #print 'induced map of wq on F: q = ' , q, signF,PF
        if signE == -1:
            return 'wq for q = ' + str(q)+'acts as -1, so P_E,F is torsion point'
        elif PE != (0,0):
            pass
        elif signF == 1 and PF != (0,0):
            pass
        else:
            wbarlist.append((wqbar,q))
            if signF == -1:
                greedypoints.add((PF[0]/2,PF[1]/2))

    #print 'greedypoints ', greedypoints
    greedypoints = list(greedypoints)

    if greedypoints == []:
        #if all wbar acts as identity, then we
        #can choose any point, al will permute the fibers
        #hence we could avoid all the fixed points
        return 2**(len(wbarlist))
    results = []


    dict1 = period_mapping_conv(F)
    v1, v2 = dict1.keys()
    A = matrix([list(v1),list(v2)])
    #print 'A = ' , A
    w1, w2 = dict1[v1], dict1[v2]
    for gp in greedypoints:
        wlist2 = [atkin_lehner(q,N) for q in al_fixing(gp,wbarlist)] #all atkin-lehner s.t. wqbar(gp) = gp.
        #print 'gp ', gp
        #print 'all wq fixing gp ', [wq.q for wq in wlist2]
        W = 1+len(wlist2)

        gp = vector(gp)

        a, b = gp*(A**-1)
        gpe  = F.elliptic_exponential(a*w1+b*w2)

        v = fixed_points_above_p(wlist2,psi,gpe)
        #print 'bad points', v
        if len(v) > 0:
           results.append(W/2)
        else:
            #free action
            results.append(W)
    return min(results)


def al_fixing(gp,wbarlist):
    """
    given a greedy point on F, give back all the possible q's such that
    wqbar fixes gp.
    """
    result  = []
    for wbar,q in wbarlist:
        #print gp, wbar(gp)
        if wbar(gp) == gp:
            result.append(q)
    return result



def fixed_points_above_p(wlist,psi,P,eps = 1e-7):
   """
   Given a point P on Elliptic Curve F,
   returns all atkin-lehner fixed points on X0(N)
   that maps to P.
   Input:
       wlist -- a list of atkin-lehner operators
       that induces a map on F that fixes P
       psi -- modular parametrization X0(N) -> F
       P -- a point on F, represented as (x:y:z)
   """
   F = psi.curve()
   v = []
   dict
   for wq in wlist:
       for Q in wq.fixed_points_h():
           #print psi(Q)[0], P[0]
           #psi(Q) is of the form (x:y:z)
           if abs(psi(Q)[0] - P[0])<eps:
               v.append((wq.q, Q))
   return v

##################################
# part I: period mapping methods #
##################################

def element_gamma0N(N,p):
    """
    returns an element
    [a,b,N,p] in Gamma0(N)
    """
    g, x, y = xgcd(p, N)
    if g != 1:
        raise ValueError, "p must be coprime to N"
        # Now p*x + N*y = 1
    return matrix([[x,-y],[N,p]])

def periodg(f,gg,terms):
    """
    given a newform f attached to elliptice
    curve E, an element gg in Gamma0(N),
    computes the period
    2*pi*i*\int_{a}^{gg(a)} f(z) dz
    by truncating the q-expansion of f to a polynomial
    with degree = terms, a is the optimal point in [Cremona97].
    """
    f = f.q_expansion(terms)
    ff = f.shift(-1).integral()
    a=gg[0,0]
    b=gg[0,1]
    c=gg[1,0]
    d=gg[1,1]
    tau=-(d/c)+(CC(I)/abs(c))
    gammatau=(a*tau+b)/(c*tau+d)
    ffp = ff.polynomial()
    return ffp.substitute(exp(2*CC(pi)*CC(I)*gammatau)) - ffp.substitute(exp(2*CC(pi)*CC(I)*tau))

def period_mapping_conv(E):
    """
    returns a dictionary {v1 = (a1,b1): w1,v2 = (a2,b2): w2}
    where the vectors vi are images of modular symbols
    gamma_i under integral_period_mapping, and
    wi are the period attached to gamma_i wrt E.
    """
    N = E.conductor()
    phi = E.modular_symbol_space(0).integral_period_mapping()
    p = 1
    f = E.modular_form()
    result = {}
    while True:
        if gcd(N,p) == 1:
            print 'p=', p
            g = element_gamma0N(N,p)
            c = [0,g[0,1]/g[1,1]]
            print 'c = ', c
            v = tuple(phi(c))
            print 'v = ', v
            w = periodg(f,g,5000) # here the number of terms is chosen to be 5000, otherwise it can get
            # so slow that it's not practical to compute. If it's 500 then it's not precise enough.
            if not close_to_0(w):
                print 'w= ', w
                if result == {}:
                    if v != (0,0):
                        result[v] = w
                else:
                    v0, w0 = result.items()[0]
                    if v0[0]*v[1]-v0[1]*v[0] != 0:
                        result[v] = w
                        return result
        p = p+1
    print 'Failed to Compute'
    return None


def period_mapping_conv_fixed(E):
    """
    returns a dictionary {v1 = (a1,b1): w1,v2 = (a2,b2): w2}
    where the vectors vi are images of modular symbols
    gamma_i under integral_period_mapping, and
    wi are the period attached to gamma_i wrt E.
    """
    N = E.conductor()
    phi = E.modular_symbol_space(0).integral_period_mapping()
    p = 1
    f = E.modular_form()
    result = {}
    for g in Gamma0(N).gens():
        #print g, g[0,1],g[1,1]
        c = [0,g[0,1]/g[1,1]]
        #print 'c = ', c
        v = tuple(phi(c))
        #print 'v = ', v
        if g[1,0] != 0:
            w = periodg(f,g,5000) # here the number of terms is chosen to be 5000, otherwise it can get
            # so slow that it's not practical to compute. If it's 500 then it's not precise enough.
            if not close_to_0(w):
                #print 'w= ', w
                if result == {}:
                    if v != (0,0):
                        result[v] = w
                else:
                    v0, w0 = result.items()[0]
                    if v0[0]*v[1]-v0[1]*v[0] != 0:
                        result[v] = w
                        return result
    print 'Failed to Compute'
    return None



def close_to_0(z,eps = 1e-7):
    return abs(z) < eps


##################################
# part II: elliptic curve methods#
##################################


def image_of_cusp(E,c):
	"""
	returns the image of c under the
	modular parametrization associated to E
	and its order?
	"""
        N = E.conductor()
        M = ModularSymbols(N)
        phi = E.modular_symbol_space(0).integral_period_mapping()

	return canonicalize_Z2(phi(M([Infinity,c])))

def canonicalize_Z2(v):
	a,b = v
	return (QQ(a - a.floor()), QQ(b-b.floor()))


def cusp_over_torsion(E,P):
    """
    P is a torsion point on elliptic curve E, return
    all the inequivalent cusps that maps to P
    """
    clist = Gamma0(E.conductor()).cusps()
    v = [(c,image_of_cusp(E,c)) for c in clist]
    result = []
    for c,Q in v:
        if are_same(P,Q):
            result.append(c)
    return result

def are_same(p,q):
    return p[0] - q[0] in ZZ and p[1] - q[1] in ZZ


####################################
# part II: supporting functions    #
####################################

def canonical_rep(c,N):
    """
    returns the canonical representation of a cusp c
    """
    if not isinstance(c, Cusp):
        c = Cusp(c)
    for other in Gamma0(N).cusps():
        if c.is_gamma0_equiv(other,N):
            return other

def find_min(a,b,c,d,q,N):
    """
    |c+ d*tau|, where tau = a+b*sqrt(-q)

    finds Lc, Ld (mod N) so that
    minimizes|Lc+ Ld( a+b(sqrt(-q)))|
    and (Lc,Ld,N) = 1
    """
    def leng(c1,d1):
        return (a*d1+c1)**2 + ((b*d1)**2)*q
    v =  []
    for L in range(1,N):
        if gcd(gcd(L*c,L*d),N) == 1:
            c1, d1  = ZZ(Mod(L*c,N)), ZZ(Mod(L*d,N))
            v.append((c1,d1,leng(c1,d1)))
    f = lambda x: x[2]
    return min(v,key = f)[:2]


def cyclic_subgroup(d,kernel):
    """
    Input:
        kernel: a tuple representation of a
        finite abelian group G that's either
        cyclic or of form Z/2 + Z/2k,

        d: an integer dividing |G|

    Output:
        a list of generators for all cyclic subgroups
        of G with order d.
    """
    #print "computing cyclic subgroup of order", d , " for a group of form ", kernel
    if len(kernel) == 1:
        gen, order = kernel[0]
        return [gen*order/d]
    else:
        #print kernel[0], kernel[1]
        gen1,d1 = kernel[0]
        gen2,d2 = kernel[1]
        if Mod(d,2) != 0:
            return [gen2*d2/d]
        elif Mod(d,4) != 0:
            return [gen2*d2/d, gen1+gen2*d2/d, gen1+2*gen2*d2/d]
        else:
            if Mod(d2,d) == 0:
                return [gen2*d2/d, gen1+gen2*d2/d]
            else:
                #no solution: the power of 2 in d is too big
                return []


def solve_congruence(q,N):
    """
    return all x mod(N/q) s.t.
    N/q | x^2 +q
    """
    d = ZZ(N/q)
    R = IntegerModRing(d)
    x = PolynomialRing(R,'x').gen()
    f = x**2 + q
    return f.roots(multiplicities = False)


def kernel(A,tau):
    """
    return a basis for the kernel of alpha
    acting on [1,tau] by A
    here we get tau from the binary qf
    and alpha = either sqrt(-q) or sqrt(-q)-q
    """
    A = A.change_ring(ZZ)
    D,P,Q = A.smith_form()
    # Now  PAQ  =  D
    d1 = D[0][0]
    d2 = D[1][1]
    v1 = P[0][0]+P[0][1]*tau
    v2 = P[1][0]+P[1][1]*tau
    if d1 == 1:
        return [(v2/d2,d2)]
    else:
        return [(v1/d1, d1),(v2/d2,d2)]

###################################
# part III: Atkin-lehner class    #
###################################

# This is the main thing: Atkin-Lehner class
class atkin_lehner(object):
    def __init__(self, q,N):
        """
        returns the atkin-lehner matrix w_q of level N
        """
        g, x, y = xgcd(q, -N//q)
        if g != 1:
            raise ValueError, "q must exactly divide N"
        # Now q*x - (N//q)*y = 1
        else:
            self.q = q
            self.N = N
            self.matrix = matrix([[q*x, y],[N,q]])

    def __call__(self,c):
        row1,row2 = self.matrix
        if isinstance(c,Cusp):
            N = self.N
            if c != Infinity:
                return canonical_rep( Cusp((row1[0]*QQ(c)+row1[1])/(row2[0]*QQ(c)+row2[1])),N )
            else:
                return canonical_rep( Cusp(row1[0]/row2[0]), N )
        else:
            return (row1[0]*c+row1[1])/(row2[0]*c+row2[1])

    def sign(self,E):
        """
        return the 1 or -1 eigenvalue of this operator acting
        on the elliptic curve E
        """
        if E.conductor() != self.N:
            raise NotImplementedError, "conductor must equal to level"
        return E.modular_form().atkin_lehner_eigenvalue(self.q)

    def induced_map(self,E):
        """
        return the induced involution bar{w_q} on the
        elliptic curve E. It's given by this formula

        bar{w_q}(x) = phi(w_q(Infinity)) + eps_q * x
        wher eps is the output of the sign function above.

        """
        c = self(Cusp(Infinity))
        P = image_of_cusp(E,c)
        def _f(Q):
            """
            return the image of Q under the induced map
            everything represented in QQ^2.
            """
            return canonicalize_Z2((P[0]+ self.sign(E)*Q[0],P[1]+ self.sign(E)*Q[1]))
        return (_f, self.sign(E), P)


    def number_of_fixed_points(self):
        #needs to be fixed
        q = self.q
        N = self.N
        M = ModularSymbols(N,sign=1)
        S = M.cuspidal_submodule()
        w = S.atkin_lehner_operator(q).matrix()
        g1 = (w-1).kernel().dimension()
        return 2*(Gamma0(N).genus())-4*g1+2

    def fixed_points_h(self):
        """
        return the fixed points of self acting on Y_0(N)
        """
        q = self.q
        N = self.N
        return fixed_points_w(q,N)




######################################
# part IV: atkin-lehner fixed points #
######################################


class AtkinLehnerFixedPoint():
    """
    The class of atkin-lehner fixed point
    """
    def __init__(self, N, q, D, z):
        """
        INPUT:

            - `N` -- (positive integer) the level

            - `q` -- positive integer such that (q,N/q) = 1

            - `D` -- negative integer, the discriminant of the corresponding CM point. Either -4q or -q.

            - `z` -- an imaginary quadratic point in the upper half plane

        """
        self._N = N
        self._q = q
        self._D = D
        self._z = z

    def z(self):
        return self._z

    def q(self):
        return self._q

    def disc(self):
        return self._D

    def N(self):
        return self._N

    def _repr_(self):
        q = self.q()
        N = self.N()
        alpha = 'sqrt(-'+str(q) + ')'

        tau = repr(self.z()).replace('alpha',alpha)
        return 'atkin-lehner fixed point  %s  of w%s on X0(%s)'%(tau,q,N)

    def str(self):
        return self._repr_()

def fixed_points_w(q,N):
    """
    returns all fixed points of w_q on Y0(N)
    of the form  a + b*alpha, a,b in Q, alpha
    lies in some imaginary quadratic field
    """
    result = []
    xlist = solve_congruence(q,N)
    if xlist  == []:
        return []
    else:
        result += find_points(xlist,q,N,2)
        if q == 2:
            result += find_points_2(N)
        elif q == 3:
            result += find_points_3_fixed(N)
        elif Mod(q,4) == 3:
            result += find_points(xlist,q,N,1)
    return result

# The find_points function is the only original algorithm written by Hao Chen.


def find_points(xlist,q,N,cond):
    """
    returns all the possible fixed point of atkin-lehner operators
    w_q acting on X_0(N) with CM ring of conductor cond.

    Input:
        -xlist: a list of integers x s.t. N/q | x^2 +q

        -q,N: for w_q and X_0(N)

        -cond: The curves with  CM ring ZZ[sqrt(-q)] when cond = 2
          CM ring ZZ[sqrt(-q)+1 /2] when cond = 1 and q = 3 mod 4

    Output:
        the full list of fixed points of w_q
        represented as elements in the imaginary quadratic field K = Q(sqrt(-q)) and
        in the upper half plane. Use CDF() to convert to complex points.

    Argument following [Kenku]
    """
    try:
        formlist =  BinaryQF_reduced_representatives(-q*(cond**2),primitive_only = True)
    except:
        return []
    result = []
    d =  ZZ(N/q)
    K.<alpha> = QuadraticField(-q);
    """
    if cond == 2:
        print ' alpha  := sqrt(-'+str(q) + ')'
    elif cond == 1:
        print ' alpha  := (sqrt(-'+str(q) + ') + 1)/2'
    """

    A = alpha.matrix()
    for f in formlist:
        a = f[0]
        b= f[1]
        tau = (-b+cond*alpha)/(2*a)
        #print 'lattice: [1, ' + str(tau) + ']'
        B_tau = matrix([[1,0],[-b/(2*a),cond/(2*a)]])
        A_tau = B_tau * A * B_tau^(-1)
        A_tau = A_tau.change_ring(ZZ)
        #print 'matrix of alpha acting on [1, tau]: ', A_tau
        #print 'kernel of alpha: ', kernel(A_tau,tau)
        g1  = kernel(A_tau,tau)[0][0]
        #print 'g1 = ', g1
        for x in xlist:
            #print 'x = ', x
            # compute the kernel of (alpha - x)
            A_x = (alpha-ZZ(x)).matrix()
            A_xtau = B_tau * A_x * B_tau^(-1)
            A_xtau = A_xtau.change_ring(ZZ)
            #print A_xtau
            #print 'kernel of alpha - ' + str(x), kernel(A_xtau,tau)
            #gen,order = kernel(A_xtau,tau)[0]
            #print gen, order
            for g2 in cyclic_subgroup(d,kernel(A_xtau,tau)):
            #g2 = gen*(order/d) #g2 generates a cyclic group of order d.
                #print 'g2 = ', g2
                g = g1+g2
                #print 'g = ', g
                #print 'Ng = , ', N*g
                coord = vector([N*g.matrix()[0][0], N*g.matrix()[0][1]])
                #print coord
                #Now the point on X_0(N) we want is (tau,<g>)
                #just need to put it in standard form.
                #suppose Ng = m+n*tau, (c,d) = 1.
                #
                #Then choose a,b so that a+btau, c+dtau is another basis for [1,tau]
                #then we have ([a+btau, c+dtau], (c+dtau) /N)
                #so the point should be [a+btau/c+dtau]
                m,n = coord*B_tau^(-1)
                m, n = m/gcd(m,n), n/gcd(m,n)
                #print 'm, n = ', m, n
                #print 'a,b,cond = ', a,b,cond
                m,n = find_min(-b/(2*a),cond/(2*a),m,n,q,N)
                #print 'new m, n = ', m, n
                a1,b1,c1,d1 = lift_to_sl2z(m,n,0)
                #print a,b,c1,d1
                point = -(a1+b1*tau)/(c1+d1*tau)
                #print 'point ', point
                result.append(AtkinLehnerFixedPoint(N,q,-q*cond**2,point))
    return result


def find_points_2(N):
    """
    special case: N = 2d, d odd
    find the extra fixed points
    given by CM ring Z[i]
    the only form is x^2+y^2,
    corresponding to the point i
    """
    result = []
    d =  ZZ(N/2)
    K.<i> = QuadraticField(-1)
    alpha =  1 + i
    # print ' alpha := 1+i '
    A = alpha.matrix()
    gen,order = kernel(A,i)[0]
    #print gen, order
    g1 = gen
    #print ' g1 ', g1
    #print A
    xlist = solve_congruence(1,N)
    #Note this gives all x such that  N | x^2+1,
    #but what I want is N | deg(alpha-x) = (x-1)^2 +1
    xlist = [x + 1 for x in xlist]
    for x in xlist:
        A_x = (alpha - ZZ(x)). matrix()
        #print'x = ', x
        gen,order = kernel(A_x,i)[0]
        #print gen,order
        g2 = gen*(order/d)
        #print g2
        g = g1 + g2
        #print 'Ng', N*g
        m, n = N*g.matrix()[0][0], N*g.matrix()[0][1]
        m, n = m/gcd(m,n), n/gcd(m,n)
        #print 'm, n = ', m, n
        #m, n =  Mod(m,N)* Mod(n,N)^-1, 1
        m,n = find_min(0,1,m,n,2,N)
        #print 'new m, n = ', m, n
        a,b,c1,d1 = lift_to_sl2z(m,n,0)
        #print a,b,c1,d1
        point = -(a+b*i)/(c1+d1*i)
        #print 'point ', point
        result.append(point)
    return result

###################################################
#part IV: Gamma0(N)-equivalence in imag quad field#
###################################################

##############
# This part modifies codes written by William Stein in chow-heegner.py
# so that it works for imaginary quadratic field with provable
# correctedness (no precision issues).
########

def sl2z_rep_in_fundom_q(z):
    """
    z is in QQ(-q),
    """
    a,b = z.matrix()[0]
    # now a = a+b*sqrt(-q)
    gamma = SL2Z(1)
    S, T = SL2Z.gens()
    change = True
    while change:
        change = False
        t = z.parts()[0]
        if abs(t) > 1/2:
            change = True
            # |t - n| <= 1/2
            # -1/2 <= t-n <= 1/2
            # n - 1/2 <= t < = 1/2+n
            # n <= t + 1/2 <= n + 1, so n = floor(t+1/2)
            n = (t + 1/2).floor()  # avoid rounding error with 0.5
            z -= n
            gamma *= T**n
        if QQ(abs(z)**2) < 1:
            change = True
            z = -1/z
            gamma *= S
    return z, gamma**(-1)

def canonicalize_sl2z_q(a, g=None):
    """
    Assume that a = g(z) is in the fundamental domain for SL2Z.
    Adjust a by applying T^(-1) or S so that a is the canonical
    representative in the fundamental domain, so a is not on the right
    edge, and if a is on the unit circle, then it is on the left hand
    side.  Also, modify g so that the relation a = g(z) continues to
    hold.

    Special case: if g is None, just adjust a, ignoring g.
    """
    if g is not None:
        S, T = g.parent().gens()
    if a.parts()[0] == 1/2:
        a -= 1
        if g is not None: g = T**(-1)*g
    elif QQ(abs(a)**2) == 1 and a.parts()[0] > 0:
    # points are sl2z equivalent on boundary of unit circle
        a = -1/a
        if g is not None: g = S*g
    return a, g


def is_sl2z_equivalent_q(z1, z2):
    w1, _ = sl2z_rep_in_fundom(z1)
    w2, _ = sl2z_rep_in_fundom(z2)
    a1, _ = canonicalize_sl2z_q(w1)
    a2, _ = canonicalize_sl2z_q(w2)
    return a1 == a2


def is_gamma0N_equiv(z1,z2,N):
    """
    Given two CM points z1, z2 on H, determine
    whether they are equivalent mod Gamma0(N).
    """
    w1, g1 = sl2z_rep_in_fundom_q(z1)  # g1(z1) = w1 = canonical rep
    w2, g2 = sl2z_rep_in_fundom_q(z2)  # g2(z2) = w2 = canonical rep
    a1, g1 = canonicalize_sl2z_q(w1, g1)
    a2, g2 = canonicalize_sl2z_q(w2, g2)
    assert g1.acton(z1) == a1
    assert g2.acton(z2) == a2
    if a1 != a2:
        # The points are not even sl2z-equivalent, so they can't be
        # Gamma_0(N) equivalent
        return False
    # So now we have z := g1(z1) = g2(z2), both in the standard
    # fundamental domain.
    #
    # The nontrivial elements of PSL2Z with a fixed point z in the
    # standard fundamental domain for the upper half plane are
    # Stab(z), where
    #
    #     * z = i, so Stab(z) generated by S (where S has order 2)
    #     * z = rho = exp(2*pi*i/3) so Stab(z) generated by S*T
    #     * z = -rhobar = exp(pi*i/3) so Stab(z) generated by T*S
    # The elements in PSL2Z that send z1 to z2 are the elements
    # g2^(-1)*A*g1 for A in Stab(z), so we just check if any are in
    # Gamma0(N).
    #print 'a1, g1: ',a1,g1, 'a2, g2: ',a2,g2
    g2i = g2**(-1)
    g = g2i*g1
    # g satisfies: g(z1) == z2.

    if g[1,0]%N == 0:
        return True
    S, T = g1.parent().gens()
    #I = matrix([1,0,0,1],2)
    if a1**2 == -1:
        #lst = [I,S,-S,-I]
        return (g2i*S*g1)[1,0]%N == 0
    elif a1**3 == 1:
        #lst = [(S*T),(S*T)^2,]
        return (g2i*S*T*g1)[1,0]%N == 0 or (g2i*S*T*S*T*g1)[1,0]%N == 0
    return False



#############################################################
#part VI: what does atkin-lehner map to on rank 2 curves ?  #
#############################################################

def qlist(N):
    """
    returns all q > 1 , q||N.
    """
    qlist0 = []
    for d in N.divisors():
        if d != 1 and gcd(d,N/d) == 1:
            qlist0.append(d)
    return qlist0



def al_point(E,q):
    """
    Let E be an elliptic curve over Q
    with conductor N, w_q be an atkin
    -lehner operator acting on X0(N).
    {z_1, ..., z_k} being its set of
    fixed points. with z1,...,zi of disc -4Q,
    zi+1,...,zk of disc -Q
    phi: X0(N) -> E.
    Then we define the atkin-lehner point
    on E to be al(E,q) =  sum(phi(z_i)).
    The conjecture is al(E,q) is always
    torsion when r(E) >=2.
    If a counter example is found,
    then it's a new way to find points on curve of rank >=2.

    """
    N = E.conductor()
    if N % q != 0 or gcd(q,N/q) != 1:
        raise ValueError, 'q must exactly divide N'
    result = []
    phi = ModularParametrization(E)
    xlist = solve_congruence(q,N)
    v1 = find_points(xlist,q,N,2)
    X = CDF(sum([phi(CDF(a)) for a in v1]))
    result.append((E.elliptic_exponential(X),-4*q))
    if q == 2:
        v2 = find_points_2(N)
        Y = CDF(sum([phi(CDF(a)) for a in v2]))
        result.append((E.elliptic_exponential(Y),-4))
    elif q == 3:
        v2= find_points_3_fixed(N)
        Y = CDF(sum([phi(CDF(a)) for a in v2]))
        result.append((E.elliptic_exponential(Y),-3))
    elif Mod(q,4) == 3:
        v2 = find_points(xlist,q,N,1)
        Y = CDF(sum([phi(CDF(a)) for a in v2]))
        result.append((E.elliptic_exponential(Y),-q))
    return result


def close_points2(P, search_bound, eps=1e-3):
    """
    Return points close to self, sorted by distance (with closed
    point first).
    INPUT:

    - ``search_bound`` -- we search through points of the form n*P + t,
      where P is a generator for the Mordell-Weil, t is any torsion point,
      and -search_bound<=n<=search_bound.

    - ``eps`` -- (default: 1e-3)
    """
    E = rational_curve(P)
    if E.rank() != 2:
        raise ValueError, 'curve must have rank 2'
    g1,g2 = E.gens()
    T = E.torsion_points()
    v = []
    from sage.groups.generic import multiples
    for m in range(0, search_bound):
        Q = m*g1
        for n in range(-search_bound,search_bound):
            R = n*g2
            for t in T:
                P1 = Q + R + t
            if abs(P1[0] - P[0]) < eps and abs(P1[1] - P[1]) < eps:
                # record pair (distance, point)
                v.append(  ((P1[0] - P[0])**2 + (P1[1] - P[1])**2, P1)  )
    v.sort()
    return [R for _, R in v]


def rational_curve(P):
    """
    Return the elliptic curve that this point is on, but over the
    rational numbers.
    """
    return EllipticCurve(QQ,[int(a.real()) for a in P.curve().a_invariants()])



def find_points_3_fixed(N):
    """
    special case: N = 3d
    find the extra fixed points
    given by CM ring Z[w], w^3 = 1
    the only form is x^2+xy + y^2
    Kenku claims that associates permute the
    subgroups.
    alpha = sqrt(-3)*(a unit in Q(w))
    """
    result = []
    d =  ZZ(N/3)
    x = var('x')
    K.<alpha>  = QuadraticField(-3)
    print ' alpha := sqrt(-3) '
    w = (alpha-1)/2
    A = matrix([[1,2],[-2,-1]])
    g1,order = kernel(A,w)[0]
    #print gen, order
    #print ' g1 ', g1
    xlist = solve_congruence(3,N)
    #Note this gives all x such that  N | x^2+3,
    #print ' xlist ', xlist
    for x in xlist:
        A_x = matrix([[1-ZZ(x),2],[-2,-1-ZZ(x)]])
        #print'x = ', x
        #print 'A_x = ', A_x
        #print 'kernel of alpha -x ', kernel(A_x,w)
        if Mod(d,2) == 1:
            if len(kernel(A_x,w)) == 1:
                gen, order = kernel(A_x,w)[0]
            else:
                gen, order = kernel(A_x,w)[1]
            #print 'gen, order:', gen, order
            g2 = gen*(order/d)
            #print 'g2 = ', g2
        elif Mod(d,4) == 2:
            #this is where the main argument
            #I got 3 g2's but they are permuted via
            #the autormophism w, so it suffices to consider 1 of them.
            gen, order = kernel(A_x,w)[1]
            g2 = gen*(order/d)
        else:
            return []
        g = g1 + g2
        #print 'g = ', g
        #print 'Ng, ', N*g.matrix()
        m, n = N*g.matrix()[0][0], N*g.matrix()[0][1]
        #print 'm, n = ', m, n
        # Ng = m +　n*(sqrt(-3)) = m +n(2w+1) = (m + n) + (2*n)*w
        m,n = m+n, 2*n
        #print 'm, n = ', m, n
        # now Ng = m+nw
        m, n = m/gcd(m,n), n/gcd(m,n)
        m,n = find_min(-0.5,0.5,m,n,3,N)
        #print 'new m, n = ', m, n
        a,b,c1,d1 = lift_to_sl2z(m,n,0)
        #print a,b,c1,d1
        point = -(a+b*w)/(c1+d1*w)
        #print 'point ', point
        result.append(point)

    return result


#############################################################
#part VII: some test functions for correctness check        #
#############################################################

def equal_to_zero(Q):
    return Q == 0 or max(abs(Q[0]),abs(Q[1])) >= 1e15

def is_close_ell(P,Q,eps):
    if equal_to_zero(P) and equal_to_zero(Q):
        return True
    elif abs(P[0] - Q[0]) < eps and abs(P[1] - Q[1]) < eps and abs(P[2] - Q[2]) < eps:
        return True
    else:
        return False

def al_dict(F,d):
    '''
    Input: elliptic curve F of conductor N, integer d
    so that there exists atkin-lehner wd.
    Output: a pre-image dictionary with keys the image of fix(wd), the
    item associated with Q is a list consisting of points in fix(wd) that maps to q.
    eps- the small positive used to compare two points, can be changed.
    '''
    eps = 1e-5
    result = {}
    N = F.conductor()
    phi = ModularParametrization(F)
    wd = atkin_lehner(d,N)
    fix = wd.fixed_points_h()
    for z in fix:
        P = F.elliptic_exponential(phi(CDF(z)))
        if result == {}:
            result[P] = [z];
        else:
            is_new_point = True
            for Q in result.keys():
                if is_close_ell(Q, P, eps):
                    result[Q] += [z]
                    is_new_point = False
                    break
            if is_new_point:
                result[P] = [z];

    return result

def img_of_wdoo(F,d):
    N = F.conductor()
    wd = atkin_lehner(d,N)
    vd = image_of_cusp(F,wd(Cusp(Infinity)))
    vd = vector(vd)
    dict1 = period_mapping_conv_fixed(F)
    v1, v2 = dict1.keys()
    w1, w2 = dict1[v1],dict1[v2]
    A = matrix([list(v1),list(v2)])
    a, b = vd*(A**-1)
    Qd = F.elliptic_exponential(a*w1+b*w2)
    print Qd
    Qd = NumericalPoint(Qd, 1e-4)
    try:
        P = Qd.identify()
        return P
    except:
        print 'failed to compute'
        return None


def find_bad_pair(p):
    """
    try to find curves in the worst case to test for the
    even index conjecture : Conductor is prime,
    F has rank 0, 0 maps to a torsion point which is not in 2 F(Q).
    to-do: to find a way that avoids calling img_of_wdoo, which
    can be quite slow and sometimes inaccurate.
    """
    flist = []
    for f in cremona_optimal_curves(p):
        wp = atkin_lehner(p,p)
        if f.rank() == 0:
            good = False
            P = img_of_wdoo(f,p)
            for Q in f.torsion_points():
                if 2*Q == P:
                    good = True
                    break
            if not good:
                flist.append(f.label())
    return flist


def test_img(F,eps):
    """
    test whether the img_of_cusp function works
    by consistency with img_of_fixed_wd.
    Basically, when wd has sign -1 on F,
    it computes the point phi_F(wd(oo))
    and phi_F(p) for p in a list of fixed points of w_d
    If everything's correct, 2*phi_F(p) = phi_F(wd(oo))
    """
    N = F.conductor()
    large = 1e15 # a large number
    dlist = []
    for d in N.divisors():
        if d != 1 and gcd(d,N/d) == 1:
            dlist.append(d)

    for d in dlist:
        wd = atkin_lehner(d,N)
        fix_wd = wd.fixed_points_h()
        if fix_wd == [] or wd.sign(F) == 1: #nothing to check
            pass
        else:
            phi  = ModularParametrization(F)
            vd = image_of_cusp(F,wd(Cusp(Infinity)))
            vd = vector(vd)
            dict1 = period_mapping_conv(F)
            v1, v2 = dict1.keys()
            w1, w2 = dict1[v1],dict1[v2]
            A = matrix([list(v1),list(v2)])
            a, b = vd*(A**-1)
            Qd = F.elliptic_exponential(a*w1+b*w2)
            print Qd
            for f in fix_wd:
                P  = F.elliptic_exponential(phi(CDF(f)))
                P = 2*P
                print P
                # then check if P is close to Qd
                if equal_to_zero(P) and equal_to_zero(Qd):
                    pass
                elif abs(P[0] - Qd[0]) < eps and abs(P[1] - Qd[1]) < eps and abs(P[2] - Qd[2]) < eps:
                    pass
                else:
                    print (F.label(),d,Qd,P)
                    return False
    return True

def img_of_cusp_complex(F,c):
    if not isinstance(c,Cusp):
        c = Cusp(c)
    vd = image_of_cusp(F,c)
    vd = vector(vd)
    dict1 = period_mapping_conv_fixed(F)
    v1, v2 = dict1.keys()
    w1, w2 = dict1[v1],dict1[v2]
    A = matrix([list(v1),list(v2)])
    a, b = vd*(A**-1)
    return a*w1+b*w2











