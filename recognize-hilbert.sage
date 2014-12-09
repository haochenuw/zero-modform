from itertools import combinations

exceptional_discs = [4*a for a in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 18, 21, 22, 24, 25, 28, 30, 33, 37, 40, 42, 45, 48, 57, 58, 60, 70, 72, 78, 85, 88, 93, 102, 105, 112, 120, 130, 133, 165, 168, 177, 190, 210, 232, 240, 253, 273, 280, 312, 330, 345, 357, 385, 408, 462, 520, 760, 840, 1320, 1365, 1848]]


exceptional_discs +=[3,7,11,15,19,27,35,43,51,67,75,91,99,115,123,147,163,187,195,235,267,315,403,427,435,483,555,595,627,715,795,1155,1435,1995,3003,3315]


## generators for negative discriminants.
def neg_discs(n):
    D = -3
    while D > n:
        yield D
        D = next_disc(D)

def next_disc(D):
    while True:
        D = D-1
        if D%4 ==0:
            m = D.divide_knowing_divisible_by(4)
            if m.is_squarefree() and m%4 in [2,3]:
                return D
        elif D%4 == 1:
            if D.is_squarefree():
                return D

class HilbertSearch():
    def __init__(self,d,m,v):
        self.disc = d
        self.m = m
        self.v = v


    def search(self,state):
        if state[0] == 1 and self.is_goal(state):
            result = [state]
            for newState in self.successors(state):
                result += self._search(newState)
            return result
        else:
            return self._search(state)


    def _search(self,state):
        if self.is_deadend(state):
            return []

        elif self.is_goal(state):
            return [state]
        else:
            result = []
            for newState in self.successors(state):
                result += self._search(newState)
            return result

    def is_goal(self,state):
        return state[-1] == self.m

    def is_deadend(self,state):
        return Mod(self.m, state[-1]) != 0

    def successors(self,state):
        from sage.schemes.elliptic_curves.heegner import RingClassField
        result =[]
        d = self.disc
        n,_ = state
        for p in self.v:
            result.append((n*p, RingClassField(d,p*n).degree_over_K()))
        return result

def all_sublists(n):
    v = list(range(n))
    result = []
    for i in range(1,n+1):
        for a in combinations(v,i):
            result.append(tuple(a))
    return result

def sub_product(v,sub):
    return prod([v[i] for i in range(len(v)) if i in sub])


def all_products(v):
    result = []
    for i in range(len(v)+1):
        for a in combinations(v,i):
            result.append(prod(a))
    return result


def _field_disc(f):
    """
    the field discriminant.
    """
    D = f.disc()
    prod = 1
    for p, mult in list(D.factor()):
        if Mod(mult,2) == 1:
            prod= p
    return -prod

from sage.schemes.elliptic_curves.heegner import RingClassField

def is_hilbert(f):
    """
    if f is a Hilbert class poly, then return the True, disc D.
    else; return False, None
    """

    try:
        f = QQ[x](f)
    except:
        pass

    if not f.is_irreducible():
        return False, None



    deg = ZZ(f.degree())
    if deg.prime_factors() == [2] or deg == 1:
        ### check for the exceptional discriminants, one class per genus.
            for d1 in exceptional_discs:
                if QQ[x](hilbert_class_polynomial(-d1)) == f:
                    verbose("case1: one class per genus.")
                    return True, -d1



    #what are the ramified primes? Contained in the ramification

    F.<a> = NumberField(f)
    D = F.disc()

    v = D.prime_divisors()
    d = prod(v)

    verbose('ramified primes are %s'%v)


    K.<b>  = QuadraticField(-d)

    fields = [K]

    if mod(d,2) == 0:  # if 2 is ramified. Then we have two options:
        # K = QQ(sqrt(d)) or K = QQ(sqrt(2d))
        L.<c> = QuadraticField(-2*d)
        fields.append(L)


    for field in fields:
        flag,disc = _is_hilbert_givenfield(f,field)
        if flag: return flag,disc


    return False, None



def _is_hilbert_givenfield(f,K):
    """
    K is an imag quad field. f is a polynomial...
    We want to check all orders contained in OK if f is the
    ring class field of that order.
    """
    dK = K.disc()


    verbose('field disc = %s'%dK)


    clK = K.class_number()

    deg = ZZ(f.degree())

    verbose('deg = %s'%deg)

    clKprime = QQ(clK*2) / K.number_of_roots_of_unity()

    quotient = deg/clKprime

    if quotient.denominator() > 1:
        v = []

    else:
        quotient = ZZ(quotient)
        verbose('quotient = %s'%quotient)
        v = quotient.prime_divisors()


    verbose('possible prime factors of conductor = %s'%v)



    H = HilbertSearch(dK,deg,v)

    candidates = [m for m,_ in H.search((1,clK))]

    verbose('candidates = %s'%candidates)

    for c in candidates:
        if QQ[x](hilbert_class_polynomial(dK*c**2)) == f:
            return True, dK*c**2

    return False, None





def rec_hilbert(f):
    d = f.disc()
    v = []
    prod0 = 1
    for p, mult in list(d.factor()):
        if Mod(mult,2) == 1:
            prod0*= p
        else:
            v.append(p)
    N = len(v)
    verbose('N = %s'%N)
    verbose('prod0 = %s'%prod0)
    cands = []
    for sub in all_sublists(N):
        n  = prod0*sub_product(v,sub)
        if check(n,v):
            cands.append(n)
    verbose('len = %s'%len(cands))
    for i in range(1,10):
        newcands =[]
        newcands = filtering(cands,v,f,i)
        cands = newcands
        verbose('len = %s'%len(cands))
        if len(cands) <= 1:
            break

    return cands

def check(n,v):
    for b in v:
        if kronecker(-n,b) == 1:
            return False
    return True


def filtering(cands,v,f,i):

    result = []
    verbose('i = %s'%i)
    w = ZZ(f(i)).prime_divisors()
    verbose('w= %s'%w)
    for n in cands:
        for p in w:
            passed = True
            print 'n,p =', n,p
            if (kronecker(-n,p) == 1 and p < n):
                passed = False
                break
        if passed:
            result.append(n)
    return result
