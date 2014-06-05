def weight_index(p,k):
    """
    return the tuple (A,B,C) in the writeup in misc math project/zero-polynomial
    of modform. Here p = level, k = weight.
    """
    B = (1 + kronecker(-3,p))*k
    C = (1 + kronecker(-1,p))*k/2
    A = (k*(p+1)-4*B-6*C)/12
    try:
        return (ZZ(A),ZZ(B),ZZ(C))
    except:
        return None


def prec(p,k):
    """
    returns the precisions needed for computing the zero
    polynomial of a modular form f.
    return a tuple (k(g-1)+2,kg+1,p(kg+1)). Ordered from
    small to large, and they are the precisions needed
    for (Delta,E4,E6), (f), (f(q^{1/p})), respectively.
    We are going to add 1 for safety
    """
    g = Gamma0(p).genus()
    return (k*(g-1)+ 2 + 1 , k*g+1 + 1, p*(k*g+1 + 1))