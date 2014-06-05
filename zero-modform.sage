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