load('~/critical-point/modified-typespace.sage')


def ram_index(E,p):
    """
    Input:

    E = Elliptic curve over QQ of conductor N.
    p = a prime dividing the conductor of E.

    Output:

    the ramification index of the modular param X0(N) -> E at the cusps of
    denominator p^r. where
        r = ord_p(N)/2 if p = 2 or 3 and ord_p(N) is even. (If odd, then all cusps of denominator a power of p is unramified).
        r = 1 if p >= 5.
    """
    if isinstance(E,str):
        E = EllipticCurve(E)
    if not p.is_prime():
        raise ValueError('p must be a prime')
    N = E.conductor()
    if p >= 5:
        return 1
    elif N.valuation(p) % 2:
        return 1
    else:
        r = N.valuation(p) //2

    f = E.modular_form()
    T = TypeSpaceModified(f,p)
    N1 = T.maximal_twist_level()
    assert N % N1 == 0
    adjustment = 1
    if N1.valuation(p) < 2:
        adjustment = p**(2-N1.valuation(p))
    return N//(N1*adjustment)






