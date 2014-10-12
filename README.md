zero-modform
============

Description of the code:

zero-modform.sage: Computes j-polynomial of atkin-lehner eigenform of square free conductor N. (**Norm** method)

exact-point.sage: Take as input the j-polynomial of a modular form, computes upper half plane representatives [z] for its zeros.

exceptional-point.sage: computes the set of points on X_0(p) where the map z \mapsto (j(z), j(pz)) is not injective, i.e., when
[z] maps to a singular point on the planar model of X_0(p).

916c/916c.sage: compute a polynomial relation between two modular functions on X0(916) and use them to obtain a critical polynomial
for E = 916c.

atkin-lehner.sage: Computes set of fixed points of any atkin_lehner operator $w_d$ on $X_0(N)$.

multimod.sage: Computes an u-polynomial of atkin-lehner eigenforms of any conductor, with a chosen u satisfying some integrality assumptions( which u = j always satisfy).(**Multimodular** method). N does not have to be square free.

polyrel-ZZ.sage: Computes u-polynomial of atkin-lehner eigenforms with r and u having concentrated poles. N does not have to be square free. (**Yang-product** method)

recognize-hilbert.sage: Given an irreducible polynomial f(x) \in ZZ[x], determine if there exists a negative discriminant D
such that $f = H_D(x)$. In the first case it finds such a $D$. 


results/:

critpolys.json: contains an incomplete list critical polynomials of elliptic curves with conductor \<1000. (Now it contains
all rank 2 curves of conductor <= 1000)

crit-poly.txt: contains all critical polynomials of elliptic curves with prime conductor p \< 1000 such that genus(X_0(p)) > 1.

389-zero-pol.txt: contains zero polynomials for traces of the 5 newforms of weight 2 and level 389.


