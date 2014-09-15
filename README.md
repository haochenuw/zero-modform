zero-modform
============

Description of the code:

zero-modform.sage: computes the polynomial satisfied by j-invariants of zeros of atkin-lehner eigenforms, based
on the assumption that the level of the form is square free.

exact-point.sage: computes approximately the upper half plane representatives [z] for zeros of atkin-lehner eigenforms.
(still in progress)

exceptional-point.sage: computes the set of points on X_0(p) where the map z \mapsto (j(z), j(pz)) is not injective, i.e., when
[z] maps to a singular point on the planar model of X_0(p).

916c/916c.sage: compute a polynomial relation between two modular functions on X0(916) and use them to obtain a critical polynomial
for E = 916c.

(to-do: add the code for computing 664a1 and 944e1)

in results:

critpolys.json: contains an incomplete list critical polynomials of elliptic curves with conductor <1000. (Now it contains
all rank 2 curves of conductor <= 1000)

crit-poly.txt: contains all critical polynomials of elliptic curves with prime conductor p <1000 such that
genus(X_0(p)) > 1.

389-zero-pol.txt: contains zero polynomials for traces of the 5 newforms of weight 2 and level 389.


