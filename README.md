zero-modform
============

There are 3 main files containing code: 

zero-modform: computes the polynomial satisfied by j-invariants of zeros of atkin-lehner eigenforms, based 
on the assumption that level = square free and f is nonzero at cusps.

exact-point: computes approximately the upper half plane representatives [z] for zeros of atkin-lehner eigenforms.
(not finished. Need to use excpetional-point)

exceptional-point: computes the set of points on X_0(p) where the map z \mapsto (j(z), j(pz)) is not injective, i.e., when 
[z] maps to a singular point on the planar model of X_0(p). 

The folder 'results' contains two files: 

crit-poly.txt: contains all critical polynomials of elliptic curves with prime conductor p <1000 such that 
genus(X_0(p)) > 1. 

389


