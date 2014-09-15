* in the file critpolys.json the data is formatted as curve: factorization over QQ of the critical polynomial of that curve.

* When the level is square free. The critpoly is defined by
    $\prod_{z : f_E(z) = 0 }(x - j(z))$.

* When the level is not square free. The critpoly is defined by
    $\prod_{z : f_E(z) = 0 }(x - u(z))$
Where u is a rational function on $X_0(N)$, usually u != j. 

For curve = 664a1, $u = EtaProduct(8,{8:8,2:4,4:-12})^{-1}$. 

For curve = 944e1, $u = EtaProduct(16,{8:-6,4:2,16:4})^{-1}$ and one should igonre the small factors before x^224+...
since those are extra zeros introduced in the computation. 

For curve = 916c1. $u = Eta product of level 916 : (eta_1)^-3 (eta_2)^7 (eta_4)^-2 (eta_229)^-1 (eta_458)^5 (eta_916)^-6$, and the factor x^2 should be ignored.
