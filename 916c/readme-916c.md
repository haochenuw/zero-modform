
###We find the polynomial relation between two modular functions r and u on X0(916),
where div(r) = Z_f + [0] + [1/2] + [1/229] + [1/458] - 228(oo) and div(u) =  86[1/2] + 29[1/458] + 28[1/4] - 143(oo).

r is got from multiplying $fj/dj$ by a suitable eta product of level 4 and then apply the Fricke involution $z \mapsto 1/916z$
u is a eta product of level 916. We chose r,u so that they have only one pole and the valuations at this pole are coprime.

###By YfYang's argument one sees that there's a relation in r,u of form
    $p(r,u) = -r^143 + u^228 + \sum_{143a+228b \leq 143+228} c_{a,b}r^au^b = 0$

(By the way, since r and u has leading coefficient 1 and integer coefficients, we know $p(r,u) \in ZZ[r,u]$, so we could
also do a multimodular algorithm.)

### We use the code in 916c-new.sage to cancel the poles until we get a holomorphic function which must be a constant.

###Then we define F(x) = p(0,x). The zeros of F(x) will be the u-values at the zeros of r. It turns out that
        F(x)  = x^2h(x), h irreducible over $\mathbb{Q}[x]$, deg(h) = 228
The $x^2$ term is 'cause the two shared zeros of r and u at the cusps [1/2] and [1/458].

###Hence $h(x) =  \prod (x - u(z_i))$ where z_i are the zeros of the modular form $f_916c1$.

