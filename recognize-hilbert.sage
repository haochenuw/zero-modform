from itertools import *

def all_sublists(n):
    v = list(range(n))
    result = []
    for i in range(1,n+1):
        for a in combinations(v,i):
            result.append(list(a))
    return result

def sub_product(v,sub):
    return prod([v[i] for i in range(len(v)) if i in sub])


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
