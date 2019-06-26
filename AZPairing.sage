

x = polygen(ZZ)
#computes the nth iterate of f
def iter(f, n):
    out = f
    for i in range(n - 1):
        out = f(out)
    return out

#computes the height of a root of f
#f must be irreducible over Q
#returns height
def irreducible_root_height(f):
    K.<y> = NumberField(f)
    h = y.global_height()
    return h

#computes the height of all roots of f
#returns a list containing the height of all roots
#doesn't work with multiplicity
def root_height(f):
    list = f.factor()
    heights = []
    for poly in list:
        d = poly[1]*(poly[0].degree())
        heights.append((irreducible_root_height(poly[0]), d))
    return heights

#computes the value of the nth iterate of f at a
def iter_value(f, a, n):
    out = a
    for i in range(n):
        out = f(out)
    return out

#computes the canonical height at a via a 10th iterate
#assumes a is rational, rational coefficients
def canonical_height(f, a):
    value = iter_value(f, a, 15)
    height = value.global_height()
    c = height/(f.degree())^15
    return c

    

#computes the average height of the roots
#only assumes rational coefficients
#needs fractions, not RDF, as coefficients
def average_height_roots(f):
    #hack to get content to work, test if it has integer coefficients
    g = x^2 + 1
    s = g.parent()
    if (f.parent() == s):
        scalar = f.leading_coefficient()/f.content()
    else:
        f = f*f.denominator()
        scalar = f.leading_coefficient()
    roots = f.complex_roots()
    height = 0
    height += float(log(scalar))
    for num in roots:
        val = num.abs()
        if (val > 1):
            height += log(val)
    return height/f.degree()



#computers average height of the kth preimages up to n
#computes canonical height
#allows non-zero basepoint
#this is equal to pairing with x^2
def average_compare_basepoint(f, n, b):
    print "Canonical height of 0 is ", canonical_height(f, 0)
    for i in range(n):
        g = iter(f, i + 1) - b
        print "average height with basepoint at step", i + 1, " is ", average_height_roots(g)






#f should be gives as a polynomial, g should be gives as a homogenous thing
#computes the nth term in the limit theorem for AZ pairing
#my algorithm uses the formula as average over points of small height
#chooses a point p, looks at preimages of p under f, computer canonical height with respect to g
#p should be rational, and should be chosen so that f^n - p is irreducible
#example: AZPairing(x^2, [c^2 + 15*d^2, d^2], i + 1, 7)
def AZPairing(f, g, n, p):
    poly = iter(f, n) - p
    K.<z> = NumberField(poly)
    Q = ProjectiveSpace(K, 1)
    H = End(Q)
    a = H(g)
    b = a.as_dynamical_system()
    point = Q([z, 1])
    return b.canonical_height(point)
