import numpy as np
from numpy.polynomial import Polynomial
from LLL import Basis
from rsa import RSA


#Opere la transformation polynôme-vecteur exposee dans le dossier
def poly_to_vect(poly, X, size=-1):
    d = poly.degree()
    n = max(size,d+1)
    vect = np.zeros(n)
    for i in range(n):
        if i<d+1:
            vect[i] = poly.coef[i] * X ** i
        
    return vect

# Opere la transformation polynome-vecteur exposee dans le dossier
def vect_to_poly(vect, X):
    n = len(vect)
    poly = Polynomial([0])
    for i in range(n):
        monomial = Polynomial([0] * i + [1])
        poly += (vect[i] / (X ** i)) * monomial

    return poly


# Construit la base de vecteurs pour la methode de la section 4.1
def poly_to_basis2(poly, N, X, h):
    d = poly.degree()
    a = np.zeros((d * h, d * h))
    for j in range(h):
        for i in range(d):
            monomial = Polynomial([0] * i + [1])
            G = (N**(h-1-j))*(poly**j)*monomial
            vect = poly_to_vect(G, X, size=h*d)
            a[i + d * j] += vect
    return Basis(a)

# Construit la base de vecteurs pour la methode de Coppersmith
def poly_to_basis1(poly, N, X):
    d = poly.degree()
    a = np.zeros((d+1 , d+1 ))
    
    for i in range(d):
        monomial = Polynomial([0] * i + [1])
        G = N*monomial
        vect = poly_to_vect(G, X, d+1)
        a[i] = vect
    a[d] = poly_to_vect(poly, X) 
    return Basis(a)




poly = Polynomial([-222, 5000,  10, 1])
print(poly)
n = 10001
X = 10
h = 3

print('\n')

base1 = poly_to_basis1(poly, n, X)
base1.lll()
p1 = vect_to_poly(base1.basis[0],X)
print(p1)
print(p1.roots())

print('\n')

base2 = poly_to_basis2(poly, n, X, h)
base2.lll()
p2 = vect_to_poly(base2.basis[0],X)
print(p2)
print(p2.roots())


