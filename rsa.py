import random


# ALGORITHMES D'ARITHMETIQUE USUELS 

#Calcule m^e mod n (expodentiation rapide)
def mod_exp(m, e, n):
    result = 1
    base = m % n
    while e > 0:
        if e % 2 == 1:
            result = (result * base) % n
        base = (base * base) % n
        e //= 2
    return result

#Calcule le r = pgcd(a,b) ainsi que des entiers u et v tels que au+bv = r  
def euclide_etendu(a, b):
    u1, u2, u3 = 1, 0, a
    v1, v2, v3 = 0, 1, b

    while v3 != 0:
        q = u3 // v3
        v1, v2, v3, u1, u2, u3 = (
            u1 - q * v1,
            u2 - q * v2,
            u3 - q * v3,
            v1,
            v2,
            v3,
        )
    return u3, u1, u2


# GENERATION DE NOMBRES PREMIERS

# Test de temoin
def test_temoin(a, n):
    k = 0
    d = n - 1
    while d % 2 == 0:
        k += 1
        d //= 2
    x = mod_exp(a, d, n)
    if x == 1 or x == n - 1:
        return False
    for _ in range(k - 1):
        x = (x * x) % n
        if x == n - 1:
            return False
    return True

# Test de Miller-Rabin
def miller_rabin(n, k):
    if n <= 3:
        return n == 2 or n == 3
    if n % 2 == 0:
        return False

    for _ in range(k):
        a = random.randint(2, n - 2)
        if test_temoin(a, n):
            return False
    return True



# Generateur de nombres premiers
def generate_prime(bits, k=100):
    while True:
        n = random.getrandbits(bits)
        n |= (1 << (bits - 1)) | 1
        if miller_rabin(n, k):
            return n

#Genere deux nombre premiers distincts
def generate_primes(bits = 1024):
    p = generate_prime(bits)
    q = generate_prime(bits)
    while p==q :
        q = generate_primes
    return p,q

# GENERATION CLES RSA 


class RSA :

    def __init__(self) :
        self.generate_keys()

    #Cle publique : (e,n)
    #Cle privee : d 
    def generate_keys(self, bits = 1024):
        p,q = generate_primes()
        self.n = p*q
        phi_n = (p-1)*(q-1)
        while True :
            e = random.randint(2, phi_n)
            r,u,v = euclide_etendu(e,phi_n)
            if r == 1:
                self.e = e
                self.d = u % phi_n
                break


# ENCODAGE D'UN MESSAGE 
   
    def encode(self, m):
        return mod_exp(m, self.e, self.n)
    
# DECODAGE D'UN MESSAGE 
    def decode(self, c):
        return mod_exp(c, self.d, self.n)


if __name__== "__main__":
    rsa = RSA()
    print('Clé publique n  : ', rsa.n)
    print("\n")
    print('Clé publique e  : ', rsa.e)
    print("\n")
    print('Clé privée d : ', rsa.d)
    print("\n")

    #premieres decimales de pi
    m = 314159265358979323846264338327950288419716939937510 
    print('Message : ', m)
    print("\n")
    
    c = rsa.encode(m)
    print('Message codé : ', c)
    print("\n")
    
    m2 = rsa.decode(c)
    print('Message décodé : ', m2)