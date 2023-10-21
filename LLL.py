import numpy as np


class Basis: 

    # Les constructeurs de cet objet sont : une matrice représentant la base, une matrice contenant les mu_{i,j} et un vecteur contenant les B_i
    def __init__(self, basis):
        self.basis = basis
        self.gram_schmidt()

    #Calcule, à partir d'une base, les coefficients de Gram-Schmidt
    def gram_schmidt(self):
        basis = self.basis
        n = len(basis)
        m = len(basis[0])
        orthogonal_basis = np.copy(basis)
        
        #Matrices contenant les coefficients de Gram-Schmidt
        mu = np.zeros((n,m))
        B = np.zeros(n)
        B[0] = np.dot(basis[0],basis[0])

        for i in range(1,n):
            for j in range(i):
                mu[i,j] = np.dot(basis[i], orthogonal_basis[j])/np.dot(orthogonal_basis[j], orthogonal_basis[j])
                orthogonal_basis[i] -= mu[i,j]*orthogonal_basis[j]
            B[i] = np.dot(orthogonal_basis[i], orthogonal_basis[i])
        self.mu = mu
        self.B = B


    # Echange b_k et b_{k-1} et calcule les nouveaux coefficients de Gram-Schmidt
    def lovasz (self, k):
        n = len(self.basis)
        mu0 = self.mu[k,k-1]
        b = self.B[k] + mu0*mu0*self.B[k-1]
        self.mu[k,k-1] = mu0*self.B[k-1]/b
        self.B[k] = self.B[k-1]*self.B[k]/b
        self.B[k-1] = b
        self.basis[[k, k-1]] = self.basis[[k-1, k]] # Echange b_k et b_{k-1}

        # Actualisation des coefficients de Gram-Schmidt 
        for j in range(k-1):
            self.mu[k,j], self.mu[k-1,j] = self.mu[k-1,j], self.mu[k,j]
        for i in range(k+1, n):
            mu1 = self.mu[i,k-1]
            self.mu[i,k-1] = self.mu[k,k-1]*self.mu[i,k-1]+(1-self.mu[k,k-1]*mu0)*self.mu[i,k]
            self.mu[i,k] = mu1 - mu0*self.mu[i,k]
        
        return max(k-1,1)
            

    # Opère la transformation b_k = b_k - [mu_{k,l}]b_l et actualise les coefficiens de Gram-Schmidt
    def taille(self, k, l):
        if abs(self.mu[k,l])> 0.5+ 10**(-12) : 
            # Le + 10^(-12) permet de résoudre un disfonctionnement dans le cas d'égalité
            self.basis[k] = self.basis[k] - round(self.mu[k,l])*self.basis[l]

            # Actualisation des coefficients de Gram-Schmidt 
            for j in range(l):
                self.mu[k,j] = self.mu[k,j] - round(self.mu[k,l])*self.mu[l,j]
            self.mu[k,l] = self.mu[k,l] - round(self.mu[k,l])
            
            


    # Applique l'algorithme LLL à la base
    def lll(self):
        n = len(self.basis)
        self.gram_schmidt()
        k = 1
        while (k < n):
            self.taille(k,k-1)
            if 0.75*self.B[k-1] > self.B[k] + self.mu[k,k-1]*self.mu[k,k-1]*self.B[k-1]:
                k = self.lovasz(k)
                continue
            for l in range(k-2, -1, -1):
                self.taille(k,l)
            k+=1


if __name__ == "__main__":
    
    '''
    a = np.array([[47.,215.],[95.,460.]])
    b = Basis(a)
    print("Base : \n ", b.basis)

    b.lll()
    print("Base réduite : \n", b.basis)
    '''
    
    a = np.array([[1.,1.,1.],[-1.,0.,2.],[3.,5.,6.]])
    b = Basis(a)
    print("Base: \n", b.basis)
    b.lll()
    print("Base réduite : \n", b.basis)
    