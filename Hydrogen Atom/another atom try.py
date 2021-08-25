# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 19:51:06 2021

@author: smaya
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import linalg

#Constants
e=3.795
r1=0.53
nth=5
epsilonpi=0.159#9*10**9
hpi=1973
m= (0.511*10**6)
n=1000
#function for nth radius
def orbital(n):
    return r1*n*n

d=(orbital(nth)-r1)/n


#making kinetic energy operator
k=(hpi*hpi)/(m*d)
A=1
B=-2
C=1

X=np.ones(n-1)*A*k
Y=np.ones(n)*B*k
Z=np.ones(n-1)*C*k    
K=np.diag(X,-1)+np.diag(Y,0)+np.diag(Z,1)

#making laplace operator
l=0 
L=np.zeros((n,n))
for i in range(0,n):
    r=r1+(i+1)*d
    L[i,i]= (l*(l+1))/(r)**2
#print(L)

#making potential operator
def potential(r):
    return (e*e*epsilonpi)/(r)
V=np.zeros((n,n))
for i in range(0,n):
    r=r1+(i+1)*d
    V[i,i]= potential(r) 
#print(V)

H=K+V
side=k*np.ones((n-1))
main=-2*k*np.ones((n))


#energyvalues,vectors= linalg.eigh_tridiagonal(main,side)
energyvalues,vectors= sp.sparse.linalg.eigs(H,k=n)
#energyvalues=np.linalg.eigvalh(H)

#print(energyvalues)
x=np.linspace(r1,orbital(nth),n)
#plt.plot(x,energyvalues)
plt.plot(x, vectors[2].T**2)

for i in range(1,nth+1):
    plt.axvline(x=orbital(i),color='k',ls='--',alpha=0.7)



plt.show()