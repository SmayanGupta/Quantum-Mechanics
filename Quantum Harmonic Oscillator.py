# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 14:41:33 2021

@author: smaya
"""

'''Smayan Gupta
19/17067
Section-A
Bsc Physics Hons'''
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh_tridiagonal


#Constants
e=3.795
r1=0.53

k1=100
hpi=1973
m= (0.511*10**6)
n=1000
a=5
ri=-a
rf=a
x=np.linspace(ri,rf,n)
d=((rf-ri)/n)#x[1]-x[0]
#print(d)

#making kinetic energy operator
k=(hpi*hpi)/(2*m*d*d)
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
    
    L[i,i]= (l*(l+1))/(x[i])**2
#print(L)

#making potential operator
def potential(r):
    return (1/2)*k1*(r)**2
V=np.zeros((n,n))
for i in range(0,n):
    r=ri+ (i+1)*d
    V[i,i]= potential(r)#x[i]) 
#print(V)

#print(K)
 
H=-(K+V-L)
#print(H)

#energyvalues,vectors=sp.linalg.eig(H)
energyvalues,vectors=np.linalg.eig(H)
#energyvalues,vectors=eigsh(H,k=10)
#energyvalues,vectors=eigh_tridiagonal(H)
#print(vectors)

#plt.plot(x,np.diag(V))


#energyvalues = np.sort(energyvalues)
print(energyvalues[0:3])
plt.legend(fancybox=True)


for i in range(4):
    plt.plot(x,(vectors.T[i]**2),label='{} Energy level'.format(i))

    
plt.legend()
plt.show()