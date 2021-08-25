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

#Constants
e=3.795
r1=0.53
nth=5
epsilonpi=1#0.159#9*10**9
hpi=1973
m= (0.511*10**6)
n=1000
#function for nth radius
def orbital(n):
    return r1*n*n
ri=1e-5
rf=10
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
l=2
L=np.zeros((n,n))
for i in range(0,n):
    
    L[i,i]= (l*(l+1))/(x[i])**2
#print(L)

#making potential operator
def potential(r):
    return (e*e*epsilonpi)/(r)
V=np.zeros((n,n))
for i in range(0,n):
    r=ri+ (i+1)*d
    V[i,i]= potential(r)#x[i]) 
#print(V)
print('----')
#print(K)
 
H=-(K+V-L)
#print(H)

#energyvalues,vectors=sp.linalg.eig(H)
energyvalues,vectors=np.linalg.eigh(H)
#energyvalues,vectors=eigsh(H,k=10)
#print(vectors)

x=np.linspace(ri,rf,n)
#plt.plot(x,np.diag(V))


#energyvalues = np.sort(energyvalues)
print(energyvalues[0:3])
plt.legend(fancybox=True)
#plt.plot(x,(vectors[4].T)**2)
#plt.plot(x,energyvalues*(-13.6))
#print(vectors[i].T)

#for i in range(0,nth+1):
 #   plt.axvline(x=orbital(i),color='k',ls='--',alpha=0.7)
for i in range(3):
    plt.plot(x,(vectors.T[i]**2),label='{} Energy level'.format(i))
plt.legend()
plt.show()


