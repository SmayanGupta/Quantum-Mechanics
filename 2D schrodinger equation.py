# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 17:27:15 2021

@author: smaya
"""
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import scipy as sp
from scipy.sparse.linalg import eigsh

n=100
n2=n**2
a=np.ones(n-1)
b=np.ones(n)*(-2)
dx=np.diag(a,-1)+np.diag(b,0)+np.diag(a,1)
X, Y = np.meshgrid(np.linspace(0,1,n, dtype=float),
                   np.linspace(0,1,n, dtype=float))

def potential(x,y):
    return 0*x
V=potential(X,Y)
V.reshape(n2)
U=np.diag(V,0)
K=-sparse.kronsum(dx,dx)
H=K+U

eigenvalues, eigenvectors = eigsh(H, k=10, which='SM')
def get_e(i):
    return eigenvectors.T[i].reshape((n,n))

plt.figure(figsize=(9,9))
plt.contourf(X, Y, get_e(3)**2, 20)

