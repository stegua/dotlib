# -*- coding: utf-8 -*-
"""
Created on Thu May 25 16:39:14 2017

@author: gualandi
"""
import matplotlib.pyplot as plt
import numpy as np

def DrawMatrix(A):
    plt.figure(figsize=(10,10))
    plt.imshow(A, cmap='gray')
    plt.show()
    
def cart2pol(x, y):
    return np.arctan2(y, x)

n = 256
A = np.ones((2*n+1,2*n+1))
D = set()
K = {}
Ls = [i for i in range(n)]
for i in Ls:
    for j in Ls:
        if (i,j) != (0,0):
            f = int(100*cart2pol(i,j))
            if f not in D:
                D.add(f)
                K[f] = 1
                A[n+i, n+j] = 0
                A[n-i, n-j] = 0
                A[n+i, n-j] = 0
                A[n-i, n+j] = 0
            else:
                K[f] += 1
            

DrawMatrix(A)
print(K)
print(len(D), n*n, n*np.log(n), sum([K[k] for k in K]))