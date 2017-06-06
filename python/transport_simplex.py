# -*- coding: utf-8 -*-
"""
Created on Sat May 27 08:48:23 2017

@author: gualandi
"""

import numpy as np

class Base(object):
    def __init__(self):
        self.i = []
        self.j = []
        self.X = {}

    def add(self, a, b, f):
        self.i.append(a)
        self.j.append(b)
        self.X[a,b] = f

    def remove(self, a, b):
        self.i.remove(a)
        self.j.remove(b)
        del self.X[a,b]


def Dist(x, y):
    return (x[0]-y[0])**2 + (x[1]-y[1])**2

def TransSimplex(A, B):
    pass

def InitBasis(X, Y, C):
    A = np.copy(X)
    B = np.copy(Y)
    i = 0
    j = 0
    n = len(A)
    m = len(B)
    
    X = {}
    base = Base()
    cost = 0
    while i < n and j < m:
        #print(A[i], B[j])
        #print(B)
        X[i,j] = round(min(A[i], B[j]))
        base.add(i,j, X[i,j])
        cost += C[i,j]*X[i,j]
        print('step:', i+1, j+1, C[i,j])
        if A[i] < B[j]:
            B[j] = B[j] - A[i]
            i += 1
        elif A[i] > B[j]:
            A[i] = A[i] - B[j]
            j += 1
        else:
            A[i] = A[i] - B[j]
            print('equal:', i+1, j+1, A[i], B[j])
            j += 1
        
    return base, X, cost
    
# Find dual variables for the given basic solution
def FindDual(A, B, C, X):
    n = len(A)
    m = len(B)
    
    u = np.zeros(n)
    v = np.zeros(m)
    
    i,j = next(iter(X))
    u[i] = C[i,j]
    v[j] = 0
    
    cc = set()
    cc.add((i,j))
    xx = set(list(X))
    while len(cc) > 0:
        ii, jj = next(iter(cc))
        print(ii,jj)
        cc.remove((ii,jj))
        xx.remove((ii,jj))
        for a,b in xx:
            if a == ii:
                v[b] = v[jj] + C[a,b] - C[ii,jj]
                cc.add((a,b))
            if b == jj:
                u[a] = u[ii] + C[a,b] - C[ii,jj]
                cc.add((a,b))

    # Return the found duals
    return u, v     

def MostNegative(u, v, X):
    M = 0
    for i in range(len(u)):
        for j in range(len(v)):
            print(u[i] + v[j] - C[i,j], end=' ')
            if u[i] + v[j] - C[i,j] > M:
                a,b = i,j
                M = u[i] + v[j] - C[i,j]
        print()
        
    return M, a, b
    
def OutOfBase(X, base):
    xmin = X[0,0]
    a = 0
    b = 0
    for i, j in X:
        if len(base.i) > 1 and len(base.j) > 1:
            if X[i,j] < xmin:
                xmin = X[i,j]
                a = i
                b = j
    
    return a,b, xmin
    
#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    #I1 = [[2,3,1,4], [1,2,3,4], [4,3,1,2], [2,3,1,4]]
    #I2 = [[2,3,4,1], [1,3,2,4], [4,3,1,2], [4,3,1,2]]

    I1 = [1.0, 5.0, 7.0]
    I2 = [3.0, 3.0, 3.0, 2.0, 2.0]
    C = np.array([[3, 2, 1, 2, 3], [5, 4, 3, -1, 1], [0, 2, 3, 4, 5]])
    print(C)
    
    A = np.array(I1).flatten()
    B = np.array(I2).flatten()
    
    n = len(A)
    eps = 0.0001
    A = A+eps
    B[n-1] = B[n-1] + n*eps
    print(B[0]+n*eps)
    
    base, X, cost = InitBasis(A, B, C)
    print(len(X), cost)
    
    u, v = FindDual(A, B, C, X)
    M, k, l = MostNegative(u, v, X)
    
    i,j, theta = OutOfBase(X, base)
    print(i,j,theta)