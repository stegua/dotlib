# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 18:30:41 2018

@author: gualandi
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 16:08:59 2018

@author: gualandi
"""
import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial.distance import euclidean, correlation, russellrao

from math import sqrt, exp
from numpy import genfromtxt
from gurobipy import Model, GRB, quicksum, tuplelist, QuadExpr
from numpy.linalg import norm

def Learn(images, L):
    """ Find the Wasserstein distance using a cutting plane on the dual """
    K = range(len(images))
    n = len(images[0])
    s = int(np.sqrt(n))
    
    KP = []
    for k1 in K:
        for k2 in K:
            if k1 < k2:
                KP.append((k1,k2))
                
    def ID(x,y):
        return x*s+y
        
    # Build model
    m = Model()
    m.setParam(GRB.Param.Method, 2)
    #m.setParam(GRB.Param.Threads, 7)
    m.setParam(GRB.Param.NumericFocus, 1)
    m.setParam(GRB.Param.Crossover, 0)
    m.setAttr(GRB.Attr.ModelSense, 1)         
    
    print('1. Start building model')
    # Create variables
    
    X = {}
    for k in K:
        for i in range(s):
            for j in range(s):
                X[k, ID(i,j)] = {}
    X[k, n] = {}
        
    P = []
    for k1,k2 in KP:
        for i in range(s):
            for j in range(s-1):
                X[k1,ID(i,j)][k2,ID(i,j+1)] = m.addVar()
                X[k1,ID(i,j+1)][k2,ID(i,j)] = m.addVar()
    
                if k1 == 0:
                    P.append((ID(i,j), ID(i,j+1)))
                    P.append((ID(i,j+1), ID(i,j)))

        for i in range(s-1):
            for j in range(s):
                X[k1,ID(i,j)][k2,ID(i+1,j)] = m.addVar()
                X[k1,ID(i+1,j)][k2,ID(i,j)] = m.addVar()
           
                if k == 0:
                    P.append((ID(i,j), ID(i+1,j)))
                    P.append((ID(i+1,j), ID(i,j)))

        for i in range(s-1):
            for j in range(s-1):
                X[k1,ID(i,j)][k2,ID(i+1,j+1)] = m.addVar()
                X[k1,ID(i+1,j+1)][k2,ID(i,j)] = m.addVar()
                X[k1,ID(i,j+1)][k2,ID(i+1,j)] = m.addVar()
                X[k1,ID(i+1,j)][k2,ID(i,j+1)] = m.addVar()

                if k1 == 0:
                    P.append((ID(i,j), ID(i+1,j+1)))
                    P.append((ID(i+1,j+1), ID(i,j)))
                    P.append((ID(i,j+1), ID(i+1,j)))
                    P.append((ID(i+1,j), ID(i,j+1)))
                 
    P = tuplelist(P)      

    Z = {}
    for i in range(n):
        Fs = P.select(i,'*')
        for i,j in Fs:
            Z[i,j] = m.addVar(lb=1)
    
    m.update()
    
    obj = QuadExpr()
    for i in range(n):
        Fs = P.select(i,'*')
        for i,j in Fs:
            for k1,k2 in KP:
                if L[k1] == L[k2]:
                    obj += Z[i,j]*X[k1,i][k2,j]
                else:
                    obj -= 0.00005*Z[i,j]*X[k1,i][k2,j]
                    
    m.setObjective(obj)
    m.update()

    
    # Flow variables
    print('2. Add initial constraint sets')
    for i in range(n):
        Fs = P.select(i,'*')
        Bs = P.select('*',i)
        for k1,k2 in KP:
            m.addConstr(quicksum(X[k1,i][k2,j] for _,j in Fs) - quicksum(X[k1,j][k2,i] for j,_ in Bs) == images[k1][i] - images[k2][i])

    for i in range(n):
        Fs = P.select(i,'*')
        for i,j in Fs:
            Gs = P.select(j,'*')
            for j,h in Gs:
                if (i,h) in Z:
                    m.addConstr(Z[i,j] + Z[j,h] >= Z[i,h])
    
    print('3. Start solution phase')
    
    # Solve the model
    m.write('learn.mps')
    m.optimize()
    
    return
    

#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    
    n = 5
    testset = genfromtxt('E:\\GitHub\\dotlib\\data\\test_classify_knn.csv', delimiter=',', skip_header=1, max_rows=n+1)

    L = []
    images = []
    for k,t in enumerate(testset):
        L.append(t[0])
        A = np.array(t[1:])
        A = A/sum(A)
        images.append(A)
    
    Learn(images, L)