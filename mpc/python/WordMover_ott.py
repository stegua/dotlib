# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 13:01:19 2023

@author: Gualandi
"""

import numpy as np

def GetData():
    X = np.load('X.npy')
    H = np.load('H.npy')
    
    print('X', X.shape)
    # X (15340, 50)
    
    print('H', H.shape)
    # (653, 15340)
    
    for i in range(len(H)):
        C = H[i]/np.min(H[i][np.nonzero(H[i])])
        if C.sum() > 4000:
            print(C.sum())
        
    indices = np.argsort(np.sum(H, axis=0))
    nmax = 4000
    H = H[:,indices[-nmax:]]
    H = np.array(H)
    # Renormalize so that lines sum to 1
    H = H / np.sum(H, axis=1)[:, None]
    # Drop embeddings not used.
    X = X[indices[-nmax:],:]
    
    print(f'{H.shape[0]} texts supported on up to {H.shape[1]} words of dimension {X.shape[1]}')
    
    return H, X

from scipy.spatial import distance_matrix
def ComputeMatrix(X):
    C = distance_matrix(X, X, 2)
    return C
          

from gurobipy import Model, GRB, quicksum, tuplelist
def DualOT(A, B, C):

    if sum(A) > sum(B):
        A, B = B, A
    
    # Build model
    mod = Model()
    mod.setParam(GRB.Param.OutputFlag,  1)
    # mod.setAttr(GRB.Attr.ModelSense,   -1)
    mod.setParam(GRB.Param.Method,      0)
    mod.setParam(GRB.Param.Threads,      1)
    # mod.setParam(GRB.Param.OutputFlag,  1)
    # mod.setParam(GRB.Param.Crossover,   0)
    #mod.setParam(GRB.Param.NumericFocus, 3)
    
    # Create variables
    x = {}
    for i, a in enumerate(A):
        if a > 0:
            for j, b in enumerate(B):            
                if b > 0:
                    x[i,j] = mod.addVar(obj=C[i,j])
    
    mod.update()
    
    E = tuplelist([(i,j) for (i,j) in x])
    
    # Objective Function
    for i, a in enumerate(A):
        if a > 0:
            mod.addConstr(quicksum(x[i,j] for _,j in E.select(i,'*')) >= a)

    for j, b in enumerate(B):
        if b > 0:
            mod.addConstr(quicksum(x[i,j] for i,_ in E.select('*', j)) <= b)
    
    # Train: Solve the model
    mod.optimize()

    return mod.getAttr(GRB.Attr.ObjVal), mod.Runtime

#-----------------------------------------------
# MAIN function
#-----------------------------------------------
if __name__ == "__main__":   
    H, X = GetData()
    
    C = ComputeMatrix(X)
    
    from time import time
    
    start = time()
    
    for i in range(len(H)):
        for j in range(i+1, len(H)):
            print(DualOT(H[0], H[1], C))
            
    end = time()

    print(f'It took {end - start} seconds!')

    
    
    
    
    