# -*- coding: utf-8 -*-
"""
Created on Mon May  7 14:02:08 2018

@author: gualandi
"""

import numpy as np
from gurobipy import Model, GRB, quicksum, tuplelist

def ComputeDistanceMatrix(n, p=2):
    """ Compute the ground distance with power p of an n*n image """
    C = {}
    for i in range(n):
        for j in range(n):
            C[i,j] = {}
            for v in range(n):
                for w in range(n):
                    C[i,j][v,w] = abs(i - v)**p + abs(j - w)**p
    
    return C

def W2_Primal(h1, h2):
    """ Find the Wasserstein distance using a cutting plane on the dual """
    n = len(h1)
    P = []
    for i in range(n):
        for j in range(n):
            P.append((i,j))                        

    # Build model
    m = Model()
    #m.setParam(GRB.Param.TimeLimit, 300)
    #m.setParam(GRB.Param.Presolve,    0)    
    #m.setParam(GRB.Param.Threads,     1)
    #  Options are: 
    #      -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 
    #                    3=concurrent, 4=deterministic concurrent.
    #m.setParam(GRB.Param.Method, 1)
                      
    print('1. Start building model')
    # Create variables
    X = {}
    Y = {}
    # Between first and second partition
    for i,j in P:
        for v in range(n):
            X[i,j,i,v] = m.addVar(lb=0, obj=abs(j - v)**2)
    # Between second and third partition    
    for i,j in P:
        for v in range(n):
            Y[i,j,v,j] = m.addVar(lb=0, obj=abs(i - v)**2)
        
    m.update()
    print('2. Add initial constraint sets')
    for i,j in P:
        m.addConstr(quicksum(X[i,j,i,v] for v in range(n)) <= h1[i,j])    
    for i,j in P:
        m.addConstr(quicksum(Y[v,j,i,j] for v in range(n)) >= h2[i,j])    
    for i,j in P:
        m.addConstr(quicksum(X[i,v,i,j] for v in range(n))-quicksum(Y[i,j,v,j] for v in range(n)) == 0)    

    # Solve plain model        
    m.optimize()
    
    return m.getAttr(GRB.Attr.ObjVal)

#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    filename1 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data128_1001.csv'
    M1 = np.loadtxt(open(filename1, "rb"), delimiter=",")

    filename2 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data128_1003.csv'
    M2 = np.loadtxt(open(filename2, "rb"), delimiter=",")

    print(W2_Primal(M1, M2))