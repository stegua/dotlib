# -*- coding: utf-8 -*-
"""
Created on Sat May  6 19:42:25 2017

@author: gualandi
"""

import numpy as np
from gurobipy import Model, GRB, quicksum, tuplelist

def ComputeDistanceMatrix(n):
    """ Compute the ground distance with power p of an n*n image """
    C = {}
    for i in range(n):
        for j in range(n):
            C[i,j] = {}
            for v in range(n):
                for w in range(n):
                    C[i,j][v,w] = abs(i - v) + abs(j - w)
    
    return C


def WassersteinMinFlowL1(h1, h2):
    """ Find the Wasserstein distance using a cutting plane on the dual """
    def ID(x,y):
        return x*s+y

    n = len(h1)
    s = int(np.sqrt(n))
    print(n,s)
    # Build model
    m = Model()
    #m.setParam(GRB.Param.TimeLimit, 300)
    #m.setParam(GRB.Param.Presolve,    0)    
    #m.setParam(GRB.Param.Threads,     1)
    #  Options are: 
    #      -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 
    #                    3=concurrent, 4=deterministic concurrent.
    #m.setParam(GRB.Param.Method, 1)
    
    # Set a maximization problem
    m.setAttr(GRB.Attr.ModelSense, 1) 
                      
    print('1. Start building model')
    # Create variables
    X = {}
    for i in range(s):
        for j in range(s):
            X[ID(i,j)] = {}
    X[n] = {}
        
    P = []
    for i in range(s):
        for j in range(s-1):
            X[ID(i,j)][ID(i,j+1)] = m.addVar(obj=1)
            X[ID(i,j+1)][ID(i,j)] = m.addVar(obj=1)

            P.append((ID(i,j), ID(i,j+1)))
            P.append((ID(i,j+1), ID(i,j)))

            
    for i in range(s-1):
        for j in range(s):
            X[ID(i,j)][ID(i+1,j)] = m.addVar(obj=1)
            X[ID(i+1,j)][ID(i,j)] = m.addVar(obj=1)
       
            P.append((ID(i,j), ID(i+1,j)))
            P.append((ID(i+1,j), ID(i,j)))                    

    P = tuplelist(P)            

    m.update()
    
    # Flow variables
    print('2. Add initial constraint sets')
    for i in range(n):
        b = h2[i] - h1[i]        
        m.addConstr(quicksum(X[i][j] for j in X[i])-quicksum(X[j][i] for j,i in P.select('*',i)) == b)

    #m.addConstr(quicksum(X[n][i] for i in range(n)) >= abs(sum(h1)-sum(h2)))
        
    print('3. Start solution phase')
    # Solve the model
    m.optimize()
    
    return m.getAttr(GRB.Attr.ObjVal)

#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    filename1 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data64_1001.csv'
    M1 = np.loadtxt(open(filename1, "rb"), delimiter=",")    
    #M1 = np.ones((4,4))
    M1 = np.array(M1.flatten())
    
    filename2 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data64_1003.csv'
    M2 = np.loadtxt(open(filename2, "rb"), delimiter=",")
    #M2 = np.ones((4,4))
    M2 = np.array(M2.flatten())    
    print('0. Reading data finished')
    
    #C = ComputeDistanceMatrix(len(M1))
    
    print(WassersteinMinFlowL1(M1, M2))
    
    