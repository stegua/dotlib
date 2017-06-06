# -*- coding: utf-8 -*-
"""
Created on Thu May  4 10:45:47 2017

@author: gualandi
"""

import numpy as np
from gurobipy import Model, GRB, quicksum, tuplelist

def ComputeDistanceMatrix(n, p=2, q=1):
    """ Compute the ground distance with power p of an n*n image """
    def m(x,y):
        return x*s+y
    
    s = int(np.sqrt(n))
    C = np.zeros((n, n))
    for i in range(s):
        for j in range(s):
            for v in range(s):
                for w in range(s):
                    C[m(i,j)][m(v,w)] = pow(pow(abs(i - v)**p + abs(j - w)**p, 1/p), q)
    return C


def Preprocessing(h1, h2):
    assert(h1.size == h2.size)
    n = len(h1)
    for i in range(n):
        offset = min(h1[i], h2[i])
        h1[i] -= offset
        h2[i] -= offset  
              
    return h1, h2
              
    
def WassersteinDualSimplex(h1, h2, M):
    """ Find the Wasserstein distance using the dual simplex """
    n = len(h1)

    # Build model
    m = Model()
    #m.setParam(GRB.Param.TimeLimit, 300)
    #m.setParam(GRB.Param.Presolve,    0)    
    #m.setParam(GRB.Param.Threads,     1)
    #  Options are: 
    #      -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 
    #                    3=concurrent, 4=deterministic concurrent.
    #m.setParam(GRB.Param.Method, 0)
    print('1. Start building model')
    # Create variables
    x = {}            
    P = set()
    D = []
    for i in range(n):
        if h1[i] > 0:
            x[i] = {}
            for j in range(n):
                if h2[j] > 0:
                    x[i][j] = m.addVar(ub=min(h1[i], h2[j]), obj=M[i][j])
                    D.append((i,j))
                    P.add(j)
    D = tuplelist(D)
    m.update()
    print('2. Add initial constraint sets')
    for i in x:
        m.addConstr(quicksum(x[i][j] for j in x[i]) <= h1[i])
    
    for j in P:
        m.addConstr(quicksum(x[i][j] for i,j in D.select('*',j)) >= h2[j])
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
    M1 = np.array(M1.flatten())

    filename2 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data64_1003.csv'
    M2 = np.loadtxt(open(filename2, "rb"), delimiter=",")
    M2 = np.array(M2.flatten())

    #M1,M2 = Preprocessing(M1,M2)


    C = ComputeDistanceMatrix(len(M1), p=2, q=1)
    print(WassersteinDualSimplex(M1, M2, C))