# -*- coding: utf-8 -*-
"""
Created on Thu May  4 11:33:01 2017

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
                    #C[i,j][v,w] = pow(abs(i - v)**p + abs(j - w)**p, 1/p)
                    C[i,j][v,w] = abs(i - v)**p + abs(j - w)**p
    
    return C

def WassersteinDualCutting(h1, h2, M):
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
    m.setParam(GRB.Param.Method, 1)
    
    # Set a maximization problem
    m.setAttr(GRB.Attr.ModelSense, -1) 
                      
    print('1. Start building model')
    # Create variables
    V = {}
    U = {}
    for i,j in P:
        # First set of dual variables
        V[i,j] = m.addVar(lb=-GRB.INFINITY, ub=0, obj=h1[i,j])
        # Second set of dual variables
        u_ub = sum([M[v,w][i,j] for v,w in P])
        U[i,j] = m.addVar(lb=0, ub=u_ub, obj=h2[i,j])
        
    m.update()
    print('2. Add initial constraint sets')
    for i,j in P:
        for v,w in P:
            #if M[i,j][v,w] <= 16.001: # Threshold for first set of constraints
            m.addConstr(V[i,j] + U[v,w], GRB.LESS_EQUAL, M[i,j][v,w])    

    # Solve plain model        
    m.optimize()

    # Solve the model
    print('3. Start Cutting planes')
    it = 0
    stime = 0
    while False: # set to true to test the method
        it += 1
        m.optimize()
        stime += m.RunTime
        flag = True        
        max_depth = 0
        for i,j in P:
            depth = -1
            a,b,c,d = -1,-1,-1,-1
            for v,w in P:
                if V[i,j].X + U[v,w].X - M[i,j][v,w] > depth:
                    a,b,c,d = i,j,v,w
                    depth = V[i,j].X + U[v,w].X - M[i,j][v,w]
            if (max_depth == 0 and depth > 0.001) or (depth >= max_depth):
                max_depth = max(max_depth, depth)
                flag = False
                m.addConstr(V[a,b] + U[c,d], GRB.LESS_EQUAL, M[a,b][c,d])
        print('ITERATION:', it,' MAX DEPTH:',round(max_depth,3),' Time:',round(stime,3))
        if flag:
            break
        #else:
        #    m.addConstr(V[a,b] + U[c,d], GRB.LESS_EQUAL, M[a,b][c,d])
    
    return m.getAttr(GRB.Attr.ObjVal)

#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    filename1 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data32_1001.csv'
    M1 = np.loadtxt(open(filename1, "rb"), delimiter=",")

    filename2 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data32_1002.csv'
    M2 = np.loadtxt(open(filename2, "rb"), delimiter=",")

    C = ComputeDistanceMatrix(len(M1),p=2)
    
    print(WassersteinDualCutting(M1, M2, C))