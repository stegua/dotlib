# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 18:24:13 2017

@author: gualandi
"""

import matplotlib.pyplot as plt
import numpy as np

def DrawMatrix(M):
    # Crea una figura, disegna l'immagine data dalla matrice, aggiungi la colorbar sulla destra
    plt.figure(figsize=(6,6))
    # Uso la colormap 'gray' per avere la schacchiera in bianco&nero
    plt.imshow(M, cmap='gray')
    plt.show()

def ComputeDistanceMatrix(n, p=2):
    """ Compute the ground distance with power p of an n*n image """
    C = {}
    for i in range(n):
        for j in range(n):
            C[i,j] = {}
            for v in range(n):
                for w in range(n):
                    C[i,j][v,w] = pow(abs(i - v)**p + abs(j - w)**p, 1/p)
    
    return C

from gurobipy import Model, GRB, quicksum

def WassersteinDualSimplex(h1, h2, M):
    """ Find the Wasserstein distance using the dual simplex """
    n = len(h1)
    P = []
    for i in range(n):
        for j in range(n):
            P.append((i,j))

    # Build model
    m = Model()
    m.setParam(GRB.Param.TimeLimit, 300)
    m.setParam(GRB.Param.Threads, 1)
    #  Options are: 
    #      -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 
    #                    3=concurrent, 4=deterministic concurrent.
    m.setParam(GRB.Param.Method, 0)
    
    # Create variables
    x = {}            
    for i,j in P:
        x[i,j] = {}
        for v,w in P:
            x[i,j][v,w] = m.addVar(obj=M[i,j][v,w])    
    m.update()
    
    for i,j in P:
        m.addConstr(quicksum(x[i,j][v,w] for v,w in P) <= h1[i][j])

    for v,w in P:
        m.addConstr(quicksum(x[i,j][v,w] for i,j in P) >= h2[v][w])

    # Solve the model
    m.optimize()
    
    return m.getAttr(GRB.Attr.ObjVal)

#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    filename1 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data64_1001.csv'
    M1 = np.loadtxt(open(filename1, "rb"), delimiter=",")
    DrawMatrix(M1)

    filename2 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data64_1003.csv'
    M2 = np.loadtxt(open(filename2, "rb"), delimiter=",")
    DrawMatrix(M2)

    C = ComputeDistanceMatrix(len(M1))
    
    print(WassersteinDualSimplex(M1, M2, C))