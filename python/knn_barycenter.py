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
from gurobipy import Model, GRB, quicksum, tuplelist
from numpy.linalg import norm

def G(x,y):
    return 73*exp(-((x-14)**2/64 + (y-14)**2/64))


def W1(H1, H2):
    """ Find the Wasserstein distance using a cutting plane on the dual """
    n = len(H1)
    s = int(np.sqrt(n))
    
    def ID(x,y):
        return x*s+y
    
    H1 = np.array(H1)
    H1 = H1/sum(H1)

    H2 = np.array(H2)
    H2 = H2/sum(H2)

    # Build model
    m = Model()
    m.setParam(GRB.Param.OutputFlag, 0)
    #m.setParam(GRB.Param.Method, 1)
    #m.setParam(GRB.Param.Threads, 7)
    #m.setParam(GRB.Param.NumericFocus, 1)
    #m.setParam(GRB.Param.Crossover, 0)
    m.setAttr(GRB.Attr.ModelSense, 1)         
    
    # Create variables
    X = {}
    for i in range(s):
        for j in range(s):
            X[ID(i,j)] = {}
        
    P = []
    for i in range(s):
        for j in range(s-1):
            X[ID(i,j)][ID(i,j+1)] = m.addVar(obj=norm(np.array([i,j,G(i,j)]) - np.array([i,j+1,G(i,j+1)])))
            X[ID(i,j+1)][ID(i,j)] = m.addVar(obj=norm(np.array([i,j+1,G(i,j+1)]) - np.array([i,j,G(i,j)])))

            P.append((ID(i,j), ID(i,j+1)))
            P.append((ID(i,j+1), ID(i,j)))

    for i in range(s-1):
        for j in range(s):
            X[ID(i,j)][ID(i+1,j)] = m.addVar(obj=norm(np.array([i,j,G(i,j)]) - np.array([i+1,j,G(i+1,j)])))
            X[ID(i+1,j)][ID(i,j)] = m.addVar(obj=norm(np.array([i+1,j,G(i+1,j)]) - np.array([i,j,G(i,j)])))
       
            P.append((ID(i,j), ID(i+1,j)))
            P.append((ID(i+1,j), ID(i,j)))

    for i in range(s-1):
        for j in range(s-1):
            X[ID(i,j)][ID(i+1,j+1)] = m.addVar(obj=norm(np.array([i,j,G(i,j)]) - np.array([i+1,j+1,G(i+1,j+1)])))
            X[ID(i+1,j+1)][ID(i,j)] = m.addVar(obj=norm(np.array([i+1,j+1,G(i+1,j+1)]) - np.array([i,j,G(i,j)])))
            X[ID(i,j+1)][ID(i+1,j)] = m.addVar(obj=norm(np.array([i,j+1,G(i,j+1)]) - np.array([i+1,j,G(i+1,j)])))
            X[ID(i+1,j)][ID(i,j+1)] = m.addVar(obj=norm(np.array([i+1,j,G(i+1,j)]) - np.array([i,j+1,G(i,j+1)])))

            P.append((ID(i,j), ID(i+1,j+1)))
            P.append((ID(i+1,j+1), ID(i,j)))
            P.append((ID(i,j+1), ID(i+1,j)))
            P.append((ID(i+1,j), ID(i,j+1)))
             
    P = tuplelist(P)      

    m.update()
    
    # Flow variables
    for i in range(n):
        Fs = P.select(i,'*')
        Bs = P.select('*',i)
        m.addConstr(quicksum(X[i][j] for _,j in Fs) - quicksum(X[j][i] for j,_ in Bs) == H1[i] - H2[i])

    # Solve the model
    m.optimize()
    
    return m.getAttr(GRB.Attr.ObjVal)

def W2(H1, H2):
    """ Find the Wasserstein distance using a cutting plane on the dual """
    H1 = np.array(H1)
    H1 = H1/sum(H1)

    H2 = np.array(H2)
    H2 = H2/sum(H2)

    h1 = np.array(H1).reshape(28,28)
    h2 = np.array(H2).reshape(28,28)
    
    n = len(h1)
    P = []
    for i in range(n):
        for j in range(n):
            P.append((i,j))                        

    # Build model
    m = Model()
    m.setParam(GRB.Param.OutputFlag, 0)
                      
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
    for i,j in P:
        m.addConstr(quicksum(X[i,j,i,v] for v in range(n)) <= h1[i,j])    
    for i,j in P:
        m.addConstr(quicksum(Y[v,j,i,j] for v in range(n)) >= h2[i,j])    
    for i,j in P:
        m.addConstr(quicksum(X[i,v,i,j] for v in range(n))-quicksum(Y[i,j,v,j] for v in range(n)) == 0)    

    # Solve plain model        
    m.optimize()
    
    return m.getAttr(GRB.Attr.ObjVal)

def Euclidean(H1, H2):
    H1 = np.array(H1)
    H1 = H1/sum(H1)

    H2 = np.array(H2)
    H2 = H2/sum(H2)

    return euclidean(H1, H2)    
    
def kNN(H, Bs):
    Ls = list(map(lambda B: W1(H[1:], B), Bs))
    best_i = 0
    best_d = Ls[0]
    for i,d in enumerate(Ls):
        if d < best_d:
            best_i = i
            best_d = d
    print('True:', int(H[0]), ' ===> Predict:', best_i)
    if int(H[0]) != best_i:
        print(Ls)
        f, axarr = plt.subplots(1, 2)
        A = np.array(H[1:]).reshape(28,28)
        axarr[0].imshow(A, cmap='coolwarm')
        B = np.array(Bs[best_i]).reshape(28,28)
        axarr[1].imshow(B, cmap='coolwarm')
        plt.subplots_adjust(wspace=0, hspace=0)        
        plt.tight_layout()
        plt.show()

    return int((int(H[0])-best_i) == 0)
    

#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    #barycenters = genfromtxt('E:\\GitHub\\dotlib\\data\\barycenter_Euclidean.csv', delimiter=',')
    barycenters = genfromtxt('E:\\GitHub\\dotlib\\data\\barycenter_L2_all.csv', delimiter=',')
    
    n = 20
    testset = genfromtxt('E:\\GitHub\\dotlib\\data\\test_classify_knn.csv', delimiter=',', skip_header=1, max_rows=n+1)
    
    accuracy = 0 
    for t in testset[:n]:
        accuracy += kNN(t, barycenters)
    print('Accurcay:', round(accuracy/n*100, 3))