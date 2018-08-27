# -*- coding: utf-8 -*-
"""
Created on Wed May  9 15:56:52 2018

@author: gualandi
"""

from gurobipy import Model, GRB, quicksum, tuplelist
from math import log
import numpy as np
import pylab as plt

def WassersteinDualSimplex(h1, h2, M,x0):
    """ Find the Wasserstein distance using the dual simplex """
    n = len(h1)
    m = len(h2)

    # Build model
    mod = Model()
    #mod.setParam(GRB.Param.NumericFocus, 3)
    #m.setParam(GRB.Param.TimeLimit, 300)
    #m.setParam(GRB.Param.Presolve,    0)    
    #m.setParam(GRB.Param.Threads,     1)
    #  Options are: 
    #      -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 
    #                    3=concurrent, 4=deterministic concurrent.
    mod.setParam(GRB.Param.Method, 0)
    print('1. Start building model')
    # Create variables
    x = {}            
    P = set()
    D = []
    for i in range(n):
        if h1[i] > 0:
            x[i] = {}
            for j in range(m):
                if h2[j] > 0:
                    x[i][j] = mod.addVar(ub=min(h1[i], h2[j]), obj=M[i][j])
                    D.append((i,j))
                    P.add(j)
    D = tuplelist(D)
    mod.update()
    print('2. Add initial constraint sets')
    for i in x:
        mod.addConstr(quicksum(x[i][j] for j in x[i]) <= h1[i])
    
    for j in P:
        mod.addConstr(quicksum(x[i][j] for i,j in D.select('*',j)) >= h2[j])
    
    #for i,j in x0:
    #    mod.addConstr(x[i][j] == x0[i,j])
        
    print('3. Start solution phase')
    # Solve the model
    mod.optimize()
    
    nnz = 0
    for i,j in D:
        if x[i][j].x > 0:
            nnz += 1
    print('nnz',nnz)
    return mod.getAttr(GRB.Attr.ObjVal)
    
def GreedyConst(h1,h2,C):
    m = len(h1)
    n = len(h2)
    
    x = {}
    
    a = [i for i in h1]
    b = [i for i in h2]
    
    M = list(range(m))
    N = list(range(n))
    
    epsilon=10
    vmax = sum(a)*sum(b)
    c = 0
    it = 0
    while len(M) > 0 and len(N) > 0:
        it += 1 
        vmin = vmax
        uhk = 0
        h,k = -1,-1
        np.random.shuffle(M)
        np.random.shuffle(N)
        for i in M:
            if vmin == vmax:
                for j in N:
                    uij = min(a[i], b[j])
                    if uij > 0:
                        v = C[i,j]*uij + epsilon*uij*(log(uij)-1)
                        if v <= vmin:
                            vmin = v
                            h, k = i, j
                            uhk = uij
            
        if h != -1:
            x[h,k] = uhk
            a[h] -= uhk
            b[k] -= uhk
            c += C[h,k]*uhk
            if it % 50 == 0:
                print(it,c)
            if a[h] == 0:
                M.remove(h)
            if b[k] == 0:
                N.remove(k)
    while True:
        As = [(h,k) for h,k in x if x[h,k] > 0.0001]
        tol = -1e-15
        best = tol
        ii,jj,hh,kk=-1,-1,-1,-1
        np.random.shuffle(As)
        flag = False
        for a in range(len(As)):
            i,j = As[a]
            if best == tol:
                for b in range(a+1,len(As)):
                    h,k = As[b]
                    delta = min(x[i,j], x[h,k])
                    
    #                d1 = -(C[i,j]*x[i,j] + epsilon*x[i,j]*(log(x[i,j])-1) + C[h,k]*x[h,k] + epsilon*x[h,k]*(log(x[h,k])-1))
    #                if x[i,j]-delta > 0:
    #                    d1 += C[i,j]*(x[i,j]-delta) + epsilon*(x[i,j]-delta)*(log(x[i,j]-delta)-1) 
    #                if x[h,k]-delta > 0:
    #                    d1 += C[h,k]*(x[h,k]-delta) + epsilon*(x[h,k]-delta)*(log(x[h,k]-delta)-1)
    #                
    #                xhj = x.get((h,j),0)
    #                xik = x.get((i,k),0)
    #                if xik > 0:
    #                    d1 += -(C[i,k]*xik + epsilon*xik*(log(xik)-1))
    #                if xhj > 0:
    #                    d1 += -(C[h,j]*xhj + epsilon*xhj*(log(xhj)-1))
    #                d1 += C[i,k]*(xik+delta) + epsilon*(xik+delta)*(log(xik+delta)-1) 
    #                d1 += C[h,j]*(xhj+delta) + epsilon*(xhj+delta)*(log(xhj+delta)-1)
    
                    d1 = (-C[h,k]-C[i,j]+C[h,j]+C[i,k])*delta
    
                    if d1 < best:
                        flag = True
                        best = d1              
                        ii,jj = i,j
                        hh,kk = h,k

        if flag:
            h,k = hh, kk
            i,j = ii,jj
            delta = min(x[i,j], x[h,k])
            x[h,k] -= delta
            x[i,j] -= delta
            x[h,j] = x.get((h,j), 0) + delta
            x[i,k] = x.get((i,k), 0) + delta
            c += delta*(-C[h,k]-C[i,j] + C[h,j] + C[i,k])
            it += 1 
            print(it, c, len(As))
        else:
            break

    return c,x

    
#-----------------------------------------------
# MAIN function
#-----------------------------------------------
if __name__ == "__main__":
    np.random.seed(13)
    N = [300,200]
    d = 2
    x = np.random.rand(2,N[0])-.5
    
    theta = 2*np.pi*np.random.rand(1,N[1])
    r = .8 + .2*np.random.rand(1,N[1])
    y = np.vstack((np.cos(theta)*r,np.sin(theta)*r))
    
#    plotp = lambda x,col: plt.scatter(x[0,:], x[1,:], s=200, edgecolors="k", c=col, linewidths=2)
#    
#    plt.figure(figsize=(10,10))
#    
#    plotp(x, 'b')
#    plotp(y, 'r')
#    
#    plt.axis("off")
#    plt.xlim(np.min(y[0,:])-.1,np.max(y[0,:])+.1)
#    plt.ylim(np.min(y[1,:])-.1,np.max(y[1,:])+.1)
#    
#    plt.show()
    
    x2 = np.sum(x**2,0)
    y2 = np.sum(y**2,0)
    C = np.tile(y2,(N[0],1)) + np.tile(x2[:,np.newaxis],(1,N[1])) - 2*np.dot(np.transpose(x),y)
                         
    h1 = np.ones(N[0])/N[0]
    h2 = np.ones(N[1])/N[1]
    
    c1,x0 = GreedyConst(h1,h2,C)
    c2 = WassersteinDualSimplex(h1,h2,C,x0)
    print('error:', (c1-c2)/c2, ' c1:',c1,' c2:',c2)