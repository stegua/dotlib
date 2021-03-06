# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 15:07:22 2017

@author: gualandi
"""

import numpy as np
import networkx as nx
from time import time

from gurobipy import Model, GRB, quicksum

def SolveCuttingPlane(h1, h2, M):
    m = Model()
    m.setParam(GRB.Param.TimeLimit, 300)
    m.setParam(GRB.Param.Method, 1)
    # Create variables
    x = {}
    n = len(h1)
    for i in range(n):
        for j in range(n):
            x[i,j] = m.addVar(obj=M[i][j], name='x'+str(i)+'_'+str(j))
    m.update()
    
    for i in range(n):
        m.addConstr(quicksum(x[i,j] for j in range(n)) == h1[i])

    for j in range(n):
        m.addConstr(quicksum(x[i,j] for i in range(n)) == h2[j])

    # Solve the model
    m.optimize()
    
    return m.getAttr(GRB.Attr.ObjVal)

def Wasserstein(h1, h2, M):
    """ Compute the Wasserstein distance between the two histograms 
        h1 and h2 using distance matrix M """
    d = len(h1)
    # Build the graph for max flow 
    G = nx.DiGraph()
    # add d nodes for each histrogam
    # (d+1) source, (d+2) target
    for i in range(d):
        G.add_node(i,   demand=-h1[i])
        G.add_node(d+i, demand=+h2[i])
    # Add all edges
    for i in range(d):
        for j in range(d):
            G.add_edge(i, d+j, weight=M[i][j], capacity=min(h1[i], h2[j]))
    
    #flowCost, flowDict = nx.capacity_scaling(G, heap=nx.utils.heaps.PairingHeap)
    flowCost, flowDict = nx.capacity_scaling(G, heap=nx.utils.heaps.BinaryHeap)     
    #flowCost, flowDict = nx.network_simplex(G)
    return flowCost
    
def MakeHistogram(d):
    """ Make a normalized random histogram on the simplex """
    hist = np.random.permutation(range(d))
    hist = [np.random.uniform(0,1) for _ in range(d)]
    hsum = sum(hist)
    hist = [h/hsum for h in hist]
    return hist
    
def MakeCostMatrix(d):
    """ Make a ransom metrix matrix as described in Cuturi 2013 (Figure 4) """
    G = nx.erdos_renyi_graph(d, 0.5)
    for u,v,w in G.edges(data=True):
        w['weight'] = np.random.uniform(0,1)
    # All pair shortest path
    M = nx.floyd_warshall(G)
    return M

#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    start = time()
    # Random graph as in 
    d = 512
    M = MakeCostMatrix(d)
    print("Build matrix time: ", time()-start)

    # Create two random histograms of dimension d
    h1 = MakeHistogram(d)
    h2 = MakeHistogram(d)

    print(SolveCuttingPlane(h1, h2, M)) 
    print("Gurobi time: ", time()-start)
    #print(Wasserstein(h1, h2, M))
    #print("Total time: ", time()-start)
