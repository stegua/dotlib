# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 16:47:42 2017

@author: gualandi
"""

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

def DrawMatrix(A):
    plt.figure(figsize=(10,10))
    plt.imshow(A, cmap='gray')
    plt.show()

def TestMatrix(FF, n):
    A = np.ones((n, n))
    cx, cy = 16, 16
    for a,b in FF:
        A[cx+a, cy+b] = 0
    
    for a,b in FF:
        A[cx-a, cy-b] = 0
    
    for a,b in FF:
        A[cx+a, cy-b] = 0
    
    for a,b in FF:
        A[cx-a, cy+b] = 0
    
    for a,b in FF:
        A[cx+b, cy+a] = 0
    
    for a,b in FF:
        A[cx-b, cy-a] = 0
    
    for a,b in FF:
        A[cx-b, cy+a] = 0
    
    for a,b in FF:
        A[cx+b, cy-a] = 0
    DrawMatrix(A)

def IsValid(a, b, n):
    if a < 0 or b < 0:
        return False
    if a >= n or b >= n:
        return False
    return True
    
def GenerateEdges(cx, cy, FF, n):
    Es = set()
    
    for a,b in FF:
        for aa,bb in filter(lambda x: IsValid(x[0], x[1], n), 
                            [(cx+a, cy+b), (cx-a, cy+b), (cx-a, cy-b), (cx+a, cy-b),
                             (cx+b, cy+a), (cx-b, cy+a), (cx-b, cy-a), (cx+b, cy-a)]):        
            Es.add((aa,bb))
    
    return Es


def MakeGraph(n, h, Fs):
    G = nx.DiGraph()
    N = set()
    E = set()
    for i in range(n):
        for j in range(n):
            N.add((i,j))
            G.add_node((i,j))
            for e in GenerateEdges(i, j, Fs, n):
                E.add(((i,j), e))

    print("n:", n, " h:", h, " Fh:", len(Fs), " nodes:", len(N), " edges:", len(E))

# conta il numero di archi per vedere se sono giusti
F1 = [(0,1), (1,1)]
F2 = [(0,1), (1,2), (1,1)]
F3 = [(0,1), (1,3), (1,2), (2,3), (1,1)]
F4 = [(0,1), (1,4), (1,3), (1,2), (2,3), (3,4), (1,1)]
F5 = [(0,1), (1,5), (1,4), (1,3), (2,5), (1,2), (3,5), (2,3), (3,4), (4,5), (1,1)]
F6 = [(0,1), (1,6), (1,5), (1,4), (1,3), (2,5), (1,2), (3,5), (2,3), (3,4), (4,5), (5,6), (1,1)]
F7 = [(0,1), (1,7), (1,6), (1,5), (1,4), (2,7), (1,3), (2,5), (3,7), (1,2), (4,7), 
      (3,5), (2,3), (5,7), (3,4), (4,5), (5,6), (6,7), (1,1)]
F8 = [(0,1), (1,8), (1,7), (1,6), (1,5), (1,4), (2,7), (1,3), (3,8), (2,5), (3,7), 
      (1,2), (4,7), (3,5), (5,8), (2,3), (5,7), (3,4), (4,5), (5,6), (6,7), (7,8), (1,1)]

FF = [F1, F2, F3, F4, F5, F6, F7, F8]

for n in [32, 64, 128, 256, 512]:
    for h, Fh in enumerate(FF):
        MakeGraph(n, h+1, Fh)

     
