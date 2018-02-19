# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 09:48:53 2017

@author: gualandi
"""

import networkx as nx


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

def farey(n):
    def gcd(a,b):
        while b: a,b = b,a%b
        return a

    def simplify(a,b):
        g = gcd(a,b)
        return (int(a/g),int(b/g))

    fs = dict()
    for i in range(1,n+1):
        for i2 in range(1,i+1):
            if i2 < n and i != i2:
                r = simplify(i2,i)
                fs[float(i2)/i] = r

    return [(0,1)] + [fs[k] for k in sorted(fs.keys())] + [(1,1)]

#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    # Farey sequence up to H=8
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
    
    #MakeGraph(8,8,F8)
    #print(len(farey(512)))
    #F32 = farey(32)
    #MakeGraph(32, 32, F32)
    
    FareySeq = farey(8)
    print(FareySeq)
    MakeGraph(16, 16, FareySeq)
    
    #for n in [32, 64, 128, 256, 512]:
     #   for h, Fh in enumerate(FF):
      #      MakeGraph(n, h+1, Fh)