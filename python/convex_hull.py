# -*- coding: utf-8 -*-
"""
Created on Wed May 27 08:26:05 2020

@author: Gualandi
"""

import numpy as np
from numpy.random import randint

def RandomGridPoint(n, m, N):
    Xs = randint(0, n, N)
    Ys = randint(0, m, N)
    return [(x,y) for x,y in zip(Xs, Ys)]


def PlotPoints(Ps, Bs):
    n = max([p[0] for p in Ps])+1
    m = max([p[1] for p in Ps])+1

    # Report solution value
    import matplotlib.pyplot as plt
    import numpy as np
    import pylab as pl
    from matplotlib import collections  as mc

    fig, ax = pl.subplots(2)
        
    # ax[0].set_xlim(-1,n+1)
    # ax[0].set_ylim(-1,m+1)
    # ax[1].set_xlim(-1,n+1)
    # ax[1].set_ylim(-1,m+1)
    
    ax[0].plot([p[0] for p in Ps], [p[1] for p in Ps], 
             'ro', color='red')
    

    Rs = []
    for a,b in zip(Bs[:-1], Bs[1:]):
        tmp = WalkGrid(a,b)
        Rs += tmp
        
    Rs += WalkGrid(Bs[-1], Bs[0])
    Ts = FillHull(Rs)
    ax[1].plot([p[0] for p in Ts], [p[1] for p in Ts], 
             'ro', color='blue')

    ax[0].autoscale()
    ax[0].margins(0.1)
    ax[1].autoscale()
    ax[1].margins(0.1)

    plt.show()
    
def FilterAxis(Ps):
    # First filter along an axis
    n = max([p[1] for p in Ps])+1
    
    Xs = [[] for _ in range(n)]
    for x,y in Ps:
        Xs[y].append((x,y))

    Bs = []
    for l in Xs:
        if l:
            l.sort(key=lambda p: p[0])
            Bs.append(l[0])
            if l[-1] != l[0]:
                Bs.append(l[-1])
                
    # Then filter according to the second axis
    m = max([p[1] for p in Bs])+1
    
    Ys = [[] for _ in range(m)]
    for x,y in Bs:
        Ys[x].append((x,y))

    Bs = []
    for l in Ys:
        if l:
            l.sort(key=lambda p: p[1])
            Bs.append(l[0])
            if l[-1] != l[0]:
                Bs.append(l[-1])

    return Bs


from math import atan2    
def PolarAngle(a, b=None):
    if b == None: b = anchor 
    x_span = a[0] - b[0]
    y_span = a[1] - b[1]
    return atan2(y_span, x_span)

def Distance(a, b=None):
    if b == None: b = anchor 
    x_span = a[0] - b[0]
    y_span = a[1] - b[1]
    return y_span**2 +  x_span**2

def Det(a, b, c):
    return  (b[0] - a[0]) * (c[1] - a[1]) \
           -(b[1] - a[1]) * (c[0] - a[0])
    
def PolarQuickSort(Ls):
    if len(Ls) <= 1: return Ls
    smaller, equal, larger= [], [], []
    
    pivot_ang = PolarAngle(Ls[randint(0,len(Ls)-1)])
    
    for p in Ls:
        p_ang = PolarAngle(p)
        if p_ang < pivot_ang: 
            smaller.append(p)
        elif p_ang == pivot_ang:
            equal.append(p)
        else:
            larger.append(p)
            
    return PolarQuickSort(smaller) \
        + sorted(equal, key=Distance) \
        + PolarQuickSort(larger)
    
    
def ConvexHull(Ps):
    global anchor

    # Preprocessing    
    Cs = FilterAxis(Ps)
    
    # Find anchor point
    min_idx = None
    for i, (x,y) in enumerate(Cs):
        if min_idx == None or y < Cs[min_idx][1]:
            min_idx = i
        if y == Cs[min_idx][1] and x < Cs[min_idx][0]:
            min_idx = i
            
    anchor = Cs[min_idx]
    Ss = PolarQuickSort(Cs)
    del Ss[Ss.index(anchor)]
    
    Hull = [anchor, Ss[0]]
    for s in Ss[1:]:
        while Det(Hull[-2], Hull[-1], s) <= 0:
            del Hull[-1]
            if len(Hull) < 2: break
        Hull.append(s)
        
    return Hull

def WalkGrid(p0, p1):
    """ https://www.redblobgames.com/grids/line-drawing.html """
    print(p0, p1)
    dx = p1[0] - p0[0]
    dy = p1[1] - p0[1]
    
    if dx == 0:
        ps = []
        for y in range(min(p0[1], p1[1]), max(p0[1], p1[1])):
            ps.append((p1[0], y))
        return ps

    if dy == 0:
        ps = []
        for x in range(min(p0[0], p1[0]), max(p0[0], p1[0])):
            ps.append((x, p1[1]))
        return ps

    
    nx = abs(dx)
    ny = abs(dy)

    sign_x = 1 if dx > 0 else -1
    sign_y = 1 if dy > 0 else -1

    p = [p0[0], p0[1]]
    points = [(p[0], p[1])]
    
    ix, iy = 0, 0
    while ix < nx or iy < ny:
        if (0.5+ix)/nx < (0.5+iy)/ny:
            p[0] += sign_x
            ix += 1
        else:
            p[1] += sign_y
            iy += 1

        points.append((p[0], p[1]))
    
    return points

def Min(a, b):
    return a if a < b else b

def Max(a, b):
    return a if a > b else b

def FillHull(Ps):        
    X = {}
    for p in Ps:
        x, y = p[0], p[1]
        if x not in X:          
            X[x] = (y,y)
        X[x] = (Min(X[x][0], y), Max(X[x][1], y))
    
    Rs = []
    for x in X:
        print(x)
        ymin, ymax = X[x]
        for y in range(ymin, ymax+1):
            Rs.append((x,y))
    
    return Rs
        
        
# -----------------------------------------------
#   MAIN function
# -----------------------------------------------
if __name__ == "__main__":
    np.random.seed(13)
    #Ps = [(int(x), int(y)) for x,y in np.random.multivariate_normal([20,20], [[10,5],[5,10]], 2000)]
    #RandomGridPoint(20, 20, 100)
    Ps = [(2,3), (3,7), (1,4)]
    Hull = ConvexHull(Ps)
    
    PlotPoints(Ps, Hull)
    
    # https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm
    # https://www.redblobgames.com/grids/line-drawing.html
    
    print(Hull)
    