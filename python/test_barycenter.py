# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 16:32:08 2018

@author: gualandi
"""

import time
import csv
import numpy as np
from math import sqrt

from numpy import genfromtxt
import matplotlib.pyplot as plt
from gurobipy import Model, GRB, quicksum, tuplelist

def DrawDigit(A, label=''):
    """ Draw single digit as a greyscale matrix"""
    fig = plt.figure(figsize=(6,6))
    # Uso la colormap 'gray' per avere la schacchiera in bianco&nero
    img = plt.imshow(A, cmap='gray_r')
    plt.xlabel(label)
    plt.show()
    
      
def BarycenterL1(images):
    """ Find the Wasserstein distance using a cutting plane on the dual """
    K = len(images)
    n = len(images[0])
    s = int(np.sqrt(n))
    
    def ID(x,y):
        return x*s+y
    
    print(K,n,s)
    # Build model
    m = Model()
    m.setParam(GRB.Param.Method, 2)
    #m.setParam(GRB.Param.Threads, 7)
    #m.setParam(GRB.Param.NumericFocus, 1)
    m.setParam(GRB.Param.Crossover, 0)
    m.setAttr(GRB.Attr.ModelSense, 1)         
    
    print('1. Start building model')
    # Create variables
    Z = {}
    for i in range(n):
        Z[i] = m.addVar(obj=0)#,name='z'+str(i), ub=1.0)
    
    X = {}
    for k in range(K):
        for i in range(s):
            for j in range(s):
                X[k, ID(i,j)] = {}
    X[k, n] = {}
        
    P = []
    for k in range(K):
        for i in range(s):
            for j in range(s-1):
                X[k,ID(i,j)][k,ID(i,j+1)] = m.addVar(obj=1)#, name='x_'+str(k)+'_'+str(i)+'_'+str(j)+'_'+str(i)+'_'+str(j+1))
                X[k,ID(i,j+1)][k,ID(i,j)] = m.addVar(obj=1)#, name='x_'+str(k)+'_'+str(i)+'_'+str(j+1)+'_'+str(i)+'_'+str(j))
    
                if k == 0:
                    P.append((ID(i,j), ID(i,j+1)))
                    P.append((ID(i,j+1), ID(i,j)))

    for k in range(K):            
        for i in range(s-1):
            for j in range(s):
                X[k,ID(i,j)][k,ID(i+1,j)] = m.addVar(obj=1)#, name='x_'+str(k)+'_'+str(i)+'_'+str(j)+'_'+str(i+1)+'_'+str(j))
                X[k,ID(i+1,j)][k,ID(i,j)] = m.addVar(obj=1)#, name='x_'+str(k)+'_'+str(i+1)+'_'+str(j)+'_'+str(i)+'_'+str(j))
           
                if k == 0:
                    P.append((ID(i,j), ID(i+1,j)))
                    P.append((ID(i+1,j), ID(i,j)))

    P = tuplelist(P)      

    m.update()
    
    # Flow variables
    print('2. Add initial constraint sets')
    for i in range(n):
        Fs = P.select(i,'*')
        Bs = P.select('*',i)
        for k in range(K):
            m.addConstr(quicksum(X[k,i][k,j] for _,j in Fs) - quicksum(X[k,j][k,i] for j,_ in Bs) == images[k][i] - Z[i])

    m.addConstr(quicksum(Z[i] for i in range(n)) == 1.0)
    print('3. Start solution phase')
    # Solve the model
    m.optimize()
    
    return [Z[i].X for i in range(n)]


def BarycenterLinf(images):
    """ Find the Wasserstein distance using a cutting plane on the dual """
    K = len(images)
    n = len(images[0])
    s = int(np.sqrt(n))
    
    def ID(x,y):
        return x*s+y
    
    print(K,n,s)
    # Build model
    m = Model()
    m.setParam(GRB.Param.Method, 2)
    #m.setParam(GRB.Param.Threads, 7)
    m.setParam(GRB.Param.NumericFocus, 1)
    m.setParam(GRB.Param.Crossover, 0)
    m.setAttr(GRB.Attr.ModelSense, 1)         
    
    print('1. Start building model')
    # Create variables
    Z = {}
    for i in range(n):
        Z[i] = m.addVar(obj=0)
    
    X = {}
    for k in range(K):
        for i in range(s):
            for j in range(s):
                X[k, ID(i,j)] = {}
    X[k, n] = {}
        
    P = []
    for k in range(K):
        for i in range(s):
            for j in range(s-1):
                X[k,ID(i,j)][k,ID(i,j+1)] = m.addVar(obj=1)
                X[k,ID(i,j+1)][k,ID(i,j)] = m.addVar(obj=1)
    
                if k == 0:
                    P.append((ID(i,j), ID(i,j+1)))
                    P.append((ID(i,j+1), ID(i,j)))

    #for k in range(K):            
        for i in range(s-1):
            for j in range(s):
                X[k,ID(i,j)][k,ID(i+1,j)] = m.addVar(obj=1)
                X[k,ID(i+1,j)][k,ID(i,j)] = m.addVar(obj=1)
           
                if k == 0:
                    P.append((ID(i,j), ID(i+1,j)))
                    P.append((ID(i+1,j), ID(i,j)))


    #for k in range(K):            
        for i in range(s-1):
            for j in range(s-1):
                X[k,ID(i,j)][k,ID(i+1,j+1)] = m.addVar(obj=1)
                X[k,ID(i+1,j+1)][k,ID(i,j)] = m.addVar(obj=1)
                X[k,ID(i,j+1)][k,ID(i+1,j)] = m.addVar(obj=1)
                X[k,ID(i+1,j)][k,ID(i,j+1)] = m.addVar(obj=1)

                if k == 0:
                    P.append((ID(i,j), ID(i+1,j+1)))
                    P.append((ID(i+1,j+1), ID(i,j)))
                    P.append((ID(i,j+1), ID(i+1,j)))
                    P.append((ID(i+1,j), ID(i,j+1)))
                 
    P = tuplelist(P)      

    m.update()
    
    # Flow variables
    print('2. Add initial constraint sets')
    for i in range(n):
        Fs = P.select(i,'*')
        Bs = P.select('*',i)
        for k in range(K):
            m.addConstr(quicksum(X[k,i][k,j] for _,j in Fs) - quicksum(X[k,j][k,i] for j,_ in Bs) == images[k][i] - Z[i])

    m.addConstr(quicksum(Z[i] for i in range(n)) == 1.0)
    print('3. Start solution phase')
    # Solve the model
    m.optimize()
    
    return [Z[i].X for i in range(n)]


def BarycenterL2(images):
    """ Find the Wasserstein distance using a cutting plane on the dual """
    K = len(images)
    n = len(images[0])
    s = int(np.sqrt(n))
    
    def ID(x,y):
        return x*s+y
    
    print(K,n,s)
    # Build model
    m = Model()
    m.setParam(GRB.Param.Method, 2)
    #m.setParam(GRB.Param.Threads, 7)
    m.setParam(GRB.Param.NumericFocus, 1)
    m.setParam(GRB.Param.Crossover, 0)
    m.setAttr(GRB.Attr.ModelSense, 1)         
    
    print('1. Start building model')
    # Create variables
    Z = {}
    for i in range(n):
        Z[i] = m.addVar(obj=0)
    
    X = {}
    for k in range(K):
        for i in range(s):
            for j in range(s):
                X[k, ID(i,j)] = {}
    X[k, n] = {}
        
    P = []
    for k in range(K):
        for i in range(s):
            for j in range(s-1):
                X[k,ID(i,j)][k,ID(i,j+1)] = m.addVar(obj=1)
                X[k,ID(i,j+1)][k,ID(i,j)] = m.addVar(obj=1)
    
                if k == 0:
                    P.append((ID(i,j), ID(i,j+1)))
                    P.append((ID(i,j+1), ID(i,j)))

    #for k in range(K):            
        for i in range(s-1):
            for j in range(s):
                X[k,ID(i,j)][k,ID(i+1,j)] = m.addVar(obj=1)#, name='x_'+str(k)+'_'+str(i)+'_'+str(j)+'_'+str(i+1)+'_'+str(j))
                X[k,ID(i+1,j)][k,ID(i,j)] = m.addVar(obj=1)#, name='x_'+str(k)+'_'+str(i+1)+'_'+str(j)+'_'+str(i)+'_'+str(j))
           
                if k == 0:
                    P.append((ID(i,j), ID(i+1,j)))
                    P.append((ID(i+1,j), ID(i,j)))


    #for k in range(K):            
        for i in range(s-1):
            for j in range(s-1):
                X[k,ID(i,j)][k,ID(i+1,j+1)] = m.addVar(obj=sqrt(2))
                X[k,ID(i+1,j+1)][k,ID(i,j)] = m.addVar(obj=sqrt(2))
                X[k,ID(i,j+1)][k,ID(i+1,j)] = m.addVar(obj=sqrt(2))
                X[k,ID(i+1,j)][k,ID(i,j+1)] = m.addVar(obj=sqrt(2))

                if k == 0:
                    P.append((ID(i,j), ID(i+1,j+1)))
                    P.append((ID(i+1,j+1), ID(i,j)))
                    P.append((ID(i,j+1), ID(i+1,j)))
                    P.append((ID(i+1,j), ID(i,j+1)))
                 
    P = tuplelist(P)      

    m.update()
    
    # Flow variables
    print('2. Add initial constraint sets')
    for i in range(n):
        Fs = P.select(i,'*')
        Bs = P.select('*',i)
        for k in range(K):
            m.addConstr(quicksum(X[k,i][k,j] for _,j in Fs) - quicksum(X[k,j][k,i] for j,_ in Bs) >= images[k][i] - Z[i])

    m.addConstr(quicksum(Z[i] for i in range(n)) <= 1.0001)
    print('3. Start solution phase')
    m.write('mnist_ineq_6.mps')
    # Solve the model
    m.optimize()
    
    return [Z[i].X for i in range(n)]

def EuclideanBarycenter(images):
    n = len(images[0])
    N = range(n)
    avg = [0 for _ in N]
    for row in images:
        for i in N:
            avg[i] += row[i] 
    
    for i in N:
        avg[i] = avg[i]/n
        
    return avg
    
#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    # Misura il tempo per le operazioni principali
    start = time.time()
    DIR = 'E:\\GitHub\\dotlib\\data\\barycenter_'
    SFX = '.csv'
    NUM = ['6']#str(i) for i in range(10)]
    
    for n in NUM:
        FILEIN = DIR+n+SFX
        my_data = genfromtxt(FILEIN, delimiter=',', skip_header=1)
        print('Reading time:', time.time()-start)
        start = time.time()
        
        # Normalize pixels
        images = []
        for row in my_data:
            # Documentation for function 'reshape':
            # https://docs.scipy.org/doc/numpy/reference/generated/numpy.reshape.html
            A = np.array(row[1:])
            A = A/sum(A)
            images.append(A)
            #DrawDigit(A, 'Digit: ' + str(int(row[0])))
            #print(sum(sum(A)))
        
        A = BarycenterL2(images)
        #A = EuclideanBarycenter(images)
        print('Solution:', sum(A))
        #A = np.array(A).reshape(28,28)        
        #DrawDigit(A, 'Barycenter')
        #A = np.array(A)
        
        # Dump barycenter
        #FILEOUT = DIR+'Euclidean.csv'
        #with open(FILEOUT, "a") as output:
        #    writer = csv.writer(output, lineterminator='\n')
        #    writer.writerow(A)
