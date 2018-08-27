# -*- coding: utf-8 -*-
"""
Created on Sat Jun  2 07:45:40 2018

@author: gualandi
"""

import time
import csv
import numpy as np

from numpy import genfromtxt
import matplotlib.pyplot as plt

def DrawDigit(A, label=''):
    """ Draw single digit as a greyscale matrix"""
    # Uso la colormap 'gray' per avere la schacchiera in bianco&nero
    img = plt.imshow(A, cmap='gray_r')
    plt.xlabel(label)
    plt.show()


#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    my_data = genfromtxt('E:\\GitHub\\dotlib\\data\\barycenter_Euclidean.csv', delimiter=',')

    fig = plt.figure(figsize=(20,5))
    
    f, axarr = plt.subplots(2, 5)
    i = 0
    j = 0
    for row in my_data:
        A = np.array(row).reshape(28,28)
        axarr[j,i].imshow(A, cmap='coolwarm')
        axarr[j,i].set_xticklabels([])
        axarr[j,i].set_yticklabels([])
        axarr[j,i].axis('off')
        axarr[j,i].set_aspect('equal')
        print(i,j)
        i+=1        
        if i == 5:
            i = 0
            j = 1

    plt.subplots_adjust(wspace=0, hspace=0)        
    plt.tight_layout()
    plt.show()
        
