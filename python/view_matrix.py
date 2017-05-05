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

#------------------------------------------
#              MAIN ENTRY POINT
#------------------------------------------
if __name__ == "__main__":
    filename1 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data32_1001.csv'
    M1 = np.loadtxt(open(filename1, "rb"), delimiter=",")
    DrawMatrix(M1)

    filename2 = 'D:\Ricerca\DOTA\data\DOTmark_1.0\Data\ClassicImages\data32_1003.csv'
    M2 = np.loadtxt(open(filename2, "rb"), delimiter=",")
    DrawMatrix(M2)