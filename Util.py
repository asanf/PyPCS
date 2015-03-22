# -*- coding: utf-8 -*-
#
# Scritto da Sanfelice Antonio
# (c) 2010

import numpy as np

class PyNurbsError(Exception):
    def __init__(self,description):
        self.description = description
    def __str__(self):
        return self.description


def upoly(u,k):
    """
        Metodo che restituisce un array di k elementi
        contenente le potenze di u da 0 a k-1
        nella forma array([pow(u,k-1), pow(u,k-2),...,pow(u,2),u,1])
    """
    try:
        u = float(u)
        k = int(k)
    except ValueError as detail:
        print "Errore di tipo: ",detail
        return
    # restituisco un array contenente pow(u,i) per ogni i da k-1 a 2,
    # concatenato a u,1. (evito di fargli calcolare le potenze per i = 1,0)
    return np.asarray([pow(u,i) for i in xrange(k-1,1,-1)]+[u,1],np.double)

def deg2rad(angle):
    return np.pi * angle / 180.

def rad2deg(angle):
    return angle * 180./np.pi

def htm(rot_x = 0, rot_y = 0, rot_z = 0, t_x = 0, t_y = 0, t_z = 0, px = 1, py = 1, pz = 1, scale = 0):
    
    """
        Metodo che restituisce una matrice di trasformazione omogenea
        @param rot_x angolo di rotazione intorno all'asse x
        @param rot_y angolo di rotazione intorno all'asse y
        @param rot_z angolo di rotazione intorno all'asse z
        @param t_x unità traslazione lungo l'asse x
        @param t_y unità traslazione lungo l'asse y
        @param t_z unità traslazione lungo l'asse z
        @param px cambio di prospettiva rispetto a x
        @param py cambio di prospettiva rispetto a y
        @param pz cambio di prospettiva rispetto a z
        @param scale fattore di scalatura
    """
    mat = np.eye(4)
    mat[3,3] = scale
    mat[3,:3] = np.asarray(px, py, pz)
    mat[:3,3] = np.asarray(t_x, t_y, t_z)

    mat[1,1] = mat[2,2] = np.cos(rot_x)
    sn = np.sin(rot_x)
    mat[1,2] = -sn
    mat[2,1] = sn

    cs = np.cos(rot_y)
    mat[0,0] = cs
    mat[2,2] += cs
    sn = np.sin(rot_y)
    mat[0,2] = sn
    mat[2,0] = -sn

    cs = np.cos(rot_z)
    mat[0,0] += cs
    mat[1,1] += cs
    sn = np.sin(rot_z)
    mat[0,1] = -sn 
    mat[1,0] = sn

    return mat


