#-*- coding: utf-8 -*-
#
# Scritto da Sanfelice Antonio
# (c) 2010

import numpy as np
from pylab import plot,show
from Curve import Points,BSpline


class CurveFit(BSpline):
    '''
       Classe che offre i metodi necessari
       a determinare i vertici di un poligono di controllo
       tale da generare una curva BSpline che interpoli
       i punti della curva dati in input.
    '''

    def __init__(self,curvePoints,numCpoints,npts,k):
        '''
            Metodo costruttore
            INPUT:
            @param curvePoints lista contenente i punti della curva noti, nella forma [[x0,y0],[x1,y1],...,[xn,yn]]
                   oppure nome del file contenente i punti.
            @param numCpoints numero di vertici desiderato per il poligono di controllo
            @param npts numero di punti della curva da calcolare
            @param k ordine della curva
            @param x vettore knot o stringa indicante il tipo di vettore knot desiderato, 'periodic' o 'open'
        '''


        BSpline.__init__(self,np.zeros((numCpoints,2)),npts,k,'open')
        # se curvePoints è una stringa
        if isinstance(curvePoints,str):
            try:
                # provo a leggere i punti dal file il cui nome è indicato da curvePoints
                self.__dict__['curvePoints'] = np.loadtxt(curvePoints).view(Points)
            except IOError as detail:
                print "Errore nel leggere dal file {0}: {1}".format(curvePoints,detail)
                raise
        else:
            self.__dict__['curvePoints'] = np.asarray(curvePoints).view(Points)

        self.__dict__['cntrl'] = self.getControlPoints()
   

    def getParametersValue(self,pvec):
        ''' Metodo che calcola il valori del parametro associati
            ai punti noti della curva.
            @param pvec vettore contente i punti della curva
        '''     
        # determino il numero di punti noti della curva
        n = len(pvec)

        # inizializzo il vettore che alla fine dell'esecuzione conterrà i valori del parametro
        t = np.zeros(n,np.double)

        # dato che t0 = 0, inizio il ciclo da 1
        for i in range(1,n):

            # il valore del parametro associato a t[i] è uguale al valore del parametro
            # associato a p[i-1] più la distanza fra p[i] e p[i-1]
            t[i] = t[i-1] + pvec[i].distance(pvec[i-1])

        # divido il vettore per l'ultimo valore che contiene (che è il valore più grande) 
        # in modo da normalizzare il vettore nell'intervallo [0,1]
        t = t/t[-1]

        return t

    def getNMatrix(self,x,t,n):
        """
            Metodo che calcola la matrice N contenente
            le funzioni di base Ni,k(t)
            @param x vettore knot
            @param t valore del parametro
            @param n numero di punti di controllo
        """

        #scalo il vettore knot nell'intervallo [0,1]
        x = np.asarray(x,np.double) 
        x = x / x[-1]

        # determino il numero di punti
        numPoints = len(t)
        
        # inizializzo la matrice come matrice numPoints x n di zeri
        N = np.zeros((numPoints,n))
        
        # per i che va da 0 a numPoints-1
        for i in range(numPoints):
            # l'i-esima riga di N è uguale alle funzioni di base associate al parametro t[i]
            # si ricordi che computeBasis restituisce un vettore contenete le funzioni base 
            # Ni,k(u) per ogni i che va da 0 a n
            N[i] = self.computeBasis(n,self.k,x,t[i])
        return N
        

    def getControlPoints(self):
        '''
            Metodo che determina i vertici del poligono di controllo
        ''' 
        
        # calcolo i valori del parametro associati ai punti della curva noti
        t = self.getParametersValue(self.curvePoints)

        # determino la matrice N contenente le funzioni di base
        N = self.getNMatrix(self.x, t, self.n)

        # Pongo D = al vettore dei punti della curva noti.
        # assegnazione fatta solo per una quetione di leggibilià
        D = self.curvePoints

        # determino la dimensione della matrice N
        r,c = N.shape

        # se la matrice non è quadrata
        if(r != c):
            # moltiplico la matrice dei termini noti e la matrice N per la trasposta di N,
            # in modo da ottenere una matrice dei coefficienti quadrata
            D = np.dot(N.T,D)
            N = np.dot(N.T,N)
        
        # inizializzo la matrice B, che conterrà i punti di controllo
        B = np.zeros((self.n,2))

        # provo a risolvere il sistema di equazioni
        try:
            B = np.linalg.solve(N,D)

        # catturo l'eccezione dovuta ad un eventuale matrice singolare
        except np.linalg.LinAlgError as detail:
            print "Errore: la matrice N è singolare: {0}".format(detail)
        
        return B

    def plot(self,**args):
        plot(self.curvePoints[:,0],self.curvePoints[:,1],'go')
        BSpline.plot(self,**args)
        return

    def __setattr__(self,name,value):
        if name is 'curvePoints':
            if not isinstance(value,(list,tuple,np.ndarray)):
                raise TypeError("{0} deve essere una lista,tupla o istanza di numpy.ndarray, non {1}".format(name,type(value)))
            else:
                self.__dict__[name] = value
                self.__dict__['cntrl'] = self.getControlPoints()
                self.calculate()
        elif name is 'n':
            if not isinstance(value,int):
                raise TypeError("L'attributo n deve essere un intero, non un {0}".format(type(value)))
            else:
                self.__dict__[name] = value
                self.__dict__['cntrl'] = self.getControlPoints()
                self.calculate()
        else:
            BSpline.__setattr__(self,name,value)

