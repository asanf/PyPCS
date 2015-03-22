#-*- coding: utf-8 -*-
#
# Scritto da Sanfelice Antonio
# (c) 2010

try:
    import numpy as np
    from Curve import Points
    from Surface import BSplineSurf
    from CurveFitting import CurveFit
    import pylab as pl
    from mpl_toolkits.mplot3d.axes3d import Axes3D
except ImportError as detail:
    print "Errore nell'importare delle librerie: {0}".format(detail)

class SurfaceFit(CurveFit,BSplineSurf):
    """
        Classe che consente di determinare
        una superfice BSpline che approssima
        dei punti dati in input
    """   
    
    def __init__(self, surfacePoints, n, m, npts, mpts, k, l):

        '''
            Metodo Costruttore:
            INPUT:
            @param surfacePoints punti noti appartenenti alla superfice o nome del file che li contiene
            @param n numero di punti di controllo da calcolare lungo la direzione u
            @param m numero di punti di controllo da calcolare lungo la direzione v
            @param npts numero di punti della superfice da calcolare lungo la direzione u
            @param mpts numero di punti della superfice da calcolare lungo la direzione v
            @param k ordine della BSpline per la direzione u
            @param l ordine della BSpline per la direzione v
        '''    
        BSplineSurf.__init__(self, np.zeros((n, m, 3)), npts, mpts, k, l, 'open','open')
        if isinstance(surfacePoints, str):
            try:
                self.__dict__['surfacePoints'] = self.loadFromFile(surfacePoints)
            except IOError as detail:
                print "Errore nel leggere i dati dal file {0}: {1}".format(surfacePoints,detail)
                raise
        elif not isinstance(surfacePoints, (list, tuple, np.ndarray)):
            raise TypeError("Errore, il parametro surfacePoints deve essere una lista,\
                             una tupla o un istanza di np.ndarray, non di tipo {0}".format(type(surfacePoints)))
        else:
            self.__dict__['surfacePoints'] = np.asarray(surfacePoints).view(Points)

        self.__dict__['cntrl'] = self.getControlPoints()

    def loadFromFile(self, nomeFile):
        BSplineSurf.loadFromFile(self, nomeFile)

    def getParametersMatrix(self):
        '''
           Metodo che calcola i valori del parametro 
           associati a tutti i punti contenuti in self.surfacePoints.
           La matrice ha dimensione r x s, dove r è il numero di punti
           della superfice noti lungo la direzione u, e s è il numero di 
           punti della superfice noti lungo la direzione w.
        '''
        # determino la dimensione della matrice dei punti noti della superfice
        r,s = self.surfacePoints.shape[:2]
        
        # inizializzo uw come matrice r x s x 2 di zeri
        # la dimensione 2 è dovuta al fatto che l'ingresso i,j della matrice 
        # conterrà il valore del parametro ui e del parametro wj
        uv = np.zeros((r, s, 2),np.double)

        # per ogni colona della matrice self.surfacePoints
        for i in range(s):

            # determino il valore del parametro u
            uv[:,i,0] = self.getParametersValue(self.surfacePoints[:,i])
        
        # per ogni riga della matrice self.surfacePoints
        for i in range(r):

            # determino il valore del parametro w
            uv[i,:,1] = self.getParametersValue(self.surfacePoints[i,:])

        return uv

    def getNMatrix(self):
        ''' 
            Metodo che calcola la matrice contenente i prodotti delle funzioni di base
            Nik*Mjl, che verrà utilizzata come matrice dei coefficienti nella risoluzione del sistema
            di equazioni N cntrl = D           
        '''

        # scalo i vettori knot nell'intervallo [0,1]
        x = np.asarray(self.x,np.double) 
        x = x / x[-1]

        y = np.asarray(self.y,np.double)
        y = y / y[-1]

        
        # determino righe e colonne della matrice self.surfacePoints
        r,s = self.surfacePoints.shape[:2]

        # inizializzo N come matrice di zeri r*s x n*m
        N = np.zeros((r * s, self.n * self.m), np.double)

        # determino i valori dei parametri u e w per i punti noti della superfice
        uv = self.getParametersMatrix()

        # inizializzo un contatore riga per la matrice C
        rowCount = 0
        
        # per uparam che va da 0 a r-1
        for uparam in range(r):

            # per wparam che va da 0 a s-1
            for wparam in range(s):

                # Qui calcolo una riga della matrice C, contenente tutti i prodotti Ni,k(uparam)*Mj,l(wparam)

                # estraggo dalla matrice uw i valori u e w che mi servono, più per una questione di leggibilità
                # che per altro
                u = uv[uparam, wparam, 0]
                v = uv[uparam, wparam, 1]
        
                #determino le funzioni di base Nk e Ml
                Nbasis = self.computeBasis(self.n, self.k, x, u)
                Mbasis = self.computeBasis(self.m, self.l, y, v)
                
                # inizializzo un contatore colonna per la matrice C
                colCount = 0
                
                # per ogni funzione di base N
                for Nik in Nbasis:

                    # e per ogni funzione di base M
                    for Mjl in Mbasis:

                        # calcolo l'entrata della matrice N
                        N[rowCount, colCount] = Nik * Mjl

                        colCount+=1

                rowCount+=1

        return N


    def getControlPoints(self):
        '''
            Metodo che calcola i punti di controllo necessari per calcolare la superfice
        '''
        r, s = self.surfacePoints.shape[:2]

        # nel seguente frammento di codice, opero dei redimensionamenti alle matrici
        # che compongono il sistema di equazioni. Dato che N è una matrice r*s x n*m
        # D deve essere una matrice r*s x 3, e B deve essere una matrice n*m x 3
        # ricordiamo che r x s è la dimensione della matrice contenente i punti noti 
        # della superfice, mentre n e m sono il numero di punti di controllo
        # lungo le direzioni u e v
        
        D = self.surfacePoints.reshape(r * s, 3)
        cntrl = np.zeros((self.n * self.m, 3))
        
        # determino la matrice C        
        N = self.getNMatrix()

        # recupero il numero di righe e colonne della matrice C
        rows,columns = N.shape
    
        #se la matrice non è quadrata        
        if(rows != columns):
            print "La matrice N non è quadrata"
            # Aggiusto le matrici C e D moltiplicandole per C trasposto
            D = np.dot(N.T,D)
            N = np.dot(N.T,N)
        
        # determino B risolvendo il sistema N cntrl = D
        try:
            cntrl = np.linalg.solve(N,D)
        except np.linalg.LinAlgError as detail:
            print "Matrice Singolare: {0}".format(detail)

        # sistemo le dimensioni di cntrl, dato che nel metodo calculate() cntrl deve essere una matrice n x m x 3
        cntrl = cntrl.reshape((self.n, self.m, 3))
        return cntrl

    def calculate(self):
        BSplineSurf.calculate(self)

    def calculateWithMatrixNotation(self):
        BSplineSurf.calculateWithMatrixNotation(self)

    def plot(self, **args):
        BSplineSurf.plot(self, **args)

