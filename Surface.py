#-*- coding: utf-8 -*-
#
# Scritto da Sanfelice Antonio
# (c) 2010

"""
@package PyNurbs
@brief Insieme di classi che consentono di calcolare diversi tipi di superfici parametriche.

Questa parte della libreria consente di calcolare i seguenti tipi di superfici parametriche:
    - Bézier
    - BSpline
    - NURBS
"""
try:
    import numpy as np
    from Util import *
    from Curve import Points,Bezier,Spline,BSpline,Nurbs
    from mpl_toolkits.mplot3d.axes3d import Axes3D
    import pylab as pl
except ImportError as detail:
    print "Errore nell'importare le librerie:\n",detail


class Surface:
    """
        Classe interfaccia che definisce i parametri e i metodi
        di base che devono essere implementati per modellare
        una superfice.
    """

    def __init__(self, cntrl, npts, mpts):
        """
                Metodo Costruttore
                
                INPUT:
                @param cntrl vettore dei punti di controllo da cui generare la superfice
                @param npts numero di linee parametriche lungo la direzione u
                @param mpts numero di linee parametriche lungo la direzione w
        """    
        
        self.__dict__['mpts'] = mpts
        
        #self.points ora ha 3 dimensioni e dovrà contenere npts * mpts punti
        self.__dict__['points'] = np.zeros((self.npts * self.mpts, 3)).view(Points)

        if isinstance(cntrl,str):
            self.loadFromFile(cntrl)

        elif not(self.cntrl.ndim == 3):
            raise PyNurbsError("I punti di controllo della superfice devono essere tridimensionali")


    def __setattr__(self, name, value):
        """
            Metodo che sovrascrive l'operatore di assegnazione
            per gli attributi, permettendo di eseguire funzioni
            in caso di cambiamenti
            @param name il nome dell'attributo che deve essere settato
            @param value il valore che si vuole assegnare al parametro 
        """

        
        if name is 'mpts':
            self.__dict__[name] = value
            self.calculate()

    def loadFromFile(self, nomeFile):
        if not isinstance(nomeFile,str):
            raise TypeError('nomeFile non è una stringa ma un {0}'.format(type(nomeFile)))

        try:
            tmp = np.loadtxt(nomeFile).view(Points)
        except IOError as detail:
            print "Errore nel caricare i dati dal file {0}: {1}".format(nomeFile,detail)

        n, m, k = tmp[-1]
        return tmp[:-1].reshape((n,m,k))
                            

        
    
    def plot(self, plotcp=False, eps=1, dpi=100):
        '''
            Metodo che effettua il plot della superfice
            INPUT:
            @param dpi dot per inch desiderati per l'immagine del plot
            @param eps spaziatura fra una linea e la successiva della superfice
        
        
        se plotcp == True mostra nel plot anche il poligono di controllo
        che ha generato la superfice
        eps indica quanto preciso si vuole la superfice, più eps si avvicina ad 1
        più strette saranno le mesh
        dpi indica la grandezza che si desidera per il plot, in punti per pollice
        '''
        

        X = self.points[:, 0].reshape(self.npts, self.mpts)
        Y = self.points[:, 1].reshape(self.npts, self.mpts)
        Z = self.points[:, 2].reshape(self.npts, self.mpts)

        fig = pl.figure(1, dpi = dpi)
        ax = Axes3D(fig)

        ax.plot_surface(X, Y, Z,cstride = eps, rstride = eps)
        if(plotcp):
            ax.plot_wireframe(self.cntrl[:, :, 0], self.cntrl[:, :, 1], self.cntrl[:, :, 2], color="#cc0000")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        pl.axis('equal')
        pl.show()
        return


class BezSurf(Bezier, Surface):
    """
        Classe che modella una superfice di Bézier
        
    """
    def __init__(self, cntrl, npts, mpts):
        '''
            Metodo costruttore
            @param cntrl matrice contenente le coordinate dei vertici del poligono di controllo
            @param npts numero di curve parametriche da calcolare lungo la direzione u
            @param mpts numero di curve parametriche da calcolare lungo la direzione w
        '''
        Bezier.__init__(self, cntrl, npts)
        Surface.__init__(self, cntrl, npts, mpts)
        
        self.__dict__['n'], self.__dict__['m'] = self.cntrl.shape[:2]
        
        self.__dict__['u_degree'] = self.n - 1
        self.__dict__['w_degree'] = self.m - 1

    def loadFromFile(self, nomeFile):
        Surface.loadFromFile(self,nomeFile)


    def calcWithBernstein(self):
        '''
            Metodo che calcola la superfice di Bézier utilizzando la definizione formale
        '''

        self.points = np.zeros((self.npts * self.mpts, 3)) 

        # determino i valori dei parametri
        upts = np.linspace(0, 1, self.npts)
        wpts = np.linspace(0, 1, self.mpts)

        # inizializzo un contatore di riga
        iCount = 0
        
        # per ogni valore del parametro u
        for u in upts:
            # per ogni valore del parametro w
            for w in wpts:
                # per ogni i da 0 a n-1
                for i in xrange(self.n):
                    # calcolo la funzione di base per la direzione u
                    jni = self.bernstein(self.u_degree, i, u)

                    # per ogni j da 0 a m-1
                    for j in xrange(self.m):
                        # calcolo la funzione di base per la direzione w
                        kmj = self.bernstein(self.w_degree, j, w)
                        # l'iCount-esimo punto della superfice è uguale al prodotto tra il punto di controllo
                        # i,j per la funzione di base i nella direzione u e la funzione di base j nella direzione w
                        self.points[iCount] += self.cntrl[i, j] * jni * kmj
                iCount+=1

    def __setattr__(self, name, value):
        if name in('u_degree', 'w_degree'):
            print "Impossibile impostare manualmente il valore di {0}".format(name)
        else:
            Bezier.__setattr__(self, name, value)

    def plot(self, plotcp = False):
        Surface.plot(self, plotcp)
        return


class BSplineSurf(Surface, BSpline):
    """
        Classe che modella una superfice Spline
    """
    def __init__(self, cntrl, npts, mpts, k, l, x, y):
        """
            Costruttore
            @param cntrl matrice dei punti di controllo
            @param npts numero di punti della superfice da calcolare lungo la direzione u
            @param mpts numero di punti della superfice da calcolare lungo la direzione v
            @param k ordine della bspline lungo la direzione u
            @param l ordine della bspline lungo la direzione v
            @param x vettore knot per la direzione u o stringa di valore 'open' oppure 'periodic'
            @param y vettore knot per la direzione w o stringa di valore 'open' oppure 'periodic'
        """
        BSpline.__init__(self, cntrl, npts, k, x)
        Surface.__init__(self, cntrl, npts, mpts)
        self.__dict__['n'], self.__dict__['m'] = self.cntrl.shape[:2]
        self.__dict__['l'] = l
        self.__dict__['y'] = self.knot(y, self.l, self.m)
    
    def calculate(self):

        self.points = np.zeros((self.npts * self.mpts, 3)).view(Points)

        # determino i valori per i parametri u e v
        uvalues = np.linspace(self.x[self.k - 1], self.x[self.n], self.npts)
        vvalues = np.linspace(self.y[self.l - 1], self.y[self.m], self.mpts)

        ptCount = 0
        for u in uvalues:
            # calcolo le basi per la direzione u, Ni,k(u)
            u_basis = self.computeBasis(self.n, self.k, self.x, u)

            for v in vvalues:
                # calcolo la funzione di base per la direzione v, Mj,l(v)
                v_basis = self.computeBasis(self.m, self.l, self.y, v)

                # calcolo le coordinate del punto, moltiplicando la matrice  
                # dei punti di controllo per i vettori delle funzioni di base
                x = np.dot(u_basis, np.dot(self.cntrl[:,:,0], v_basis))
                y = np.dot(u_basis, np.dot(self.cntrl[:,:,1], v_basis))
                z = np.dot(u_basis, np.dot(self.cntrl[:,:,2], v_basis))
                self.points[ptCount] = Points([[x,y,z]])
                ptCount+=1


    def calculateWithMatrixNotation(self):
        """
            Metodo che calcola i punti della superfice
            tramite il prodotto matriciale
            $[u][N][cntrl][M]^T[v]$
        """

        # dato che il calcolo viene effettuato generando delle "sottosuperfici" k*l
        # devo determinare il numero di tratti in entrambe le direzioni u e v
        # in modo da determinare il numero di punti da calcolare per ogni singola
        # sotto superfice
        num_tratti_u = self.n - self.k + 1
        num_tratti_v = self.m - self.l + 1

        # divido il numero di punti da calcolare per il numero di tratti lungo le due direzioni
        scaled_npts = self.npts / num_tratti_u
        scaled_mpts = self.mpts / num_tratti_v

        # aggiorno npts e mpts, dato che nel determinare il numero di punti
        # da calcolare per singolo tratto, ci sono stati degli arrotondamenti
        # e quindi non verranno calcolati esattamente npts e mpts punti lungo le due direzioni
        # lasciare invariati npts e mpts porterebbe alla visualizzazione di artefatti durante
        # plotting
        self.__dict__['npts'] = scaled_npts * num_tratti_u
        self.__dict__['mpts'] = scaled_mpts * num_tratti_v
    
        # aggiorno la grandezza del vettore points
        self.__dict__['points'] = np.zeros((self.npts * self.mpts,3))

        # determimo le matrici di base per entrambe le direzioni, N e M
        N = self.computePeriodicMatrix(self.k)
        M = self.computePeriodicMatrix(self.l)

        # determino i valori dei parametri per cui calcolare la curva
        u = np.arange(0, 1, 1/float(scaled_npts))
        v = np.arange(0, 1, 1/float(scaled_mpts))

        ptCount = 0
        for s in xrange(num_tratti_u):
        
            # determino la sottomatrice di controllo contenente le righe dalla s alla s+k
            cntrl_subnet_rows = self.cntrl[s:s + self.k, :, :]
    
            # per ogni u
            for i in u:
                # determino il polinomio in u
                ipoly = upoly(i, self.k)
                left = np.dot(ipoly, N)

                for t in xrange(num_tratti_v):
                    # determino la sottomatrice di cntrl_subnet_rows contenente solo le colonne dalla t alla t+l
                    cntrl_subnet = cntrl_subnet_rows[:, t:t + self.l, :]

                    # calcolo la prima parte del prodotto matriciale
                    left_x = np.dot(left, np.dot(cntrl_subnet[:, :, 0], M.T))        
                    left_y = np.dot(left, np.dot(cntrl_subnet[:, :, 1], M.T))
                    left_z = np.dot(left, np.dot(cntrl_subnet[:, :, 2], M.T))
                    
                    # for j in v
                    for j in v:

                        # determino il polinomio in v
                        jpoly = upoly(j, self.l)

                        # completo il calcolo moltiplicando la parte sinistra per la parte mancante,
                        # ovvero il polinomio in v
                        x = np.dot(left_x, jpoly)
                        y = np.dot(left_y, jpoly)
                        z = np.dot(left_z, jpoly)

                        # aggiungo le coordinate x,y,z trovate ai punti della superfice
                        self.points[ptCount] = Points([[x, y, z]])
                        ptCount += 1
                  

    def __setattr__(self, name, value):
        if name in ('k','l'):
            maxOrder = {'k': self.n, 'l': self.m}[name]
            if not isinstance(value, int):
                raise TypeError("{0} deve essere di tipo int, non {1}".format(name, type(value)))
            elif value > maxOrder:
                raise ValueError("{0} può essere al più {1} con l'attuale poligono di controllo".format(name, maxOrder))
        else:
            BSpline.__setattr__(self, name, value)
            Surface.__setattr__(self, name, value)


class NurbsSurf(BSplineSurf):
    """
        Classe che modella una superfice NURBs
    """
    def __init__(self,cntrl,npts,mpts,k,l,x,y):
        '''
            Metodo costruttore
            @param cntrl array nxmx4 contenente le coordinate omogenee dei punti di controllo, pesi compresi
            @param npts numero di linee parametriche lungo la direzione u
            @param npts numero di linee parametriche lungo la direzione w
            @param k ordine della b-spline razionale lungo la direzione u
            @param l ordine della b-spline razionale lungo la direzione w
            @param x vettore knot per la direzione u, si possono passare i valori 'open', 'periodic' o direttamente il vettore
            @param y vettore knot per la direzione w, si possono passare i valori 'open', 'periodic' o direttamente il vettore
        '''
        BSplineSurf.__init__(self,cntrl,npts,mpts,k,l,x,y)
        
        # variabile test che serve a controllare se uno dei parametri della superfice è stato cambiato
        # utile nell'algoritmo che ricava la superfice nurbs per evitare calcoli inutili 
        self.__dict__['itest'] = -1
        
        # array che conterranno le funzioni di base per ogni punto della superfice lungo
        # le due direzioni
        self.__dict__['niku'] = []
        self.__dict__['mjlw'] = []
        
        # array che conterrà la i valori della funzione Sum(u,w) per ogni valore di u e w
        # la funzione sum(u,w) è la funzione che compare al denominatore nel calcolo della base
        # razionale delle nurbs
        self.__dict__['sumuw'] = []

        # variabile booleana che controlla se c'è stato un cambiamento nei pesi
        self.__dict__['weightsAreChanged'] = True

        # se la matrice dei punti di controllo non ha 4 coordinate per punto
        if not(self.cntrl.shape[2] == 4):
            raise PyNurbsError("I punti di controllo della superfice devono essere quadridimensionali")

    def changeWeight(self,i,j,h):
        '''
            Metodo che consente di variare il singolo peso di un punto di controllo
            Sebbene sia possibile effettuare questo cambio accedendo direttamente 
            all'attributo cntrl, ciò porterebbe a dei problemi con l'algoritmo
            utilizzato dal metodo calculate. Questo metodo infatti segnala se 
            c'è stato un cambiamento nei pesi.
        '''
        self.cntrl[i,j,3] = h
        self.__dict__['weightsAreChanged'] = True

    def calculate(self):
        # se uno fra n, m, k, l, npts, mpts cambia, ricalcolo le funzioni di base 
        if self.itest != (self.n + self.m + self.npts + self.mpts + self.k + self.l):
            
            upts = np.linspace(self.x[self.k - 1], self.x[self.n], self.npts)
            wpts = np.linspace(self.y[self.l - 1], self.y[self.m], self.mpts)
            
            # calcolo le funzioni di base per ogni valore di u e le salvo in un array
            for u in upts:
                base = self.computeBasis(self.n, self.k, self.x, u)
                self.niku.append(base)

            # calcolo le funzioni di base per ogni valore di w e le salvo in un array
            for w in wpts:
                base = self.computeBasis(self.m, self.l, self.y, w)
                self.mjlw.append(base)

            self.itest = self.n + self.m + self.k + self.l + self.npts + self.mpts

        if self.weightsAreChanged:
            # calcolo la funzione somma di prodotti di funzioni di base per pesi dei punti di controlli
            for nbasis in self.niku:
                for mbasis in self.mjlw:
                    s = 0
                    for i in xrange(self.n):
                        for j in xrange(self.m):
                            s += (self.cntrl[i,j,3] * nbasis[i] * mbasis[j])
                    # conservo il reciproco in modo da non effettuare divisoni durante il loop principale
                    # per il calcolo della superfice
                    self.sumuw.append(1./s)
            self.__dict__['weightsAreChanged'] = False
        
        ptCount = 0
        for nbasis in self.niku:
            for mbasis in self.mjlw:    
                for i in range(self.n):
                    # se la funzione di base i-esima è 0 evito il calcolo
                    if nbasis[i] != 0:
                        for j in range(self.m):
                            # se la funzione di base j-esima è 0 evito il calcolo
                            if mbasis[j] != 0:
                                # calcolo la base razionale (self.cntrl[i,j,3] è il peso del punto di controllo i,j)
                                rbasis = self.cntrl[i,j,3] * nbasis[i] * mbasis[j] * self.sumuw[ptCount]
                                    
                                # aggiorno il punto della superfice che si sta calcolando
                                # si noti che self.cntrl[i,j,:3] contiene solo le coordinate x,y,z e non il peso

                                self.points[ptCount] += Points(self.cntrl[i,j,:3] * rbasis)
                                    
                ptCount += 1

        
