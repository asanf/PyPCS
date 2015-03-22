
#-*- coding: utf-8 -*-
#
# Scritto da Sanfelice Antonio
# (c) 2010
"""
@package PyNurbs
@brief Insieme di classi che fornisce metodi per calcolare curve parametriche. 

Questa parte della libreria consente di calcolare diversi tipi di curve
parametriche, in particolare:
    - Natural Cubic Spline
    - Hermite Spline
    - Cardinal Spline
    - Bézier Curve
    - BSpline
    - NURBS
"""
try:
    import numpy as np
    from Util import *
    from scipy.misc import comb
    from scipy import factorial as fact
    import pylab as pl
except ImportError as detail:
    print "Errore nell'importare le librerie:\n",detail
    
    


class Points(np.ndarray):
    """
        Classe che modella un array di punti
        Aggiunge a numpy.ndarray i metodi "distance","chordlenghts" e "convexComb"
    """
    def __new__(subclass,data,dtype = np.double):
        obj = np.asarray(data,dtype).view(subclass)
        return obj
    

    def distance(self,p):
        """
            Metodo che calcola la distanza euclidea fra due punti.
            Si noti che con questa notazione il numero di dimensioni
            è irrilevante.
            @param p il punto con cui calcolare la distanza
        """
        return np.sqrt(sum(pow(p-self,2)))


    def chordLenght(self,i=0,j=None):
        """
            Metodo che calcola la somma delle distanze fra
            l'i-esimo e il j-esimo punto dell'array
            Chiamato senza argomenti calcola la somma delle distanze
            dal primo all'ultimo punto
            @param i l'indice del primo punto della sequenza
            @param j l'indice dell'ultimo punto della sequenza
        """
        return sum([self[k].distance(self[k+1]) for k in xrange(len(self[i:j])-1)])



    def convexComb(self,p,u):
        """
            Metodo che restituisce la combinazione convessa u
            con un punto p
            @param p il punto con cui si vuole calcolare la combinazione convessa
            @param u il valore della combinazione, compreso fra 0 e 1
        """
        if u < 0 or u > 1:
            raise ValueError("il parametro u deve essere compreso fra 0 e 1")
        return (1-u)*self + u*p



class Curve:
    """
        Classe interfaccia che definisce i parametri e i metodi
        di base che devono essere implementati per modellare
        una curva.
    """

    def __init__(self,cntrl,npts):
    
        """
                Costruttore
                
                INPUT:
                @param cntrl vettore dei punti di controllo o nome del file che contiene i punti
                @param npts numero di punti della curva da calcolare
                
                ATTRIBUTI:
                @var points punti che costituiscono la curva
        """
        
        if isinstance(cntrl,str):
            self.__dict__['cntrl'] = self.loadFromFile(cntrl)
        else:
            try:
                self.__dict__['cntrl'] = Points(cntrl)
            except Exception as detail:
                raise PyNurbsError("Errore formato punti di controllo:{0}".format(detail))
                        
        self.__dict__['npts'] = npts
        #inzializzo il vettore dei punti come vettore di zeri npts x 2
        self.__dict__['points'] = np.zeros((self.npts, 2)).view(Points)
        

    def calculate(self):
        """
            Metodo che effettua i calcoli necessari a determinare
            self.points
        """
        pass


    def loadFromFile(self,nomeFile):
        """
            Metodo che carica da file i punti di controllo
        """
        if not isinstance(nomeFile,str):
            raise TypeError("Il nome del file deve essere una stringa, non di tipo {0}".format(type(nomeFile)))
        try:
            return np.loadtxt(nomeFile).view(Points)
        except IOError as detail:
            print "Errore nel caricare i dati dal file {0}: {1}".format(nomeFile,detail)
            raise
        
        

    def plot(self, plotcp = True):
        """
            Metodo che esegue il plot dei valori contenuti in self.points
            @param plotcp flag booleano, se True allora visualizza anche le linee di collegamento fra i punti di controllo
            @param show flag booleano, se True visualizza il plot(se false, utile per overwriting,show() deve essere l'ultima operazone di plot)
        """
        
        pl.plot(self.points[:, 0], self.points[:, 1])

        if plotcp:
            pl.plot(self.cntrl[:, 0],self.cntrl[:, 1],'ro')
            pl.plot(self.cntrl[:, 0], self.cntrl[:, 1],'r--')

        pl.axis('equal')
        pl.show()

        return

    def __add__(self, c):
        """
            Metodo che sovrascrive l'operatore di somma
        """

        if not(isinstance(c,type(self))):
            raise TypeError("Il secondo operando deve essere un istanza di {0}".format(type(self)))

        # copio i punti di controllo della seconda curva, dato che devo modificarli per
        # "avvicinare" la seconda curva alla prima
        other_curve = c.cntrl.copy()      
        
        # calcolo la differenza di posizione fra l'ultimo punto della prima curva e il primo della seconda
        diff = self.cntrl[-1] - other_curve[0]

        # traslo i singoli punti della seconda curva 
        for pt in other_curve:
            pt += diff

        # creo un nuovo insieme di punti di controllo unendo i punti di controllo delle due curve
        new_cntrl_points = np.append(self.cntrl,other_curve[1:],0).view(Points)

        # deepcopy è un metodo di python (contenuto nel pacchetto copy) che permette di duplicare un oggetto
        from copy import deepcopy        

        new_curve = deepcopy(self)
        new_curve.cntrl = new_cntrl_points
    
        return new_curve
            
    
    def __setattr__(self,name,value):
        """
            Metodo che sovrascrive l'operatore di assegnazione
            per gli attributi, permettendo di eseguire funzioni
            in caso di cambiamenti
            @param name il nome dell'attributo che deve essere settato
            @param value il valore che si vuole assegnare al parametro 
        """

        # se il nome dell'attributo non è npts o cntrl, ma è comunque un
        # attributo esistente, setta attributo = value
        if name not in ('npts','cntrl') and self.__dict__.has_key(name):
            self.__dict__[name]=value
            return
        # se l'attributo da modificare è npts
        elif name is 'npts':
            # controllo se value è un intero
            if not isinstance(value,int):
                raise TypeError("npts deve essere un intero")
            else:
                self.__dict__[name] = value
        elif name is 'cntrl':
            if not (isinstance(value,Points) and value.ndim == self.cntrl.ndim):
                raise PyNurbsError("cntrl deve essere un istanza di numpy.ndarray di %d dimensioni"%(self.cntrl.ndim))
            else:
                self.__dict__[name] = value
        elif name is 'points':
            print "mannaggia"
        
                


class Bezier(Curve):
    """
        Classe che modella una curva di Bézier
    """
    def __init__(self,cntrl,npts):
        """
            Metodo costruttore: cntrl va immesso nella forma [(x0,y0),(x1,y1),...,(xn,yn)]
            INPUT:
            @param cntrl punti di controllo da cui ottenere la curva
            @param npts numero di punti da calcolare

            ATTRIBUTI:
            @var deg grado della curva (num_punti_di_controllo -1)
        """

        Curve.__init__(self,cntrl,npts)
        self.__dict__['deg'] = len(self.cntrl)-1
    
    def bernstein(self,n,i,t):
        """
            Metodo che calcola la funzione di base di bernstein
            @param n il grado della curva (ordine della curva di Bézier)
            @param i l'indice della funzione di base da calcolare
            @param t il valore del parametro
        """

        return comb(n, i)*pow(t, i)*pow(1-t, n-i)

    def calcWithBernstein(self):

        """
            Metodo che calcola la curva come da definizione formale
        """

        tvalues = np.linspace(0,1,self.npts)
        for i,t in enumerate(tvalues):
            # l'i-esimo punto è uguale alla sommatoria per i da 0 a deg di cntrl[i]*bernstein(deg,i,t)
            self.points[i] = sum([self.cntrl[j]*self.bernstein(self.deg,j,t) for j in xrange(self.deg+1)])

    def deCasteljau(self,cntrl,deg,t):

        """
            Metodo che calcola il punto della curva corrispondente al valore di parametro t
            @param t il valore del parametro corrispondente al punto desiderato
        """
            
        #inizializzo un vettore points contenente una prima combinazione convessa in t dei punti di controllo
        tmp = Points([cntrl[i].convexComb(cntrl[i+1],t) for i in xrange(deg)])
        
        #itero il calcolo delle combinazioni lineari fino ad ottenere un unico punto
        for i in xrange(deg -1):
            for j in xrange(deg -1 -i):
                tmp[j] = tmp[j].convexComb(tmp[j+1],t)

        #points[0] contiene il punto appartenente alla curva di Bezier
        return tmp[0]

    def calculate(self):
        #points è un array contenente gli output di deCasteljau per ogni valore di t
        self.__dict__['points'] = Points([self.deCasteljau(self.cntrl,self.deg,t) for t in np.linspace(0,1,self.npts)])

    # sovrascrivo l'operatore di assegnamento in modo da proteggere deg da modifiche
    def __setattr__(self,name,value):
        if name is 'deg':
            raise PyNurbsError("Il grado della curva non può essere modificato manualmente")
        elif name is 'cntrl':
            Curve.__setattr__(self,name,value)
            self.__dict__['deg'] = len(self.cntrl)-1
        Curve.__setattr__(self,name,value)
            


class Spline(Curve):
    """
        Classe che modella una generica spline
        utilizza la notazione matriciale per effettuare
        i calcoli
    """

    def __init__(self,cntrl,npts,rpp,k):
        """
            Metodo costruttore
            INPUT:
            @param cntrl vettore dei punti di controllo da cui ottenere la curva
            @param npts numero di punti della curva da calcolare
            @param rpp numero di righe della matrice di base dedicate ad un singolo\
                   punto di controllo
            @param k ordine della curva desiderato (grado desiderato + 1)

            ATTRIBUTI:
            @var B matrice di base che identifica il tipo di spline
            @var coeff coefficienti del polinomio in u
        """

        Curve.__init__(self,cntrl,npts)
        self.__dict__['k'] = k
        self.__dict__['rpp'] = rpp
        self.__dict__['coeff'] = np.array([[]],np.double).reshape(0,2)
        #inizializzo la matrice di base come matrice di base k x k
        self.__dict__['B'] = np.eye(self.k)

    def calculate(self):
        """
            Metodo che effettua il calcolo.
            L'algoritmo determina i coefficienti 'a' del polinomio
            ak-1*u^k-1 + ak-2*u^k-2 + ... + a2*u^2 + a1*u + a0
            fatto ciò valuta i polinomi ottenendo i punti della curva
        """


        self.points = np.zeros((self.npts,2))

        # Il calcolo deve fermarsi k punti prima dell'ultimo
        # in modo che il dot product non vada oltre la grandezza dell'array
        # numTratti è un contatore che conterà il numero di tratti di cui
        # la spline è composta
        numTratti = 0
        for i in xrange(0,len(self.cntrl)-self.k+1,self.rpp):
            numTratti+=1
            self.coeff = np.append(self.coeff,np.dot(self.B,self.cntrl[i:i+self.k]),0)


        # determino l'eps nello spazio del parametro u
        eps = float(numTratti)/self.npts


        # determino il valore di u. Utilizzo "arange" in luogo di linspace perché
        # arange calcola un intervallo [a,b[ mentre linspace calcola un intervallo chiuso
        # che comprende gli estremi. dato che la fine di un tratto coincide con l'inizio del
        # successivo, utilizzare arange evita di calcolare due volte lo stesso punto
        upts = np.arange(0,1,eps)

        pts = np.array([[]],np.double).reshape(0,2)
        for i in xrange(numTratti):
            c = i*self.k
            # Calcolo i punti dell'intero tratto effettuando il dot product fra il polinomio in u
            # e la sotto matrice costituita dalle righe dalla i alla i+self.k della matrice dei
            # coefficienti e lo aggiungo alla lista dei punti
            pts = np.append(pts,[np.dot(upoly(u,self.k),self.coeff[c:c+self.k]) for u in upts],0)
        self.points = pts
        return



class NaturalCubicSpline(Spline):
    """
        Classe che modella una spline cubica naturale
    """
    def __init__(self,cntrl,npts):
        """
            Costruttore
            @param cntrl punti di controllo indicanti punto posizione, vettore della derivata prima e 
                         della derivata seconda nella forma [(x0,y0),(d'0x,d'0y),(d''0x,d''0y),(x1,y1),...]
            @param npts il numero di punti della curva che si vuole calcolare
        """

        Spline.__init__(self,cntrl,npts,3,4)
        self.__dict__['B'] = np.array([[-1,-1,-.5,1],[0,0,.5,0],[0,1,0,0],[1,0,0,0]],np.double)

    def plot(self):
        """
            Metodo che effettua il plot della curva, mostrando le condizioni di controllo
            (derivata prima e seconda) come frecce
        """

        positionPoints = self.cntrl[::self.rpp]
        firstDerivatives = self.cntrl[1::self.rpp]
        secondDerivatives = self.cntrl[2::self.rpp]

        xmin = positionPoints[:,0].min()
        xmax = positionPoints[:,0].max()
        ymin = positionPoints[:,1].min()
        ymax = positionPoints[:,1].max()

        offset = max(abs(np.append(firstDerivatives,secondDerivatives)))
        
        pl.plot(self.points[:,0],self.points[:,1])

        ax = pl.gca()
        
        for pt,d1,d2 in zip(positionPoints,firstDerivatives,secondDerivatives):
            ax.add_patch(pl.Arrow(pt[0],pt[1],d1[0],d1[1],label='d1',color='r',width=0.2,alpha=0.3))
            ax.add_patch(pl.Arrow(pt[0],pt[1],d2[0],d2[1],label='d2',color='g',width=0.2,alpha=0.3))

        axis_limits = [xmin-offset,xmax+offset,ymin-offset,ymax+offset]

        pl.axis(axis_limits)
        return


class HermiteSpline(Spline):
    """
        Classe che modella una spline di Hermite
    """
    def __init__(self,cntrl,npts):
        """
            Costruttore
            @param cntrl punti di controllo indicanti punto posizione e vettore derivata prima nella forma:\
                         [(x0,y0),(d0x,d0y),(x1,y1),(d1x,d1y),...]
            @param npts numero di punti della curva da calcolare
        """

        Spline.__init__(self,cntrl,npts,2,4)
        self.__dict__['B'] = np.array([[2,1,-2,1],[-3,-2,3,-1],[0,1,0,0],[1,0,0,0]],np.double)

    def plot(self):
        """
            Metodo che effettua il plot della curva, mostrando le condizioni di controllo
            (derivata prima) come frecce
        """

        positionPoints = self.cntrl[::self.rpp]
        firstDerivatives = self.cntrl[1::self.rpp]

        xmin = positionPoints[:,0].min()
        xmax = positionPoints[:,0].max()
        ymin = positionPoints[:,1].min()
        ymax = positionPoints[:,1].max()

        offset = max(abs(np.ravel(firstDerivatives)))
        
        pl.plot(self.points[:,0],self.points[:,1])

        ax = pl.gca()
        
        for pt,d1 in zip(positionPoints,firstDerivatives):
            ax.add_patch(pl.Arrow(pt[0],pt[1],d1[0],d1[1],label='d1',color='r',width=0.2,alpha=0.3))

        axis_limits = [xmin-offset,xmax+offset,ymin-offset,ymax+offset]

        pl.axis(axis_limits)

	pl.show()

        return


class CardinalSpline(Spline):
    """
        Classe che modella una Cardinal spline
    """

    def __init__(self,cntrl,npts,t):
        """
            Costruttore
            @param cntrl vettore indicante la posizione dei punti di controllo
            @param npts numero di punti della curva da calcolare
            @param t parametro di \'tensione\' della curva
        """

        Spline.__init__(self,cntrl,npts,1,4)
        s = (1-t)/2
        self.__dict__['B'] = np.array([[-s,(2-s), (s-2), s],[2*s,(s-3),(3-2*s),-s],[-s,0,s,0],[0,1,0,0,]],np.double)


class BSpline(Spline):
    """
        Classe che modella una BSpline
    """
    def __init__(self,cntrl,npts,k,x):
        """
            Costruttore
            INPUT:
            @param cntrl vettore dei vertici del poligono di controllo
            @param npts numero di punti della curva da calcolare
            @param k ordine della bspline (grado polinomio +1)
            @param x può essere una stringa di valore "periodic" o "open": in questo caso
                     viene generato un vettore knot periodico o aperto, con molteplicità k agli estremi.
                     Altrimenti x può essere un vettore knot deciso dall'utente, passato sotto forma di
                     lista, tupla o numpy.ndarray

            ATTRIBUTI:
            @var n numero di punti di controllo
            @var B matrice di base nel caso si scelga il calcolo matriciale
        """

        Spline.__init__(self,cntrl,npts,1,k)
        self.__dict__['n'] = len(self.cntrl)
        self.__dict__['x'] = self.knot(x,self.k,self.n)
        self.__dict__['B'] = None
        

    def knot(self,x,k,n):
        """
            Metodo che genera un vettore knot
            il metodo genera un nuovo vettore se il parametro x è una stringa
            di valore 'periodic' o 'open'. Se invece x è una lista, una tupla
            o un istanza di numpy.ndarray che definisce un vettore knot,
            il metodo controlla che tale vettore sia valido, ordinandolo
            a priori e verificando se ci sono elementi diversi da 
            int o float all'interno di esso.
            @param x stringa indicante il tipo di knot vector desiderato, o un knot vector
            @param k l'ordine della bspline
            @param n il numero di punti di controllo
        """
        knotlen = k + n
        if x is 'periodic':
            return np.array(range(knotlen+1),np.double)
        elif x is 'open':
            left = [0]*k
            right = [n-k+1]*k
            center = range(1,right[0])
            return np.array(left+center+right)
        elif not isinstance(x,(list,tuple,np.ndarray)):
            raise TypeError("Il vettore dei knot deve essere passato sotto forma di lista, tupla o numpy.ndarray")
        else:
            x.sort()
            for element in x:
                if not isinstance(element,(int,float)):
                    raise TypeError("Il vettore knot deve essere composto di interi o float")
            return x
    
    def computeBasis(self,n,k,x,t):
        """
            Metodo che restituisce un vettore delle funzioni di base Ni,k(t)
            per ogni i da 0 a numero punti di controllo -1, calcolati tramite
            la ricorsione di Cox - De Boor
            INPUT:
            @param n numero di punti di controllo
            @param k ordine della bspline
            @param x vettore knot
            @param t valore del parametro

            Variabili:
            @var numKnots grandezza del vettore knot
            @var temp vettore temporaneo su cui verrà effettuata la ricorsione di Cox-De Boor
        """

        numKnots = len(x)
        
        # inizializzo temp con la base della ricorsione, Ni,1 = 1 se x[i] <= t < x[i+1], 0 altrimenti
        temp = np.array([1 if t>=x[i] and t < x[i+1] else 0 for i in xrange(numKnots - 1)],np.double)
    
        for k in xrange(2,k+1):
            for i in xrange(numKnots-k):
                # se Ni,k-1 è diverso da 0
                if temp[i]:
                    # calcolo il primo termine della ricorsione di Cox - De Boor
                    term1 = ( (t - x[i]) / (x[i+k-1] - x[i])) * temp[i]
                else:
                    # altrimenti lo pongo direttamente a 0
                    term1 = 0
                # se Ni+1,k-1 è diverso da 0
                if temp[i+1]:
                    # calcolo il secondo termine della ricorsione di Cox - De Boor
                    term2 = ((x[i+k] - t) / (x[i+k] - x[i+1])) * temp[i+1]
                else:
                    # altrimenti lo pongo direttamente a 0
                    term2 = 0
                
                # calcolo Ni,k sommando i due termini della ricorsione di Cox - De Boor
                temp[i] = term1 + term2

            # se il parametro coincide con l'ultimo valore knot
        if t == x[-1]:
            # recupero l'ultimo punto
            temp[n-1] = 1
            
        # restituisco i primi n elementi di temp
        return temp[:n]

    # metodo che calcola la matrice di base per una bspline periodica,
    # utilizzando la formula descritta a pag. 76 di Introduction to Nurbs
    # B[i,j] = 1/(k-1)! * comb(k-1,i) * sum per l=j a k di (k-(l+1))^i * -1^(l-j) * comb(k,l-j)
    # per ogni i,j appartenenti a [0,k-1]
    # il metodo inoltre, salva in un file denominato "order_k_periodic_BSpline_Base_Matrix.npy"
    # dove k è l'ordine della bspline, la matrice calcolata. Il primo tentativo che il metodo
    # effettua è di caricare da file la matrice se esiste, altrimenti la calcola e la salva,
    # creando, utilizzo dopo utilizzo, un archivio di matrici di base per BSpline periodiche.
    def computePeriodicMatrix(self,k):
        # imposto il nome del file
        fileName = 'order_'+str(k)+'_periodic_BSpline_Basis_Matrix.npy'
        
        # provo a leggere il file 
        try:
            B = np.load(fileName)

        except IOError as detail:
            # se il file non esiste, stampo un messaggio di avviso e la matrice verrà calcolata
            print "La matrice di ordine %d non è presente fra quelle salvate, la calcolo..."%(k)
        else:
            # se nessuna eccezione ha avuto luogo, restituisco la matrice letta da file
            #print "La matrice di ordine %d è presente fra quelle salvate"%(k)
            return B

        # inizializzo una matrice k x k di zeri
        B = np.zeros((k,k))

        # calcolo per primo il fattore moltiplicativo 1/(k-1)! che è fisso
        fact_coeff = 1/fact(k-1)
        
        # inizializzo gli iteratori di indici per righe e colonne
        rows = cols = range(k)
        
        for i in rows:
            # calcolo le combinazioni di k-1 su i
            temp = comb(k-1,i)
            # e le moltiplico per il coefficienti moltiplicativo 1/(k-1)!
            # ottenendo il valore della parte della formula che moltiplica la sommatoria
            a = fact_coeff * temp
            for j in cols:  
                # b è un accumulatore che conterrà il valore della sommatoria                
                b=0
                # il ciclo seguente calcola la sommatoria
                for l in range(j,self.k):
                    b += pow((self.k - (l+1)),i)*pow(-1,l-j)*comb(self.k,l-j)
                B[i,j] = a*b

        # salvo la matrice
        np.save(fileName,B)
        return B

    def calculateWithMatrixNotation(self):
        # con il seguente controllo verifico se la matrice self.B esiste e ne recupero il numero di righe e colonne
        # se non esiste pongo il numero di righe e di colonne pari a -1
        r,c = self.B.shape if isinstance(self.B,np.ndarray) else (-1,-1)
        # se le dimensioni della matrice non corrispondono all'attuale ordine della bspline o,
        # implicitamente, la matrice B ancora non esiste, calcolo la matrice di base
        if (r,c) != (self.k,self.k):
            self.__dict__['B'] = self.computePeriodicMatrix(self.k)
        # e invoco il metodo calculate della classe Spline
        Spline.calculate(self)
        return
        

    def calculate(self):

        self.points = np.zeros((self.npts, 2))
    
        #in caso contrario, determino i valori del parametro per cui calcolare la curva
        tvalues = np.linspace(self.x[self.k - 1],self.x[self.n],self.npts)
        
        for i,t in enumerate(tvalues):
            N = self.computeBasis(self.n, self.k, self.x, t)
            self.points[i] = np.dot(N, self.cntrl)

                
    # sovrascrivo l'operatore di assegnazione, in modo da proteggere gli attributi B e n
    # e da effettuare i dovuti controlli per i rimanenti
    def __setattr__(self, name, value):
        if name in ('B', 'n'):
            print "L'attributo {0} non può essere modificato dall'utente direttamente.".format(name)
        elif name is 'x':
            self.__dict__['x'] = self.knot(value, self.k, self.n)
            self.calculate()
        elif name is 'k':
            if not isinstance(value,int):
                raise TypeError("Il grado k deve essere di tipo int,non {0}".format(type(value)))
            elif value > len(self.cntrl):
                raise ValueError("Il grado k può essere al più uguale al numero di control points, che è {0}".format(len(self.cntrl)))
            else:
                self.__dict__['k'] = value
                self.calculate()
        else:
            Spline.__setattr__(self,name,value)
        if name is 'cntrl':
            self.__dict__['n'] = len(self.cntrl)
            self.__dict__['x'] = self.knot('periodic',self.k,self.n)


class Nurbs(BSpline):
    """
        Classe che modella una curva NURBS
    """
    def __init__(self,cntrl,npts,k,x,h):
        """
            Costruttore
            @param h lista, tupla o numpy.ndarray contenente i pesi dei punti di controllo
        """

        BSpline.__init__(self,cntrl,npts,k,x)
        
        # se h non è né una lista, né una tupla, né un numpy.ndarray
        if not isinstance(h,(list,tuple,np.ndarray)):
            raise TypeError("h deve essere una lista, una tupla o un numpy.ndarray, non di tipo {0}\n".format(type(h)))
        elif len(h) != len(self.cntrl):
            raise ValueError("h ha lunghezza {0} mentre cntrl ha lunghezza {1}".format(len(h),len(self.cntrl)))
        self.__dict__['h'] = h

    def computeBasis(self,n,k,x,h,t):
        
        # calcolo le funzioni di base non razionali
        temp = BSpline.computeBasis(self,n,k,x,t)
        # moltiplico le funzioni di base per i pesi
        temp = np.asarray([h[i]*basis for i,basis in enumerate(temp)])
        # divido tutto per la somma degli hi*Ni,k(t), ottenendo la base razionale
        temp = temp/sum(temp)
        return temp

    def calculate(self):
        tvalues = np.linspace(self.x[self.k - 1],self.x[self.n],self.npts)
        
        for i,t in enumerate(tvalues):
            N = self.computeBasis(self.n,self.k,self.x,self.h,t)
            self.points[i] = np.dot(N,self.cntrl)

    def __setattr__(self,name,value):
        if name is 'h':
            if not isinstance(value,(list,tuple,np.ndarray)):
                raise TypeError("h deve essere una lista, una tupla o un istanza di numpy.ndarray, non di tipo {0}".format(type(h)))
            elif len(value) != len(self.cntrl):
                raise ValueError("h ha lunghezza {0} mentre cntrl ha lunghezza {1}".format(len(value),len(self.cntrl)))
            self.__dict__[name] = value
            self.calculate()
        else:
            BSpline.__setattr__(self,name,value)

    def __add__(self,c):
        new_curve = Curve.__add__(self,c)
        new_curve.h = np.append(self.h,c.h[1:])
        return new_curve
