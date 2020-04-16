#!/usr/bin/python3
# -*-coding:Utf-8 -*
import numpy as np
from numpy import sqrt,arctan,pi,cos,sin,exp
from scipy.optimize import brentq
from matplotlib import pyplot as plt
np.set_printoptions(precision=4)  # For compact display.

class Guide:
    """ Classe définissant un guide d'onde achiral:
        Usage:
      import Achiral as ga
      g=ga.Guide(ns=1.46,ng=1.61,d=3,n0=1,wl=0.633) # default 
      g=ga.guide(ns=1.46,ng=1.697,d=0.25,n0=1,wl=0.633,anisotropie=1.042) #heli3
      g=Guide()
      g  # affiche la l objet
      help(g) #cette aide
      g.GrapheDisp(dmin=0.,dmax=5.,Npts=20) # dispersion 
      g.tabdisp #le tableau de la disp
      #sauvegarde de ce tableau
      import numpy as np
      np.savetxt('test.out',g.tabdisp,fmt='%1.4e')
      #  x et mode 0..2
      np.savetxt('test.out',g.tabdisp[:,0:4],fmt='%1.4e')
      g.GrapheTE(0) #graphique du mode 0    
      g.AfficheProfilTE() #Affecte et affiche les profils de champ
                      # de tous les modes TE pour gnuplot
        """
    def __init__(self,ns=1.46,ng=1.61,d=3,n0=1,wl=0.633,anisotropie=1.):
        """Constructeur 
        calcule et affecte tous les tableaux sauf champ
            """
        self.ns=ns
        self.ng=ng
        self.an=anisotropie*anisotropie
        self.d=d
        self.n0=n0
        self.rs=ns*ns/ng/ng
        self.r0=n0*n0/ng/ng
        self.wl=wl
        self.k0=2*pi/wl
        #coupures:
        self.dcTE=[] 
        self.dcTM=[]
        #neff
        self.NeffTE=[]
        self.NeffTM=[]
        # parametres de continuite pour les champ
        self.c1=[]
        self.c2=[]
        # listes des profils de champ, en gnuplot
        self.champs=[]
        # courbe de dispertion
        self.tabdisp=[]
        #nombre de mode
        self.NTE=0
        self.NTM=0
        self.CalculeLesNeff()
    def __repr__(self):
        """ Affichage cool """
        reponse="-"*20+"\n"
        interf="\n"+reponse
        reponse+=" n0="+str(self.n0)
        reponse+=interf
        reponse+="ng="+str(self.ng)+" d="+str(self.d)
        reponse+=interf
        reponse+="ns="+str(self.ns)
        reponse+=interf
        reponse+=str(self.wl)+"nm ,"
        reponse+=str(self.NTE)+" TE et "+str(self.NTM)+" TM\n"
        reponse+="Coupure(TE)="+str(self.dcTE)+"\n"
        reponse+="Coupure(TM)="+str(self.dcTM)+"\n"
        reponse+="Neff(TE)="+str(self.NeffTE)+"\n"
        reponse+="Neff(TM)="+str(self.NeffTM)+"\n"
        diff=np.asarray(self.NeffTE)-np.asarray(self.NeffTM)
        reponse+="nTE-nTM="+str(diff)+"\n"
        reponse+="champs:"+str(self.champs)+"\n"
        return reponse

    def CalculeLesNeff(self):
        """ Determine et affecte dans l'ordre:
          dcTE/TM et  NTE/TM  puis NeffTE/TM 
         on calcule d'abord le nombre de modes 
         à partir des coupures puis on résoud 
         avec Brentq l'équ. de dispertion """
        na=self.ns
        nb=self.ng
        self.NeffTE=[]
        self.NeffTM=[]
        self.dcTE=[] 
        self.dcTM=[]
        self.NbreModesTE()
        # avec range(0) ça ne fait rien 
        # -> ça s arrete avant "0"
        for i in range(self.NTE):
            sol=brentq(self.ArcTE,na,nb,args=(i))
            self.NeffTE.append(sol)                              
        self.NbreModesTM()
        for i in range(self.NTM):
            sol=brentq(self.ArcTM,na,nb*np.sqrt(self.an),args=(i))
            self.NeffTM.append(sol) 
        if (self.NTE>0):
            self.Diametre()
            
    def CalculeDisp(self,dmin=0.,dmax=5.,Npts=50):
        """ Calcul la courbe de dispersiion
        entre dmin et dmax sur Npts.
        le resultat est affecte au tableau 2D tabdisp
        constitué de la colonne "d" puis des neff """
        # nombre max de modes:
        self.d=dmax
        self.NbreModesTE()
        Nmax=self.NTE
        self.NbreModesTM()
        if self.NTM> Nmax:
            Nmax=self.NTM
        #un premier passage pour la ligne 1
        self.d=dmin
        self.CalculeLesNeff()
        completion=self.ns*np.ones((Nmax-len(self.NeffTE)))
        completionTM=self.ns*np.ones((Nmax-len(self.NeffTM)))
        self.tabdisp=np.hstack((self.d,self.NeffTE,completion,self.NeffTM,completionTM))
        for di in np.linspace(dmin,dmax,num=Npts):
            self.d=di
            self.CalculeLesNeff()
            completion=self.ns*np.ones((Nmax-len(self.NeffTE)))
            completionTM=self.ns*np.ones((Nmax-len(self.NeffTM)))
            ligne=np.hstack((self.d,self.NeffTE,completion,self.NeffTM,completionTM))
            self.tabdisp=np.vstack((self.tabdisp,ligne))

    def CalculeLambdaDisp(self,lbmin=0.4,lbmax=0.8,Npts=50):
        """ Calcul la courbe de dispersiion
        entre lbmin et lbmax sur Npts.
        le resultat est affecte au tableau 2D tabdisp
        constitué de la colonne "lambda" puis des neff """
        # nombre max de modes:
        self.wl=lbmin
        self.k0=2*pi/self.wl
        self.CalculeLesNeff()
        self.NbreModesTE()
        Nmax=self.NTE
        #un premier passage pour la ligne 1
        #completion=self.ns*np.ones((Nmax-len(self.NeffTE)))
        #completionTM=self.ns*np.ones((Nmax-len(self.NeffTM)))
        self.tabdisp=np.hstack((self.wl,self.NeffTE,self.NeffTM))
        for i in range(Npts-1):
            self.wl=lbmin+(i+1)*lbmax/Npts
            self.k0=2*pi/self.wl
            self.CalculeLesNeff()
            completion=self.ns*np.ones((Nmax-len(self.NeffTE)))
            completionTM=self.ns*np.ones((Nmax-len(self.NeffTM)))
            ligne=np.hstack((self.wl,self.NeffTE,completion,self.NeffTM,completionTM))
            self.tabdisp=np.vstack((self.tabdisp,ligne))

    def ArcTE(self,neff,Mode):
        """ Relation de dispertion TE """
        v=sqrt(neff*neff-self.n0*self.n0)
        u=sqrt(self.ng*self.ng-neff*neff)
        w=sqrt(neff*neff-self.ns*self.ns)
        wd=self.k0*u*self.d
        if u>0:
            retour=wd-arctan(v/u)-arctan(w/u)-Mode*pi
        else:
            retour=wd-3.141592653589-Mode*pi
        return retour

    def CoupureTE(self,Mode):
        """ coupure du mode  en TE"""
        v=sqrt(self.ns*self.ns-self.n0*self.n0)
        u=sqrt(self.ng*self.ng-self.ns*self.ns)
        return (arctan(v/u)+Mode*pi)/u/self.k0

    def ArcTM(self,neff,Mode):
        """ Relation de dispertion TM """
        v=sqrt(neff*neff-self.n0*self.n0)/self.r0
        u=sqrt(self.ng*self.ng*self.an-neff*neff)
        w=sqrt(neff*neff-self.ns*self.ns)/self.rs
        wd=self.k0*u*self.d
        if u>0:
            retour=wd-arctan(v/u)-arctan(w/u)-Mode*pi
        else:
            retour=wd-3.141592653589-Mode*pi
        return retour
        
    def CoupureTM(self,Mode):
        """ coupure du mode  TM """
        v=sqrt(self.ns*self.ns-self.n0*self.n0)/self.r0
        u=sqrt(self.ng*self.ng*self.an-self.ns*self.ns)
        return (arctan(v/u)+Mode*pi)/self.k0/u
    

    def NbreModesTE(self):
        Nmodes=0
        while self.CoupureTE(Nmodes)<self.d:
            self.dcTE.append(self.CoupureTE(Nmodes))
            Nmodes+=1
        self.NTE=Nmodes

    def NbreModesTM(self):
        """ Calcul et affecte coupureTM NTM """
        Nmodes=0
        while self.CoupureTM(Nmodes)<self.d:
            self.dcTM.append(self.CoupureTM(Nmodes))
            Nmodes+=1
        self.NTM=Nmodes

    def AfficheProfilTE(self):
        """ Function profil TE pour gnuplot 
        mises dans champ sous la forme
        champ[0]="E0(x)=...."
        champ[1]="E1(x)=..." """
        self.champs=[]
        ne=np.array(self.NeffTE)
        u=sqrt(self.ng**2-ne**2)*self.k0
        w=sqrt(ne**2-self.ns**2)*self.k0
        v=sqrt(ne**2-self.n0**2)*self.k0
        self.c1=-v/u
        self.c2=cos(u*self.d)-self.c1*sin(u*self.d)
        for i in range(self.NTE):
            E0="(x>0)*exp(-{0}*x)".format(v[i])
            E1="(x<0)*(x>-{0})*(cos({1}*x)+{2}*sin({1}*x))".\
                    format(self.d,u[i],self.c1[i])
            E2="(x<-{1})*(cos({0}*{1})-{2}*sin({0}*{1}))*exp({3}*(x+{1}))".\
                    format(u[i],self.d,self.c1[i],w[i])
            Ex="E{0}(x)=".format(i)+E0+"+"+E1+"+"+E2
            print(Ex)
            self.champs.append(Ex)

    def Diametre(self):
        """ Calcule le waist du champ du mode 0 """
        Npts=2000
        ne=self.NeffTE[0]
        v=sqrt(ne**2-self.n0**2)*self.k0
        Max=3./v
        w=sqrt(ne**2-self.ns**2)*self.k0
        Min=-3./w-self.d
        x=np.arange(Min,Max,(Max-Min)/Npts)
        y=np.piecewise(x,[x>0,x<-self.d],\
                   [self.EvanAir,self.EvanSub,self.Ecoeur],0)
        Ncentre=np.argmax(y)
        ymax=y[Ncentre]
       # print("Max:{}, Min:{},ymax*e-1:{}".format(Max,Min,ymax*exp(-1)))
       # print(" y1:{},y-1:{}".format(y[1],y[-1]))
        if y[1]<ymax*exp(-1):
            N0=next(data[0] for data in enumerate(y) if data[1]>ymax*exp(-1))
        else:
            N0=0
        self.xm1=x[N0]
        self.xmin=x[Ncentre]-1.5*(x[Ncentre]-x[N0])
        if y[-1]<ymax*exp(-1):
            Nplus=N0+next(data[0] for data in enumerate(y[N0:]) if
                data[1]<ymax*exp(-1))
        else:
            Nplus=-1
        self.xp1=x[Nplus]
        self.xmax=x[Ncentre]+1.5*(x[Nplus]-x[Ncentre])
        self.diam=self.xp1-self.xm1
      #  print("Diametre: {}".format(self.diam))

        

    def GrapheDisp(self,dmin=0.,dmax=5.,Npts=50):
        """ Graphique dispersion 
        ! rentrer les nombres avec les "." """
        self.CalculeDisp(dmin=dmin,dmax=dmax,Npts=Npts)
        x=self.tabdisp[:,0]
        plt.figure(1)
        for i in range(2*self.NTE):
            plt.plot(x, self.tabdisp[:,i+1])
        plt.xlabel(r'd($\mu$m)')
        plt.ylabel(r'n$_{eff}$')
        plt.grid(True)
        plt.title("n0={0},ng={1},ns={2},{3}nm"\
                  .format(self.n0,self.ng,self.ns,self.wl))
        plt.show()
        return 1 

    def GrapheTE(self,Mode=0):
        """ Graphique du  l'intesnite du mode TEi """
        Npts=100
        self.Diametre()
        gr=1.5
        x=np.arange(self.xmin*gr,abs(self.xmax)*gr,self.diam*gr/Npts)
        y=np.piecewise(x,[x>0,x<-self.d],\
                   [self.IvanAir,self.IvanSub,self.Icoeur],Mode)

        plt.figure(1)
        plt.ioff()
        plt.plot(x, y)
        plt.xlabel(r'x($\mu$m)')
        plt.ylabel('E')
        plt.grid(True)
        plt.title("n0={0},ng={1},ns={2},d={3}"\
                  .format(self.n0,self.ng,self.ns,self.d))
        plt.show()
        return 1 

    def EvanAir(self,x,Mode):
        """ Profil TE du mode i dans le superstrat"""
        ne=self.NeffTE[Mode]
        v=sqrt(ne**2-self.n0**2)*self.k0
        return exp(-v*x)

    def Ecoeur(self,x,Mode):
        """ Profil TE du mode i dans le coeur """
        ne=self.NeffTE[Mode]
        v=sqrt(ne**2-self.n0**2)*self.k0
        u=sqrt(self.ng**2-ne**2)*self.k0
        return  cos(u*x)-v/u*sin(u*x)

    def EvanSub(self,x,Mode):
        """ Profil TE du mode i dans le substrat"""
        ne=self.NeffTE[Mode]
        v=sqrt(ne**2-self.n0**2)*self.k0
        u=sqrt(self.ng**2-ne**2)*self.k0
        w=sqrt(ne**2-self.ns**2)*self.k0
        E=(cos(u*self.d)+v/u*sin(u*self.d))*exp(w*(x+self.d))
        return E
    
    def IvanAir(self,x,Mode):
        """ Profil TE du mode i dans le superstrat"""
        ne=self.NeffTE[Mode]
        v=sqrt(ne**2-self.n0**2)*self.k0
        return exp(-2*v*x)

    def Icoeur(self,x,Mode):
        """ Profil TE du mode i dans le coeur """
        ne=self.NeffTE[Mode]
        v=sqrt(ne**2-self.n0**2)*self.k0
        u=sqrt(self.ng**2-ne**2)*self.k0
        E=cos(u*x)-v/u*sin(u*x)
        return  E*E

    def IvanSub(self,x,Mode):
        """ Profil TE du mode i dans le substrat"""
        ne=self.NeffTE[Mode]
        v=sqrt(ne**2-self.n0**2)*self.k0
        u=sqrt(self.ng**2-ne**2)*self.k0
        w=sqrt(ne**2-self.ns**2)*self.k0
        E=(cos(u*self.d)+v/u*sin(u*self.d))*exp(w*(x+self.d))
        return E*E
#sg=Guide(1.62,1.68,2.8)
##print(sg.__doc__)
#sg.CalculeLesNeff()
#sg.AfficheProfilTE()
#sg.GrapheTE(0)

