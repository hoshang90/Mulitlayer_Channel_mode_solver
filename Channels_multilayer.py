#!/usr/bin/python3
"""Fully vectorial finite-difference mode solver example."""
#basé sur EMpy: https://github.com/lbolla/EMpy
#    grandeurs en microns
#    grille de pas: pas ds les deux directions
#    hauteur grille (y): Hc+H1+H2+Hsub
#    largeur grille (x): LB+LC+LB
# utilisation avec ipython
#---------------------------------------------------------------
# generation du profile: ipython3 profile create ridge
# puis editer .ipython/profile_ridge/ipython_config.py

import numpy as np
import EMpy
import matplotlib.pyplot as plt
import time

import Achiral
class Channel():
    """ modes propres dans guide canal  
    basé sur EMpy: https://github.com/lbolla/EMpy
    grandeurs en microns
    grille de pas: pas ds les deux directions
    hauteur grille (y): Hc+H1+H2+Hsub
    largeur grille (x): LB+LC+LB_R
    ___________________________________
        nc    |  nc_2    | nc_3, Hc
    ___________________________________
    <-- LB --> <-- LC --> <-- LB_R -->
        n1B   |   n1C    | n1D , H1
    --------------------------------------
    <-- LB --> <-- LC --> <-- LB_R-->
        n2B   |   n2C    | n2D  , H2
    -------------------------------------
      nsub    |  nsub_2  | nsub_3,Hsub
    _______________________________________
    """
    def __init__(self,wl=0.64,pas=0.16,nc=1., nc_2=1., nc_3=1.,Hc=0.8,\
            LB=5, LC=3.1,LB_R=5,n1B=1.,n1C=1.620,n1D=1.,H1=2.3,\
            n2B=1.62,n2C=1.62,n2D=1.62,H2=2.,\
            nsub=1.61,nsub_2=1.61,nsub_3=1.61,Hsub=2.,OR=4):
        self.Hc=Hc
        self.nc2=nc**2;self.nc2_2 = nc_2**2 ;self.nc2_3 = nc_3**2
        self.LB=LB;self.LB_R=LB_R;self.LC=LC
        self.n1C2=n1C**2; self.n1B2=n1B**2;self.n1D2=n1D**2
        self.H1=H1
        self.n2C2=n2C**2; self.n2B2=n2B**2;self.n2D2=n2D**2
        self.H2=H2
        self.ns2=nsub**2;self.ns2_2=nsub_2**2;self.ns2_3=nsub_3**2;self.Hsub=Hsub
        self.wl=wl
        self.pas=pas
        self.OR=OR

        self.initGrille()
        print("Grille {}x{} micron2".format(self.LargGrille,self.HautGrille))
        print("pas={}, Npts={}".format(self.pas,self.x.size*self.y.size))
  
    def initGrille(self):
        ''' Initizalisation de la grille de calcul '''
        pas=self.pas
        self.HautGrille=self.Hc+self.H1+self.H2+self.Hsub
        self.LargGrille=self.LC+self.LB+self.LB_R
        self.x = np.arange(0, self.LargGrille+pas, pas)
        self.y = np.arange(0, self.HautGrille+pas, pas)

    def __repr__(self):
        reponse=self.LaStructure()
        Nmodes=len(self.sol.modes)
        neffplanMax=self.NmaxPlan()
        for i in range(Nmodes):
            neff=self.sol.modes[i].neff
            reponse+="*Mode[{}]: neff={:.3f}, ".format(i,neff)
            reponse+="%TE= {:.2f}, ".format(100*self.sol.modes[i].TEfrac())
            reponse+="%Conf= {:.2f}".format(100*self.Confinement(Nmode=i))
            if neff<neffplanMax:
                reponse+=" Leak "
            reponse+="\n"
        DN=self.sol.modes[0].neff-self.sol.modes[1].neff
        reponse+="n0-n1={:g}, ".format(DN)
        CB=self.OR*1e-3/180*self.wl
        ecc=CB/(DN+np.sqrt(CB*CB+DN*DN))
        reponse+="ecc={:.3f}\n".format(ecc)
        reponse+="S3={:.3f},{:.3f}\n".format(np.sin(2*np.arctan((np.real(ecc))**-1)),np.sin(2*np.arctan(np.real(ecc))))
        return reponse

    def NumOf_guidedModes(self,Nmodes=2): # to get a list of (all guided modes, not leaky modes).
        neigs = Nmodes
        #self.initGrille()
        tol = 1e-9
        pas =self.pas
        boundary = '0000'
        self.sol = EMpy.modesolvers.FD.VFDModeSolver(self.wl, self.x, self.y, \
                                                     self.epsfunc, boundary).solve(neigs, tol)
        Nmodes_len = len(self.sol.modes)
        neffplanMax = self.NmaxPlan()
        NumOf_guidedModes_list = [neffplanMax]
        for i in range(Nmodes_len):
            neff = self.sol.modes[i].neff
            if neff>=neffplanMax:
                NumOf_guidedModes_list.append(neff)
        return np.asarray(NumOf_guidedModes_list)

    def LaStructure(self):
        """ Dessine la structure sur le terminal """
        reponse = " \n"
        reponse += "." * 20 + " All units are in (" + '\u03BCm' + ") " + "." * 19 + "\n"
        reponse += "|<--LB={} -->| <----  LC={} ---->| <------ LB_R={} ------>|\n".format(self.LB, self.LC, self.LB_R)
        reponse += "|<--nc={:.3f},>| <---- nc_2={:.3f}-->| <--nc_3={:.3f}, Hc={}-->| \n". \
            format(np.sqrt(self.nc2_3), np.sqrt(self.nc2_2), np.sqrt(self.nc2), self.Hc)
        reponse += "|" + "-" * 60 + "|" + "\n"
        # reponse+=" "*12+"-"*14+"\n"
        reponse += "|<-n1B={:.3f},>| <----- n1C={:.3f}-->| <---n1D={:.3f}, H1={}-->| \n". \
            format(np.sqrt(self.n1D2), np.sqrt(self.n1C2), np.sqrt(self.n1B2), self.H1)
        reponse += "|" + "-" * 60 + "|" + "\n"
        reponse += "|<-n2B={:.3f},>| <----- n2C={:.3f}-->| <---n2D={:.3f}, H2={}-->| \n". \
            format(np.sqrt(self.n2D2), np.sqrt(self.n2C2), np.sqrt(self.n2B2), self.H2)
        reponse += "|" + "-" * 60 + "|" + "\n"
        reponse += "|<-nsub={:.3f},>| <--nsub_2={:.3f}-->| <nsub_3={:.3f},Hsub={}->| \n". \
            format(np.sqrt(self.ns2_3), np.sqrt(self.ns2_2), np.sqrt(self.ns2), self.Hsub)
        reponse += "." * 62 + "\n"
        reponse += " \n"
        return reponse

    def Calcule(self,Nmodes=2,verbose=True,trace=True):
        neigs=Nmodes
        self.initGrille()
        tol=1e-9
        boundary='0000'
        self.sol=EMpy.modesolvers.FD.VFDModeSolver(self.wl, self.x, self.y,\
                self.epsfunc, boundary).solve(neigs, tol)
        if verbose:
            print(self.__repr__())
        if trace:
            self.Intensite()

    def epsfunc(self,x_, y_):
        """Return a matrix describing a 2d material.
        :param x_: x values
        :param y_: y values
        :return: 2d-matrix
        """
        xx, yy = np.meshgrid(x_,y_)
        # to fill Hc x LB
        cg = np.where((yy.T >= self.H1 + self.H2 + self.Hsub)
                      * (np.abs(xx.T - self.LB-self.LC - self.LB_R) < self.LB), self.nc2,self.ns2_3 )
        #we can also write as: cg = np.where((np.abs(yy.T -self.Hc - self.H1 - self.H2 - self.Hsub) <= self.Hc)
                      #* (np.abs(xx.T ) <= self.LC+self.LB_R), self.nc2_2, cg)
        # to fill Hc x LC
        cg = np.where((yy.T >= self.H1 + self.H2 + self.Hsub)
                      * (np.abs(xx.T ) <= self.LC+self.LB_R), self.nc2_2, cg)
        # to fill Hc x LB_2
        cg = np.where((yy.T >= self.H1 + self.H2 + self.Hsub)
                      * (np.abs(xx.T ) <= self.LB_R), self.nc2_3, cg)
        # to fill H1 x LB
        cg = np.where((np.abs(yy.T - 0.5 * self.H1 - self.H2 - self.Hsub) <= 0.5 * self.H1)
                      * (np.abs(xx.T - self.LB-self.LC - self.LB_R) < self.LB), self.n1B2, cg)
        # to fill H1 x LC
        cg = np.where((np.abs(yy.T - 0.5 * self.H1 - self.H2 - self.Hsub) <= 0.5 * self.H1)
                      * (np.abs(xx.T ) <= self.LC+self.LB_R), self.n1C2, cg)
        # to fill H1 x LB_2
        cg = np.where((np.abs(yy.T - 0.5 * self.H1 - self.H2 - self.Hsub) <= 0.5 * self.H1)
                      * (np.abs(xx.T ) <= self.LB_R), self.n1D2, cg)
        # to fill H2 x LB
        cg = np.where((np.abs(yy.T - 0.5 * self.H2 - self.Hsub) <= 0.5 * self.H2)
                      * (np.abs(xx.T - self.LB - self.LC - self.LB_R) < self.LB), self.n2B2, cg)
        # to fill H2 x LC
        cg = np.where((np.abs(yy.T - 0.5 * self.H2 - self.Hsub) <= 0.5 * self.H2)
                      * (np.abs(xx.T ) <= self.LC+self.LB_R), self.n2C2, cg)
        # to fill H2 x LB_2
        cg = np.where((np.abs(yy.T - 0.5 * self.H2 - self.Hsub) <= 0.5 * self.H2)
                      * (np.abs(xx.T ) <= self.LB_R), self.n2D2, cg)
        # to fill Hsub x LB
        cg = np.where((yy.T <= self.Hsub)
                      * (np.abs(xx.T - self.LB - self.LC - self.LB_R) < self.LB), self.ns2, cg)
        # to fill Hsub x LC
        cg = np.where((yy.T <= self.Hsub)
                      * (np.abs(xx.T ) <= self.LC+self.LB_R), self.ns2_2, cg)
        # to fill Hsub x LB_2
        cg = np.where((yy.T <= self.Hsub)
                      * (np.abs(xx.T ) <= self.LB_R), self.ns2_3, cg)
        return cg

    def Traceindice(self):
        ''' plot l'indice '''
        indices=self.epsfunc(self.x,self.y)
        levels=np.sqrt([1,self.nc2, self.nc2_2, self.nc2_3,self.n1B2,self.n1C2, self.n1D2,self.n2B2,self.n2C2, self.n2D2,self.ns2, self.ns2_2, self.ns2_3])
        fig=plt.figure()
        plt.contour(self.x,self.y,np.sqrt(indices.T),np.unique(levels))
        plt.contourf(self.x,self.y,np.sqrt(indices.T),np.unique(levels))
        plt.colorbar()
        plt.show()

    def NmaxPlan(self):
        ''' indice TE0 entre n2B '''
        nsub=np.sqrt(self.ns2)
        ng=np.sqrt(self.n2B2)
        n0=np.sqrt(self.n1B2)
        plan=Achiral.Guide(ns=nsub,ng=ng,d=self.H2,n0=n0,wl=self.wl)
        neffplanMax=0
        if plan.NTE>0:
            neffplanMax=plan.NeffTE[0]
        return neffplanMax

    def Intensite(self,Nmode=0):
        ''' plot l'intensité du mode  '''
        #le champ E a un pont de moins, on recupere les bons x
        x0=EMpy.utils.centered1d(self.x)
        y0= EMpy.utils.centered1d(self.y)
        Itot=np.abs(self.sol.modes[Nmode].intensity())
        fig = plt.figure()
        plt.hot()
        levels=np.unique([1,self.nc2, self.nc2_2, self.nc2_3,self.n1B2,self.n1C2, self.n1D2,self.n2B2,self.n2C2, self.n2D2,self.ns2, self.ns2_2, self.ns2_3])
        Nlevel=15
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(x0,y0 ,Itot.T, Nlevel)
        print("%Conf: {:.2f}".format(100*self.Confinement(Nmode=Nmode)))
        print("%TE: {:.1f}".format(100*self.sol.modes[Nmode].TEfrac()))
        #return Itot

    def Confinement(self,Nmode=0):
        ''' retourne  I(cover)/Itotale '''
        #le champ E a un pont de moins, on recupere les bons x
        x0=EMpy.utils.centered1d(self.x)
        y0= EMpy.utils.centered1d(self.y)
        Itot=np.abs(self.sol.modes[Nmode].intensity())
        epsi=self.epsfunc(x0,y0)
        masque=np.where(epsi<self.ns2,True,False)
        return np.sum(Itot*masque)/np.sum(Itot)


    def TraceMode(self,Nmode=0):
        ''' trace les 6 composantes du champ '''
        Nlevel=15
        levels=np.unique([1,self.nc2, self.nc2_2, self.nc2_3,self.n1B2,self.n1C2, self.n1D2,self.n2B2,self.n2C2, self.n2D2,self.ns2, self.ns2_2, self.ns2_3])
        #le champ E a un pont de moins, on recupere les bons x
        x0=EMpy.utils.centered1d(self.x)
        y0= EMpy.utils.centered1d(self.y)

        fig = plt.figure(figsize=(12, 6))
        plt.hot()
        fig.add_subplot(2, 3, 1)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(x0,y0,abs(self.sol.modes[Nmode].Ex).T, Nlevel)
        plt.title('Ex')
        plt.colorbar()

        fig.add_subplot(2, 3, 2)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(x0,y0,abs(self.sol.modes[Nmode].Ey).T, Nlevel)
        plt.title('Ey')
        plt.colorbar()

        fig.add_subplot(2, 3, 3)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(x0,y0,abs(self.sol.modes[Nmode].Ez).T, Nlevel)
        plt.title('Ez')
        plt.colorbar()

        fig.add_subplot(2, 3, 4)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(self.x,self.y,abs(self.sol.modes[Nmode].Hx).T, Nlevel)
        plt.title('Hx')
        plt.colorbar()

        fig.add_subplot(2, 3, 5)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(self.x,self.y,abs(self.sol.modes[Nmode].Hy).T, Nlevel)
        plt.title('Hy')
        plt.colorbar()

        fig.add_subplot(2, 3, 6)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(self.x,self.y,abs(self.sol.modes[Nmode].Hz).T, Nlevel)
        plt.title('Hz')
        plt.colorbar()

        plt.show()

    def TraceEx(self,Nmode=0):
        ''' Trace le profil de Ex au milieu du guide'''
        milieu=int(self.sol.modes[Nmode].Ex.shape[0]/2.)
        Ex=self.sol.modes[Nmode].Ex[milieu,:]
        scale=np.sqrt(self.n1C2)/np.amax(np.abs(Ex))
        indices=np.sqrt(self.epsfunc(self.x,self.y))[milieu,:]
        ymax=self.Hc
        ymin=-(self.H1+self.H2+self.Hsub)
        fig=plt.figure()
        plt.plot(np.linspace(ymin,ymax,indices.size),indices,np.linspace(ymin,ymax,Ex.size),np.abs(Ex)*scale)
        plt.grid()
        plt.xlabel(r'$\mu m$')
        plt.show()

    def VariaX(self,X0=0.1, X1=5, dX=0.2,variable='LC',trace=True,fich="DataDn"):
        ''' calcule Dn et Sc en bouclant sur la variable consideree: 
            retourne les tableaux 1D Dn et Sc
            avec trace, trace la courbe
            pour sauver: np.savetxt("file",tableau,fmt="%.2g") '''
        lesX=np.arange(X0,X1,dX)
        lesDn=np.zeros((lesX.size))
        lesSc=np.zeros((lesX.size))
        for i,X in enumerate(lesX):
            self.FixeVariable(X,variable)
            self.Calcule(verbose=False)
            dn=self.sol.modes[0].neff-self.sol.modes[1].neff
            Sc=self.Confinement()
            lesDn[i]=dn
            lesSc[i]=Sc
            print(variable,"={:.2f} => 10^4 dn={:.2g},%Sc={:.3f}".format(X,1e4*dn,100*Sc))
        if trace:
            fig = plt.figure(figsize=(12, 6))
            fig.add_subplot(2,1,1)
            plt.plot(lesX,1e4*np.abs(lesDn))
            plt.ylabel("10^4*Dn")
            plt.grid(True)
            plt.xlabel(variable)
            fig.add_subplot(2,1,2)
            plt.plot(lesX,100*lesSc)
            plt.ylabel("100*Sc")
            plt.grid(True)
            plt.xlabel(variable)
            plt.show()
        return (lesDn,lesSc)

    def VariaXY(self,X0=0.1, X1=5,
            dX=0.2,vX='LC',Y0=0.1,Y1=5,dY=0.2,vY='H1',fich="dataDn"):
        """ idem VariaX mais deux boucles imbriquees """
        lesX=np.arange(X0,X1,dX)
        lesY=np.arange(Y0,Y1,dY)
        self.lesDn=np.zeros((lesX.size,lesY.size))
        self.lesSc=np.zeros((lesX.size,lesY.size))
        self.lesEcc=np.zeros((lesX.size,lesY.size))
        CB=self.OR*1e-3/180*self.wl
        for i,X in enumerate(lesX):
            self.FixeVariable(X,vX)
            t0= time.process_time()
            for j,Y in enumerate(lesY):
                self.FixeVariable(Y,vY)
                self.Calcule(verbose=False,trace=False)
                #on ne garde que les solution sans leaking
                if self.sol.modes[0].neff>self.NmaxPlan() :
                    dn=abs(self.sol.modes[0].neff-self.sol.modes[1].neff)
                else :
                    dn=-1
               # Sc=self.Confinement()
                self.lesDn[i,j]=dn
            te=time.process_time() - t0
            print("Boucle {:.2f}, reste \
                        {:.2f}min".format(te,te*(lesX.size-i)/60))
        print("---------------")
        print('Fichiers enregistres: self.lesDn.T: '+fich)
        print("f(x)=x<0?1/0:DnCB/(x+sqrt(DnCB*DnCB+x*x))")
        print("plot DnCB=7e-5,\""+fich+"\" matrix\
                u ($1*{}+{}):($2*{}+{}):(f($3)) w image".format(dX,X0,dY,Y0))
        entete=self.LaStructure()+"f(x)=cb/(x+sqrt(cb*cb+x*x))\n"
        entete+="plot cb=4e-3/180*0.64,\"+fich+\" matrix u ($1*{}+{}):($2*{}+{}):(f($3)) w image".format(dX,X0,dY,Y0)
        entete += "Ploting variables are X0={},X1={},dX={},Y0={},Y1={},dY={},".format(X0,X1,dX,Y0,Y1,dY)
        np.savetxt(fich,self.lesDn.T,fmt='%.2g',header=entete, encoding='utf-8')
#        entete=self.LaStructure()+"plot \"dataSc\" matrix u ($1*{}+{}):($2*{}+{}):($3) w image".format(dX,X0,dY,Y0)
#        np.savetxt("dataSc",self.lesSc.T,fmt='%.2g',header=entete)
#        entete=self.LaStructure()+"plot \"dataEcc\" matrix u ($1*{}+{}):($2*{}+{}):($3) w image".format(dX,X0,dY,Y0)
#        np.savetxt("dataEcc",self.lesEcc.T,fmt='%.2g',header=entete)

    def FixeVariable(self,x,var):
        if var=='LC':
            self.LC=x
        elif var=='H1':
            self.H1=x
        else :
            self.H2=x
        self.initGrille()

if __name__ == '__main__':
    g=Channel()
    # np.set_printoptions(precision=2)
    #print(g.epsfunc(g.x,g.y))
    #g.plotindice()
    g.Traceindice()
    # g.TraceMode()

