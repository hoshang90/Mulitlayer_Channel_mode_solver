#!/usr/bin/python3
"""Fully vectorial finite-difference mode solver example."""
# utilisation avec ipython
#---------------------------------------------------------------
# generation du profile: ipython3 profile create channelWG
# puis editer .ipython/profile_channelWG/ipython_config.py

# import numpy as np
# np.set_printoptions(linewidth=180,threshold=5000)
#'import RidgeAchiraux as ra',
#'g=ra.Ridge()',
#'dt=np.arange(1,10,1)',
#'%matplotlib'

# lancer ipython3 avec:
# ipython3 --profile=channelWG
#---------------------------------------------------------------
# g=ra.Ridge(nsub=2,nguide=3,ncover=1,wl=1,largeur=2,hauteur=6,souscouche=2,ex=1,ey=1,pas=0.1)
# g.Traceindice()
# g.Calcule(2) => calcule des modes [0..1]
# g.Infos()
# g.TraceModes(0)

#---------------------------------------------------------------
# import matplotlib.pyplot as plt  #=> config
#neff=np.arange(0.900,0.990,0.005);dneff=g.TraceDetneff(neff)
#plt.plot(neff,dneff,'o');plt.show()
#-----------------------------------------

import numpy as np
import EMpy
import matplotlib.pyplot as plt
import time
class Ridge():
    """ modes propres dans guide canal de type ridge 
    basé sur EMpy: https://github.com/lbolla/EMpy
    grandeurs en microns
    grille de pas: pas ds les deux directions
    hauteur grille (y): Hc+H1+H2+Hsub
    largeur grille (x): LB+LC+LB
    ___________________________________
        nc, Hc
    ___________________________________
    <-- LB --> <-- LC --> <-- LB -->
        n1B   |   n1C    | n1B , H1
    --------------------------------------
    <-- LB --> <-- LC --> <-- LB -->
        n2B   |   n2C    | n2B  , H2
    -------------------------------------
        nsub,Hsub
    _______________________________________
    """
    def __init__(self,wl=0.64,pas=0.16,\
            nc=1.33,Hc=0.5,\
            LB=3,LC=2,\
            n1B=1.33,n1C=1.620,H1=0.5,\
            n2B=1.620,n2C=1.620,H2=1,\
            nsub=1.615,Hsub=1.5,OR=4,ng=0
            ):
        # --------- ridge
        self.nc2=nc*nc;self.Hc=Hc
        self.LB=LB;self.LC=LC
        # les bords dans le cladding
        self.n1B2=self.nc2
        # n1c=n2
        self.n1C2=ng*ng
        self.n2B2=self.n1C2
        self.n2C2=self.n1C2
        #self.n1C2=n1C*n1C; self.n1B2=n1B*n1B;
        self.H1=H1
        #self.n2C2=n2C*n2C; self.n2B2=n2B*n2B;
        self.H2=H2
        self.ns2=nsub*nsub;self.Hsub=Hsub
        self.wl=wl
        self.pas=pas
        self.OR=OR/360.

        self.initGrille()
        print("Grille {}x{} micron2".format(self.LargGrille,self.HautGrille))
        print("pas={}, Npts={}".format(self.pas,self.x.size*self.y.size))

    def initGrille(self):
        ''' Initizalisation de la grille de calcul '''
        pas=self.pas
        self.HautGrille=self.Hc+self.H1+self.H2+self.Hsub
        self.LargGrille=self.LC+2*self.LB
        self.x = np.arange(0, self.LargGrille+pas, pas)
        self.y = np.arange(0, self.HautGrille+pas, pas)



    def __repr__(self):
        reponse="-"*20+"\n"
        interf="\n"+reponse
        reponse+=" nc={:.3f}, Hc={}\n".format(np.sqrt(self.nc2),self.Hc)
        reponse+="<- LB={} -->| <- LC={} ->| <- LB ->\n".format(self.LB,self.LC)
        reponse+="_"*30+"\n"
        reponse+=" n1B={:.3f}  |   n1C={:.3f}  | n1B , H1={}\n".\
                format(np.sqrt(self.n1B2),np.sqrt(self.n1C2),self.H1)
        reponse+="-"*20+"\n"
        reponse+="n2B={:.3f}   |   n2C={:.3f}  | n2B , H2={:.2f}\n".\
                format(np.sqrt(self.n2B2),np.sqrt(self.n2C2),self.H2)
        reponse+="-"*30+"\n"
        reponse+=" nsub={:.3f}, Hsub={}\n".format(np.sqrt(self.ns2),self.Hsub)
        reponse+="-"*30+"\n"
        Nmodes=len(self.sol.modes)
        for i in range(Nmodes):
            neff=self.sol.modes[i].neff
            reponse+="*Mode[{}]: neff={:.3f}, ".format(i,neff)
            reponse+="%TE= {:.2f}, ".format(100*self.sol.modes[i].TEfrac())
            reponse+="%Conf= {:.2f}\n".format(100*self.Confinement(num=i))
        DN=self.sol.modes[0].neff-self.sol.modes[1].neff
        reponse+="n0-n1={:g}\n".format(DN)
        ORN=self.OR*1e-3*self.wl
        ecc=ORN/(DN+np.sqrt(ORN*ORN+DN*DN))
        reponse+="ecc={:.3f}\n".format(ecc)
        return reponse

    def Calcule(self,Nmodes=2,verbose=False,trace=True):
        neigs=Nmodes
        self.initGrille()
        tol=1e-8
        boundary='0000'
        self.sol=EMpy.modesolvers.FD.VFDModeSolver(self.wl, self.x, self.y,\
                self.epsfunc, boundary).solve(neigs, tol)
        print(self.__repr__())
        if verbose:
            print("-------------------------")
            print("g.Infos()")
            print("g.TraceMode(num=0): champ du mode num")
            print("I=g.Intensite(num=0): intensité du mode num")
            print("g.TraceEx(num=0): profil au centre")
        if trace:
            self.Intensite()
    

    def epsfunc(self,x_, y_):
        """Return a matrix describing a 2d material.
        :param x_: x values
        :param y_: y values
        :return: 2d-matrix
        """
        xx, yy = np.meshgrid(x_,y_)
        #on remplit  en partant de la couche 1
        cg=np.where( (np.abs(yy.T-0.5*self.H1-self.H2-self.Hsub)<=0.5*self.H1)
            *(np.abs(xx.T-self.LargGrille/2.)<self.LC/2),self.n1C2,self.n1B2)
        # cover: 
        cg=np.where(yy.T>=self.H1+self.H2+self.Hsub,self.nc2,cg)
        #couche 2 /bords
        cg=np.where( (np.abs(yy.T-0.5*self.H2-self.Hsub)<=0.5*self.H2),self.n2B2,cg)
        #couche 2 /centre
        cg=np.where( (np.abs(yy.T-0.5*self.H2-self.Hsub)<=0.5*self.H2)
            *(np.abs(xx.T-self.LargGrille/2.)<self.LC/2),self.n2C2,cg)
        # substrat
        return np.where(yy.T<= self.Hsub, self.ns2,cg)

    def Traceindice(self):
        ''' plot l'indice '''
        indices=self.epsfunc(self.x,self.y)
        levels=np.sqrt([1,self.nc2,self.n1B2,self.n1C2,self.n2B2,self.n2C2,self.ns2])
        fig=plt.figure()
        plt.contour(self.x,self.y,np.sqrt(indices.T),np.unique(levels))
        plt.contourf(self.x,self.y,np.sqrt(indices.T),np.unique(levels))
        plt.colorbar()
        plt.show()


    def Intensite(self,num=0):
        ''' plot l'intensité '''
        #le champ E a un pont de moins, on recupere les bons x
        x0=EMpy.utils.centered1d(self.x)
        y0= EMpy.utils.centered1d(self.y)
        Itot=np.abs(self.sol.modes[num].intensity())
        fig = plt.figure()
        plt.hot()
        levels=np.unique([1,self.nc2,self.n1B2,self.n1C2,self.n2B2,self.n2C2,self.ns2])
        Nlevel=15
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(x0,y0 ,Itot.T, Nlevel)
        print("%Conf: {:.2f}".format(100*self.Confinement(num=num)))
        print("%TE: {:.1f}".format(100*self.sol.modes[num].TEfrac()))
        #return Itot

    def Confinement(self,num=0):
        ''' retourne  I(cover)/Itotale '''
        #le champ E a un pont de moins, on recupere les bons x
        x0=EMpy.utils.centered1d(self.x)
        y0= EMpy.utils.centered1d(self.y)
        Itot=np.abs(self.sol.modes[num].intensity())
        epsi=self.epsfunc(x0,y0)
        masque=np.where(epsi<self.ns2,True,False)
        return np.sum(Itot*masque)/np.sum(Itot)


    def TraceMode(self,num=0):
        ''' trace les 6 composantes du champ '''
        Nlevel=15
        levels=np.unique([1,self.nc2,self.n1B2,self.n1C2,self.n2B2,self.n2C2,self.ns2])
        #le champ E a un pont de moins, on recupere les bons x
        x0=EMpy.utils.centered1d(self.x)
        y0= EMpy.utils.centered1d(self.y)

        fig = plt.figure(figsize=(12, 6))
        plt.hot()
        fig.add_subplot(2, 3, 1)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(x0,y0,abs(self.sol.modes[num].Ex).T, Nlevel)
        plt.title('Ex')
        plt.colorbar()

        fig.add_subplot(2, 3, 2)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(x0,y0,abs(self.sol.modes[num].Ey).T, Nlevel)
        plt.title('Ey')
        plt.colorbar()

        fig.add_subplot(2, 3, 3)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(x0,y0,abs(self.sol.modes[num].Ez).T, Nlevel)
        plt.title('Ez')
        plt.colorbar()

        fig.add_subplot(2, 3, 4)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(self.x,self.y,abs(self.sol.modes[num].Hx).T, Nlevel)
        plt.title('Hx')
        plt.colorbar()

        fig.add_subplot(2, 3, 5)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(self.x,self.y,abs(self.sol.modes[num].Hy).T, Nlevel)
        plt.title('Hy')
        plt.colorbar()

        fig.add_subplot(2, 3, 6)
        plt.contour(self.x,self.y,self.epsfunc(self.x,self.y).T,levels,colors='white')
        plt.contourf(self.x,self.y,abs(self.sol.modes[num].Hz).T, Nlevel)
        plt.title('Hz')
        plt.colorbar()

        plt.show()

    def TraceEx(self,num=0):
        ''' Trace le profil de Ex au milieu du guide'''
        milieu=int(self.sol.modes[num].Ex.shape[0]/2.)
        Ex=self.sol.modes[num].Ex[milieu,:]
        scale=np.sqrt(self.n1C2)/np.amax(np.abs(Ex))
        indices=np.sqrt(self.epsfunc(self.x,self.y))[milieu,:]
        ymax=self.Hc
        ymin=-(self.H1+self.H2+self.Hsub)
        fig=plt.figure()
        plt.plot(np.linspace(ymin,ymax,num=indices.size),indices,np.linspace(ymin,ymax,num=Ex.size),np.abs(Ex)*scale)
        plt.grid()
        plt.xlabel(r'$\mu m$')
        plt.show()

    def VariaX(self,X0=0.1, X1=5, dX=0.2,variable='LC',trace=True):
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
            dX=0.2,vX='LC',Y0=0.1,Y1=5,dY=0.2,vY='H1'):
        lesX=np.arange(X0,X1,dX)
        lesY=np.arange(Y0,Y1,dY)
        self.lesDn=np.zeros((lesX.size,lesY.size))
        self.lesSc=np.zeros((lesX.size,lesY.size))
        for i,X in enumerate(lesX):
            self.FixeVariable(X,vX)
            t0= time.clock()
            for j,Y in enumerate(lesY):
                self.FixeVariable(Y,vY)
                self.Calcule(verbose=False)
                dn=abs(self.sol.modes[0].neff-self.sol.modes[1].neff)
                Sc=self.Confinement()
                self.lesDn[i,j]=dn
                self.lesSc[i,j]=Sc
                #print(j,end=",")
                #print("\n")
            te=time.clock() - t0
            if i==0:
                print("Boucle {:.2f}, total  {:.2f}s".format(te,te*(lesX.size-i)))
        print("---------------")
        print('np.savetxt("dataDn",g.lesDn.T,fmt=\'%.2g\');',end=" ")
        print('np.savetxt("dataSc",g.lesSc.T,fmt=\'%.2g\')')
        print("plot \"dataDn\" matrix u ($1*{}+{}):($2*{}+{}):(-log10($3)) w image".format(dX,X0,dY,Y0))
        print("plot \"dataSc\" matrix u ($1*{}+{}):($2*{}+{}):(-log10($3)) w image".format(dX,X0,dY,Y0))



    def FixeVariable(self,x,var):
        if var=='LC':
            self.LC=x
        elif var=='H1':
            self.H1=x
        else :
            self.H2=x
        self.initGrille()


            

        
        


if __name__ == '__main__':
    g=Ridge()
    np.set_printoptions(precision=2)
    #print(g.epsfunc(g.x,g.y))
    #g.plotindice()
    g.Calcule()
    g.TraceMode()

