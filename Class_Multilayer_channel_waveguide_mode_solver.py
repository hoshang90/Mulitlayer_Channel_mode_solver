
from tkinter import *
from tkinter import messagebox
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from Channels_multilayer import Channel
import matplotlib.pyplot as plt
from tkinter import filedialog
from scipy.stats import multivariate_normal # these are for the line profile profile
import scipy

class Mul_Ch_Wav_Mod_Sol(Frame):
    def __init__(self, master):
        super().__init__(master)
        self.master = master
        self.grid()
        self.place()
        # creates variable to be updated by the user. values are further obtained by the .get() method
        self.nC_1 = DoubleVar();self.nC_1.set(1.); self.nC_2 = DoubleVar();self.nC_2.set(1.); self.nC_3 = DoubleVar();self.nC_3.set(1.);
        self.nL1_1 = DoubleVar();self.nL1_1.set(1.);self.nL1_2 = DoubleVar();self.nL1_2.set(1.62);self.nL1_3 = DoubleVar();self.nL1_3.set(1.)
        self.nL2_1 = DoubleVar();self.nL2_1.set(1.62);self.nL2_2 = DoubleVar();self.nL2_2.set(1.62);self.nL2_3 = DoubleVar();self.nL2_3.set(1.62);
        self.nSub_1 = DoubleVar();self.nSub_1.set(1.61);self.nSub_2 = DoubleVar();self.nSub_2.set(1.61);self.nSub_3 = DoubleVar();self.nSub_3.set(1.61);
        self.TL_1 = DoubleVar();self.TL_1.set(2.3);self.TL_2 = DoubleVar();self.TL_2.set(2.);self.TChan = DoubleVar();self.TChan.set(3.1)
        self.WL = DoubleVar();self.WL.set(0.64);self.OPR = DoubleVar();self.OPR.set(4.0);self.vVaryH1LC_H1 = StringVar();self.vVaryH1LC_H1.set('1.0 to 4.0');
        self.vVaryH1LC_LC = StringVar();self.vVaryH1LC_LC.set('1.5 to 5.0');self.vVaryH1LC_H2 = DoubleVar();self.vVaryH1LC_H2.set(1.0);
        self.vVaryH1LC_SaveAs = StringVar();self.vVaryH1LC_SaveAs.set('Filename');self.vSimulat = DoubleVar();self.vSimulat.set(np.zeros((1,1)))
        self.LC_1= DoubleVar();self.LC_1.set(1.5); self.LC_2= DoubleVar();self.LC_2.set(5.0)
        self.H1_1=DoubleVar();self.H1_1.set(1.); self.H1_2=DoubleVar();self.H1_2.set(4.)
        self.H2value=DoubleVar();self.H2value.set(1.)
        self.grating_period = DoubleVar();self.grating_period.set(0.83);self.diffraction_mode = IntVar();self.diffraction_mode.set(1);
        self.Nmodes = IntVar();self.Nmodes.set(2)
        self.FileName=StringVar();self.FileName.set("test")
        #for the profile
        self.x_cut = DoubleVar();self.x_cut.set(3.1);self.y_cut = DoubleVar();self.y_cut.set(2.3)
        #self.vVaryH1LC_SaveAs_loc = StringVar();self.vVaryH1LC_SaveAs_loc.set('C:/Users/Home')#location to save
        self.create_labels()
        self.create_Entrys()
        self.create_buttons()

    def create_Entrys(self):
        Entry(self, width=10, borderwidth=5, textvariable=self.nC_1).grid(row=1, column=3, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nC_2).grid(row=1, column=2, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nC_3).grid(row=1, column=1, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nL1_1).grid(row=2, column=3, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nL1_2).grid(row=2, column=2, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nL1_3).grid(row=2, column=1, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nL2_1).grid(row=3, column=3, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nL2_2).grid(row=3, column=2, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nL2_3).grid(row=3, column=1, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nSub_1).grid(row=4, column=3, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nSub_2).grid(row=4, column=2, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.nSub_3).grid(row=4, column=1, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.TL_1).grid(row=2, column=4, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.TL_2).grid(row=3, column=4, padx=5, pady=5)
        Entry(self, width=10, borderwidth=5, textvariable=self.TChan).grid(row=6, column=2, padx=5, pady=5)
        Entry(self, width=8, borderwidth=5, textvariable=self.WL).grid(row=6, column=5, padx=5, pady=5)
        Entry(self, width=8, borderwidth=5, textvariable=self.OPR).grid(row=8, column=5, padx=5, pady=5)
        Entry(self, width=8, borderwidth=5, textvariable=self.H1_1).grid(row=13, column=0, padx=5, pady=5)
        Entry(self, width=8, borderwidth=5, textvariable=self.H1_2).grid(row=13, column=1, padx=5, pady=5)
        Entry(self, width=8, borderwidth=5, textvariable=self.LC_1).grid(row=13, column=2, padx=5, pady=5)
        Entry(self, width=8, borderwidth=5, textvariable=self.LC_2).grid(row=13, column=3, padx=5, pady=5)
        Entry(self, width=5, borderwidth=5, textvariable=self.H2value).grid(row=13, column=6, padx=5, pady=5)
        Entry(self.labelFrame1, width=10, borderwidth=5, textvariable=self.FileName).pack()
        # for the grating coupler
        Entry(self, width=5, borderwidth=5, textvariable=self.diffraction_mode).grid(row=16, column=5, padx=5, pady=5)
        Entry(self, width=5, borderwidth=5, textvariable=self.grating_period).grid(row=17, column=5, padx=5, pady=5)
        Entry(self, width=5, borderwidth=5, textvariable=self.Nmodes).grid(row=11, column=5, padx=5, pady=5)
        #self.Simulat_file = Entry(self, width=10, borderwidth=5, textvariable=self.vSimulat_file);self.Simulat_file.grid(row=14, column=6, padx=5, pady=5)
        ########
        Entry(self, width=5, borderwidth=5, textvariable=self.x_cut).grid(row=15, column=3, padx=5, pady=5)
        Entry(self, width=5, borderwidth=5, textvariable=self.y_cut).grid(row=15, column=4, padx=5, pady=5)

    def create_labels(self):
        Label(self, text = "Refractive indices").grid(row = 0, column= 2)
        Label(self, text = "    Thickness").grid(row = 0, column= 4)
        Label(self, text = '\u03BCm'+" (H1)",font = 'Verdana 10').place(x=415,y=63)
        Label(self, text = '\u03BCm'+" (H2)",font = 'Verdana 10').place(x=415,y=100)
        Label(self, text = "Cover").grid(row = 1)
        Label(self, text = "Layer 2").grid(row = 2)
        Label(self, text = "Layer 1").grid(row = 3)
        Label(self, text = "Substrate").grid(row = 4)
        Label(self, text = "     Channel  (LC)").grid(row = 5, column= 2)
        Label(self, text = '\u03BCm',font = 'Verdana 10').place(x=240,y=195) # unicode careters
        Label(self, text = "").grid(row = 7);Label(self, text = "").grid(row = 8);Label(self, text = "").grid(row = 9)#to creat space
        Label(self, text = "").grid(row = 10);Label(self, text = "").grid(row = 11);Label(self, text = "").grid(row = 12)
        Label(self, text = "").grid(row = 13);Label(self, text = "").grid(row = 14);Label(self, text = "").grid(row = 15)
        Label(self, text = "").grid(row = 16);Label(self, text = "").grid(row = 17);Label(self, text = "").grid(row = 18)
        Label(self, text = '  \u03BB',font = 'Times 12').place(x=435,y=170) # unicode careters
        Label(self, text ='\u03BCm',font = 'Times 12').place(x=490,y=195)
        Label(self, text = "OR",font = 'Times 12').place(x=430,y=225) # unicode careters
        Label(self, text =u'\u2070/mm',font = 'Verdana 10').place(x=490,y=255) # unicode careters
        self.labelFrame = LabelFrame(self, text="Choose file\nto plot");self.labelFrame.grid(column=5, row=14)#######################################################
        self.labelFrame1 = LabelFrame(self, text="Save as");self.labelFrame1.grid(column=2, row=14)
        #create the structure
        Label(self, text="<---- LC ---->").place(x=81, y=235)
        Label(self, text="-" * 15).place(x=80, y=250)
        Label(self, text="|").place(x=80, y=258)
        Label(self, text="|").place(x=157, y=258)
        Label(self, text="|").place(x=80, y=275)
        Label(self, text="|").place(x=157, y=275)
        Label(self, text="-" * 15).place(x=170, y=285)
        Label(self, text="\u2191", font=10).place(x=165, y=253)
        Label(self, text="\u2193", font=10).place(x=165, y=271)
        Label(self, text="H1").place(x=165, y=268.5)  # "H\u2081" H Subscript 1, Superscript is \uu207x(1,2,3...)
        Label(self, text="-" * 15).place(x=1, y=282)
        Label(self, text="-" * 50).place(x=1, y=320)
        Label(self, text="\u2191", font=10).place(x=245, y=288)
        Label(self, text="\u2193", font=10).place(x=245, y=306)
        Label(self, text="H2").place(x=245, y=303.5)
        Label(self, text="Vary H1 (\u03BCm)").place(x=20, y=415)
        Label(self, text="Start").place(x=5, y=370);Label(self, text="Stop").place(x=80, y=370)
        Label(self, text="Vary LC (\u03BCm) ").place(x=200, y=415)
        Label(self, text="Start").place(x=170, y=370);Label(self, text="Stop").place(x=260, y=370)
        Label(self, text="@ constant H2(\u03BCm) =").place(x=382, y=390)
        #Label(self, text="Save as:", font="Calibri 15").place(x=165, y=405)
        # Label(self, text="Choose file to plot", font="Calibri 12").place(x=365, y=405)
        Label(self, text="Diffraction order # =").place(x=320, y=605)
        Label(self, text="Grating period (\u039B) =").place(x=320, y=640)
        Label(self, text="# of modes =").place(x=355, y=332)
        Label(self, text="Horizontal & vertical lines intersecting @ X =").place(x=25, y=565)
        Label(self, text="& Y =").place(x=315, y=565)
    def create_buttons(self):
        Button(self, text="Plot_indices", fg='green', command=self.Plot_indices).grid(row=18, column=3)
        Button(self, text="Solve", fg='green', command=self.Solve).grid(row=18, column=2)
        Button(self, text="Close", fg='red', command=self.close_all).grid(row=18, column=0)
        Button(self, text="Close plots", fg='red', command= self.close_plots).grid(row=18, column=1)
        Button(self, text="Simulate", fg='Green', command=self.Simulate).grid(row=14, column=1)
        # Button(self, text="Plot image map\nof the simulation", fg='Green', command="").place(x=365, y=435)
        Button(self.labelFrame, text="Browse A File", command=self.fileDialog).grid(row=14, column=6)#######################################################
        Button(self.labelFrame, text="Plot Ecc.", command=self.plot_simulate_ecc).grid(row=15, column=6)
        Button(self.labelFrame, text="Plot \u0394n", command=self.plot_simulate_dn).grid(row=16, column=6)
        Button(self.labelFrame, text="Plot profile",fg='green', command=self.plot_line_profile).grid(row=17, column=6)
        Button(self, text="Coupling \u03B8", fg='green', command=self.NumOf_guided_modes).grid(row=18, column=5)#grating
        #Button(self.labelFrame1, text="Browse\n directory", command=self.fileDialog_saveas).grid(row=14, column=3)
    def Plot_indices(self):
        self.g = Channel(wl=self.WL.get(), pas=0.1, \
                       nc=self.nC_1.get(), nc_2=self.nC_2.get(), nc_3=self.nC_3.get(), Hc=0.8, \
                       LC=self.TChan.get(), LB=5., LB_R=5., \
                       n1B=self.nL1_1.get(), n1C=self.nL1_2.get(), n1D=self.nL1_3.get(), H1=self.TL_1.get(), \
                       n2B=self.nL2_1.get(), n2C=self.nL2_2.get(), n2D=self.nL2_3.get(), H2=self.TL_2.get(), \
                       nsub=self.nSub_1.get(), nsub_2=self.nSub_2.get(), nsub_3=self.nSub_3.get(), Hsub=2.,
                       OR=self.OPR.get())
        self.g.Traceindice()
    def fileDialog_saveas(self): # to browse locaation of a file
        self.vVaryH1LC_SaveAs_loc = filedialog.askdirectory(initialdir="", title="Select the location")#(("all files", "*.*"), ("csv files", "*.csv")))
    def fileDialog(self): # to brows the file
        self.vSimulat = filedialog.askopenfilename(initialdir="", title="Select A File", filetype=
        (("all files", "*.*"), ("csv files", "*.csv")))#(("all files", "*.*"), ("csv files", "*.csv")))
        self.df = pd.read_csv(self.vSimulat, sep=' ', comment='#', header=None)
        self.dr=self.df.to_numpy()
        #print(self.filename)
    def plot_line_profile(self):
        fig, main_ax = plt.subplots(figsize=(6, 6))
        divider = make_axes_locatable(main_ax)
        top_ax = divider.append_axes("top", 1.05, pad=0.1, sharex=main_ax)
        right_ax = divider.append_axes("right", 1.05, pad=0.1, sharey=main_ax)
        self.curX = self.x_cut.get()  # position of the vertical line  They should be always a float with one decimal like 1.1 or 1.2 etc...
        self.curY = self.y_cut.get()  # position of the horizontal line
        self.w_array = np.arange(self.LC_1.get(), self.LC_2.get(),0.1)  # introduce the x axis scale (xmin,xmax,step) we should know all these three parameters from the file we introduce in the  next step
        self.h_array = np.arange(self.H1_1.get(),self.H1_2.get(), 0.1)  # introduce the y axis scale (ymin,ymax,step)
        self.DnCB=((0.001*(float(self.WL.get()))*float(self.OPR.get()))/180)
        ecc = (self.DnCB / (self.dr + np.sqrt(self.DnCB ** 2 + self.dr ** 2)))  # define eccentricity matrix
        # make some labels invisible
        top_ax.xaxis.set_tick_params(labelbottom=False)
        right_ax.yaxis.set_tick_params(labelleft=False)

        main_ax.set_xlabel('W (\u03BCm)')
        main_ax.set_ylabel('H (\u03BCm)')
        top_ax.set_ylabel(r'E$_{cc}$')
        right_ax.set_xlabel(r'E$_{cc}$')
        z_max = 1  # z.max()
        self.curX = np.around(float(self.curX), 2)
        self.curY = np.around(float(self.curY), 2)
        # print((ecc[(np.argmax(np.where(np.around(self.w_array,2)==self.curY,self.w_array,0))),:]))############################
        im = main_ax.imshow(self.DnCB / (self.dr + np.sqrt(self.DnCB**2 + self.dr**2)),cmap="hot", extent=[self.LC_1.get(), self.LC_2.get(), self.H1_1.get(),self.H1_2.get()], origin='lower')
        main_ax.autoscale(enable=False)
        right_ax.autoscale(enable=False)
        top_ax.autoscale(enable=False)
        right_ax.set_xlim(right=z_max)
        top_ax.set_ylim(top=z_max)
        self.v_line = main_ax.axvline(self.curX, color='b')
        self.h_line = main_ax.axhline(self.curY, color='g')
        # print(ecc[:,(np.argmax(np.where(np.around(self.w_array,2)==self.curY,self.w_array,0)))])#############################
        self.v_prof, = right_ax.plot(ecc[:, (np.argmax(np.where(np.around(self.w_array, 2) == self.curX, self.w_array, 0)))], self.h_array,'b-')  # (np.argmax(np.where(np.around(self.h_array,2)==self.curY,self.h_array,0)))
        self.h_prof, = top_ax.plot(self.w_array, ecc[(np.argmax(np.where(np.around(self.h_array, 2) == self.curY, self.h_array, 0))), :],'g-')  # (np.argmax(np.where(np.around(self.w_array,2)==self.curX,self.w_array,0)))
        # define the colorbar##################################
        cax = divider.new_vertical(size="5%", pad=0.4, title="Eccentricity")
        fig.add_axes(cax)
        fig.colorbar(im, cax=cax, orientation="horizontal")
        cax.set_xlim(0, 1)
        cax.set
        # plt.savefig('colorbar_positioning_03.png', format='png', bbox_inches='tight')##################
        plt.show()
    def plot_simulate_ecc(self):
        fig, ax = plt.subplots(figsize=(8, 6))
        #plt.title("Eccentricy as a function of channel dimensions")
        print("Optcal rotaion is: "+ str(self.OPR.get()))
        self.DnCB=((0.001*(float(self.WL.get()))*float(self.OPR.get()))/180)
        print("Circular bireferengence (CB) is: "+ str(self.DnCB))
        im = plt.imshow(self.DnCB / (self.dr + np.sqrt(self.DnCB**2 + self.dr**2)),cmap="hot", extent=[self.LC_1.get(), self.LC_2.get(), self.H1_1.get(),self.H1_2.get()], origin='lower')
        plt.xlabel("H1 (\u03BCm)")
        plt.ylabel("LC (\u03BCm)")
        divider = make_axes_locatable(ax)
        cax = divider.new_vertical(size="5%", pad=0.4, title="Eccentricity")
        fig.add_axes(cax)
        fig.colorbar(im, cax=cax, orientation="horizontal")
####################################################################################################### mshu x y labela wash karw
        # plt.savefig('colorbar_positioning_03.png', format='png', bbox_inches='tight')
        plt.show()
        plt.close()
    def printt(self):
        print(self.filename)
    def plot_simulate_dn(self):

        fig, ax = plt.subplots(figsize=(8, 6))
        #plt.title("Modal bireferengence as a function of channel dimensions")
        im = plt.imshow(self.dr,cmap="hot",  extent=[self.LC_1.get(), self.LC_2.get(), self.H1_1.get(),self.H1_2.get()], origin='lower')
        plt.xlabel("H1 (\u03BCm)")
        plt.ylabel("LC (\u03BCm)")
        divider = make_axes_locatable(ax)
        cax = divider.new_vertical(size="5%", pad=0.4, title="Modal bireferengence")
        fig.add_axes(cax)
        fig.colorbar(im, cax=cax,orientation="horizontal")

        # plt.savefig('colorbar_positioning_03.png', format='png', bbox_inches='tight')
        plt.show()
        plt.close()
    def close_all(self):
        # Button(self.root,text = 'Click Me', command=lambda:[self.funcA(), self.funcB(), self.funcC()])
        self.master.destroy()
        plt.close('all')
    def Solve(self):
        self.g = Channel(wl=self.WL.get(), pas=0.1, \
                       nc=self.nC_1.get(), nc_2=self.nC_2.get(), nc_3=self.nC_3.get(), Hc=0.8, \
                       LC=self.TChan.get(), LB=5., LB_R=5., \
                       n1B=self.nL1_1.get(), n1C=self.nL1_2.get(), n1D=self.nL1_3.get(), H1=self.TL_1.get(), \
                       n2B=self.nL2_1.get(), n2C=self.nL2_2.get(), n2D=self.nL2_3.get(), H2=self.TL_2.get(), \
                       nsub=self.nSub_1.get(), nsub_2=self.nSub_2.get(), nsub_3=self.nSub_3.get(), Hsub=2.,
                       OR=self.OPR.get())
        self.g.Calcule(Nmodes=self.Nmodes.get())
        #self.g.Intensite()
        plt.show()

    def NumOf_guided_modes(self): # This will print angles which can give coupled light.
        print("Coupling neff is calculated with equation\n n_eff=n_top.sin(\u03B8) + (m.(\u03BB/\u039B) \n where "
              "n_eff is effective refractive index, n_top is refractive index of the cover\n"
              "\u03B8 is coupling angle,m is diffraction mode, \u03BB is wavelength, and \u039B is grating period")
        self.g = Channel(wl=self.WL.get(), pas=0.1, \
                       nc=self.nC_1.get(), nc_2=self.nC_2.get(), nc_3=self.nC_3.get(), Hc=0.8, \
                       LC=self.TChan.get(), LB=5., LB_R=5., \
                       n1B=self.nL1_1.get(), n1C=self.nL1_2.get(), n1D=self.nL1_3.get(), H1=self.TL_1.get(), \
                       n2B=self.nL2_1.get(), n2C=self.nL2_2.get(), n2D=self.nL2_3.get(), H2=self.TL_2.get(), \
                       nsub=self.nSub_1.get(), nsub_2=self.nSub_2.get(), nsub_3=self.nSub_3.get(), Hsub=2.,
                       OR=self.OPR.get())
        if self.g.NmaxPlan()==0: # blow H2==0.9 NmaxPlan is 0 and it does not give good value.
            messagebox.showinfo("Hey", 'It looks like NmaxPlan is {}'.format(self.g.NmaxPlan()))
            pass
        else:
            for i in np.real((self.g.NumOf_guidedModes(Nmodes=self.Nmodes.get())[1:])):
                print("For the {:.4f} mode there is one coupling angle at {:.4f}".format(i,self.Nmode_effective(i)))
    
    def Nmode_effective(self,n_eff): # to find n_e that is effective refractive index of the grating
        '''This method returns guided modes of the structure (not leaky modes)'''
        #print(self.grating_period,self.diffraction_mode)
        n_top=np.average(np.array([self.nC_1.get(),self.nC_2.get(),self.nC_3.get()]))
        m=self.diffraction_mode.get(); WL=self.WL.get(); GP=self.grating_period.get()
        TheSin=( n_eff-m*WL/GP )/n_top
        if np.abs(TheSin) >1:
            theta = 0
        else:
            theta=np.rad2deg(np.arcsin(TheSin) )
        return theta


    def Simulate(self):
        if messagebox.askyesno("Warning! ", 'This may take several minutes\n Are you sure you want to continue?'):
            messagebox.showinfo("Simulation",'You can see the see the remaining time in the run window')
            self.g = Channel(wl=self.WL.get(), pas=0.1, \
                       nc=self.nC_1.get(), nc_2=self.nC_2.get(), nc_3=self.nC_3.get(), Hc=0.8, \
                       LC=self.TChan.get(), LB=5., LB_R=5., \
                       n1B=self.nL1_1.get(), n1C=self.nL1_2.get(), n1D=self.nL1_3.get(), H1=self.TL_1.get(), \
                       n2B=self.nL2_1.get(), n2C=self.nL2_2.get(), n2D=self.nL2_3.get(), H2=self.TL_2.get(), \
                       nsub=self.nSub_1.get(), nsub_2=self.nSub_2.get(), nsub_3=self.nSub_3.get(), Hsub=2.,
                       OR=self.OPR.get())
            self.g.H2 = self.H2value.get()
            self.g.VariaXY(X0=self.LC_1.get(), X1=self.LC_2.get(), dX=0.1, Y0=self.H1_1.get(), Y1=self.H1_2.get(), dY=0.1, fich=self.FileName.get())
        else:
            messagebox.showinfo("Canceled",'You have successfully canceled the simulation')
    def close_plots(self):
        plt.close('all')
if __name__ == '__main__':
    root = Tk()  # we car write any name instead of root
    root.title('Multilayer channel waveguide mode solver')
    # root.iconbitmap('rib_waveguide.ico')
    root.geometry("550x700")
    app = Mul_Ch_Wav_Mod_Sol(root)
    app.mainloop()
