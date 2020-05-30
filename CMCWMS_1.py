
from tkinter import *
from tkinter import messagebox
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import numpy as np
from Channels_multilayer import Channel
import matplotlib.pyplot as plt
from tkinter import filedialog

class Mul_Ch_Wav_Mod_Sol(Frame):
    '''Help on Class multilayer channel waveguide mode solver:
    Usage:
    It can be used to solve light confinement modes inside a channel dielectric chiral structure.
    '''
    def __init__(self, master):
        super().__init__(master)
        self.master = master
        self.grid()
        self.place()
        # creat variable to be updated by the user
        self.vnc1 = DoubleVar();self.vnc1.set(1.); self.vnc2 = DoubleVar();self.vnc2.set(1.); self.vnc3 = DoubleVar();self.vnc3.set(1.);
        self.vnL1_1 = DoubleVar();self.vnL1_1.set(1.);self.vnL1_2 = DoubleVar();self.vnL1_2.set(1.62);self.vnL1_3 = DoubleVar();self.vnL1_3.set(1.)
        self.vnL2_1 = DoubleVar();self.vnL2_1.set(1.62);self.vnL2_2 = DoubleVar();self.vnL2_2.set(1.62);self.vnL2_3 = DoubleVar();self.vnL2_3.set(1.62);
        self.vnSub1 = DoubleVar();self.vnSub1.set(1.61);self.vnSub2 = DoubleVar();self.vnSub2.set(1.61);self.vnSub3 = DoubleVar();self.vnSub3.set(1.61);
        self.vTL1 = DoubleVar();self.vTL1.set(2.3);self.vTL2 = DoubleVar();self.vTL2.set(2.);self.vTc = DoubleVar();self.vTc.set(3.1)
        self.vWL = DoubleVar();self.vWL.set(0.64);self.vOR = DoubleVar();self.vOR.set(4.0);self.vVaryH1LC_H1 = StringVar();self.vVaryH1LC_H1.set('1.0 to 4.0');
        self.vVaryH1LC_LC = StringVar();self.vVaryH1LC_LC.set('1.5 to 5.0');self.vVaryH1LC_H2 = DoubleVar();self.vVaryH1LC_H2.set(1.0);
        self.vVaryH1LC_SaveAs = StringVar();self.vVaryH1LC_SaveAs.set('Filename');self.vSimulat = DoubleVar();self.vSimulat.set(np.zeros((1,1)))

        self.vgrating_period = DoubleVar();self.vgrating_period.set(0.83);self.vdiffraction_mode = DoubleVar();self.vdiffraction_mode.set(1.0);
        self.vNmodes = IntVar();self.vNmodes.set(2)
        #self.vVaryH1LC_SaveAs_loc = StringVar();self.vVaryH1LC_SaveAs_loc.set('C:/Users/Home')#location to save
        self.create_labels()
        self.create_Entrys()
        self.create_buttons()
    def create_Entrys(self):
        '''Here all the entries are defined'''
        self.nCover_1 = Entry(self, width=10, borderwidth=5, textvariable=self.vnc1);self.nCover_1.grid(row=1, column=3, padx=5, pady=5)
        self.nCover_2 = Entry(self, width=10, borderwidth=5, textvariable=self.vnc2);self.nCover_2.grid(row=1, column=2, padx=5, pady=5)
        self.nCover_3 = Entry(self, width=10, borderwidth=5, textvariable=self.vnc3);self.nCover_3.grid(row=1, column=1, padx=5, pady=5)
        self.nLayer1_1 = Entry(self, width=10, borderwidth=5, textvariable=self.vnL1_1);self.nLayer1_1.grid(row=2, column=3, padx=5, pady=5)
        self.nLayer1_2 = Entry(self, width=10, borderwidth=5, textvariable=self.vnL1_2);self.nLayer1_2.grid(row=2, column=2, padx=5, pady=5)
        self.nLayer1_3 = Entry(self, width=10, borderwidth=5, textvariable=self.vnL1_3);self.nLayer1_3.grid(row=2, column=1, padx=5, pady=5)
        self.nLayer2_1 = Entry(self, width=10, borderwidth=5, textvariable=self.vnL2_1);self.nLayer2_1.grid(row=3, column=3, padx=5, pady=5)
        self.nLayer2_2 = Entry(self, width=10, borderwidth=5, textvariable=self.vnL2_2);self.nLayer2_2.grid(row=3, column=2, padx=5, pady=5)
        self.nLayer2_3 = Entry(self, width=10, borderwidth=5, textvariable=self.vnL2_3);self.nLayer2_3.grid(row=3, column=1, padx=5, pady=5)
        self.nSubstrate_1 = Entry(self, width=10, borderwidth=5, textvariable=self.vnSub1);self.nSubstrate_1.grid(row=4, column=3, padx=5, pady=5)
        self.nSubstrate_2 = Entry(self, width=10, borderwidth=5, textvariable=self.vnSub2);self.nSubstrate_2.grid(row=4, column=2, padx=5, pady=5)
        self.nSubstrate_3 = Entry(self, width=10, borderwidth=5, textvariable=self.vnSub3);self.nSubstrate_3.grid(row=4, column=1, padx=5, pady=5)
        self.Thickness_Layer_1 = Entry(self, width=10, borderwidth=5, textvariable=self.vTL1);self.Thickness_Layer_1.grid(row=2, column=4, padx=5, pady=5)
        self.Thickness_Layer_2 = Entry(self, width=10, borderwidth=5, textvariable=self.vTL2);self.Thickness_Layer_2.grid(row=3, column=4, padx=5, pady=5)
        self.Thickness_channel = Entry(self, width=10, borderwidth=5, textvariable=self.vTc);self.Thickness_channel.grid(row=6, column=2, padx=5, pady=5)
        self.Wavelength = Entry(self, width=10, borderwidth=5, textvariable=self.vWL);self.Wavelength.grid(row=6, column=5, padx=5, pady=5)
        self.Optical_rotation = Entry(self, width=10, borderwidth=5, textvariable=self.vOR);self.Optical_rotation.grid(row=8, column=5, padx=5, pady=5)
        self.VaryH1LC_H1 = Entry(self, width=10, borderwidth=5, textvariable=self.vVaryH1LC_H1);self.VaryH1LC_H1.grid(row=13, column=1, padx=5, pady=5)
        self.VaryH1LC_LC = Entry(self, width=10, borderwidth=5, textvariable=self.vVaryH1LC_LC);self.VaryH1LC_LC.grid(row=13, column=3, padx=5, pady=5)
        self.VaryH1LC_H2 = Entry(self, width=5, borderwidth=5, textvariable=self.vVaryH1LC_H2);self.VaryH1LC_H2.grid(row=13, column=6, padx=5, pady=5)
        self.VaryH1LC_SaveAs = Entry(self.labelFrame1, width=10, borderwidth=5, textvariable=self.vVaryH1LC_SaveAs);self.VaryH1LC_SaveAs.grid(row=15, column=3, padx=5, pady=5)

        self.diffraction_mode_E = Entry(self, width=5, borderwidth=5, textvariable=self.vdiffraction_mode);self.diffraction_mode_E.grid(row=16, column=5, padx=5, pady=5)
        self.grating_period_E = Entry(self, width=5, borderwidth=5, textvariable=self.vgrating_period);self.grating_period_E.grid(row=17, column=5, padx=5, pady=5)
        self.NModes_E = Entry(self, width=5, borderwidth=5, textvariable=self.vNmodes);self.NModes_E.grid(row=10, column=5, padx=5, pady=5)
        #self.Simulat_file = Entry(self, width=10, borderwidth=5, textvariable=self.vSimulat_file);self.Simulat_file.grid(row=14, column=6, padx=5, pady=5)
        ########
        self.nC_1 = float(self.nCover_1.get());self.nC_2 = float(self.nCover_2.get());self.nC_3 = float(self.nCover_3.get())
        self.nL1_1 = float(self.nLayer1_1.get());self.nL1_2 = float(self.nLayer1_2.get());self.nL1_3 = float(self.nLayer1_3.get())
        self.nL2_1 = float(self.nLayer2_1.get());self.nL2_2 = float(self.nLayer2_2.get());self.nL2_3 = float(self.nLayer2_3.get())
        self.nSub_1 = float(self.nSubstrate_1.get());self.nSub_2 = float(self.nSubstrate_2.get());self.nSub_3 = float(self.nSubstrate_3.get())
        self.TL_2 = float(self.Thickness_Layer_2.get());self.TL_1 = float(self.Thickness_Layer_1.get());self.TChan = float(self.Thickness_channel.get())
        self.WL = float(self.Wavelength.get());self.OPR = float(self.Optical_rotation.get());self.VaryH1LC_H1_array = self.VaryH1LC_H1.get()
        self.VaryH1LC_LC_array = self.VaryH1LC_LC.get();self.VaryH1LC_H2_Value = float(self.VaryH1LC_H2.get())
        # print(self.VaryH1LC_H1_array)
        self.H1_1 = np.fromstring(self.VaryH1LC_H1_array, dtype=float, sep='to')[0]
        self.H1_2 = np.fromstring(self.VaryH1LC_H1_array, dtype=float, sep='to')[1]
        self.LC_1 = np.fromstring(self.VaryH1LC_LC_array, dtype=float, sep='to')[0]
        self.LC_2 = np.fromstring(self.VaryH1LC_LC_array, dtype=float, sep='to')[1]
        #for the grating coupler
        self.filename = self.VaryH1LC_SaveAs.get();self.grating_period=self.grating_period_E.get();
        self.diffraction_mode=self.diffraction_mode_E.get();self.Nmodes_value=int(self.NModes_E.get())
    def create_labels(self):
        '''Here all the labels are defined'''
        Label(self, text = "Refractive indices").grid(row = 0, column= 2)
        Label(self, text = "    Thickness").grid(row = 0, column= 4)
        Label(self, text = '\u03BCm'+" (H1)",font = 'Verdana 10').place(x=400,y=63)
        Label(self, text = '\u03BCm'+" (H2)",font = 'Verdana 10').place(x=400,y=100)
        Label(self, text = "Cover").grid(row = 1)
        Label(self, text = "Layer 2").grid(row = 2)
        Label(self, text = "Layer 1").grid(row = 3)
        Label(self, text = "Substrate").grid(row = 4)
        Label(self, text = "     Channel  (LC)").grid(row = 5, column= 2)
        Label(self, text = '\u03BCm',font = 'Verdana 10').place(x=225,y=195) # unicode careters
        Label(self, text = "").grid(row = 7);Label(self, text = "").grid(row = 8);Label(self, text = "").grid(row = 9)#to creat space
        Label(self, text = "").grid(row = 10);Label(self, text = "").grid(row = 11);Label(self, text = "").grid(row = 12)
        Label(self, text = "").grid(row = 13);Label(self, text = "").grid(row = 14);Label(self, text = "").grid(row = 15)
        Label(self, text = "").grid(row = 16);Label(self, text = "").grid(row = 17);Label(self, text = "").grid(row = 18)
        Label(self, text = '  \u03BB',font = 'Verdana 10').place(x=420,y=172) # unicode careters
        Label(self, text ='\u03BCm',font = 'Times 12').place(x=480,y=195)
        Label(self, text = " OR",font = 'Times 12').place(x=415,y=230) # unicode careters
        Label(self, text =u'\u2070/mm',font = 'Verdana 10').place(x=480,y=255) # unicode careters
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
        Label(self, text="Vary\nH1 (\u03BCm) =").place(x=0, y=375)
        Label(self, text="Vary\nLC (\u03BCm) =").place(x=180, y=375)
        Label(self, text="@ constant H2 (\u03BCm) =").place(x=360, y=390)
        #Label(self, text="Save as:", font="Calibri 15").place(x=165, y=405)
        # Label(self, text="Choose file to plot", font="Calibri 12").place(x=365, y=405)
        Label(self, text="Diffraction mode # =").place(x=300, y=560)
        Label(self, text="Grating period (\u039B) =").place(x=305, y=595)
        Label(self, text="# of modes =").place(x=340, y=310)
    def create_buttons(self):
        '''Here all the buttons are defined'''
        Button(self, text="Plot_indices", fg='green', command=self.Plot_indices).grid(row=18, column=3)
        Button(self, text="Solve", fg='green', command=self.Solve).grid(row=18, column=2)
        Button(self, text="Close", fg='red', command=self.close_all).grid(row=18, column=0)
        Button(self, text="Close plots", fg='red', command= self.close_plots).grid(row=18, column=1)
        Button(self, text="Simulate", fg='Green', command=self.Simulate).grid(row=14, column=1)
        # Button(self, text="Plot image map\nof the simulation", fg='Green', command="").place(x=365, y=435)
        Button(self.labelFrame, text="Browse A File", command=self.fileDialog).grid(row=14, column=6)#######################################################
        Button(self.labelFrame, text="Plot Ecc.", command=self.plot_simulate_ecc).grid(row=15, column=6)
        Button(self.labelFrame, text="Plot \u0394n", command=self.plot_simulate_dn).grid(row=16, column=6)
        Button(self, text="Coupling \u03B8", fg='green', command=self.NumOf_guided_modes).grid(row=18, column=5)#grating
        #Button(self.labelFrame1, text="Browse\n directory", command=self.fileDialog_saveas).grid(row=14, column=3)
    def Plot_indices(self):
        '''This method plots refractive index of the structure'''
        self.create_Entrys()
        self.g = Channel(wl=self.WL, pas=0.1, \
                       nc=self.nC_1, nc_2=self.nC_2, nc_3=self.nC_3, Hc=0.8, \
                       LC=self.TChan, LB=5., LB_R=5., \
                       n1B=self.nL1_1, n1C=self.nL1_2, n1D=self.nL1_3, H1=self.TL_1, \
                       n2B=self.nL2_1, n2C=self.nL2_2, n2D=self.nL2_3, H2=self.TL_2, \
                       nsub=self.nSub_1, nsub_2=self.nSub_2, nsub_3=self.nSub_3, Hsub=2.,
                       OR=self.OPR)
        self.g.Traceindice()
    def fileDialog_saveas(self): # to browse locaation of a file
        self.vVaryH1LC_SaveAs_loc = filedialog.askdirectory(initialdir="", title="Select the location")#(("all files", "*.*"), ("csv files", "*.csv")))
    def fileDialog(self): # to brows the file
        '''This method gets the calls the location of a file from user call'''
        self.vSimulat = filedialog.askopenfilename(initialdir="", title="Select A File", filetype=
        (("all files", "*.*"), ("csv files", "*.csv")))#(("all files", "*.*"), ("csv files", "*.csv")))
        self.df = pd.read_csv(self.vSimulat, sep=' ', comment='#', header=None)
        self.dr=self.df.to_numpy()
        #print(self.filename)
    def plot_simulate_ecc(self):
        '''This method plots eccentricity of the simulated data in the file called by the method fileDialog'''
        self.create_Entrys()
        fig, ax = plt.subplots(figsize=(8, 6))
        #plt.title("Eccentricy as a function of channel dimensions")
        print("Optcal rotaion is: "+ str(self.OPR))
        self.DnCB=((0.001*(self.WL)*self.OPR)/180) #CB=self.OPR*1e-3/180*self.WL
        print("Circular bireferengence (CB) is: "+ str(self.DnCB))
        im = plt.imshow(self.DnCB / (self.dr + np.sqrt(self.DnCB**2 + self.dr**2)), extent=[self.LC_1, self.LC_2, self.H1_1,self.H1_2], origin='lower')
        plt.xlabel("H1 (\u03BCm)")
        plt.ylabel("LC (\u03BCm)")
        divider = make_axes_locatable(ax)
        cax = divider.new_vertical(size="5%", pad=0.4, title="Eccentricity")
        fig.add_axes(cax)
        fig.colorbar(im, cax=cax, orientation="horizontal")
####################################################################################################### mshu x y labela wash karw
        # plt.savefig('colorbar_positioning_03.png', format='png', bbox_inches='tight')
        plt.hot()
        plt.show()
        plt.close()
    def printt(self):
        print(self.filename)
    def plot_simulate_dn(self):
        '''This method plots modal birefringence of the simulated data in the file called by the method fileDialog'''
        self.create_Entrys()
        fig, ax = plt.subplots(figsize=(8, 6))
        #plt.title("Modal bireferengence as a function of channel dimensions")
        im = plt.imshow(self.dr,  extent=[self.LC_1, self.LC_2, self.H1_1,self.H1_2], origin='lower')
        plt.xlabel("H1 (\u03BCm)")
        plt.ylabel("LC (\u03BCm)")
        divider = make_axes_locatable(ax)
        cax = divider.new_vertical(size="5%", pad=0.4, title="Modal bireferengence")
        fig.add_axes(cax)
        fig.colorbar(im, cax=cax,orientation="horizontal")

        # plt.savefig('colorbar_positioning_03.png', format='png', bbox_inches='tight')
        plt.hot()
        plt.show()
        plt.close()
    def close_all(self):
        '''This method closes everything it just has self.master.destroy() and plt.close('all')'''
        self.master.destroy()
        plt.close('all')
    def Solve(self):
        self.create_Entrys()
        self.g = Channel(wl=self.WL, pas=0.1, \
                       nc=self.nC_1, nc_2=self.nC_2, nc_3=self.nC_3, Hc=0.8, \
                       LC=self.TChan, LB=5., LB_R=5., \
                       n1B=self.nL1_1, n1C=self.nL1_2, n1D=self.nL1_3, H1=self.TL_1, \
                       n2B=self.nL2_1, n2C=self.nL2_2, n2D=self.nL2_3, H2=self.TL_2, \
                       nsub=self.nSub_1, nsub_2=self.nSub_2, nsub_3=self.nSub_3, Hsub=2.,
                       OR=self.OPR)
        self.g.Calcule(Nmodes=self.Nmodes_value)
        #self.g.Intensite()
        plt.show()
    def NumOf_guided_modes(self): # This will print angles which can give coupled light.
        '''This method finds any differaction mode of a grating that matches
         with the guided modes of the structure and returns the angle'''
        self.create_Entrys()
        print("Coupling neff is calculated with equation\n n_eff=n_top.sin(\u03B8) + (m.(\u03BB/\u039B) \n where "
              "n_eff is effective refractive index, n_top is refractive index of the cover\n"
              "\u03B8 is coupling angle,m is diffraction mode, \u03BB is wavelength, and \u039B is grating period")
        self.g = Channel(wl=self.WL, pas=0.1, \
                       nc=self.nC_1, nc_2=self.nC_2, nc_3=self.nC_3, Hc=0.8, \
                       LC=self.TChan, LB=5., LB_R=5., \
                       n1B=self.nL1_1, n1C=self.nL1_2, n1D=self.nL1_3, H1=self.TL_1, \
                       n2B=self.nL2_1, n2C=self.nL2_2, n2D=self.nL2_3, H2=self.TL_2, \
                       nsub=self.nSub_1, nsub_2=self.nSub_2, nsub_3=self.nSub_3, Hsub=2.,
                       OR=self.OPR)
        if self.g.NmaxPlan()==0: # blow H2==0.9 NmaxPlan is 0 and it does not give good value.
            messagebox.showinfo("Hey", 'It looks like NmaxPlan is {}'.format(self.g.NmaxPlan()))
            pass
        else:
            print("neff of the guided modes of the structure are: {}".format(self.g.NumOf_guidedModes(Nmodes=self.Nmodes_value)[1:]))
            just_Modes = ((self.g.NumOf_guidedModes(Nmodes=self.Nmodes_value)[1:]))
            theta = np.arange(0, 90, 0.001)# array of angles in degrees
            All_n_effectives_are = np.around(self.Nmode_effective(theta), decimals=4)#around function round to 4 decimals and neglect the rest.
            guided_modes_n_effectives_are_array = np.where((np.where(
                All_n_effectives_are <= (np.amax(just_Modes)), All_n_effectives_are, 0)) >= (np.amin(just_Modes)),
                                                           All_n_effectives_are, 0)
            guided_modes_n_effectives_are_index = np.where(guided_modes_n_effectives_are_array > 0)
            guided_modes_n_effectives_are = All_n_effectives_are[guided_modes_n_effectives_are_index]
            guided_modes_n_effectives_theta_are = np.around(theta[guided_modes_n_effectives_are_index], decimals=1)
            s = guided_modes_n_effectives_are.size
            if s>0:
                for i in range(0, s):
                    print("There is one guided mode at angle {} that has neff = {}".format(
                    (guided_modes_n_effectives_theta_are[i])
                    , guided_modes_n_effectives_are[i]))
            else:
                messagebox.showinfo("Sorry",'The is no coupling angle that match with guided modes')
    def Nmode_effective(self,one_theta=5.2): # to find n_e that is effective refractive index of the grating
        '''This method returns guided modes of the structure (not leaky modes)'''
        self.create_Entrys()
        #print(self.grating_period,self.diffraction_mode)
        n_top=np.average(np.array([float(self.nC_1),float(self.nC_2),float(self.nC_3)]))
        n_e = ((n_top * (np.sin(np.rad2deg(one_theta)))) + (float(self.diffraction_mode) * (float(self.WL)/float(self.grating_period))))#n_top= np.average(np.arange(float(self.nC_1),float(self.nC_2),float(self.nC_3)))
        return n_e

    def Simulate(self):
        '''This method finds modal birefringence of the structure
         and simulates this for the defined channel dimentions'''
        self.create_Entrys()
        if messagebox.askyesno("Warning! ", 'This may take several minutes\n Are you sure you want to continue?'):
            messagebox.showinfo("Simulation",'You can see the see the remaining time in the run window')
            self.g = Channel(wl=0.64, pas=0.1, \
                       nc=1, nc_2=1, nc_3=1, Hc=0.8, \
                       LC=3.1, LB=5., LB_R=5., \
                       n1B=1, n1C=1.62, n1D=1, H1=2.3, \
                       n2B=1.62, n2C=1.62, n2D=1.62, H2=2, \
                       nsub=1.61, nsub_2=1.61, nsub_3=1.61, Hsub=2.,
                       OR=4)
            self.g.H2 = self.VaryH1LC_H2_Value
            self.g.VariaXY(X0=self.LC_1, X1=self.LC_2, dX=0.1, Y0=self.H1_1, Y1=self.H1_2, dY=0.1, fich=str(self.filename))
        else:
            messagebox.showinfo("Canceled",'You have successfully canceled the simulation')
    def close_plots(self):
        plt.close('all')
if __name__ == '__main__':
    root = Tk()  # we car write any name instead of root
    root.title('Multilayer channel waveguide mode solver')
    root.iconbitmap('rib_waveguide.ico')
    root.geometry("600x660")
    app = Mul_Ch_Wav_Mod_Sol(root)
    app.mainloop()