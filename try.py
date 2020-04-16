import numpy as np
diffraction_mode=2
WL=0.64
grating_period=1.1
def Nmode_effective(n_eff):  # to find n_e that is effective refractive index of the grating
    '''This method returns guided modes of the structure (not leaky modes)'''
    #create_Entrys()
    # print(self.grating_period,self.diffraction_mode)
    n_top = 1# np.average(np.array([nC_1, nC_2, nC_3]))
    if ((n_eff - (diffraction_mode * (WL / grating_period))) / (n_top)) > 1:
        theta = 0
    elif ((n_eff - (diffraction_mode * (WL / grating_period))) / (n_top)) < 0:
        theta = 0
    else:
        theta = np.rad2deg(np.arcsin((n_eff - (diffraction_mode * (WL / grating_period))) / (n_top)))
    # n_e = ((n_top * (np.sin(np.rad2deg(one_theta)))) + (float(self.diffraction_mode) * (float(self.WL)/float(self.grating_period))))#n_top= np.average(np.arange(float(self.nC_1),float(self.nC_2),float(self.nC_3)))
    return theta
print(Nmode_effective(1.619))