import numpy as np
from tkinter import messagebox
from Channels_multilayer import Channel
g=Channel(wl=0.64,pas=0.1,nc=1., nc_2=1., nc_3=1.,Hc=0.8,\
            LB=5, LC=3.2,LB_R=5,n1B=1.,n1C=1.620,n1D=1.,H1=2.3,\
            n2B=1.62,n2C=1.62,n2D=1.62,H2=0.9,\
            nsub=1.61,nsub_2=1.61,nsub_3=1.61,Hsub=2.,OR=4)
print ("All units are in \u03BCm")
m = 1
lam = 0.64
theta=np.arange(0,90,0.01)
#theta = theta.tolist()
# print (theta)
# print((g.NumOf_guidedModes(Nmodes=10)))
# print((g.NumOf_guidedModes(Nmodes=10)[1:]))
just_Modes=((g.NumOf_guidedModes(Nmodes=5)[1:]))
# print(np.amax(g.NumOf_guidedModes(Nmodes=10)[1:]))
# print(np.amin(g.NumOf_guidedModes(Nmodes=10)[1:]))
# print ("The angles are is "+str(theta))
n_top = 1.0# float(input("Please inter refractive index of cover layer: "))
n_guiding = 1.62# float(input("Please inter refractive index of guding layer: "))
n_bottom = 1.61# float(input("Please inter refractive index of guding layer: "))
grating_period = 1000/1200#np.arange(0.2,2,0.01)#float(input("Please inter grating period "))
def Nmode_effective (one_theta):
    n_e = ((n_top *(np.sin(np.rad2deg(one_theta))))+(m*(lam/grating_period)))
    return n_e


if g.NmaxPlan() == 0:  # blow H2==0.9 NmaxPlan is 0 and it does not give good value.
    messagebox.showinfo("Hey", 'It looks like NmaxPlan is {}'.format(g.NmaxPlan()))
    pass
else:
    just_Modes = ((g.NumOf_guidedModes(Nmodes=10)[1:]))
    theta = np.arange(0, 90, 0.01)  # array of angles in degrees
    All_n_effectives_are = np.around(Nmode_effective(theta),
                                     decimals=4)  # around function round to 4 decimals and neglect the rest.
    guided_modes_n_effectives_are_array = np.where((np.where(
        All_n_effectives_are <= (np.amax(just_Modes)), All_n_effectives_are, 0)) >= (np.amin(just_Modes)),
                                                   All_n_effectives_are, 0)
    guided_modes_n_effectives_are_index = np.where(guided_modes_n_effectives_are_array > 0)
    guided_modes_n_effectives_are = np.around(All_n_effectives_are[guided_modes_n_effectives_are_index],decimals=3)
    guided_modes_n_effectives_theta_are = np.around(theta[guided_modes_n_effectives_are_index], decimals=1)
    # print(g.NumOf_guidedModes(Nmodes=10)[1:])
    # for i in g.NumOf_guidedModes(Nmodes=10)[1:]:
    #     guided_modes_n_effectives_are_true=np.where(guided_modes_n_effectives_are==i,guided_modes_n_effectives_are,0)
    guided_modes_n_effectives_are_true=guided_modes_n_effectives_are
    s = guided_modes_n_effectives_are_true.size
    if s > 0:
        for i in range(0, s):
            print("There is one guided mode at angle {} that has neff = {}".format(
                (guided_modes_n_effectives_theta_are[i])
                , guided_modes_n_effectives_are_true[i]))
    else:
        messagebox.showinfo("Sorry", 'The is no coupling angle that match with guided modes')