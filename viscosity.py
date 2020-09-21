import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix' # math fontの設定
plt.rcParams['font.size'] = 25

# input files
t1, h1, velo1, m1 = np.loadtxt("./01/out.dat", skiprows=1, unpack=True)
h1, rho1, temp1, mu1 = np.loadtxt("./01/atmos_model.dat", skiprows=1, unpack=True)

h1 /= 1000

mu_sutherland = np.array(len(h1))
mu_power_law = np.array(len(h1))

#for i in range(len(h1)):
#    mu_sutherland[i] = 1.458e-6*(temp1[i]**1.5)/(temp1[i] + 110.4) # Sutherland
#    mu_power_law[i] = 1.716e-5*((temp1[i]/273.0)**(2./3.)) # power law by Maxwell, Rayleigh (for rarefied gas)
mu_sutherland = 1.458e-6*(temp1**1.5)/(temp1 + 110.4) # Sutherland
mu_power_law = 1.716e-5*((temp1/273.0)**(2./3.)) # power law by Maxwell, Rayleigh (for rarefied gas)

plt.figure(figsize=(12.8,7.8))
plt.plot(mu_sutherland,h1,color='r',linewidth=4,label='Sutherland Formula')
plt.plot(mu_power_law,h1,color='#3180b6',linewidth=4,linestyle='dashed',label='Power Law')
plt.xlabel('Air Viscosity [Pa$\cdot$s]')
plt.ylabel('Altitude [km]')
plt.legend(frameon=False)
plt.ylim(0,400)
plt.savefig('./figs/viscosity.pdf',format='pdf',dpi=1200,transparent=True,bbox_inches='tight',pad_inches=0)
plt.show()
