import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 25

# input files
t, h, velo, m, dis, th = np.loadtxt("./out.dat", skiprows=1, unpack=True)
h, rho, temp, mu = np.loadtxt("./atmos_model.dat", skiprows=1, unpack=True)
h, mach, re_num, Cd, Ch, kn = np.loadtxt("./coeff.dat", skiprows=1, unpack=True)
h, bright, tau, magni, lumi = np.loadtxt("./brightness.dat", skiprows=1, unpack=True)


# [kg] -> [g]
m  *= 1000

# [m] -> [km]
h  /= 1000.0

# single figure
plt.figure(figsize=(12.8,7.8),dpi=100)
plt.plot(velo,h,linewidth=4,color='red')
#plt.plot(kn5,h5,linewidth=2,color='red')
plt.xlabel('Velocity [m/s]')
#plt.xscale('log')
plt.xlim(0,8000)
plt.ylabel('Altitude [km]')
plt.ylim(0,400)
#plt.legend(loc='upper right')
#plt.legend()
#plt.savefig('./figs/velo.png')
#plt.savefig('./figs/velo.svg', format='svg', dpi=1200, transparent=True, bbox_inches='tight', pad_inches=0) #format = svg
plt.savefig('.velo.pdf', format='pdf', dpi=1200, transparent=True, bbox_inches='tight', pad_inches=0) #format = pdf

plt.show()
