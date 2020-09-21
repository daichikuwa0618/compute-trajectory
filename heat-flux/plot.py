import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 25

t0, h0, velo0, m0, dis0, th0 = np.loadtxt("./100%/out.dat", skiprows=1, unpack=True)
t1, h1, velo1, m1, dis1, th1 = np.loadtxt("./86.81%/out.dat", skiprows=1, unpack=True)
t2, h2, velo2, m2, dis2, th2 = np.loadtxt("./49.62%/out.dat", skiprows=1, unpack=True)
# input files
#t, h, velo, m, dis, th = np.loadtxt("./100%/out.dat", skiprows=1, unpack=True)
#h, rho, temp, mu = np.loadtxt("./atmos_model.dat", skiprows=1, unpack=True)
#h, mach, re_num, Cd, Ch, kn = np.loadtxt("./coeff.dat", skiprows=1, unpack=True)
#h, bright, tau, magni, lumi = np.loadtxt("./brightness.dat", skiprows=1, unpack=True)
#
#t1, h1, velo1, m1, dis1, th1 = np.loadtxt("./86.81%/out.dat", skiprows=1, unpack=True)
#h1, rho, temp, mu = np.loadtxt("./atmos_model.dat", skiprows=1, unpack=True)
#h1, mach, re_num, Cd, Ch, kn = np.loadtxt("./coeff.dat", skiprows=1, unpack=True)
#h1, bright, tau, magni, lumi = np.loadtxt("./brightness.dat", skiprows=1, unpack=True)
#
#t2, h2, velo2, m2, dis2, th2 = np.loadtxt("./49.62%/out.dat", skiprows=1, unpack=True)
#h2, rho, temp, mu = np.loadtxt("./atmos_model.dat", skiprows=1, unpack=True)
#h2, mach, re_num, Cd, Ch, kn = np.loadtxt("./coeff.dat", skiprows=1, unpack=True)
#h2, bright, tau, magni, lumi = np.loadtxt("./brightness.dat", skiprows=1, unpack=True)

# [kg] -> [g]
m0  *= 1000
m1  *= 1000
m2  *= 1000

# [m] -> [km]
h0 /= 1000.0
h1 /= 1000.0
h2 /= 1000.0

# single figure
plt.figure(figsize=(12.8,7.8),dpi=100)
#plt.style.use('dark_background')
plt.style.use('dark_background')
plt.plot(velo0,h0,label='100%',linewidth=4)
plt.plot(velo1,h1,label='87%',linewidth=4)
plt.plot(velo2,h2,label='49%',linewidth=4)
#plt.plot(kn5,h5,linewidth=2,color='red')
plt.xlabel('Velocity [m/s]')
#plt.xscale('log')
#plt.xlim(1e-5,1e10)
plt.ylabel('Altitude [km]')
plt.ylim(0,400)
plt.legend(loc='upper left')
#plt.legend()
#plt.savefig('./figs/velo.png')
#plt.savefig('./figs/velo.svg', format='svg', dpi=1200, transparent=True, bbox_inches='tight', pad_inches=0) #format = svg
plt.savefig('./velo.pdf', format='pdf', dpi=1200, transparent=True, bbox_inches='tight', pad_inches=0) #format = pdf

plt.show()
