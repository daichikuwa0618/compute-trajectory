import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 18

# input files
t1, h1, velo1, m1 = np.loadtxt("./01/out.dat", skiprows=1, unpack=True)
h1, rho1, temp1, mu1 = np.loadtxt("./01/atmos_model.dat", skiprows=1, unpack=True)
h1, mach1, re_num1, Cd1, Ch1 = np.loadtxt("./01/coeff.dat", skiprows=1, unpack=True)
h1, bright1 = np.loadtxt("./01/brightness.dat", skiprows=1, unpack=True)

t2, h2, velo2, m2 = np.loadtxt("./02/out.dat", skiprows=1, unpack=True)
h2, rho2, temp2, mu2 = np.loadtxt("./02/atmos_model.dat", skiprows=1, unpack=True)
h2, mach2, re_num2, Cd2, Ch2 = np.loadtxt("./02/coeff.dat", skiprows=1, unpack=True)
h2, bright2 = np.loadtxt("./02/brightness.dat", skiprows=1, unpack=True)

t3, h3, velo3, m3 = np.loadtxt("./03/out.dat", skiprows=1, unpack=True)
h3, rho3, temp3, mu3 = np.loadtxt("./03/atmos_model.dat", skiprows=1, unpack=True)
h3, mach3, re_num3, Cd3, Ch3 = np.loadtxt("./03/coeff.dat", skiprows=1, unpack=True)
h3, bright3 = np.loadtxt("./03/brightness.dat", skiprows=1, unpack=True)

t4, h4, velo4, m4 = np.loadtxt("./04/out.dat", skiprows=1, unpack=True)
h4, rho4, temp4, mu4 = np.loadtxt("./04/atmos_model.dat", skiprows=1, unpack=True)
h4, mach4, re_num4, Cd4, Ch4 = np.loadtxt("./04/coeff.dat", skiprows=1, unpack=True)
h4, bright4 = np.loadtxt("./04/brightness.dat", skiprows=1, unpack=True)

t5, h5, velo5, m5 = np.loadtxt("./05/out.dat", skiprows=1, unpack=True)
h5, rho5, temp5, mu5 = np.loadtxt("./05/atmos_model.dat", skiprows=1, unpack=True)
h5, mach5, re_num5, Cd5, Ch5 = np.loadtxt("./05/coeff.dat", skiprows=1, unpack=True)
h5, bright5 = np.loadtxt("./05/brightness.dat", skiprows=1, unpack=True)

# [kg] -> [g]
m1 *= 1000
m2 *= 1000
m3 *= 1000
m4 *= 1000
m5 *= 1000

# single figure
plt.figure(figsize=(12.8,7.8),dpi=100)
plt.plot(bright1,h1,linewidth=2,label='A')
plt.plot(bright2,h2,linewidth=2,label='B')
plt.plot(bright3,h3,linewidth=2,label='C')
plt.plot(bright4,h4,linewidth=2,label='D')
plt.plot(bright5,h5,linewidth=2,label='E')
plt.xlabel('Brightness [W]')
#plt.xscale('log')
plt.xlim(0,3500)
plt.ylabel('Altitude [m]')
plt.ylim(0,400000)
#plt.legend(loc='upper left')
plt.savefig('./figs/bright.png')

# some figures
#fig1, (axU, axL) = plt.subplots(nrows=2, figsize=(8,10),sharey=True)
#
#axU.plot(rho1,h1,linewidth=2,label='A')
#axU.plot(rho2,h2,linewidth=2,label='B')
#axU.plot(rho3,h3,linewidth=2,label='C')
#axU.plot(rho4,h4,linewidth=2,label='D')
#axU.plot(rho5,h5,linewidth=2,label='E')
#axU.set_title('Air Density')
#axU.set_xlabel('Density [kg/m^3]')
#axU.set_xscale('log')
#axU.set_xlim(1e-15,1e5)
#axU.set_ylabel('Altitude [m]')
#axU.set_ylim(0,400000)
#axU.legend(loc='upper right')
#
#axL.plot(temp1,h1,linewidth=2,label='A')
#axL.plot(temp2,h2,linewidth=2,label='B')
#axL.plot(temp3,h3,linewidth=2,label='C')
#axL.plot(temp4,h4,linewidth=2,label='D')
#axL.plot(temp5,h5,linewidth=2,label='E')
#axL.set_title('Air Temperature')
#axL.set_xlabel('Temperature [K]')
#axL.set_xlim(0,1000)
#axL.set_ylabel('Altitude [m]')
#axL.legend(loc='upper right')

plt.show()
