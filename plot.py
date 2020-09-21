import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix' # math fontの設定
plt.rcParams['font.size'] = 25

# input files
t1, h1, velo1, m1, d1, th1 = np.loadtxt("./cd-with-coeff/out.dat", skiprows=1, unpack=True)
h1, rho1, temp1, mu1 = np.loadtxt("./cd-with-coeff/atmos_model.dat", skiprows=1, unpack=True)
h1, mach1, re_num1, Cd1, Ch1, kn1 = np.loadtxt("./cd-with-coeff/coeff.dat", skiprows=1, unpack=True)
h1, bright1, tau1, magni1, lumi1 = np.loadtxt("./cd-with-coeff/brightness.dat", skiprows=1, unpack=True)
h1, heat1 = np.loadtxt('./cd-with-coeff/heat-flux.dat', skiprows=1, unpack=True)
# h1, bright1 = np.loadtxt("./cd-with-coeff/brightness.dat", skiprows=1, unpack=True)
#
t2, h2, velo2, m2, d2, th2 = np.loadtxt("./cd-without-coeff/out.dat", skiprows=1, unpack=True)
h2, rho2, temp2, mu2 = np.loadtxt("./cd-without-coeff/atmos_model.dat", skiprows=1, unpack=True)
h2, mach2, re_num2, Cd2, Ch2, kn2 = np.loadtxt("./cd-without-coeff/coeff.dat", skiprows=1, unpack=True)
h2, bright2, tau2, magni2, lumi2 = np.loadtxt("./cd-without-coeff/brightness.dat", skiprows=1, unpack=True)
h2, heat2 = np.loadtxt('./cd-without-coeff/heat-flux.dat', skiprows=1, unpack=True)
# h2, bright2 = np.loadtxt("./cd-without-coeff/brightness.dat", skiprows=1, unpack=True)
#
#t3, h3, velo3, m3 = np.loadtxt("./03/out.dat", skiprows=1, unpack=True)
#h3, rho3, temp3, mu3 = np.loadtxt("./03/atmos_model.dat", skiprows=1, unpack=True)
#h3, mach3, re_num3, Cd3, Ch3 = np.loadtxt("./03/coeff.dat", skiprows=1, unpack=True)
##h3, bright3, tau3, magni3, lumi3 = np.loadtxt("./03/brightness.dat", skiprows=1, unpack=True)
##h3, bright3 = np.loadtxt("./03/brightness.dat", skiprows=1, unpack=True)
#
#t4, h4, velo4, m4 = np.loadtxt("./04/out.dat", skiprows=1, unpack=True)
#h4, rho4, temp4, mu4 = np.loadtxt("./04/atmos_model.dat", skiprows=1, unpack=True)
#h4, mach4, re_num4, Cd4, Ch4 = np.loadtxt("./04/coeff.dat", skiprows=1, unpack=True)
#h4, bright4, tau4, magni4, lumi4 = np.loadtxt("./04/brightness.dat", skiprows=1, unpack=True)
##h4, bright4 = np.loadtxt("./04/brightness.dat", skiprows=1, unpack=True)
#
#t5, h5, velo5, m5 = np.loadtxt("./05/out.dat", skiprows=1, unpack=True)
#h5, rho5, temp5, mu5 = np.loadtxt("./05/atmos_model.dat", skiprows=1, unpack=True)
#h5, mach5, re_num5, Cd5, Ch5, kn5 = np.loadtxt("./05/coeff.dat", skiprows=1, unpack=True)
#h5, bright5, tau5, magni5, lumi5 = np.loadtxt("./05/brightness.dat", skiprows=1, unpack=True)
##h5, bright5 = np.loadtxt("./05/brightness.dat", skiprows=1, unpack=True)

t, h, velo, m, dis, th = np.loadtxt("./out.dat", skiprows=1, unpack=True)
h, rho, temp, mu = np.loadtxt("./atmos_model.dat", skiprows=1, unpack=True)
h, mach, re_num, Cd, Ch, kn = np.loadtxt("./coeff.dat", skiprows=1, unpack=True)
h, bright, tau, magni, lumi = np.loadtxt("./brightness.dat", skiprows=1, unpack=True)
h, heat = np.loadtxt('./heat-flux.dat', skiprows=1, unpack=True)


# [kg] -> [g]
m1 *= 1000
m2 *= 1000
#m3 *= 1000
#m4 *= 1000
#m5 *= 1000
m  *= 1000

# [m] -> [km]
h /= 1000.0
h1 /= 1000.0
h2 /= 1000.0
#h3 /= 1000.0
#h4 /= 1000.0
#h5 /= 1000.0
# h  /= 1000.0
heat /= 1e6
heat1 /= 1e6
heat2 /= 1e6

# single figure
#plt.style.use('dark_background')
plt.figure(figsize=(12.8,7.8),dpi=100)
#plt.plot(temp1,h1,linewidth=4,color='red')
#plt.plot(velo2,h2,linewidth=2,label='B')
#plt.plot(velo3,h3,linewidth=2,label='C')
#plt.plot(velo4,h4,linewidth=2,label='D')
#plt.plot(velo5,h5,linewidth=2,label='E')
#plt.plot(kn5,h5,linewidth=5)
plt.plot(velo,h,linewidth=4,color='r',label='w/ coeff')
plt.plot(velo2,h2,linewidth=4,linestyle='dashed',color='#3180b6',label='w/o coeff')
plt.xlabel('Velocity [m/s]')
# plt.xscale('log')
# plt.xlim(1.5,2.5)
# plt.xticks(np.arange(1.5,2.6,0.1))
plt.ylabel('Altitude [km]')
plt.ylim(0,400)
# plt.hlines(500,4,5,color='#3180b6',linestyle='dashed',linewidth=4,label='Kimura (2018)') # not displayed on graph, but on legend.
#plt.legend(loc='upper right')
plt.legend(frameon=False)
#plt.savefig('./figs/velo.png')
#plt.savefig('./figs/velo.svg', format='svg', dpi=1200, transparent=True, bbox_inches='tight', pad_inches=0) #format = svg
plt.savefig('./figs/heat-coeff/velo.pdf', format='pdf', dpi=1200, transparent=True, bbox_inches='tight', pad_inches=0) #format = pdf
#fig, ax1 = plt.subplots()
#ax2 = ax1.twiny()
#ax1.set_xlabel('Velocity [m/s]')
#ax2.set_xlabel('Mass [g]')
#ax1.set_xlim(0,8000)
#ax1.set_ylim(0,400)
#ax2.set_xlim(0,3)
#ax1.set_ylabel('Altitude [km]')
#ax1.plot(velo5,h5,linewidth=3,color="red",label="Velocity")
#ax2.plot(m5,h5,linewidth=3,label="Mass")
#ax1.legend(bbox_to_anchor=(0,1),loc='upper left')
#ax2.legend(bbox_to_anchor=(0,0.9),loc='upper left')
#ax1.grid()

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

# plt.show()
