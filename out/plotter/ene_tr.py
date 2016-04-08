#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pylab

origin = 'lower'
#origin = 'upper'

#data = np.genfromtxt("3D, t=8.000-(forward).dat", delimiter=',')
#data = np.genfromtxt("radial-N=10.dat", delimiter=',')
data = np.genfromtxt("out/data/Evo.dat", delimiter=',')
d2 = np.genfromtxt("out/data/e(0.07), tf=10.dat", delimiter=',')
#f, (ax1,ax2) = plt.subplots(2, sharex=True)
plt.figure(figsize=(4,6))
plt.subplots_adjust(hspace=.4)

delta = 0.025

#x = y = np.arange(-3.0, 3.01, delta)
xt = data[:100,1]
tt = data[::100,0]

#print(xt)

X, Y = np.meshgrid(xt, tt)
#Z = 10 * (Z1 - Z2)
Z = data[:,2].reshape(len(xt),len(tt))
#Z = griddata( x,y,z,X,Y)

nr, nc = Z.shape

ene = d2[:,1]
tim = d2[:,0]

plt.subplot(211)
plt.loglog(tim, ene, 'or')
plt.loglog(tim, (tim/.6)**(-1./3.), '-b')
plt.loglog(tim, 10.*(tim/.6)**(-4./3.), ':b')
#ax1.set_yscale('log')
plt.ylim([0.02,1.3])
plt.xlim([0.5,11.6])
#ax1.set_xlim(0.6,10.6)

plt.subplot(212)
CS = plt.contourf(X, Y, Z, 20,
                  #[-1, -0.1, 0, 0.1],
                  #alpha=0.8,
                  cmap=plt.cm.PuBu,
                  origin=origin)

# Note that in the following, we explicitly pass in a subset of
# the contour levels used for the filled contours.  Alternatively,
# We could pass in additional levels to provide extra resolution,
# or leave out the levels kwarg to use all of the original levels.
CS2 = plt.contour(CS, levels=CS.levels[::2],
                  colors='r',
                  origin=origin,
                  hold='on')

plt.title(r'radial expansion')
plt.clabel(CS2, inline=1,fontsize=14)
plt.xlabel(r'$r$')
plt.ylabel(r'$\tau$')
plt.xlim(-0,10);
plt.xlim(-0,10);

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('energy density')
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)
#plt.show()
sopt = {'dpi':300}
plt.savefig("test1.pdf",**sopt)


