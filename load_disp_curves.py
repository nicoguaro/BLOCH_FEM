# -*- coding: utf-8 -*-
"""
Reshape and plot dispersion curves.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
 
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16 

fname = "homog1.fre"
# This example uses a 41 x 41 grid
nvals = 6  # Number of eigenvalues retrieved
nkx = 41  # Number of wavenumbers in x
nky = 41  # Number of wavenumbers in y
nks = nkx*nky  
data = np.loadtxt(fname)
data = np.reshape(data, (nkx, nky, nvals + 2))

kx = data[:,:,0]
ky = data[:,:,1]
data[:,:,2:-1] = np.sqrt(np.abs(data[:,:,2:-1]))
w0 = data[:,:,2]
w1 = data[:,:,3]

plt.contourf(2*kx/np.pi, 2*ky/np.pi, w0, cmap='hot')
#plt.pcolor(2*kx/np.pi, 2*ky/np.pi, w0, cmap='hot')
#plt.contourf(2*kx/np.pi, 2*ky/np.pi, w1, cmap='hot')
#plt.pcolor(2*kx, 2*ky/np.pi, w1, cmap='hot')
plt.colorbar()
plt.xlabel(r"$2d_x k_x/\pi$")
plt.ylabel(r"$2d_y k_y/\pi$")
plt.axis('image')
#plt.savefig('disp_contour.png', dpi=600)
plt.show()