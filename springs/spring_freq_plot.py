from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("spring.fre")
print(data.shape)
data = data.reshape((41, 41, 8))
data[:,:,2:] = np.sqrt(data[:,:,2:])

min_freq = np.min(data[:,:,2:])
max_freq = np.max(data[:,:,2:])
levels = np.linspace(min_freq, max_freq, 20)
kx = data[:,:,0]
ky = data[:,:,1]
plt.figure(figsize=(15,10))
for cont in range(2, 8):
    plt.subplot(2, 3, cont - 1)
    plt.contourf(kx, ky, data[:,:,cont], levels, cmap="YlGnBu_r")
    plt.contour(kx, ky, data[:,:,cont],levels, colors="k")
    plt.axis("image")
plt.savefig("spring_dispersion.svg")
plt.show()
