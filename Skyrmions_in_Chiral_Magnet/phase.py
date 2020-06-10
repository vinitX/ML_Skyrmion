import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

L=24

t=time.time()
chi=np.loadtxt('D_chi.csv',delimiter=',')
mag=np.loadtxt('D_mag.csv',delimiter=',')
tmag=np.loadtxt('D_tmag.csv',delimiter=',')
ham=np.loadtxt('D_ham.csv',delimiter=',')
peak=np.loadtxt('D_peak.csv',delimiter=',')

chi=gaussian_filter(np.reshape(chi,(30,5)),2)
mag=gaussian_filter(np.reshape(mag,(30,5)),2)
tmag=gaussian_filter(np.reshape(tmag,(30,5)),2)
ham=gaussian_filter(np.reshape(ham,(30,5)),2)
peak=gaussian_filter(np.reshape(peak,(30,5)),2)

chi=chi/np.max(np.abs(chi))
mag=mag/np.max(mag)
tmag=tmag/np.max(tmag)
peak=peak/np.max(peak)

print(time.time()-t)

A=np.zeros((30,5))
for i in range(30):
    for j in range(5):
        if mag[i,j]>0.7: #Ferromagnet
            A[i,j]=3
        elif tmag[i,j]>0.3 and mag[i,j]<0.2:  #Tilted Ferromagnet
            A[i,j]=2
        elif np.abs(chi[i,j])>0.5: #Skyrmion
            A[i,j]=1
        elif mag[i,j]<0.2 and peak[i,j]>0.4: #Spiral
            A[i,j]=-3
        
        else:
            A[i,j]=-0.5

plt.imshow(A)
#np.savetxt('phase.csv',A,delimiter=',')
