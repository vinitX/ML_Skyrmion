import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

L=24

t=time.time()
chi=np.reshape(np.loadtxt('chi.csv',delimiter=','),(46,15))
mag=np.reshape(np.loadtxt('mag.csv',delimiter=','),(46,15))
tmag=np.reshape(np.loadtxt('tmag.csv',delimiter=','),(46,15))
ham=np.reshape(np.loadtxt('ham.csv',delimiter=','),(46,15))
peak=np.reshape(np.loadtxt('peak.csv',delimiter=','),(46,15))

chi2=np.reshape(np.loadtxt('chi2.csv',delimiter=','),(46,15))
mag2=np.reshape(np.loadtxt('mag2.csv',delimiter=','),(46,15))
tmag2=np.reshape(np.loadtxt('tmag2.csv',delimiter=','),(46,15))
ham2=np.reshape(np.loadtxt('ham2.csv',delimiter=','),(46,15))
peak2=np.reshape(np.loadtxt('peak2.csv',delimiter=','),(46,15))
print(time.time()-t)

a=np.abs(gaussian_filter(chi[:,8],2))
b=np.abs(gaussian_filter(chi2[:,8],2))
plt.plot(a/np.max(a))
plt.plot(b/np.max(b))
plt.show()

a=np.abs(gaussian_filter(mag[:,8],2))
b=np.abs(gaussian_filter(mag2[:,8],2))
plt.plot(a/np.max(a))
plt.plot(b/np.max(b))
plt.show()

a=np.abs(gaussian_filter(tmag[:,8],2))
b=np.abs(gaussian_filter(tmag2[:,8],2))
plt.plot(a/np.max(a))
plt.plot(b/np.max(b))
plt.show()

a=np.abs(gaussian_filter(ham[:,8],2))
b=np.abs(gaussian_filter(ham2[:,8],2))
plt.plot(a/np.max(a))
plt.plot(b/np.max(b))
plt.show()

a=np.abs(gaussian_filter(peak[:,8],2))
b=np.abs(gaussian_filter(peak2[:,8],2))
plt.plot(a/np.max(a))
plt.plot(b/np.max(b))
plt.show()
