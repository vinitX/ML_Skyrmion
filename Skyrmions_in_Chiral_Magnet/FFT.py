import numpy as np
import time
import matplotlib.pyplot as plt

L=24
m=46
n=15

D=np.zeros((L,L)) #Euclidean Distance Matrix 
for i in range(L):
    for j in range(L):
        D[i,j]=(min(i,L-i)**2+min(j,L-j)**2)**0.5
        
t=time.time()
S=np.loadtxt('H.csv',delimiter=',')
print(np.shape(S))
print(time.time()-t)
S=np.reshape(S,(m,n,L,L,3))

A=25
B=10

plt.imshow(D)
plt.show()

X=S[A,B,:,:,2]-np.mean(S[A,B,:,:,2])
plt.imshow(X)
plt.show()
W=np.abs(np.fft.fft2(X))
plt.imshow(W)
plt.show()

lambd=np.sum(W*D)/np.sum(W)
print(lambd)

print(np.max(W)/np.sum(W))
print(np.sort(W,axis=None)[:-5:-1]/np.sum(W))

sigma=(np.sum(W*(D-lambd)**2)/np.sum(W))**0.5
print(sigma)