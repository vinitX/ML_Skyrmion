import numpy as np
import time
import matplotlib.pyplot as plt

L=24

t=time.time()
S=np.loadtxt('H.csv',delimiter=',')
S=np.reshape(S,(46,15,L,L,3))
F=np.loadtxt('H_fft.csv',delimiter=',')
F=np.reshape(F,(46,15,24,24))

for K in range(0,45):
    for H in range(14,15):        
        filename=str(K)+"_"+str(H)+".png"
        plt.imshow(S[K,H,:,:,2])
        plt.show()
        #plt.savefig(filename)
        filename=str(K)+"_"+str(H)+"f.png"
        plt.imshow(F[K,H,:,:])
        #plt.savefig(filename)
        
print(time.time()-t)