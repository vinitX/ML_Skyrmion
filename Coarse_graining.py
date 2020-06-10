import numpy as np   
import time

L=96

t=time.time()
#test = np.array(list(np.loadtxt("H4.csv", delimiter=","))).astype("float32")
print(time.time()-t)
t=time.time()

test=np.reshape(test,(400,L,L,3))

res=np.zeros((400,24,24,3))

for i in range(24):
    for j in range(24):
        res[:,i,j,:]=np.sum(np.sum(test[:,(L//24)*i:(L//24)*(i+1),(L//24)*j:(L//24)*(j+1),:],axis=2),axis=1)
        norm=np.sum(res[:,i,j,:]**2,axis=1)**0.5
        if norm.any()<0.00001: norm+=0.00001
        res[:,i,j,0]=res[:,i,j,0]/norm
        res[:,i,j,1]=res[:,i,j,1]/norm
        res[:,i,j,2]=res[:,i,j,2]/norm
        
np.savetxt('H4_CG.csv',np.reshape(res,(-1,24*24*3)),delimiter=',')

print(time.time()-t)
