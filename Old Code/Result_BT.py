import numpy as np
import matplotlib.pyplot as plt
import time

L=24
#heat1 = np.array(list(np.loadtxt("heat_1.csv", delimiter=","))).astype("float32")
heat2 = np.array(list(np.loadtxt("heat_2.csv", delimiter=","))).astype("float32")

f=open('heat_3.csv','ab')

def chirality(theta,phi):
    s=np.zeros((L,L,3))
    s[:,:,0]=np.sin(theta)*np.cos(phi)
    s[:,:,1]=np.sin(theta)*np.sin(phi)
    s[:,:,2]=np.cos(theta)
    
    t=np.zeros((L+2,L+2,3))
    t[1:L+1,1:L+1,:]=s
    chi=0
    for i in range(L):
        for j in range(L):
            ti=t[i+1,j+1,:]
            tx=t[i+2,j+1,:]
            tx_=t[i,j+1,:]
            ty=t[i+1,j+2,:]
            ty_=t[i+1,j,:]
            chi=chi+np.dot(ti,np.cross(tx,ty))+np.dot(ti,np.cross(ty,tx_))+np.dot(ti,np.cross(ty_,tx))+np.dot(ti,np.cross(tx_,ty_))
    return chi

def magnet(theta):
    z=np.zeros((L,L))
    z=np.cos(theta)
    return (np.sum(z)/(L*L))

heat2=np.reshape(heat2,(-1,L*L*2))

n=np.shape(heat2)[0]
print(np.shape(heat2))

X=[]
m=[]
t=time.time()
for i in range(n):
    s=np.zeros((L,L,3))
    theta=heat2[i,:(L*L)]
    phi=heat2[i,(L*L):]
    theta=np.reshape(theta,(L,L))
    phi=np.reshape(phi,(L,L))
    s[:,:,0]=np.sin(theta)*np.cos(phi)
    s[:,:,1]=np.sin(theta)*np.sin(phi)
    s[:,:,2]=np.cos(theta)
    s=np.reshape(s,(1,-1))
    np.savetxt(f,s,delimiter=',')
    #X.append(chirality(theta,phi))
    m.append(magnet(theta))
    
np.savetxt('heat_4.csv',m,delimiter=',')
#plt.plot(X)
#plt.show()
print(time.time()-t)
#plt.plot(m)
#plt.show()
    
    