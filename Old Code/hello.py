import numpy as np

L=24

def chirality1(s):
    s=np.reshape(s,(L,L,3))
    
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

def chirality2(s):
    s=np.reshape(s,(L,L,3))
    
    sx=np.roll(s,-1,axis=0)
    sx_=np.roll(s,1,axis=0)
    sy=np.roll(s,-1,axis=1)
    sy_=np.roll(s,1,axis=1)
    
    chi1=np.sum(np.multiply(np.cross(sx,sy,axisa=2,axisb=2),s))
    chi2=np.sum(np.multiply(np.cross(sy,sx_,axisa=2,axisb=2),s))
    chi3=np.sum(np.multiply(np.cross(sy_,sx,axisa=2,axisb=2),s))
    chi4=np.sum(np.multiply(np.cross(sx_,sy_,axisa=2,axisb=2),s))
    
    chi=chi1+chi2+chi3+chi4
    
    return chi

a=np.loadtxt('thetaB14T3.txt')
theta=np.reshape(a,(L,L))
b=np.loadtxt('phiB14T3.txt')
phi=np.reshape(b,(L,L))
s=np.zeros((L,L,3))
s[:,:,0]=np.sin(theta)*np.cos(phi)
s[:,:,1]=np.sin(theta)*np.sin(phi)
s[:,:,2]=np.cos(theta)

chi1=chirality1(s)
chi2=chirality2(s)
print(chi2,chi1)