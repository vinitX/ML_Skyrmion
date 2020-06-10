import numpy as np
import random
import time

L=24
theta=np.zeros((L,L))
phi=np.zeros((L,L))
s=np.zeros((L,L,3))
n_iter=100

def monte(s,theta,phi,b,T):
    J=1
    D=6**0.5
    B=[0,0,b] #acting along z-axis
    beta=1/T
    e=0.2#1
    update=0
    step=0
    
    theta2=(theta+e*(np.random.rand(L,L)-0.5))%(np.pi) #theta has to be in range (0,pi)
    phi2=(phi+e*(np.random.rand(L,L)-0.5))%(2*np.pi)
    s2=np.zeros((L,L,3))
    s2[:,:,0]=np.sin(theta2)*np.cos(phi2)
    s2[:,:,1]=np.sin(theta2)*np.sin(phi2)
    s2[:,:,2]=np.cos(theta2)
    y=np.array([0,1,0])
    x=np.array([1,0,0])
    TF=0
    for i in range(0,L):
        for j in range(0,L):
            if np.sin(theta2[i,j])<0.00001: continue #else there will be a problem in taking log of sin(theta2)
            si=s[i,j,:]                
            if i+1==L: sx=s[0,j,:]
            else: sx=s[i+1,j,:]
            if i-1==-1: sx_=s[L-1,j,:]
            else: sx_=s[i-1,j,:]            
            if j+1==L: sy=s[i,0,:]
            else: sy=s[i,j+1,:]
            if j-1==-1: sy_=s[i,L-1,:]
            else: sy_=s[i,j-1,:]            
            del_E=-J*np.dot(s2[i,j,:]-si,sx+sx_+sy+sy_) - np.dot(B,s2[i,j,:]-si) + D*np.dot(y,np.cross(s2[i,j,:]-si,sx-sx_)) - D*np.dot(x,np.cross(s2[i,j,:]-si,sy-sy_))# - T*(np.log(np.sin(theta2[i,j]))-np.log(np.sin(theta[i,j])))
            step+=1
            if del_E<0 or np.exp(-beta*del_E)>random.random():
                TF+=del_E
                update+=1
                theta[i,j]=theta2[i,j]
                phi[i,j]=phi2[i,j]
                s[i,j,:]=s2[i,j,:]  
    print(TF)            
    return [s,theta,phi,step,update]

for B in range(1):#,9,2):
    for T in range(20,21):#,-1):
        t = time.time()
        if T==20:
            theta=np.random.rand(L,L)*(np.pi)
            phi=np.random.rand(L,L)*(2*np.pi)
        else:    
            a=np.loadtxt('thetaB'+str(B)+'T'+str(T+1)+'.txt')
            theta=np.reshape(a,(L,L))
            b=np.loadtxt('phiB'+str(B)+'T'+str(T+1)+'.txt')
            phi=np.reshape(b,(L,L))

        for m in range(n_iter):
            [s,theta,phi,step,update]=monte(s,theta,phi,B/2,T/10)

            if m%100==0:
                print(int(m/100),end=' ')
        print("\n\n\n")
        print("Magnetic Field: ",B/2)
        print("Temperature: ",T/10)
        print("Acceptance Ratio: ",update/step)

        elapsed = time.time() - t
        print("Time elapsed = :",elapsed)
        print("\n\n\n")