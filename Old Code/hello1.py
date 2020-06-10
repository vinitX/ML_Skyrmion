import numpy as np
import random
import time

L=24
theta=np.zeros((L,L))
phi=np.zeros((L,L))
s=np.zeros((L,L,3))
n_iter=10001

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
                update+=1
                theta[i,j]=theta2[i,j]
                phi[i,j]=phi2[i,j]
                s[i,j,:]=s2[i,j,:]  
    return [s,theta,phi,step,update]


for T in range(16,17):
    for B in range(0,1):
        t = time.time()
        a=np.loadtxt('thetab0T20.txt')
        theta=np.reshape(a,(L,L))
        b=np.loadtxt('phib0T20.txt')
        phi=np.reshape(b,(L,L))
        chi=0 #Net Chirality
        X=[] #Stores Net Chiralities of configurations
        mag=0 #Net Magnetisation
        for m in range(n_iter):
            [s,theta,phi,step,update]=monte(s,theta,phi,B/2,T/10)
            if m%100==0 and m>n_iter/2:
                temp=chirality(theta,phi)
                chi+=temp
                mag+=magnet(theta)
                print('\n',mag*100/(m-n_iter/2),'|',chi*100/(m-n_iter/2),end='\n')

            if m%100==0:
                print(int(m/100),end=' ')
        print("\n\n\n")
        print("Magnetic Field: ",B/2)
        print("Temperature: ",T/10)
        print("Acceptance Ratio: ",update/step)
        X=np.reshape(X,(1,-1))
        
        a=np.reshape(theta,(1,L*L))
        b=np.reshape(phi,(1,L*L))
        filename_a='thetab'+str(B)+'T'+str(T)+'.txt'
        filename_b='phib'+str(B)+'T'+str(T)+'.txt'
        np.savetxt(filename_a,a)
        np.savetxt(filename_b,b)
        elapsed = time.time() - t
        print("Time elapsed = :",elapsed)
        print("\n\n\n")