import numpy as np
import time
import matplotlib.pyplot as plt

L=24
s=np.zeros((L,L,3))
n_iter=100
J=1
Dy=6**0.5
Dx=-6**0.5
#J2=0.8

#f=open('Hre.csv','ab')

white=np.zeros((L,L),dtype=bool)
black=np.zeros((L,L),dtype=bool)
for i in range(L):
    for j in range(L):
        if (i+j)%2==0:
            white[i][j]=True
        else:
            black[i][j]=True

def Hamiltonian(s):
    sx=np.roll(s,-1,axis=0)
    sy=np.roll(s,-1,axis=1)
    
    H=-J*(np.sum(s*(sx+sy))) -B*np.sum(s[:,:,2]) -K*np.sum(s[:,:,2]*s[:,:,2]) + Dy*np.sum(np.cross(s,sx,axisa=2,axisb=2)[:,:,1]) + Dx*np.sum(np.cross(s,sy,axisa=2,axisb=2)[:,:,0])
    
    return H

def tmag(s):
    mx=np.sum(s[:,:,0])**2
    my=np.sum(s[:,:,1])**2
    mz=np.sum(s[:,:,2])**2
    return ((mx+my+mz)**0.5)/(L*L)

def magnet(s):
    return np.sum(s[:,:,2])

def chirality(s):
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

def monte(s,theta,phi,B,T,K,colour):
    beta=1/T
    e=0.2#1
    
    del_theta=e*(np.random.rand(L,L)-0.5)
    del_phi=e*(np.random.rand(L,L)-0.5)
    theta2=(theta+del_theta)%(np.pi) #theta has to be in range (0,pi)
    phi2=(phi+del_phi)%(2*np.pi)
    
    s2=np.zeros((L,L,3))
    s2[:,:,0]=np.sin(theta2)*np.cos(phi2)
    s2[:,:,1]=np.sin(theta2)*np.sin(phi2)
    s2[:,:,2]=np.cos(theta2)
    y=[0,1,0]
    x=[1,0,0]
    
    sx=np.roll(s,-1,axis=0)
    sx_=np.roll(s,1,axis=0)
    sy=np.roll(s,-1,axis=1)
    sy_=np.roll(s,1,axis=1)
    
    E_J=-J*np.sum((s2-s)*(sx+sx_+sy+sy_),axis=2)
    #E_J2=J2*np.sum((s2*sx)**2+(s2*sy)**2+(s2*sx_)**2+(s2*sy_)**2-(s*sx)**2-(s*sy)**2-(s*sx_)**2-(s*sy_)**2,axis=2)
    E_B=-B*(s2-s)[:,:,2]
    E_D=+Dy*np.dot(np.cross(s2-s,sx-sx_,axisa=2,axisb=2),y) + Dx*np.dot(np.cross(s2-s,sy-sy_,axisa=2,axisb=2),x)
    #E_T=-T*(np.log(np.sin(theta2)+0.00001)-np.log(np.sin(theta)+0.00001))
    E_A=-K*(s2[:,:,2]*s2[:,:,2]-s[:,:,2]*s[:,:,2])
    del_E=E_J+E_A+E_B +E_D #E_J2
    
    trans_prob= np.logical_or(del_E<0, np.exp(-beta*del_E)>np.random.rand(L,L))
    if colour==1:
        trans_prob=np.multiply(trans_prob,white)
    else:
        trans_prob=np.multiply(trans_prob,black)

    theta=(theta+np.multiply(trans_prob,del_theta))%(np.pi) #theta has to be in range (0,pi)
    phi=(phi+np.multiply(trans_prob,del_phi))%(2*np.pi)
    p=np.zeros((L,L,3))
    p[:,:,0]=s[:,:,0]+np.multiply(trans_prob,(s2-s)[:,:,0])
    p[:,:,1]=s[:,:,1]+np.multiply(trans_prob,(s2-s)[:,:,1])
    p[:,:,2]=s[:,:,2]+np.multiply(trans_prob,(s2-s)[:,:,2])

    return [p,theta,phi]
zeta=[]
for K in np.arange(-12,-8,10):
    for B in np.arange(1.5,3,3):
        Dx=6**0.5
        for alpha in np.arange(1,2,0.1):
            tim=time.time()
            Ham=[]

            theta=np.random.rand(L,L)*(np.pi)
            phi=np.random.rand(L,L)*(2*np.pi)
            s[:,:,0]=np.sin(theta)*np.cos(phi)
            s[:,:,1]=np.sin(theta)*np.sin(phi)
            s[:,:,2]=np.cos(theta)
            
            
            for t in np.arange(500):
                T=3*np.exp(-0.02*t)
                for m in range(n_iter):
                    [s,theta,phi]=monte(s,theta,phi,B,T,K,0)
                    [s,theta,phi]=monte(s,theta,phi,B,T,K,1)
                Ham.append(Hamiltonian(s))

            #np.savetxt(f,np.reshape(s,(1,-1)),delimiter=',')
        
            print('{',K,B,Dx,'}',end=' : ')
            print(chirality(s)/(L*L),end=', ')
            print(magnet(s)/(L*L),end=', ')
            print(tmag(s))
            
            zeta.append(Ham[:])
            plt.imshow(np.reshape(s[:,:,2],(L,L)))
            plt.show()
            elapsed = time.time() - tim
            print("Time elapsed = :",elapsed)
            print("\n")

plt.plot(zeta[:,-1])
#f.close()
