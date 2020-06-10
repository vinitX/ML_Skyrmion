import numpy as np
import time
import matplotlib.pyplot as plt

L=24
s=np.zeros((L,L,3))
n_iter=500
J=1
D=6**0.5
K=1.4

f=open('HDMZA1.csv','ab')

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
    
    H=-J*(np.sum(s*(sx+sy))) + D*np.sum(np.cross(s,sx,axisa=2,axisb=2)[:,:,1]) - D*np.sum(np.cross(s,sy,axisa=2,axisb=2)[:,:,0]) - K*np.sum(s[:,:,2]*s[:,:,2])
    return H

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

def monte(s,theta,phi,B,T,colour):
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
    
    E_J=-J*np.dot(np.multiply(s2-s,sx+sx_+sy+sy_),np.ones(3))
    E_B=-B*(s2-s)[:,:,2]
    E_D=+D*np.dot(np.cross(s2-s,sx-sx_,axisa=2,axisb=2),y) - D*np.dot(np.cross(s2-s,sy-sy_,axisa=2,axisb=2),x)
    E_T=-T*(np.log(np.sin(theta2)+0.00001)-np.log(np.sin(theta)+0.00001))
    E_A=-K*(s2[:,:,2]*s2[:,:,2]-s[:,:,2]*s[:,:,2])
    del_E=E_J+E_B+E_D+E_T+E_A
    
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

for A in np.arange(-1,3,0.2):
    tim = time.time()
    for B in np.arange(0,3,0.15):
        Ham=[]

        theta=np.random.rand(L,L)*(np.pi)
        phi=np.random.rand(L,L)*(2*np.pi)
        s[:,:,0]=np.sin(theta)*np.cos(phi)
        s[:,:,1]=np.sin(theta)*np.sin(phi)
        s[:,:,2]=np.cos(theta)
        for t in np.arange(200):
            T=3*np.exp(-0.03*t)
    
            for m in range(n_iter):
                [s,theta,phi]=monte(s,theta,phi,B,T,0)
                [s,theta,phi]=monte(s,theta,phi,B,T,1)
                Ham.append(Hamiltonian(s))

        np.savetxt(f,np.reshape(s,(1,-1)),delimiter=',')
        print(chirality(s)/(L*L))
        print("{",A,B,'}',end=' ')
        plt.plot(Ham)
        plt.show()
        plt.imshow(s[:,:,2])
        plt.show()

    print("\n")
    elapsed = time.time() - tim
    print("Time elapsed = :",elapsed)
    print("\n")

f.close()
