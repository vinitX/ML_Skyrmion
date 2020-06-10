import matplotlib.pyplot as plt
import math
import numpy as np
import random

L=18
theta=np.zeros((L,L))
phi=np.zeros((L,L))
s=np.zeros((L,L,3))
chi=np.zeros((L,L))
'''
Bz=3.5

for magnet in range(6,11):
    J=1
    D=6**0.5
    B=[0,0,Bz] #acting along z-axis
    T=0.15
    beta=1/T
    e=0.2#1
    update=0
    step=0
    E_prof=[]
    TE=0

    filename_theta="theta"+str(magnet-1)+".txt"
    filename_phi="phi"+str(magnet-1)+".txt"
    a=np.loadtxt(filename_theta)
    theta=a.reshape(L,L)
    b=np.loadtxt(filename_phi)
    phi=b.reshape(L,L)
    #initializing s
    '''
for i in range(0,L):
    for j in range(0,L):
        theta[i,j]=random.random()*math.pi
        phi[i,j]=random.random()*2*math.pi
        s[i,j,0]=math.sin(theta[i,j])*math.cos(phi[i,j])
        s[i,j,1]=math.sin(theta[i,j])*math.sin(phi[i,j])
        s[i,j,2]=math.cos(theta[i,j])
'''
    for i in range(0,L):
        for j in range(0,L):
            y=np.array([0,1,0])
            x=np.array([1,0,0])
            si=s[i,j,:]
            
            if i+1==L: sx=s[0,j,:]
            else: sx=s[i+1,j,:]        
            if j+1==L: sy=s[i,0,:]
            else: sy=s[i,j+1,:]
                
            Ei=-J*np.dot(si,sx+sy)+D*np.dot(y,np.cross(si,sx))-D*np.dot(x,np.cross(si,sy))-np.dot(B,si)
            TE=TE+Ei
            
    for m in range(25001):
        for i in range(0,L):
            for j in range(0,L):
                theta2=theta[i,j]+e*(random.random()-0.5)
                phi2=phi[i,j]+e*(random.random()-0.5)
                s2=np.zeros(3)
                s2[0]=math.sin(theta2)*math.cos(phi2)
                s2[1]=math.sin(theta2)*math.sin(phi2)
                s2[2]=math.cos(theta2)
                y=np.array([0,1,0])
                x=np.array([1,0,0])
                si=s[i,j,:]
                
                if i+1==L: sx=s[0,j,:]
                else: sx=s[i+1,j,:]
                if i-1==-1: sx_=s[L-1,j,:]
                else: sx_=s[i-1,j,:]
                
                if j+1==L: sy=s[i,0,:]
                else: sy=s[i,j+1,:]
                if j-1==-1: sy_=s[i,L-1,:]
                else: sy_=s[i,j-1,:]
                
                del_E=-J*np.dot(s2-si,sx+sx_+sy+sy_) - np.dot(B,s2-si) + D*np.dot(y,np.cross(s2-si,sx-sx_)) - D*np.dot(x,np.cross(s2-si,sy-sy_))
                step+=1
                
                chi[i,j]=np.dot(si,np.cross(sx,sy))
                
                if del_E<0 or math.exp(-beta*del_E)>random.random():
                    update+=1
                    TE=TE+del_E
                    theta[i,j]=theta2
                    phi[i,j]=phi2
                    s[i,j,:]=s2[:]
        if m%10==0:
            E_prof.append(TE)
            
        if m%5000==0:
            X, Y = np.meshgrid(np.arange(0, L, 1), np.arange(0, L, 1))
            plt.figure()
            plt.title('Heisenberg Model')
            Q = plt.quiver(X, Y, s[:,:,0], s[:,:,1], units='width')
            plt.figure()
            plt.imshow(chi)
                   
    print("Field : ",B)
    print("Acceptance Ratio: ",update/step)
    print("Energy : ",TE)
    
    
                                    
    plt.plot(E_prof)
    plt.ylabel('Energy')
    plt.xlabel('Time')
    plt.show()

    a=np.reshape(theta,(1,L*L))
    b=np.reshape(phi,(1,L*L))
    filename_theta="theta"+str(magnet)+".txt"
    filename_phi="phi"+str(magnet)+".txt"
    np.savetxt(filename_theta,a)
    np.savetxt(filename_phi,b)
    
    Bz+=0.5
'''
a=np.reshape(theta,(1,L*L))
b=np.reshape(phi,(1,L*L))
np.savetxt("theta0.txt",a)
np.savetxt("phi0.txt",b)
'''