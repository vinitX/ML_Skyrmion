import matplotlib.pyplot as plt
import numpy as np
import random

L=24
theta=np.zeros((L,L))
phi=np.zeros((L,L))
s=np.zeros((L,L,3))
chi=np.zeros((L,L))

J=1
D=6**0.5
B=[0,0,0] #acting along z-axis
T=0.15
beta=1/T
e=0.2#1
update=0
step=0
E_prof=[]
TE=0
'''
a=np.loadtxt("theta0.txt")
theta=a.reshape(L,L)
b=np.loadtxt("phi0.txt")
phi=b.reshape(L,L)
    #initializing s
'''
for i in range(0,L):
    for j in range(0,L):
        theta[i,j]=random.random()*np.pi
        phi[i,j]=random.random()*2*np.pi
        s[i,j,0]=np.sin(theta[i,j])*np.cos(phi[i,j])
        s[i,j,1]=np.sin(theta[i,j])*np.sin(phi[i,j])
        s[i,j,2]=np.cos(theta[i,j])

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
        Ei=-J*np.dot(si,sx+sy) + D*np.dot(y,np.cross(si,sx)) - D*np.dot(x,np.cross(si,sy))-np.dot(B,si) - T*(np.log(np.sin(theta[i,j])+0.001))
        TE=TE+Ei


for m in range(20001):
    for i in range(0,L):
        for j in range(0,L):
            theta2=(theta[i,j]+e*(random.random()-0.5))%(np.pi) #theta has to be in range (0,pi)
            if np.sin(theta2)<0.00001: continue #else there will be a problem in taking log of sin(theta2)
        
            phi2=phi[i,j]+e*(random.random()-0.5)
            s2=np.zeros(3)
            s2[0]=np.sin(theta2)*np.cos(phi2)
            s2[1]=np.sin(theta2)*np.sin(phi2)
            s2[2]=np.cos(theta2)
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
            
            del_E=-J*np.dot(s2-si,sx+sx_+sy+sy_) - np.dot(B,s2-si) + D*np.dot(y,np.cross(s2-si,sx-sx_)) - D*np.dot(x,np.cross(s2-si,sy-sy_)) - T*(np.log(np.sin(theta2))-np.log(np.sin(theta[i,j])))
            step+=1
            
            chi[i,j]=np.dot(si,np.cross(sx,sy))
                
            if del_E<0 or np.exp(-beta*del_E)>random.random():
                update+=1
                TE=TE+del_E
                theta[i,j]=theta2
                phi[i,j]=phi2
                s[i,j,:]=s2[:]
    print(m)
    
    
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
'''
a=np.reshape(theta,(1,L*L))
b=np.reshape(phi,(1,L*L))
np.savetxt("theta0.txt",a)
np.savetxt("phi0.txt",b)
    