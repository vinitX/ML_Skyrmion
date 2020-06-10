import numpy as np
import random
import time
from numba import vectorize

L=24
theta=np.zeros((L,L))
phi=np.zeros((L,L))
s=np.zeros((L,L,3))

@vectorize(['(i4)'])
def monte(magnet):
    t = time.time()
    J=1
    D=6**0.5
    B=[0,0,magnet/2] #acting along z-axis
    T=0.15
    beta=1/T
    e=0.2#1
    update=0
    step=0
    if magnet==0:
        filename_a="theta"+str(magnet)+".txt"
        filename_b="phi"+str(magnet)+".txt"
    else:
        filename_a="theta"+str(magnet-1)+".txt"
        filename_b="phi"+str(magnet-1)+".txt"
    a=np.loadtxt(filename_a)
    theta=a.reshape(L,L)
    b=np.loadtxt(filename_b)
    phi=b.reshape(L,L)
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
                if del_E<0 or np.exp(-beta*del_E)>random.random():
                    update+=1
                    theta[i,j]=theta2
                    phi[i,j]=phi2
                    s[i,j,:]=s2[:]     
        if m%100==0 and m>10000:
            a=np.reshape(theta,(1,L*L))
            b=np.reshape(phi,(1,L*L))
            filename_a="theta"+str(magnet)+"_"+str(int(m/100)-200)+".txt"
            filename_b="phi"+str(magnet)+"_"+str(int(m/100)-200)+".txt"
            np.savetxt(filename_a,a)
            np.savetxt(filename_b,b)
            print(int(m/100),end=' ')
    print("\n\n\n")
    print("Field : ",B)
    print("Acceptance Ratio: ",update/step)
    a=np.reshape(theta,(1,L*L))
    b=np.reshape(phi,(1,L*L))
    filename_a="theta"+str(magnet)+".txt"
    filename_b="phi"+str(magnet)+".txt"
    np.savetxt(filename_a,a)
    np.savetxt(filename_b,b)
    elapsed = time.time() - t
    print("Time elapsed for B = :",magnet/2," : ",elapsed)
    print("\n\n\n")
    
    
for magnet in range(1,7):
    monte(magnet)

