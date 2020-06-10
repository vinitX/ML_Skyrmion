import numpy as np
import random
import time

L=24
theta=np.zeros((L,L))
phi=np.zeros((L,L))
s=np.zeros((L,L,3))
h=open('test4.csv','ab')

def monte(s,theta,phi,magnet,temp):
    J=1
    D=6**0.5
    B=[0,0,magnet/2] #acting along z-axis
    T=temp/10
    beta=1/T
    e=0.2#1
    update=0
    step=0
    
    theta2=(theta+e*(np.random.rand(L,L)-0.5))%(np.pi) #theta has to be in range (0,pi)
    phi2=phi+e*(np.random.rand(L,L)-0.5)
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
            del_E=-J*np.dot(s2[i,j,:]-si,sx+sx_+sy+sy_) - np.dot(B,s2[i,j,:]-si) + D*np.dot(y,np.cross(s2[i,j,:]-si,sx-sx_)) - D*np.dot(x,np.cross(s2[i,j,:]-si,sy-sy_)) - T*(np.log(np.sin(theta2[i,j]))-np.log(np.sin(theta[i,j])))
            step+=1
            if del_E<0 or np.exp(-beta*del_E)>random.random():
                update+=1
                theta[i,j]=theta2[i,j]
                phi[i,j]=phi2[i,j]
                s[i,j,:]=s2[i,j,:]  
    return [s,theta,phi,step,update]

theta=np.random.rand(1,L*L)*(np.pi)
phi=np.random.rand(1,L*L)*(2*np.pi)
filename_a="theta"+str(-1)+".txt"
filename_b="phi"+str(-1)+".txt"
np.savetxt(filename_a,theta)
np.savetxt(filename_b,phi)

for magnet in range(0,10):
    for temp in range(1,6):
        t = time.time()
        filename_a="theta"+str(magnet-1)+".txt" #building on existing data
        filename_b="phi"+str(magnet-1)+".txt"
        a=np.loadtxt(filename_a)
        theta=a.reshape(L,L)
        b=np.loadtxt(filename_b)
        phi=b.reshape(L,L)
        #theta=np.random.rand(L,L)*(np.pi)
        #phi=np.random.rand(L,L)*(2*np.pi)
        for m in range(50001):
            [s,theta,phi,step,update]=monte(s,theta,phi,magnet,temp)
            if m%1==0:# and m>20000: 
                a=np.reshape(theta,(1,L*L))
                b=np.reshape(phi,(1,L*L))
                c=np.append(a,b)
                c=np.reshape(c,(1,2*L*L))
                np.savetxt(h,c,delimiter=',')
            if m%100==0: print(int(m/100),end=' ')
        print("\n\n\n")
        print("Field : ",magnet/2)
        print("Temperature : ",temp/10)
        print("Acceptance Ratio: ",update/step)
        filename_a="theta"+str(magnet)+".txt"
        filename_b="phi"+str(magnet)+".txt"
        np.savetxt(filename_a,a)
        np.savetxt(filename_b,b)
        elapsed = time.time() - t
        print("Time elapsed for B = :",magnet/2," : ",elapsed)
        print("\n\n\n")