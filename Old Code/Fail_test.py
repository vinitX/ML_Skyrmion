import numpy as np
import time

L=24
theta=np.zeros((L,L))
phi=np.zeros((L,L))
s=np.zeros((L,L,3))
n_iter=10001

f=open('Fine2.csv','ab')

white=np.zeros((L,L),dtype=bool)
black=np.zeros((L,L),dtype=bool)
for i in range(L):
    for j in range(L):
        if (i+j)%2==0:
            white[i][j]=True
        else:
            black[i][j]=True

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

def magnet(theta):
    z=np.zeros((L,L))
    z=np.cos(theta)
    return (np.sum(z)/(L*L))

def monte(s,theta,phi,B,T,colour):
    J=1
    D=6**0.5
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
    del_E=E_J+E_B+E_D+E_T
    
    trans_prob= np.logical_or(del_E<0, np.exp(-beta*del_E)>np.random.rand(L,L))
    if colour==1:
        trans_prob=np.multiply(trans_prob,white)
    else:
        trans_prob=np.multiply(trans_prob,black)

    #print(np.linalg.norm(E_D),np.linalg.norm(E_B),np.linalg.norm(del_E))
    #print(np.sum(np.multiply(del_E,trans_prob)))
    theta=(theta+np.multiply(trans_prob,del_theta))%(np.pi) #theta has to be in range (0,pi)
    phi=(phi+np.multiply(trans_prob,del_phi))%(2*np.pi)
    p=np.zeros((L,L,3))
    p[:,:,0]=s[:,:,0]+np.multiply(trans_prob,(s2-s)[:,:,0])
    p[:,:,1]=s[:,:,1]+np.multiply(trans_prob,(s2-s)[:,:,1])
    p[:,:,2]=s[:,:,2]+np.multiply(trans_prob,(s2-s)[:,:,2])
    
    return [p,theta,phi]


for B in range(321,341):
    t = time.time()
    for T in range(41):
        temp=2*np.exp(-0.1*T)
        if T==0:
            theta=np.random.rand(L,L)*(np.pi)
            phi=np.random.rand(L,L)*(2*np.pi)
        else:
            a=np.loadtxt('thetaB'+str(B)+'T'+str(T-1)+'.txt')
            theta=np.reshape(a,(L,L))
            b=np.loadtxt('phiB'+str(B)+'T'+str(T-1)+'.txt')
            phi=np.reshape(b,(L,L))
            
        for m in range(n_iter):
            [s,theta,phi]=monte(s,theta,phi,B/100,temp,0)
            [s,theta,phi]=monte(s,theta,phi,B/100,temp,1)
            if T>20 and m%100==0:
                s_temp=np.reshape(s,(1,-1))
                np.savetxt(f,s_temp,delimiter=',')
                
        print("{0:.2f}".format(temp),end=' ')
        a=np.reshape(theta,(1,L*L))
        b=np.reshape(phi,(1,L*L))
        filename_a='thetaB'+str(B)+'T'+str(T)+'.txt'
        filename_b='phiB'+str(B)+'T'+str(T)+'.txt'
        np.savetxt(filename_a,a)
        np.savetxt(filename_b,b)
    print("\n")
    print("Magnetic Field: ",B/100)
    elapsed = time.time() - t
    print("Time elapsed = :",elapsed)
    print("\n")
f.close()