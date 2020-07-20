import numpy as np
import time
import matplotlib.pyplot as plt

L=24
s=np.zeros((L,L,3))
n_iter=100
J=1
lambd=6
D=(2**0.5)*np.tan(2*np.pi/lambd)

np.set_printoptions(precision=2)

#f=open('PureH_test.csv','ab')

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
    
    H=-J*(np.sum(s*(sx+sy))) + D*np.sum(np.cross(s,sx,axisa=2,axisb=2)[:,:,1]) - D*np.sum(np.cross(s,sy,axisa=2,axisb=2)[:,:,0]) 
    return H

def tmag(s):
    mx=np.sum(s[:,:,0])**2
    my=np.sum(s[:,:,1])**2
    mz=np.sum(s[:,:,2])**2
    return ((mx+my+mz)**0.5)/(L*L)

def magnet(s):
    return np.sum(s[:,:,2])

def approx_chirality(s):
    s=np.reshape(s,(L,L,3))
    
    sx=np.roll(s,-1,axis=0)
    sx_=np.roll(s,1,axis=0)
    sy=np.roll(s,-1,axis=1)
    sy_=np.roll(s,1,axis=1)
    
    chi1=np.sum(s*np.cross(sx,sy,axisa=2,axisb=2))
    chi2=np.sum(s*np.cross(sy,sx_,axisa=2,axisb=2))
    chi3=np.sum(s*np.cross(sy_,sx,axisa=2,axisb=2))
    chi4=np.sum(s*np.cross(sx_,sy_,axisa=2,axisb=2))
    
    chi=chi1+chi2+chi3+chi4
    return (chi/(8*np.pi))/L

def chi(a,b,c):
    t1=np.sum(a*np.cross(b,c,axisa=2,axisb=2),axis=2)
    t2=1+np.sum(a*b+b*c+c*a,axis=2)
    print(np.sum(t2<0),end=', ')
    print(t1[np.nonzero(t2<0)],t2[np.nonzero(t2<0)],np.sum(t1[np.nonzero(t2<0)]>0)-np.sum(t1[np.nonzero(t2<0)]<0))
    chi=np.sum(np.arctan(t1/t2))#+np.pi*(np.sum(t1[np.nonzero(t2<0)]>0)-np.sum(t1[np.nonzero(t2<0)]<0))
    return chi

def chirality(s):
    s=np.reshape(s,(L,L,3))
    
    sx=np.roll(s,-1,axis=0)
    sx_=np.roll(s,1,axis=0)
    sy=np.roll(s,-1,axis=1)
    sy_=np.roll(s,1,axis=1)
    
    return (chi(s,sx,sy)+chi(s,sx_,sy_)+chi(s,sy,sx_)+chi(s,sy_,sx))/(4*np.pi)/L
'''
    t1=np.sum(s*np.cross(sx,sy,axisa=2,axisb=2),axis=2)
    t2=1+np.sum(s*sx+sx*sy+sy*s,axis=2)
    np.sum(t1[np.nonzero(t2<0)]<0)
    print(np.sum(t2<0),end=', ')
    print(t1[np.nonzero(t2<0)],t2[np.nonzero(t2<0)],np.sum(t1[np.nonzero(t2<0)]>0)-np.sum(t1[np.nonzero(t2<0)]<0))
    chi1=np.sum(np.arctan(t1/t2))+np.pi*(np.sum(t1[np.nonzero(t2<0)]>0)-np.sum(t1[np.nonzero(t2<0)]<0))
    
    t1=np.sum(s*np.cross(sx_,sy_,axisa=2,axisb=2),axis=2)
    t2=1+np.sum(s*sx_+sx_*sy_+sy_*s,axis=2)
    print(np.sum(t2<0),end=', ')
    print(t1[np.nonzero(t2<0)],t2[np.nonzero(t2<0)],np.sum(t1[np.nonzero(t2<0)]>0)-np.sum(t1[np.nonzero(t2<0)]<0))
    chi2=np.sum(np.arctan(t1/t2))
    
    t1=np.sum(s*np.cross(sy,sx_,axisa=2,axisb=2),axis=2)
    t2=1+np.sum(s*sx_+sx_*sy+sy*s,axis=2)
    print(np.sum(t2<0),end=', ')
    print(t1[np.nonzero(t2<0)],t2[np.nonzero(t2<0)],np.sum(t1[np.nonzero(t2<0)]>0)-np.sum(t1[np.nonzero(t2<0)]<0))
    chi3=np.sum(np.arctan(t1/t2))
    
    t1=np.sum(s*np.cross(sy_,sx,axisa=2,axisb=2),axis=2)
    t2=1+np.sum(s*sx+sx*sy_+sy_*s,axis=2)
    print(np.sum(t2<0),end=', ')
    print(t1[np.nonzero(t2<0)],t2[np.nonzero(t2<0)],np.sum(t1[np.nonzero(t2<0)]>0)-np.sum(t1[np.nonzero(t2<0)]<0))
    chi4=np.sum(np.arctan(t1/t2))
    
    chi=chi1+chi2+chi3+chi4
    return (chi/(4*np.pi))/L
'''

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
    
    E_J=-J*np.sum((s2-s)*(sx+sx_+sy+sy_),axis=2)
    E_B=-B*(s2-s)[:,:,2]
    E_D=D*np.dot(np.cross(s2-s,sx-sx_,axisa=2,axisb=2),y) - D*np.dot(np.cross(s2-s,sy-sy_,axisa=2,axisb=2),x)
    E_T=-T*(np.log(np.sin(theta2)+0.00001)-np.log(np.sin(theta)+0.00001))
    del_E=E_J+E_B+E_D+E_T 
    
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

step=[] #list of values required for making temperature grid
for T in np.arange(2,0,-0.1):
    step.append(int(100*np.log(3/T)))
step.append(0)

for B in np.arange(1,4,0.2):
    tim=time.time()
    Ham=[]
    C=[]

    theta=np.random.rand(L,L)*(np.pi)
    phi=np.random.rand(L,L)*(2*np.pi)
    s[:,:,0]=np.sin(theta)*np.cos(phi)
    s[:,:,1]=np.sin(theta)*np.sin(phi)
    s[:,:,2]=np.cos(theta)
    
    j=0
    for t in np.arange(5):
        T=3*np.exp(-0.01*t)
        H=0
        H2=0
        for m in range(n_iter):
            [s,theta,phi]=monte(s,theta,phi,B,T,0)
            [s,theta,phi]=monte(s,theta,phi,B,T,1)
            H+=Hamiltonian(s)
            H2+=(Hamiltonian(s)**2)
        Ham.append(H/n_iter)
        C.append((H2/n_iter-(H/n_iter)**2))
        '''    
        if t==step[j]:
            j+=1
            #np.savetxt(f,np.reshape(s,(1,-1)),delimiter=',')
            if j==20: break
        '''
    #np.savetxt(f,np.reshape(s,(1,-1)),delimiter=',')
    
    #x=3*np.exp(-0.005*np.arange(t+1))
    #plt.plot(x,Ham)
    #plt.show()
    #plt.plot(x[20:],gaussian_filter(C[20:],10))
    #plt.show()
    plt.imshow(s[:,:,2])
    plt.show()
    print('{',B,'}',end=' : ')
    print(chirality(s)*L,end=', ')
    #print(approx_chirality(s)*L,end=', ')
    print(magnet(s)/(L*L),end=', ')
    print(tmag(s))

    elapsed = time.time() - tim
    print("Time elapsed = :",elapsed)
    print("\n")

#f.close()