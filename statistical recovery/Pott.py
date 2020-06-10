import numpy as np
import time
import matplotlib.pyplot as plt

L=100
n_iter=10000

f=open('Bloom.csv','ab')

E=0
hist=[]

white=np.zeros((L,L),dtype=bool)
black=np.zeros((L,L),dtype=bool)
for i in range(L):
    for j in range(L):
        if (i+j)%2==0:
            white[i][j]=True
        else:
            black[i][j]=True

def monte(s,D,T,colour,E):
    J=1
    beta=1/T
    
    s2=np.random.choice([-1,0,1],(L,L))
    
    sx=np.roll(s,-1,axis=0)
    sx_=np.roll(s,1,axis=0)
    sy=np.roll(s,-1,axis=1)
    sy_=np.roll(s,1,axis=1)
    
    del_E=D*(s2*s2-s*s) - J*(s2-s)*(sx+sx_+sy+sy_)
    
    trans_prob= np.logical_or(del_E<0, np.exp(-beta*del_E)>np.random.rand(L,L))
    if colour==1:
        trans_prob=trans_prob*white
    else:
        trans_prob=trans_prob*black
        
    s=s+(s2-s)*trans_prob

    E=E+np.sum(del_E*trans_prob)
    hist.append(E)
    return [s,E]

T=0.1




for D in np.arange(1.9,2,0.005):
    t = time.time()
    s=np.random.choice([-1,0,1],(L,L))
    hist=[]
    E=0
    #a=np.loadtxt('thetaB'+str(B)+'T'+str(T+1)+'.txt')
    #theta=np.reshape(a,(L,L))
    for m in range(n_iter):
        [s,E]=monte(s,D,T,0,E)
        [s,E]=monte(s,D,T,1,E)
        '''if m%100==0 and m>10000:
            s_temp=np.reshape(s,(1,-1))
            np.savetxt(f,s_temp,delimiter=',')

        print("{0:.2f}".format(T/10),end=' ')
        a=np.reshape(theta,(1,L*L))
        b=np.reshape(phi,(1,L*L))
        filename_a='thetaB'+str(B)+'T'+str(T)+'.txt'
        filename_b='phiB'+str(B)+'T'+str(T)+'.txt'
        np.savetxt(filename_a,a)
        np.savetxt(filename_b,b)'''

    plt.imshow(s)
    plt.show()
    print("D :",D)
    print("m :",np.mean(s))
    plt.plot(hist)
    plt.show()
    print("Time elapsed = :",time.time()-t)

f.close()