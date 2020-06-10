import numpy as np
import time
import matplotlib.pyplot as plt
import random

L=100
n_iter=50

#f=open('Bloom.csv','ab')

E=0
hist=[]

J=np.zeros(L)
for i in range(1,L):
    if i<L//2+1:
        J[i]=1/(i**3)
    else:
        J[i]=1/((L-i)**3)

def monte(s,T,E):
    beta=1/T
    for i in range(L):
        #i=np.random.randint(0,L)
    
        del_E=-np.dot(np.not_equal(s,s[i]),np.roll(J,i))+np.dot(np.equal(s,s[i]),np.roll(J,i))
        #print(np.roll(J,i),s)
        if del_E<0 or np.exp(-beta*del_E)>random.random():
            s[i]=-s[i]
            E=E+del_E
    hist.append(E)
    return [s,E]

T=0.1
mag=[]
s=np.random.choice([-1,1],(L))
t = time.time()
for T in np.arange(2,0,-0.001):
    
    hist=[]
    E=0
    #a=np.loadtxt('thetaB'+str(B)+'T'+str(T+1)+'.txt')
    #theta=np.reshape(a,(L,L))
    mag.append(0)
    for m in range(n_iter):
        [s,E]=monte(s,T,E)
        mag[-1]+=np.mean(s)
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

    #plt.plot(s)
    #plt.show()
    #print(T,end=' ')
    #print("m :",np.mean(s))
    #plt.plot(hist)
    #plt.show()
    

print("Time elapsed = :",time.time()-t)
plt.plot(np.abs(mag))
plt.show()

#f.close()