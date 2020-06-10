import numpy as np
import time
import matplotlib.pyplot as plt
import random

L=20
n_iter=10000

f=open('2D_3spin.csv','ab')
g=open('2D_3H.csv','ab')
h=open('2D_3H2.csv','ab')

J=[1,0.7,0.5]

def monte(s,T):
    beta=1/T
    R=np.random.randint(-3,L-3,size=(L))
    for i in range(L):
        for j in range(L):
            del_E=2*s[R[i],R[j]]*((s[R[i],R[j]+1]+s[R[i]+1,R[j]])*J[0] \
                     +(s[R[i]+1,R[j]+1]+s[R[i]+1,R[j]-1])*J[1] \
                     +(s[R[i]+2,R[j]]+s[R[i],R[j]+2])*J[2])

        if del_E<0 or np.exp(-beta*del_E)>random.random():
            s[R[i],R[j]]=-s[R[i],R[j]]
    return s

H=[]
H2=[]
mag=[]

s=np.random.choice([-1,1],(L,L))
t = time.time()
for T in np.arange(5,0,-0.01):
    m_temp=0
    h_temp=0
    h2_temp=0

    for m in range(n_iter):
        s=monte(s,T)
        
        E = -J[0]*(np.sum(s*np.roll(s,1,axis=0))+np.sum(s*np.roll(s,1,axis=1))) \
                -J[1]*(np.sum(s*np.roll(np.roll(s,1,axis=0),1,axis=1))+np.sum(s*np.roll(np.roll(s,1,axis=0),-1,axis=1))) \
                -J[2]*(np.sum(s*np.roll(s,2,axis=0))+np.sum(s*np.roll(s,2,axis=1)))
        
        if m>=n_iter//2:
            m_temp+=np.abs(np.mean(s))
            h_temp+=(E/(n_iter//2))
            h2_temp+=((E*E)/(n_iter//2))
            
            if m%100==99:
                s_temp=np.reshape(s,(1,-1))
                np.savetxt(f,s_temp,delimiter=',')
    mag.append(m_temp/(n_iter//2))
    H.append(h_temp)
    H2.append(h2_temp)
    print("{0:.2f}".format(T),end=' ')
    
np.savetxt(g,H,delimiter=',')
np.savetxt(h,H2,delimiter=',')
length=500
k=1
C=[]

for i in range(length):
    C.append((np.mean(H2[i*k:(i+1)*k])-np.mean(H[i*k:(i+1)*k])**2)/(5-5*i/length))

print(time.time()-t)

plt.plot(C)
plt.show()
plt.plot(H)
plt.show()
plt.plot(H2)
plt.show()
plt.plot(mag)
plt.show()

plt.show()

f.close()
g.close()
h.close()