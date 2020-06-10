import numpy as np
import time
import matplotlib.pyplot as plt
import random

L=100
n_iter=5000

f=open('3spin.csv','ab')
g=open('3H.csv','ab')
h=open('3H2.csv','ab')

J=[1,0.5,0.33]
#J=[1,0,0]

def monte(s,T):
    beta=1/T
    R=np.random.randint(-4,L-4,size=L)
    for i in range(L):
        del_E=2*s[R[i]]*(s[R[i]+1]*J[0]+s[R[i]+2]*J[1]+s[R[i]+3]*J[2])

        if del_E<0 or np.exp(-beta*del_E)>random.random():
            s[R[i]]=-s[R[i]]
    return s

H=[]
H2=[]
mag=[]

s=np.random.choice([-1,1],(L))
t = time.time()
for T in np.arange(2,0,-0.001):
    m_temp=0
    h_temp=0
    h2_temp=0

    for m in range(n_iter):
        s=monte(s,T)
        E=0
        for i in range(3):
            E = E - J[i]*np.dot(s,np.roll(s,i+1)) #incorrect
        
        if m>=n_iter//2:
            m_temp+=np.abs(np.mean(s))
            h_temp+=E
            h2_temp+=(E*E)
            
            if m%100==99:
                s_temp=np.reshape(s,(1,-1))
                np.savetxt(f,s_temp,delimiter=',')
    mag.append(m_temp/(n_iter//2))
    H.append(h_temp/(n_iter//2))
    H2.append(h2_temp/(n_iter//2))
    print("{0:.2f}".format(T),end=' ')
    
np.savetxt(g,H,delimiter=',')
np.savetxt(h,H2,delimiter=',')
length=2000
k=1
C=[]

for i in range(length):
    C.append((np.mean(H2[i*k:(i+1)*k])-np.mean(H[i*k:(i+1)*k])**2)/(2-2*i/length))

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