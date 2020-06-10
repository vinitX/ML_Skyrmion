import numpy as np
import math
import random
import matplotlib.pyplot as plt
import time

L=20
nsteps=100000
N = L * L
nbr = {i : ((i // L) * L + (i + 1) % L, (i + L) % N,
            (i // L) * L + (i - 1) % L, (i - L) % N)
                                    for i in range(N)}

def cluster(S,T):
    p  = 1.0 - math.exp(-2.0 / T)
        
    k = random.randint(0, N - 1)
    Pocket, Cluster = [k], [k]
    while Pocket != []:
        j = random.choice(Pocket)
        for l in nbr[j]:
            if S[l] == S[j] and l not in Cluster and random.uniform(0.0, 1.0) < p:
                Pocket.append(l)
                Cluster.append(l)
        Pocket.remove(j)
    for j in Cluster:
        S[j] *= -1
    return S


f=open('2D_spin.csv','ab')
g=open('2D_H.csv','ab')
h=open('2D_H2.csv','ab')

H=[]
H2=[]
mag=[]
C=[]
Q=[]
Q2=[]
m=[]
t_hist=[]
t_global=time.time()

s = np.random.choice([1,-1],size=(L,L))

for T in np.arange(5,0,-0.02):
    if T<2.2691: nsteps=10000
    t=time.time()
    m_temp=0
    h_temp=0
    h2_temp=0
    print(T,end=' ')
    for step in range(nsteps):
        s=np.reshape(cluster(np.reshape(s,(L*L)),T),(L,L))
        E=-np.sum(s*np.roll(s,-1,axis=0))-np.sum(s*np.roll(s,-1,axis=1))
        
        if step>=nsteps//2:
            m_temp+=np.abs(np.mean(s))
            h_temp+=E
            h2_temp+=((E*E)/(nsteps//2))
            if step%(nsteps//100)==((nsteps//100)-1):
                np.savetxt(f,s,delimiter=',')

    mag.append(m_temp/(nsteps//2))
    H.append(h_temp/(nsteps//2))
    H2.append(h2_temp)
    
    t_hist.append(time.time()-t)
np.savetxt(g,H,delimiter=',')
np.savetxt(h,H2,delimiter=',')


print(np.shape(H))
length=len(H)
k=(len(H)//length)
for i in range(length):
    C.append((np.mean(H2[i*k:(i+1)*k])-np.mean(H[i*k:(i+1)*k])**2)/(2-2*i/length))
    Q.append(np.mean(H[i*k:(i+1)*k]))
    Q2.append(np.mean(H2[i*k:(i+1)*k]))
    m.append(np.mean(mag[i*k:(i+1)*k]))
print(time.time()-t_global)
plt.plot(C[:])
plt.show()
plt.plot(Q)
plt.show()
plt.plot(Q2)
plt.show()
plt.plot(m)
plt.show()
plt.plot(t_hist)
plt.show()

f.close()
g.close()
h.close()