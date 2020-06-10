import numpy as np
import random
import matplotlib.pyplot as plt
import time

L = 100
J=1
nsteps = 2000
nbr = {i : ((i + 1) % L,(i - 1) % L) for i in range(L)}

def cluster(s,T):
    p  = 1.0 - np.exp(-2.0 / T)
    
    k = random.randint(0, L-1)
    Pocket, Cluster = [k], [k]
    while Pocket != []:
        j = random.choice(Pocket)
        for l in nbr[j]:
            if s[l] == s[j] and l not in Cluster \
                and random.uniform(0.0, 1.0) < p:
                Pocket.append(l)
                Cluster.append(l)
        Pocket.remove(j)
    for j in Cluster:
        s[j] *= -1
    
    return s

f=open('spin.csv','ab')
g=open('H.csv','ab')
h=open('H2.csv','ab')


s = np.random.choice([1,-1],size=L)
H=[]
H2=[]
mag=[]
C=[]
Q=[]
Q2=[]
m=[]
t=time.time()
           
for T in np.arange(2,0,-8e-5):
    m_temp=0
    h_temp=0
    h2_temp=0

    for step in range(nsteps):
        s=cluster(s,T)
        E=-J*(np.dot(s,np.roll(s,-1)))
        
        if step>=nsteps//2:
            m_temp+=np.abs(np.mean(s))
            h_temp+=E
            h2_temp+=(E*E)
            if step%100==99:
                np.savetxt(f,s,delimiter=',')
    mag.append(m_temp/(nsteps//2))
    H.append(h_temp/(nsteps//2))
    H2.append(h2_temp/(nsteps//2))
np.savetxt(g,H,delimiter=',')
np.savetxt(h,H2,delimiter=',')

print(np.shape(H))
length=25000
k=(len(H)//length)
for i in range(length):
    C.append((np.mean(H2[i*k:(i+1)*k])-np.mean(H[i*k:(i+1)*k])**2)/(2-2*i/length))
    Q.append(np.mean(H[i*k:(i+1)*k]))
    Q2.append(np.mean(H2[i*k:(i+1)*k]))
    m.append(np.mean(mag[i*k:(i+1)*k]))
print(time.time()-t)
plt.plot(C[:])
plt.show()
plt.plot(Q)
plt.show()
plt.plot(Q2)
plt.show()
plt.plot(m)
plt.show()

f.close()
g.close()
h.close()