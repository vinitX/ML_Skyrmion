import numpy as np

L=18
theta=np.zeros((L,L))
phi=np.zeros((L,L))
s=np.zeros((L,L,5))

a=np.loadtxt("theta5.txt")
theta=a.reshape(L,L)
b=np.loadtxt("phi5.txt")
phi=b.reshape(L,L)

for i in range(0,L):
    for j in range(0,L):
        s[i,j,0]=i
        s[i,j,1]=j
        s[i,j,2]=np.sin(theta[i,j])*np.cos(phi[i,j])
        s[i,j,3]=np.sin(theta[i,j])*np.sin(phi[i,j])
        s[i,j,4]=np.cos(theta[i,j])

a=np.reshape(s,(1,5*L*L))
np.savetxt("spin5.txt",a)