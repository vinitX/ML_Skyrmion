import numpy as np
import matplotlib.pyplot as plt

fine=np.array(list(np.loadtxt("Fail.csv", delimiter=","))).astype("float32")

def entropy(s,B,T):
    s=np.reshape(s,(24,24,3))
    
    J=1
    D=6**0.5
    beta=1/T
    
    y=[0,1,0]
    x=[1,0,0]
    
    sin_theta = (s[:,:,0]**2 + s[:,:,1]**2)**(0.5)
    
    sx=np.roll(s,-1,axis=0)
    sy=np.roll(s,-1,axis=1)
    
    E_J=-J*np.dot(np.multiply(s,sx+sy),np.ones(3))
    E_B=-B*(s)[:,:,2]
    E_D=+D*np.dot(np.cross(s,sx,axisa=2,axisb=2),y) - D*np.dot(np.cross(s,sy,axisa=2,axisb=2),x)
    E_T=-T*(np.log(sin_theta+0.00001))
    E=E_J+E_B+E_D+E_T
    
    p=np.exp(-beta*E)
    Z=np.sum(p)
    
    E_exp=np.sum(E*p)/Z
    
    S = E_exp/T + np.log(Z)
    C = np.sum(((E-E_exp)**2)*p)/(Z*T*T)
    
    return [S,C]

print(np.shape(fine))

fine=np.reshape(fine,(22,2020,1728))

S=np.zeros((22,2020,2))
for i in range(22):
    for j in range(2020):
        S[i,j,:]=entropy(fine[i,j,:],i/5,0.1)
        
C=np.reshape(S[:,:,1],(22,2020))
S=np.reshape(S[:,:,0],(22,2020))

plt.plot(np.mean(C,axis=1))
plt.show()
plt.plot(np.mean(S,axis=1))
plt.show()


        
        