import numpy as np
import matplotlib.pyplot as plt

L=24

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

s=np.loadtxt('HDMZA1.csv',delimiter=',')
print(np.shape(s))
s=np.reshape(s,(-1,L,L,3))
print(np.shape(s))
X=np.zeros((20,20))
for i in range(np.shape(s)[0]):
    a=i//20
    b=i%20
    X[a,b]=chirality(s[i,:,:])
plt.imshow(X)