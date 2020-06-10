import numpy as np
import time

L=24

t=time.time()
S=np.loadtxt('L3.csv',delimiter=',')
S=np.reshape(S,(15,30,L,L,3))
print(time.time()-t)
J=1

#f=open('L2_fft.csv','ab')

def Ham(s):
    sx=np.roll(s,-1,axis=0)
    sy=np.roll(s,-1,axis=1)
    
    H=-J*(np.sum(s*(sx+sy))) -K*np.sum((s[:,:,2]*s[:,:,2])) + Dy*np.sum(np.cross(s,sx,axisa=2,axisb=2)[:,:,1]) + Dx*np.sum(np.cross(s,sy,axisa=2,axisb=2)[:,:,0]) #+ J2*np.sum((s*sx)**2+(s*sy)**2) 
    return H

def tmag(s):
    mx=np.sum(s[:,:,0])**2
    my=np.sum(s[:,:,1])**2
    mz=np.sum(s[:,:,2])**2
    return ((mx+my+mz)**0.5)/(L*L)

def magnet(s):
    return np.sum(s[:,:,2])/(L*L)

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
    return chi/(L*L)

def FFT(s):
    D=np.zeros((L,L)) #Euclidean Distance Matrix 
    for i in range(L):
        for j in range(L):
            D[i,j]=(min(i,L-i)**2+min(j,L-j)**2)**0.5
        
    X=s[:,:,2]-np.mean(s[:,:,2])
    W=np.abs(np.fft.fft2(X))
    #return np.roll(np.roll(W,L//2,axis=0),L//2,axis=1)
    lambd=np.sum(W*D)/np.sum(W)
    #peak=np.max(W)/np.sum(W)
    #peaks=np.sort(W,axis=None)[:-5:-1]/np.sum(W)
    #sigma=(np.sum(W*(D-lambd)**2)/np.sum(W))**0.5
    return lambd

A=[]
B=[]
#C=[]
#D=[]
E=[]

for K in range(15):
    for H in range(30):
        A.append(chirality(S[K,H,:,:,:]))
        B.append(magnet(S[K,H,:,:,:]))
        #np.savetxt(f,np.reshape(FFT(S[K,H,:,:,:]),(1,-1)),delimiter=',')
        #C.append(tmag(S[K,H,:,:,:]))
        #D.append(Ham(S[K,H,:,:,:]))
        E.append(FFT(S[K,H,:,:,:]))

np.savetxt('L_chi.csv',A,delimiter=',')
np.savetxt('L_mag.csv',B,delimiter=',')
#np.savetxt('D3_tmag.csv',C,delimiter=',')
#np.savetxt('D3_ham.csv',D,delimiter=',')
np.savetxt('L_peak.csv',E,delimiter=',')

#f.close()



