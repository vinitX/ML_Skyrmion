import numpy as np   
import time

L=24

def chi(a,b,c):
    t1=np.sum(a*np.cross(b,c,axisa=2,axisb=2),axis=2)
    t2=1+np.sum(a*b+b*c+c*a,axis=2)
    #print(np.sum(t2<0),end=', ')
    #print(t1[np.nonzero(t2<0)],t2[np.nonzero(t2<0)])
    #print(np.sum(t1[np.nonzero(t2<0)]>0)-np.sum(t1[np.nonzero(t2<0)]<0),end=' ')
    
    #Dealing with zero denominator
    eps=1e-20
    index=np.nonzero(abs(t2)<eps)
    t2[index]=t2[index]+2*eps

    chi=np.sum(np.arctan(t1/t2))#+np.pi*(np.sum(t1[np.nonzero(t2<0)]>0)-np.sum(t1[np.nonzero(t2<0)]<0))
    return chi

def chirality(s):
    s=np.reshape(s,(L,L,3))
    
    sx=np.roll(s,-1,axis=0)
    sx_=np.roll(s,1,axis=0)
    sy=np.roll(s,-1,axis=1)
    sy_=np.roll(s,1,axis=1)
    
    return (chi(s,sx,sy)+chi(s,sx_,sy_)+chi(s,sy,sx_)+chi(s,sy_,sx))/(4*np.pi)/L

'''
t=time.time()
heat3 = np.array(list(np.loadtxt("heat3.csv", delimiter=","))).astype("float32")
print(time.time()-t)
t=time.time()
print(np.shape(heat3))

fl=open('train_label.csv','ab')
gl=open('eval_label.csv','ab')

al=np.zeros((64000,4))
bl=np.zeros((16000,4))
heat3=np.reshape(heat3,(20,20,200,576,3))

for B in range(20):
    for T in range(20):
        for j in range(200):
            if j<160:
                al[((B*20)+T)*160+j,0]=chirality(heat3[B,T,j,:,:])
                al[((B*20)+T)*160+j,1]=sum(heat3[B,T,j,:,2])/(L*L)
                al[((B*20)+T)*160+j,2]=B/20
                al[((B*20)+T)*160+j,3]=(20-T)/20
            else:
                bl[((B*20)+T)*40+(j-160),0]=chirality(heat3[B,T,j,:,:])
                bl[((B*20)+T)*40+(j-160),1]=sum(heat3[B,T,j,:,2])/(L*L)   
                bl[((B*20)+T)*40+(j-160),2]=B/20
                bl[((B*20)+T)*40+(j-160),3]=(20-T)/20

np.savetxt(fl,al,delimiter=',')
np.savetxt(gl,bl,delimiter=',')

print(time.time()-t)

fl.close()
gl.close()
'''
t=time.time()
test = np.array(list(np.loadtxt("train.csv", delimiter=","))).astype("float32")
print(time.time()-t)
print(np.shape(test))

m=20
n=20
o=160
cl=np.zeros((m,n,o,4))
test=np.reshape(test,(m,n,o,L*L,3))

for B in range(m):
    for T in range(n):
        for j in range(o):
            cl[B,T,j,0]=chirality(test[B,T,j,:,:])
            cl[B,T,j,1]=sum(test[B,T,j,:,2])/(L*L)
            cl[B,T,j,2]=B/m
            cl[B,T,j,3]=(n-T)/n

np.savetxt('train_label.csv',np.reshape(cl,(-1,4)),delimiter=',')
