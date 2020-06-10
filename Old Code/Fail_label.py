import numpy as np   
import time
import pandas as pd

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

t=time.time()
fail = pd.read_csv("Fail2.csv", Header=None).astype("float32")
#test = pd.read_csv("test.csv", Header=None).astype("float32")
print(time.time()-t)
t=time.time()
print(np.shape(fail))
#print(np.shape(test))

fl=open('train_label.csv','ab')
#gl=open('eval_label.csv','ab')
#hl=open('test_label.csv','ab')

al=[]
#bl=[]
#heat3=np.reshape(heat3,(20,20,200,576,3))

for B in range(600):
    if B<200: label=0
    elif B<400: label=1
    else: label=2
    al.append([label]*50*11)
    #bl.append([label]*21)
'''

cl=np.zeros((40000,4))
test=np.reshape(test,(20,20,100,576,3))

for B in range(20):
    for T in range(20):
        for j in range(100):
                cl[((B*20)+T)*100+j,0]=chirality(test[B,T,j,:,:])
                cl[((B*20)+T)*100+j,1]=sum(test[B,T,j,:,2])/(L*L)
                cl[((B*20)+T)*100+j,2]=B/20
                cl[((B*20)+T)*100+j,3]=(20-T)/20
'''
np.savetxt(fl,al,delimiter=',')
#np.savetxt(gl,bl,delimiter=',')
#np.savetxt(hl,cl,delimiter=',')
print(np.shape(al))#,np.shape(bl))
#print(np.shape(cl))
print(time.time()-t)

fl.close()
#gl.close()
#hl.close()