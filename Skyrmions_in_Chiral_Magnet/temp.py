import numpy as np

L=24

S1=np.loadtxt('L3.csv',delimiter=',')
S1=np.reshape(S1,(15,20,L*L*3))
S2=np.loadtxt('L4.csv',delimiter=',')
S2=np.reshape(S2,(15,5,L*L*3))
S3=np.loadtxt('L5.csv',delimiter=',')
S3=np.reshape(S3,(15,5,L*L*3))

f=open('L_temp.csv','ab')

for i in range(15):
    np.savetxt(f,np.reshape(S1[i,:,:],(-1,L*L*3)),delimiter=',')
    np.savetxt(f,np.reshape(S2[i,:,:],(-1,L*L*3)),delimiter=',')
    np.savetxt(f,np.reshape(S3[i,:,:],(-1,L*L*3)),delimiter=',')

f.close()
