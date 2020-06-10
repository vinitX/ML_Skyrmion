import numpy as np

L=24
heat1 = np.array(list(np.loadtxt("heat1.csv", delimiter=","))).astype("float32")
#f=open("heat3.csv",'ab')
print(np.shape(heat1))
'''heat2=np.reshape(heat2,(20*9*50,24,24,2))

for i in range(np.shape(heat2)[0]):
    s=np.zeros((L,L,3))
    theta=heat2[i,:,:,0]
    print(np.max(theta))
    phi=heat2[i,:,:,1]
    print(np.max(phi))
    s[:,:,0]=np.sin(theta)*np.cos(phi)
    s[:,:,1]=np.sin(theta)*np.sin(phi)
    s[:,:,2]=np.cos(theta)
    s=np.reshape(s,(1,L*L*3))
    #np.savetxt(f,s,delimiter=',')'''