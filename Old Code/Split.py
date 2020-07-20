import numpy as np   
import time

t=time.time()
heat3 = np.array(list(np.loadtxt("heat3.csv", delimiter=","))).astype("float32")

print(np.shape(heat3))

f=open('train.csv','ab')
g=open('eval.csv','ab')

heat3=np.reshape(heat3,(400,200,576,3))

for i in range(400):
    np.savetxt(f,np.reshape(heat3[i,:160,:,:],(160,1728)),delimiter=',')
    np.savetxt(g,np.reshape(heat3[i,160:,:,:],(40,1728)),delimiter=',')

print(time.time()-t)
f.close()
g.close()