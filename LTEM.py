import numpy as np   
import time
'''
t=time.time()
heat3 = np.array(list(np.loadtxt("heat3.csv", delimiter=","))).astype("float32")

print(np.shape(heat3))

f=open('train2.csv','ab')
g=open('eval2.csv','ab')

heat3=np.reshape(heat3,(400,200,576,3))

for i in range(400):
    np.savetxt(f,np.reshape(heat3[i,:160,:,:2],(160,1152)),delimiter=',')
    np.savetxt(g,np.reshape(heat3[i,160:,:,:2],(40,1152)),delimiter=',')

print(time.time()-t)
f.close()
g.close()
'''
for K in [1]:
    for p in [0.5]:
        t=time.time()
        test = np.array(list(np.loadtxt("dis"+str(K)+"_"+str(p)+".csv", delimiter=","))).astype("float32")
        print(np.shape(test))
        f=open("dis"+str(K)+"_"+str(p)+"L.csv",'ab')

        test=np.reshape(test,(-1,200,576,3))

        for i in range(np.shape(test)[0]):
            np.savetxt(f,np.reshape(test[i,:,:,:2],(200,1152)),delimiter=',')
        
        print(time.time()-t)
        f.close()
