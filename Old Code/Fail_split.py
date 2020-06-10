import numpy as np   
import time

t=time.time()
heat3 = np.array(list(np.loadtxt("Fail.csv", delimiter=","))).astype("float32")

print(np.shape(heat3))

f=open('train.csv','ab')
g=open('eval.csv','ab')

heat3=np.reshape(heat3,(22,20,101,1728))

for B in range(22):
    for T in range(20):
        if (B>6 and B<9) or (B>13 and B<17): continue
        np.savetxt(f,heat3[B,T,:80,:],delimiter=',')
        np.savetxt(g,heat3[B,T,80:,:],delimiter=',')

print(time.time()-t)
f.close()
g.close()
