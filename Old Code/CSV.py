import numpy as np   

heat4 = np.array(list(np.loadtxt("heat_4.csv", delimiter=","))).astype("float32")
heat3 = np.array(list(np.loadtxt("heat_3.csv", delimiter=","))).astype("float32")

f=open('train3.csv','ab')
g=open('eval3.csv','ab')
#h=open('test.csv','ab')
fl=open('train_label3.csv','ab')
gl=open('eval_label3.csv','ab')
al=[]
bl=[]
print(np.shape(heat3),np.shape(heat4))
for i in range(0,180):
    np.savetxt(f,heat3[i*50:i*50+40,:],delimiter=',')
    np.savetxt(g,heat3[i*50+40:i*50+50,:],delimiter=',')
    
    al.append(heat4[i*50:i*50+40])
    bl.append(heat4[i*50+40:i*50+50])
    
np.savetxt(fl,al,delimiter=',')
np.savetxt(gl,bl,delimiter=',')
print(np.shape(al),np.shape(bl))
f.close()   
g.close()
fl.close()
gl.close()