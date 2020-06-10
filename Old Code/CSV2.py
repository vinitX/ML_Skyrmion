import numpy as np   
f3=open('train3.csv','ab')
g3=open('eval3.csv','ab')
h3=open('test3.csv','ab')

train_data = np.array(list(np.loadtxt("train2.csv", delimiter=","))).astype("float32")
eval_data = np.array(list(np.loadtxt("eval2.csv", delimiter=","))).astype("float32")

x_train = train_data.reshape(-1, 24, 24, 2)
theta=x_train[:,:,:,0]
print(np.amax(theta))
phi=x_train[:,:,:,1]
print(np.amax(phi))
x=np.multiply(np.sin(theta),np.cos(phi))
y=np.multiply(np.sin(theta),np.sin(phi))
z=np.cos(theta)
'''
x_val = eval_data.reshape(-1, 24, 24, 2)


print(np.shape(c))
f.close()   
g.close()
h.close()
'''