import numpy as np   
f3=open('train1.csv','ab')
g3=open('eval1.csv','ab')
h3=open('test1.csv','ab')

train_data = np.array(list(np.loadtxt("train2.csv", delimiter=","))).astype("float32")
eval_data = np.array(list(np.loadtxt("eval2.csv", delimiter=","))).astype("float32")
test_data = np.array(list(np.loadtxt("test4.csv", delimiter=","))).astype("float32")

x_train = train_data.reshape(-1, 24, 24, 2)
theta=x_train[:,:,:,0]
print(np.amax(theta))
z=np.cos(theta)
z=np.reshape(z,(-1,576))
np.savetxt(f3,z,delimeter=',')

x_val = eval_data.reshape(-1, 24, 24, 2)
theta=x_val[:,:,:,0]
print(np.amax(theta))
z=np.cos(theta)
z=np.reshape(z,(-1,576))
np.savetxt(g3,z,delimeter=',')

x_test = test_data.reshape(-1, 24, 24, 2)
theta=x_test[:,:,:,0]
print(np.amax(theta))
z=np.cos(theta)
z=np.reshape(z,(-1,576))
np.savetxt(h3,z,delimeter=',')

print(np.shape(z))
f3.close()   
g3.close()
h3.close()