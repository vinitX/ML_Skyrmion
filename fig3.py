import numpy as np
import matplotlib.pyplot as plt
'''
a = np.array(list(np.loadtxt("y_act.csv", delimiter=","))).astype("float32")
p = np.array(list(np.loadtxt("y_pred.csv", delimiter=","))).astype("float32")
p1 = np.array(list(np.loadtxt("y_pred1.csv", delimiter=","))).astype("float32")
p2 = np.array(list(np.loadtxt("y_pred2.csv", delimiter=","))).astype("float32")

a=np.reshape(a,(20,20,100,4))
p=np.reshape(p,(20,20,100,4))
p1=np.reshape(p1,(20,20,100,4))
p2=np.reshape(p2,(20,20,100,4))
'''

plt.plot(a[:,1,50,0],'b')
plt.plot(p[:,1,50,0],'g')
plt.plot(p1[:,1,50,0],'r')
plt.plot(p2[:,1,50,0],'y')

plt.plot(0.2+a[:,5,50,0],'b')
plt.plot(0.2+p[:,5,50,0],'g')
plt.plot(0.2+p1[:,5,50,0],'r')
plt.plot(0.2+p2[:,5,50,0],'y')

plt.plot(0.5+a[:,10,50,0],'b')
plt.plot(0.5+p[:,10,50,0],'g')
plt.plot(0.5+p1[:,10,50,0],'r')
plt.plot(0.5+p2[:,10,50,0],'y')

plt.plot(1+a[:,15,50,0],'b')
plt.plot(1+p[:,15,50,0],'g')
plt.plot(1+p1[:,15,50,0],'r')
plt.plot(1+p2[:,15,50,0],'y')

plt.xlabel('B')
plt.ylabel('Chirality')
plt.show()