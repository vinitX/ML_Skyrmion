import numpy as np
import matplotlib.pyplot as plt

heat = np.array(list(np.loadtxt("H_label.csv", delimiter=","))).astype("float32")
plt.imshow(np.reshape(heat,(20,20,4))[:,:,0])
plt.show()
plt.imshow(np.reshape(heat,(20,20,4))[:,:,1])
plt.show()

heat = np.array(list(np.loadtxt("H2CG_label.csv", delimiter=","))).astype("float32")
plt.imshow(np.reshape(heat,(20,20,4))[:,:,0])
plt.show()
plt.imshow(np.reshape(heat,(20,20,4))[:,:,1])
plt.show()

heat = np.array(list(np.loadtxt("H3CG_label.csv", delimiter=","))).astype("float32")
print(np.shape(heat))
plt.imshow(np.reshape(heat,(20,20,4))[:,:,0])
plt.show()
plt.imshow(np.reshape(heat,(20,20,4))[:,:,1])
plt.show()