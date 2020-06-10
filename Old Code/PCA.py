import numpy as np
import matplotlib.pyplot as plt

X = np.array(list(np.loadtxt("train.csv", delimiter=","))).astype("float")
X_=np.matrix.transpose(X)

lambd,W=np.linalg.eig(np.matmul(X_,X))

plt.plot(np.sort(W[0,:]))
plt.show()