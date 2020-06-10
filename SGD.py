import numpy as np
import matplotlib.pyplot as plt
import random
import scipy.stats as sp

a=np.loadtxt('y_hat2_1.csv',delimiter=',')
a=np.reshape(a[:,1],(20,20))
x=np.arange(0,1,0.05)
y=a[:10,19]
x=x[:10]

m=np.sum(x*y)/np.sum(x*x)
plt.scatter(x,y)
plt.plot(x,m*x)
plt.show()
print(m)


x=np.array([0,0.5,1,1.5,2])
y=np.array([1.038,1.125,1.130,1.2217,1.349])
plt.scatter(x,y)
m,c,a,b,d=sp.linregress(x,y)
print(m,c)
plt.plot(x,m*x+c)

x=np.array([0,0.5,1,1.5,2])
y=np.array([1.038,1.199,1.323,1.495,1.8095])
plt.scatter(x,y)
m,c,a,b,d=sp.linregress(x,y)
print(m,c)
plt.plot(x,m*x+c)

x=np.array([0,0.5,1,1.5,2])
y=np.array([1.038,1.117,1.144,1.249,1.403])
plt.scatter(x,y)
m,c,a,b,d=sp.linregress(x,y)
print(m,c)
plt.plot(x,m*x+c)

x=np.array([0,0.5,1,1.5,2])
y=np.array([1.038,1.199,1.339,1.524,1.859])
plt.scatter(x,y)
m,c,a,b,d=sp.linregress(x,y)
print(m,c)
plt.plot(x,m*x+c)

'''
def SGD(m,alpha):
    hist_m=[]
    hist_L=[]
    for i in range(100):
        d=(1+m*m)**0.5
        L=np.sum(np.abs(y-m*x))/(1+m*m)
        
        dL=0
        for j in range(20):
            if (y[j]-m*x[j])>0:
                dL=-d*x[j]-np.abs(y[j]-m*x[j])*m/d
            else:
                dL=d*x[j]-np.abs(y[j]-m*x[j])*m/d
        dL=dL/(1+m*m)
        m=m-alpha*dL
        hist_L.append(L)
        hist_m.append(m)
        
    plt.plot(hist_m)
    plt.show()
    plt.plot(hist_L)
    plt.show()
    return m
        
def Newton(m):
    for i in range(10):
        d=(1+m*m)**0.5
        L=np.sum((np.abs(y-m*x)))/(1+m*m)
        
        dL=0
        ddL=0
        for j in range(20):
            if (y[j]-m*x[j])>0:
                dL=-d*x[j]-np.abs(y[j]-m*x[j])*m/d
                ddL=(d**3)*y[j]-(x[j]+m*y[j])*(3/2)*d
            else:
                dL=d*x[j]-np.abs(y[j]-m*x[j])*m/d
                ddL=-(d**3)*y[j]-(x[j]+m*y[j])*(3/2)*d
        dL=dL/(1+m*m)
        ddL=ddL/(d**6)
        
        if ddL==0:
            return m
        
        m=m-dL/ddL
        print("{0:.2f}".format(L),end=' ')
        print(m)
    return m

def monte(m,e):
    hist_m=[]
    hist_L=[]
    beta=20
    for T in range(41):
        temp=2*np.exp(-0.1*T)
        beta=1/temp
        for i in range(1000):
            m_=m+e*(random.random()-0.5)
            L=np.sum(np.abs(y-m*x))/(1+m*m)
            L_=np.sum(np.abs(y-m_*x))/(1+m_*m_)
            dL=L_-L
        
            if dL<0 or np.exp(-beta*dL)>random.random(): m=m_
            hist_L.append(L)
            hist_m.append(m)
    plt.plot(hist_m)
    plt.show()
    plt.plot(hist_L)
    plt.show()
    return m

x=10
alpha=0.1
hist_x=[]
hist_L=[]
for i in range(100):
    L=(x-2.5)**2
    dL=2*(x-2.5)
    x=x-alpha*L
    hist_L.append(L)
    hist_x.append(x)
plt.plot(hist_x)
plt.show()
plt.plot(np.log10(hist_L))
plt.show()
'''

