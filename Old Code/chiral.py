import numpy as np   

L=24

def chirality(theta,phi,k,l):
    s=np.zeros((L,L,3))
    s[:,:,0]=np.sin(theta)*np.cos(phi)
    s[:,:,1]=np.sin(theta)*np.sin(phi)
    s[:,:,2]=np.cos(theta)
    
    t=np.zeros((L+2,L+2,3))
    t[1:L+1,1:L+1,:]=s
    chiral=np.zeros((L,L))
    for i in range(L):
        for j in range(L):
            ti=t[i+1,j+1,:]
            tx=t[i+2,j+1,:]
            tx_=t[i,j+1,:]
            ty=t[i+1,j+2,:]
            ty_=t[i+1,j,:]
            chiral[i,j]=np.dot(ti,np.cross(tx,ty))+np.dot(ti,np.cross(tx_,ty))+np.dot(ti,np.cross(tx,ty_))+np.dot(ti,np.cross(tx_,ty_))
    chi=np.reshape(chiral,(1,L*L))
    filename_chi="chi"+str(k)+"_"+str(l)+".txt"
    np.savetxt(filename_chi,chi)
            
for i in range(10):
    #if i==3 or i==7: continue
    for j in range(0,501):
        filename_a="theta"+str(i)+"_"+str(j)+".txt"
        filename_b="phi"+str(i)+"_"+str(j)+".txt"
        a=np.loadtxt(filename_a)
        b=np.loadtxt(filename_b)
        a=np.reshape(a,(L,L))
        b=np.reshape(b,(L,L))
        chirality(a,b,i,j)
