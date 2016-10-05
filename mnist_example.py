from byd import *
from keras.datasets import mnist
import numpy as np

mnist_data = mnist.load_data()
train_images = mnist_data[0][0]
train_labels = mnist_data[0][1]
test_images = mnist_data[1][0]
test_labels = mnist_data[1][1]

#View training digits
#change frames using "," and "."
imagesc(train_images)

#View all training digits as 150 frames of 20 by 20,
#change colors using "C", use "Return" to check value of strange parts of the data
def create_20_20(x):
    return np.concatenate([np.concatenate(x[20*i:20*i+20],1) for i in range(20)],0)
imagesc([create_20_20(train_images[400*i:400*i+400]) for i in range(150)])

#Change the width using "[" and "]" to identify runs "1,2,...,9"
#Looks like the runs occur on a width of 109.
imagesc(blockup(train_labels[57000:],100))

#Do PCA on the 1s
ones = train_images[train_labels == 1]
ones_R = ones.reshape(len(ones),28*28)
[U,D,V] = np.linalg.svd(ones_R,0)

#Spin it around
#color some middle guys and bring back to python
cl = np.array(ggobi(U[:,:20])  )
from Collections import Counter
Counter(cl)

#View middle guys
imagesc(ones[cl==1])

#Look at right singular vectors
imagesc([V[i].reshape(28,28) for i in range(20)])

#View k dimensional Reconstructions as frames
imagesc([np.dot(np.dot(U[0,:k],np.diag(D[:k])),V[:k]).reshape(28,28) for k in range(32)])

#Find a funny one
imagesc(ones) #---> 5

#View k dimensional Reconstructions as frames for the funny one
imagesc([np.dot(np.dot(U[5,:k],np.diag(D[:k])),V[:k]).reshape(28,28) for k in range(32)])

#Plot each one in 3d using first three left singular vectors
def create_frame(L,x,y):
    x = (x - x.min()) / (x.max() - x.min())
    y = (y - y.min()) / (y.max() - y.min())
    Z = np.zeros([1028,1028])
    for i in range(len(L)):
        ii = int(x[i]*1000)        
        jj = int(y[i]*1000)
        Z[ii:ii+28,jj:jj+28] += L[i]
    return Z

from math import sin, cos, pi
ONES_3D = [create_frame(ones[:1000],U[:,1][:1000],cos(2*pi*t/100)*U[:,0][:1000] + sin(2*pi*t/100)*U[:,2][:1000]) for t in range(100)]
imagesc(ONES_3D)

#Compare PCA to simple autoencoder from keras blog (hey don't judge me!)
import keras
autoencoder = keras.models.load_model('mnist_ae.h5')
W = autoencoder.get_weights()
IN_W = W[0].T.reshape(32,28,28)
OUT_W = W[2].reshape(32,28,28)

#First we look at the weights, you can see aritfacts of the ReLU v.s. the sigmoid here
imagesc([np.concatenate([IN_W[i],OUT_W[i]],0) for i in range(32)])

#How about PCA, they look much for Fourier like, (PCA does that)
R = train_images.reshape(len(train_images),28*28)
[U,D,V] = np.linalg.svd(R,0)
imagesc([blockup(V[i],28) for i in range(32)])

#See how the neural nets energy is spread in PCA space
#We actually see a dropoff at 32
imagesc(np.corrcoef(V[:32],OUT_W.reshape(32,28*28)))
imagesc(np.corrcoef(V[:128],OUT_W.reshape(32,28*28)))





