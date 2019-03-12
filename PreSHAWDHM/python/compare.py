import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

extVar=np.loadtxt("extVar.txt")
extVar1=np.loadtxt("extVar1.txt")
#extVar2=np.loadtxt("extVar2.txt")

maskfile=open("./dongbei/dongbei5km.asc",'r')
for i in range(6):
	maskfile.readline()
mask=np.loadtxt(maskfile)

extVar[extVar==-9999]=np.float("NaN")
extVar[mask==-9999]=np.float("NaN")
extVar1[mask==-9999]=np.float("NaN")
extVar1[extVar1<=-9999]=np.float("NaN")
#extVar2[mask==-9999]=np.float("NaN")
dif=extVar-extVar1
#dif[mask==-9999]=np.float("NaN")
print np.max(dif[~np.isnan(dif)])
fig = plt.figure(1, figsize=(10,10))
imgplot = plt.imshow(dif[100:150,100:150],interpolation='nearest',cmap=cm.hot)
plt.colorbar(imgplot,shrink=0.4)
fig.savefig("./dif.png",dpi=100,bbox_inches='tight')
plt.close()

fig = plt.figure(1, figsize=(10,10))
imgplot = plt.imshow(extVar,interpolation='nearest',cmap=cm.hot)
plt.colorbar(imgplot,shrink=0.4)
fig.savefig("./extVar.png",dpi=100,bbox_inches='tight')
plt.close()

fig = plt.figure(1, figsize=(10,10))
imgplot = plt.imshow(extVar1-extVar,interpolation='nearest',cmap=cm.hot)
plt.colorbar(imgplot,shrink=0.4)
fig.savefig("./dif.png",dpi=100,bbox_inches='tight')
plt.close()

fig = plt.figure(1, figsize=(10,10))
imgplot = plt.imshow(extVar1,interpolation='nearest',cmap=cm.hot)
plt.colorbar(imgplot,shrink=0.4)
fig.savefig("./extVar1.png",dpi=100,bbox_inches='tight')
plt.close()



