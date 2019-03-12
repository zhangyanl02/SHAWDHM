from netCDF4 import Dataset
import numpy as np

def Extract_LAI(SourceNCname,lat,lon):
      NCgroup = Dataset(SourceNCname, "r")
      nlat=len(NCgroup.dimensions["lat"])
      nlon=len(NCgroup.dimensions["lon"])
      ntime=len(NCgroup.dimensions["time"])
      lat0=NCgroup.variables["lat"][0]
      latn=NCgroup.variables["lat"][nlat-1]
      lon0=NCgroup.variables["lon"][0]
      lonn=NCgroup.variables["lon"][nlon-1]
      time_in_source=NCgroup.variables["time"][:]
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      ExtLAI=np.zeros((ntime,dims[0],dims[1]))
      tmpLAI=np.zeros((366,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      print "lat0,lon0",lat0,lon0
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      #print "minrow,maxrow",minrow,maxrow
      tempdata=NCgroup.variables["lai"][:,minrow:maxrow+1,mincol:maxcol+1]
      NCgroup.close()
      
      for i in range(dims[0]):
            for j in range(dims[1]):
                  ExtLAI[:,i,j]=tempdata[:,row[i,j]-minrow,col[i,j]-mincol]      
      for j in range(ntime-1):
            firstday=time_in_source[j]
            secondday=time_in_source[j+1]
            for day in range(firstday,secondday):
                  tmpLAI[day-1,:,:]=ExtLAI[j,:,:]+(ExtLAI[j+1,:,:]-ExtLAI[j,:,:])/(secondday-firstday)*(day-firstday)
      for day in range(361,367):
            tmpLAI[day-1,:,:]=ExtLAI[ntime-1,:,:]+(ExtLAI[ntime-1,:,:]-ExtLAI[ntime-2,:,:])/(time_in_source[ntime-1]-time_in_source[ntime-2])*(day-time_in_source[ntime-1])
      return tmpLAI
