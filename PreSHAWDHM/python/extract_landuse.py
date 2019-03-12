from netCDF4 import Dataset
import numpy as np

def Extract_Landuse(SourceNCname,lat,lon):
      NCgroup = Dataset(SourceNCname, "r")
      nlat=len(NCgroup.dimensions["lat"])
      nlon=len(NCgroup.dimensions["lon"])
      lat0=NCgroup.variables["lat"][0]
      latn=NCgroup.variables["lat"][nlat-1]
      lon0=NCgroup.variables["lon"][0]
      lonn=NCgroup.variables["lon"][nlon-1]
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      print nlat,nlon 
      print lat0,latn
      print lon0,lonn
      print  latcellsize ,loncellsize  
      dims=np.shape(lon)
      ExtLanduse=np.zeros((dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)      
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      print minrow,maxrow,mincol,maxcol
      tempdata=NCgroup.variables["landuse"][minrow:maxrow+1,mincol:maxcol+1]
      NCgroup.close()
      
      for i in range(dims[0]):
            for j in range(dims[1]):
                  ExtLanduse[i,j]=tempdata[row[i,j]-minrow,col[i,j]-mincol]
      return ExtLanduse
