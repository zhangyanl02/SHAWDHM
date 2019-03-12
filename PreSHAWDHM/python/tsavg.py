from netCDF4 import Dataset
import numpy as np
from gisutil import *
import os

def ExtTsavg(SourceDir,StartYear,EndYear,nr,nc):
      tmp=np.zeros((24,nr,nc))
      tsavg=np.zeros((nr,nc))
      nd=0
      for y in range(StartYear,EndYear+1):
            for jd in range(1,366):
                  ncfilename=SourceDir+"/"+str(y)+"/"+str(y)+"-"+str(y*1000+jd)[-3:]+".nc"
                  ncfile = Dataset(ncfilename, "r", format="NETCDF4")
                  tmp=tmp+ncfile.variables["temp"][:,:,:]
                  nd=nd+1
      tmp=tmp/nd
      for r in range(nr):
            for c in range(nc):
                  tsavg[r,c]=np.mean(tmp[:,r,c])
      return tsavg
      
      
def ExtTsavgFromCalculatedSoil(SourceDir,StartYear,EndYear,nr,nc):
      tmp=np.zeros((24,nr,nc))
      tsavg=np.zeros((nr,nc))
      nd=0
      for y in range(StartYear,EndYear+1):
            for jd in range(1,366):
                  ncfilename=SourceDir+"/"+str(y)+"-"+str(y*1000+jd)[-3:]+".nc"
                  ncfile = Dataset(ncfilename, "r", format="NETCDF4")
                  var=ncfile.variables["TSDT2d"][:,9,:,:]
                  if(np.mean(var)<100.0):
                        tmp=tmp+var
                        nd=nd+1
      tmp=tmp/nd
      for r in range(nr):
            for c in range(nc):
                  tsavg[r,c]=np.mean(tmp[:,r,c])
      return tsavg     

