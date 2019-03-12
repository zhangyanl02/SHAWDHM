from netCDF4 import Dataset
from gisutil import *
import numpy as np

def GetRowCols(sourceDir,lon,lat):
      dims=np.shape(lon)
      BDfile=sourceDir+"BD.nc"
      BDgroup = Dataset(BDfile, "r")
      nlat=len(BDgroup.dimensions["lat"])
      nlon=len(BDgroup.dimensions["lon"])
      lat0=BDgroup.variables["lat"][0]
      latn=BDgroup.variables["lat"][nlat-1]
      lon0=BDgroup.variables["lon"][0]
      lonn=BDgroup.variables["lon"][nlon-1]
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      row=((lat[:,:]-lat0)//latcellsize).astype(int)
      col=((lon[:,:]-lon0)//loncellsize).astype(int)
      row[row[:,:]>nlat-1]=nlat-1
      col[col[:,:]>nlon-1]=nlon-1
      return row,col

def DetermineLayers(Depths,depth_in_org_dataset):
      depth_in_org_dataset[:]=depth_in_org_dataset[:]/100.0
      outputlayer=np.shape(Depths)[0]
      sourcelayer=np.shape(depth_in_org_dataset)[0]
      layer=np.zeros((outputlayer,2*sourcelayer+1),dtype=np.int)
      upper=np.zeros((outputlayer,1),dtype=np.float)
      lower=np.zeros((outputlayer,1),dtype=np.float)
      upper1=np.zeros((sourcelayer,1),dtype=np.float)
      lower1=np.zeros((sourcelayer,1),dtype=np.float)
      dz=np.zeros((outputlayer,sourcelayer),dtype=np.float)
      upper[0]=0.0
      for l in range(1,outputlayer):
            upper[l]=0.5*(Depths[l-1]+Depths[l])
            lower[l-1]=upper[l]
      lower[outputlayer-1]=2*Depths[outputlayer-1]-upper[outputlayer-1]
      upper1[0]=0.0
      lower1[0]=depth_in_org_dataset[0]
      for l in range(1,sourcelayer):
            upper1[l]=depth_in_org_dataset[l-1]
            lower1[l]=depth_in_org_dataset[l]
      for l in range(outputlayer):
            for d in range(sourcelayer):
                  if(upper1[d]>=upper[l] and upper1[d]<=lower[l]):
                        dz[l,d]=np.min([lower1[d],lower[l]])-upper1[d]
                  if(lower1[d]>=upper[l] and lower1[d]<=lower[l]):
                        dz[l,d]=lower1[d]-np.max([upper1[d],upper[l]])
                  if(upper1[d]<=upper[l] and lower1[d]>=lower[l]):
                        dz[l,d]=lower[l]-upper[l]
                  if(upper1[sourcelayer-2]<=upper[l] and lower1[sourcelayer-1]>=upper[l]):
                        if(lower[l]>=lower1[sourcelayer-1]):
                              dz[l,sourcelayer-1]=lower[l]-upper[l]
                  if(upper1[sourcelayer-1]<=upper[l] and upper1[sourcelayer-1]<lower[l]):
                        if(lower[l]>=lower1[sourcelayer-1]):
                              dz[l,sourcelayer-1]=lower[l]-upper[l]                        
            dz[l,:]=dz[l,:]/(lower[l]-upper[l])
            if(np.sum(dz[l,:])==0):
                  dz[l,sourcelayer-1]=1.0
            if(np.sum(dz[l,:])!=1):
                  print "error in determine layer in layer:",l
      return dz


def ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,varname,baseDepth):
      #var names in : BD,
      NCfile=sourceDir+varname+".nc"
      print varname
      NCgroup = Dataset(NCfile, "r")
      ndepth=len(NCgroup.dimensions["depth"])
      depth_in_org_dataset=NCgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      dims=np.shape(lon)
      extVar=np.zeros((ndepth,dims[0],dims[1]))
      tmpVar=np.zeros((Nlayer,dims[0],dims[1]))
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata1=NCgroup.variables[varname][:,minrow:maxrow+1,mincol:maxcol+1]
      tempdata=tempdata1[:,:].astype(float)
      NCgroup.close()
      for d in range(ndepth):
          tmp=tempdata[d,row[:,:]-minrow,col[:,:]-mincol]
          tmp[np.isnan(tmp[:,:])]=-9999
          tmp[tmp[:,:]==-999]=-9999
          extVar[d,:,:]=tmp[:,:]
          if d > 0:
             extVar[d,tmp[:,:]==-9999]=extVar[d-1,tmp[:,:]==-9999]
      if(varname=="BD"):
         extVar[extVar[:,:,:]!=-9999]=extVar[extVar[:,:,:]!=-9999]*1000.0
         extVar[extVar[:,:,:]<0.0]=1100.0
      if(varname=="SA"):
         extVar[extVar[:,:,:]<0.0]=33.0
      if(varname=="CL"):
         extVar[extVar[:,:,:]<0.0]=33.0
      if(varname=="K_SCH"):
         extVar[extVar[:,:,:]!=-9999]=extVar[extVar[:,:]!=-9999]/100.0/24.0/3600.0
         extVar[extVar[:,:,:]<0.0]=20.0/100.0/24.0/3600.0
      if(varname=="LAMBDA"):
         extVar[extVar[:,:]<0.0]=1/15.0
      if(varname=="POR"):
         extVar[extVar[:,:,:]<0.0]=0.35
      if(varname=="PSI_S"):
         extVar[extVar[:,:,:]!=-9999]=extVar[extVar[:,:]!=-9999]/100.0
         extVar[extVar[:,:,:]>=0.0]=-0.2
         extVar[extVar[:,:,:]==-9999]=-0.2
      if(varname=="SI"):
         extVar[extVar[:,:,:]<0.0]=34.0
      if(varname=="SOM"):
         extVar[extVar[:,:,:]<0.0]=0.0
      if(varname=="THSCH"):
         extVar[extVar[:,:,:]<=0.0]=0.38
      if(varname=="THR"):
         extVar[extVar[:,:,:]<=0.0]=0.05
      if(varname=="ALPHA"):
         extVar[extVar[:,:,:]<0.0]=0.03
      if(varname=="N"):
         extVar[extVar[:,:,:]<0.0]=1.2
      if(varname=="TH33"):
         extVar[extVar[:,:,:]<0.0]=0.28
      for i in range(dims[0]):
            for j in range(dims[1]):
                for nl in range(Nlayer):
                    tmpVar[nl,i,j]=np.inner(layerfactor[nl,:],extVar[:,i,j])
      if(varname=="POR"):
         tmpVar[np.asarray(Depths)>baseDepth,:,:]=0.08
      if(varname=="K_SCH"):
         tmpVar[np.asarray(Depths)>baseDepth,:,:]=2.0/100.0/24.0/3600.0
      if(varname=="THSCH"):
         tmpVar[np.asarray(Depths)>baseDepth,:,:]=0.06
      if(varname=="THR"):
         tmpVar[np.asarray(Depths)>baseDepth,:,:]=0.02
      if(varname=="BD"):
         tmpVar[np.asarray(Depths)>baseDepth,:,:]=1500.0
      if(varname=="FC"):
         tmpVar[np.asarray(Depths)>baseDepth,:,:]=0.05
      return tmpVar




      
def ExtractBDFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      BDfile=sourceDir+"BD.nc"
      BDgroup = Dataset(BDfile, "r")
      nlat=len(BDgroup.dimensions["lat"])
      nlon=len(BDgroup.dimensions["lon"])
      ndepth=len(BDgroup.dimensions["depth"])
      lat0=BDgroup.variables["lat"][0]
      latn=BDgroup.variables["lat"][nlat-1]
      lon0=BDgroup.variables["lon"][0]
      lonn=BDgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=BDgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extBD=np.zeros((ndepth,dims[0],dims[1]))
      tmpBD=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=BDgroup.variables["BD"][:,minrow:maxrow+1,mincol:maxcol+1]
      BDgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extBD[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extBD[d,i,j]) or extBD[d,i,j]==-999):
                              extBD[d,i,j]=-9999
                        if(extBD[d,i,j]!=-9999):
                              extBD[d,i,j]=extBD[d,i,j]*1000.0
                        else:
                              extBD[d,i,j]=extBD[d-1,i,j]
                        if(extBD[d,i,j]<=0.0):
                              extBD[d,i,j]=1100.00
                  for nl in range(Nlayer):
                        tmpBD[nl,i,j]=np.inner(layerfactor[nl,:],extBD[:,i,j])
      return tmpBD
      
def ExtractSAFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      SAfile=sourceDir+"SA.nc"
      SAgroup = Dataset(SAfile, "r")
      nlat=len(SAgroup.dimensions["lat"])
      nlon=len(SAgroup.dimensions["lon"])
      ndepth=len(SAgroup.dimensions["depth"])
      lat0=SAgroup.variables["lat"][0]
      latn=SAgroup.variables["lat"][nlat-1]
      lon0=SAgroup.variables["lon"][0]
      lonn=SAgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=SAgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extSA=np.zeros((ndepth,dims[0],dims[1]))
      tmpSA=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)    
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1  
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=SAgroup.variables["SA"][:,minrow:maxrow+1,mincol:maxcol+1]
      SAgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extSA[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extSA[d,i,j]) or extSA[d,i,j]==-999):
                              extSA[d,i,j]=-9999
                        if(extSA[d,i,j]==-9999):
                              extSA[d,i,j]=extSA[d-1,i,j]
                        if(extSA[d,i,j]<=0.0):
                              extSA[d,i,j]=33.00
                  for nl in range(Nlayer):
                        tmpSA[nl,i,j]=np.inner(layerfactor[nl,:],extSA[:,i,j])
      return tmpSA

def ExtractCLFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      CLfile=sourceDir+"CL.nc"
      CLgroup = Dataset(CLfile, "r")
      nlat=len(CLgroup.dimensions["lat"])
      nlon=len(CLgroup.dimensions["lon"])
      ndepth=len(CLgroup.dimensions["depth"])
      lat0=CLgroup.variables["lat"][0]
      latn=CLgroup.variables["lat"][nlat-1]
      lon0=CLgroup.variables["lon"][0]
      lonn=CLgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=CLgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extCL=np.zeros((ndepth,dims[0],dims[1]))
      tmpCL=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)   
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1   
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=CLgroup.variables["CL"][:,minrow:maxrow+1,mincol:maxcol+1]
      CLgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extCL[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extCL[d,i,j]) or extCL[d,i,j]==-999):
                              extCL[d,i,j]=-9999
                        if(extCL[d,i,j]==-9999):
                              extCL[d,i,j]=extCL[d-1,i,j]
                        if(extCL[d,i,j]<=0.0):
                              extCL[d,i,j]=33.00
                  for nl in range(Nlayer):
                        tmpCL[nl,i,j]=np.inner(layerfactor[nl,:],extCL[:,i,j])
      return tmpCL

def ExtractK_SCHFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      K_SCHfile=sourceDir+"K_SCH.nc"
      K_SCHgroup = Dataset(K_SCHfile, "r")
      nlat=len(K_SCHgroup.dimensions["lat"])
      nlon=len(K_SCHgroup.dimensions["lon"])
      ndepth=len(K_SCHgroup.dimensions["depth"])
      lat0=K_SCHgroup.variables["lat"][0]
      latn=K_SCHgroup.variables["lat"][nlat-1]
      lon0=K_SCHgroup.variables["lon"][0]
      lonn=K_SCHgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=K_SCHgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extK_SCH=np.zeros((ndepth,dims[0],dims[1]))
      tmpK_SCH=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)  
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1    
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=K_SCHgroup.variables["K_SCH"][:,minrow:maxrow+1,mincol:maxcol+1]
      K_SCHgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extK_SCH[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extK_SCH[d,i,j]) or extK_SCH[d,i,j]<=0.0):
                              extK_SCH[d,i,j]=-9999
                        if(extK_SCH[d,i,j]!=-9999):
                              extK_SCH[d,i,j]=extK_SCH[d,i,j]/100.0/24.0/3600.0
                        else:
                              extK_SCH[d,i,j]=extK_SCH[d-1,i,j]
                        if(extK_SCH[d,i,j]<=0.0):
                              extK_SCH[d,i,j]=20/100.0/24.0/3600.0
                  for nl in range(Nlayer):
                        tmpK_SCH[nl,i,j]=np.inner(layerfactor[nl,:],extK_SCH[:,i,j])
      return tmpK_SCH



def ExtractLAMBDFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      LAMBDfile=sourceDir+"LAMBDA.nc"
      LAMBDgroup = Dataset(LAMBDfile, "r")
      nlat=len(LAMBDgroup.dimensions["lat"])
      nlon=len(LAMBDgroup.dimensions["lon"])
      ndepth=len(LAMBDgroup.dimensions["depth"])
      lat0=LAMBDgroup.variables["lat"][0]
      latn=LAMBDgroup.variables["lat"][nlat-1]
      lon0=LAMBDgroup.variables["lon"][0]
      lonn=LAMBDgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=LAMBDgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extLAMBD=np.zeros((ndepth,dims[0],dims[1]))
      tmpLAMBD=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)   
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1   
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=LAMBDgroup.variables["LAMBDA"][:,minrow:maxrow+1,mincol:maxcol+1]
      LAMBDgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extLAMBD[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extLAMBD[d,i,j]) or extLAMBD[d,i,j]==-999):
                              extLAMBD[d,i,j]=-9999
                        if(extLAMBD[d,i,j]==-9999):
                              extLAMBD[d,i,j]=extLAMBD[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpLAMBD[nl,i,j]=np.inner(layerfactor[nl,:],extLAMBD[:,i,j])
                        if(tmpLAMBD[nl,i,j]==0):
                                tmpLAMBD[nl,i,j]=1/15.0
      return tmpLAMBD

def ExtractPORFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      PORfile=sourceDir+"POR.nc"
      PORgroup = Dataset(PORfile, "r")
      nlat=len(PORgroup.dimensions["lat"])
      nlon=len(PORgroup.dimensions["lon"])
      ndepth=len(PORgroup.dimensions["depth"])
      lat0=PORgroup.variables["lat"][0]
      latn=PORgroup.variables["lat"][nlat-1]
      lon0=PORgroup.variables["lon"][0]
      lonn=PORgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=PORgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extPOR=np.zeros((ndepth,dims[0],dims[1]))
      tmpPOR=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)  
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1    
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=PORgroup.variables["POR"][:,minrow:maxrow+1,mincol:maxcol+1]
      PORgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extPOR[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extPOR[d,i,j]) or extPOR[d,i,j]==-999):
                              extPOR[d,i,j]=-9999
                        if(extPOR[d,i,j]==-9999):
                              extPOR[d,i,j]=extPOR[d-1,i,j]
                        if(extPOR[d,i,j]<=0.0):
                              extPOR[d,i,j]=0.35
                        tmpPOR[nl,i,j]=np.inner(layerfactor[nl,:],extPOR[:,i,j])
      return tmpPOR

def ExtractPSI_SFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      PSI_Sfile=sourceDir+"PSI_S.nc"
      PSI_Sgroup = Dataset(PSI_Sfile, "r")
      nlat=len(PSI_Sgroup.dimensions["lat"])
      nlon=len(PSI_Sgroup.dimensions["lon"])
      ndepth=len(PSI_Sgroup.dimensions["depth"])
      lat0=PSI_Sgroup.variables["lat"][0]
      latn=PSI_Sgroup.variables["lat"][nlat-1]
      lon0=PSI_Sgroup.variables["lon"][0]
      lonn=PSI_Sgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=PSI_Sgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extPSI_S=np.zeros((ndepth,dims[0],dims[1]))
      tmpPSI_S=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)    
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1  
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=PSI_Sgroup.variables["PSI_S"][:,minrow:maxrow+1,mincol:maxcol+1]
      PSI_Sgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extPSI_S[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extPSI_S[d,i,j]) or extPSI_S[d,i,j]==-999):
                              extPSI_S[d,i,j]=-9999
                        if(extPSI_S[d,i,j]!=-9999):
                              extPSI_S[d,i,j]=extPSI_S[d,i,j]/100.0
                        else:
                              extPSI_S[d,i,j]=extPSI_S[d-1,i,j]
                        if(extPSI_S[d,i,j]==0.0):
                              extPSI_S[d,i,j]=-0.2
                        if(extPSI_S[d,i,j]==-9999):
                              extPSI_S[d,i,j]=-0.2
                  for nl in range(Nlayer):
                        tmpPSI_S[nl,i,j]=np.inner(layerfactor[nl,:],extPSI_S[:,i,j])
      return tmpPSI_S


def ExtractSIFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      SIfile=sourceDir+"SI.nc"
      SIgroup = Dataset(SIfile, "r")
      nlat=len(SIgroup.dimensions["lat"])
      nlon=len(SIgroup.dimensions["lon"])
      ndepth=len(SIgroup.dimensions["depth"])
      lat0=SIgroup.variables["lat"][0]
      latn=SIgroup.variables["lat"][nlat-1]
      lon0=SIgroup.variables["lon"][0]
      lonn=SIgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=SIgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extSI=np.zeros((ndepth,dims[0],dims[1]))
      tmpSI=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)  
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1    
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=SIgroup.variables["SI"][:,minrow:maxrow+1,mincol:maxcol+1]
      SIgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extSI[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extSI[d,i,j]) or extSI[d,i,j]==-999):
                              extSI[d,i,j]=-9999
                        if(extSI[d,i,j]==-9999):
                              extSI[d,i,j]=extSI[d-1,i,j]
                        if(extSI[d,i,j]<=0.0):
                              extSI[d,i,j]=34
                  for nl in range(Nlayer):
                        tmpSI[nl,i,j]=np.inner(layerfactor[nl,:],extSI[:,i,j])
      return tmpSI

def ExtractSOMFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      SOMfile=sourceDir+"SOM.nc"
      SOMgroup = Dataset(SOMfile, "r")
      nlat=len(SOMgroup.dimensions["lat"])
      nlon=len(SOMgroup.dimensions["lon"])
      ndepth=len(SOMgroup.dimensions["depth"])
      lat0=SOMgroup.variables["lat"][0]
      latn=SOMgroup.variables["lat"][nlat-1]
      lon0=SOMgroup.variables["lon"][0]
      lonn=SOMgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=SOMgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extSOM=np.zeros((ndepth,dims[0],dims[1]))
      tmpSOM=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)   
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1   
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=SOMgroup.variables["SOM"][:,minrow:maxrow+1,mincol:maxcol+1]
      SOMgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extSOM[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extSOM[d,i,j]) or extSOM[d,i,j]==-999):
                              extSOM[d,i,j]=-9999
                        if(extSOM[d,i,j]==-9999):
                              extSOM[d,i,j]=extSOM[d-1,i,j]
                  for nl in range(Nlayer):
                        tmpSOM[nl,i,j]=np.inner(layerfactor[nl,:],extSOM[:,i,j])
      return tmpSOM


def ExtractTHETASFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      THSCHfile=sourceDir+"THSCH.nc"
      THSCHgroup = Dataset(THSCHfile, "r")
      nlat=len(THSCHgroup.dimensions["lat"])
      nlon=len(THSCHgroup.dimensions["lon"])
      ndepth=len(THSCHgroup.dimensions["depth"])
      lat0=THSCHgroup.variables["lat"][0]
      latn=THSCHgroup.variables["lat"][nlat-1]
      lon0=THSCHgroup.variables["lon"][0]
      lonn=THSCHgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=THSCHgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extTHSCH=np.zeros((ndepth,dims[0],dims[1]))
      tmpTHSCH=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)   
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1   
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=THSCHgroup.variables["THSCH"][:,minrow:maxrow+1,mincol:maxcol+1]
      THSCHgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extTHSCH[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extTHSCH[d,i,j]) or extTHSCH[d,i,j]==-999):
                              extTHSCH[d,i,j]=-9999
                        if(extTHSCH[d,i,j]==-9999):
                              extTHSCH[d,i,j]=extTHSCH[d-1,i,j]
                        if(extTHSCH[d,i,j]<=0.0):
                              extTHSCH[d,i,j]=0.38
                  for nl in range(Nlayer):
                        tmpTHSCH[nl,i,j]=np.inner(layerfactor[nl,:],extTHSCH[:,i,j])
      return tmpTHSCH

def ExtractTHETARFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      THRfile=sourceDir+"THR.nc"
      THRgroup = Dataset(THRfile, "r")
      nlat=len(THRgroup.dimensions["lat"])
      nlon=len(THRgroup.dimensions["lon"])
      ndepth=len(THRgroup.dimensions["depth"])
      lat0=THRgroup.variables["lat"][0]
      latn=THRgroup.variables["lat"][nlat-1]
      lon0=THRgroup.variables["lon"][0]
      lonn=THRgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=THRgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extTHR=np.zeros((ndepth,dims[0],dims[1]))
      tmpTHR=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)    
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1  
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=THRgroup.variables["THR"][:,minrow:maxrow+1,mincol:maxcol+1]
      THRgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extTHR[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extTHR[d,i,j]) or extTHR[d,i,j]==-999):
                              extTHR[d,i,j]=-9999
                        if(extTHR[d,i,j]==-9999):
                              extTHR[d,i,j]=extTHR[d-1,i,j]
                        if(extTHR[d,i,j]<=0.0):
                              extTHR[d,i,j]=0.05
                  for nl in range(Nlayer):
                        tmpTHR[nl,i,j]=np.inner(layerfactor[nl,:],extTHR[:,i,j])
      return tmpTHR

def ExtractFCFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      FCfile=sourceDir+"TH33.nc"
      FCgroup = Dataset(FCfile, "r")
      nlat=len(FCgroup.dimensions["lat"])
      nlon=len(FCgroup.dimensions["lon"])
      ndepth=len(FCgroup.dimensions["depth"])
      lat0=FCgroup.variables["lat"][0]
      latn=FCgroup.variables["lat"][nlat-1]
      lon0=FCgroup.variables["lon"][0]
      lonn=FCgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=FCgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extFC=np.zeros((ndepth,dims[0],dims[1]))
      tmpFC=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize) 
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1     
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=FCgroup.variables["TH33"][:,minrow:maxrow+1,mincol:maxcol+1]
      FCgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extFC[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extFC[d,i,j]) or extFC[d,i,j]==-999):
                              extFC[d,i,j]=-9999
                        if(extFC[d,i,j]==-9999):
                              extFC[d,i,j]=extFC[d-1,i,j]
                        if(extFC[d,i,j]<=0.0):
                              extFC[d,i,j]=0.28
                  for nl in range(Nlayer):
                        tmpFC[nl,i,j]=np.inner(layerfactor[nl,:],extFC[:,i,j])
      return tmpFC
  
  
def ExtractAlphaFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      Alphafile=sourceDir+"ALPHA.nc"
      Alphagroup = Dataset(Alphafile, "r")
      nlat=len(Alphagroup.dimensions["lat"])
      nlon=len(Alphagroup.dimensions["lon"])
      ndepth=len(Alphagroup.dimensions["depth"])
      lat0=Alphagroup.variables["lat"][0]
      latn=Alphagroup.variables["lat"][nlat-1]
      lon0=Alphagroup.variables["lon"][0]
      lonn=Alphagroup.variables["lon"][nlon-1]
      depth_in_org_dataset=Alphagroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extAlpha=np.zeros((ndepth,dims[0],dims[1]))
      tmpAlpha=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)  
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1    
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=Alphagroup.variables["ALPHA"][:,minrow:maxrow+1,mincol:maxcol+1]
      Alphagroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extAlpha[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extAlpha[d,i,j]) or extAlpha[d,i,j]==-999):
                              extAlpha[d,i,j]=-9999
                        if(extAlpha[d,i,j]==-9999):
                              extAlpha[d,i,j]=extAlpha[d-1,i,j]
                        if(extAlpha[d,i,j]<=0.0):
                              extAlpha[d,i,j]=0.03
                  for nl in range(Nlayer):
                        tmpAlpha[nl,i,j]=np.inner(layerfactor[nl,:],extAlpha[:,i,j])
      return tmpAlpha

def ExtractSoilNFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon):
      SoilNfile=sourceDir+"N.nc"
      SoilNgroup = Dataset(SoilNfile, "r")
      nlat=len(SoilNgroup.dimensions["lat"])
      nlon=len(SoilNgroup.dimensions["lon"])
      ndepth=len(SoilNgroup.dimensions["depth"])
      lat0=SoilNgroup.variables["lat"][0]
      latn=SoilNgroup.variables["lat"][nlat-1]
      lon0=SoilNgroup.variables["lon"][0]
      lonn=SoilNgroup.variables["lon"][nlon-1]
      depth_in_org_dataset=SoilNgroup.variables["depth"][:]
      layerfactor=DetermineLayers(Depths,depth_in_org_dataset)
      latcellsize=(latn-lat0)/(nlat-1)
      loncellsize=(lonn-lon0)/(nlon-1)
      dims=np.shape(lon)
      extSoilN=np.zeros((ndepth,dims[0],dims[1]))
      tmpSoilN=np.zeros((Nlayer,dims[0],dims[1]))
      row=np.zeros((dims[0],dims[1]),dtype=np.int)
      col=np.zeros((dims[0],dims[1]),dtype=np.int)
      for i in range(dims[0]):
            for j in range(dims[1]):
                  row[i,j]=int((lat[i,j]-lat0)/latcellsize)
                  col[i,j]=int((lon[i,j]-lon0)/loncellsize)   
                  if(row[i,j]>nlat-1):
                     row[i,j]=nlat-1
                  if(col[i,j]>nlon-1):
                     col[i,j]=nlon-1   
      minrow=np.min(row)
      maxrow=np.max(row)
      mincol=np.min(col)
      maxcol=np.max(col)
      tempdata=SoilNgroup.variables["N"][:,minrow:maxrow+1,mincol:maxcol+1]
      SoilNgroup.close()
      for i in range(dims[0]):
            for j in range(dims[1]):
                  for d in range(ndepth):
                        extSoilN[d,i,j]=tempdata[d,row[i,j]-minrow,col[i,j]-mincol]
                        if(np.isnan(extSoilN[d,i,j]) or extSoilN[d,i,j]==-999):
                              extSoilN[d,i,j]=-9999
                        if(extSoilN[d,i,j]==-9999):
                              extSoilN[d,i,j]=extSoilN[d-1,i,j]
                        if(extSoilN[d,i,j]<=0.0):
                              extSoilN[d,i,j]=1.2
                  for nl in range(Nlayer):
                        tmpSoilN[nl,i,j]=np.inner(layerfactor[nl,:],extSoilN[:,i,j])
      return tmpSoilN

def ExtractSoilDataFromDaiYongJiu(sourceDir,OutputFile,Depths,lat,lon,baseDepth):
      dims=np.shape(lon)
      Nlayer=len(Depths)
      row,col=GetRowCols(sourceDir,lon,lat)
      BD=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"BD",baseDepth)
      SI=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"SI",baseDepth)
      SA=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"SA",baseDepth)
      CL=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"CL",baseDepth)
      SOM=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"SOM",baseDepth)
      PSI_S=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"PSI_S",baseDepth)
      LAMBD=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"LAMBDA",baseDepth)
      THETAS=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"THSCH",baseDepth)
      KS=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"K_SCH",baseDepth)
      THR=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"THR",baseDepth)
      FC=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"TH33",baseDepth)
      Alpha=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"ALPHA",baseDepth)
      soilN=ExtractSoilVarFromDaiYongjiu(sourceDir,Nlayer,Depths,lat,lon,col,row,"N",baseDepth)

      print np.shape(lon)
      soilpara = Dataset(OutputFile,"w",format="NETCDF3_CLASSIC")
      lat   = soilpara.createDimension("NY", dims[0])
      lon   = soilpara.createDimension("NX", dims[1])
      depth = soilpara.createDimension("depth", Nlayer)
      #dims =("depth","NY","NX",)
      dims =("depth","NY","NX",)

      
      THETASs= soilpara.createVariable("thetas","f4",dims)
      #THETASs.missing_value="-9999"
      THETASs.units="cm3 cm-3"
      THETASs.long_name="Saturated water content"
      THETASs[:,:,:]=THETAS[:,:,:]
      
      THRs= soilpara.createVariable("thr","f4",dims)
      #THRs.missing_value="-9999"
      THRs.units="cm3 cm-3"
      THRs.long_name="Residual water content"
      THRs[:,:,:]=THR[:,:,:]
      
      FCs= soilpara.createVariable("fc","f4",dims)
      #FCs.missing_value="-9999"
      FCs.units="cm3 cm-3"
      FCs.long_name="Field water capacity"
      FCs[:,:,:]=FC[:,:,:]
      
      Alphas= soilpara.createVariable("alpha","f4",dims)
      #Alphas.missing_value="-9999"
      Alphas.units="-"
      Alphas.long_name="Alpha"
      Alphas[:,:,:]=Alpha[:,:,:]
      
      SoilNs= soilpara.createVariable("n","f4",dims)
      #SoilNs.missing_value="-9999"
      SoilNs.units="-"
      SoilNs.long_name="N"
      SoilNs[:,:,:]=soilN[:,:,:]
      
      KSs= soilpara.createVariable("ks","f4",dims)
      #KSs.missing_value="-9999"
      KSs.units="m/s"
      KSs.long_name="Saturate hydraulic conductivity"
      KSs[:,:,:]=KS[:,:,:]    
      
      BDs= soilpara.createVariable("bd","f4",dims)
      #BDs.missing_value="-9999"
      BDs.units="kg/m3"
      BDs.long_name="Bulk Density"
      BDs[:,:,:]=BD[:,:,:]
      
      SAs= soilpara.createVariable("sa","f4",dims)
      #SAs.missing_value="-9999"
      SAs.units="0-100 %"
      SAs.long_name="Sand content"
      SAs[:,:,:]=SA[:,:,:]
      
      SIs= soilpara.createVariable("si","f4",dims)
      #SIs.missing_value="-9999"
      SIs.units="0-100 %"
      SIs.long_name="Silt content"
      SIs[:,:,:]=SI[:,:,:]   
      
      CLs= soilpara.createVariable("cl","f4",dims)
      #CLs.missing_value="-9999"
      CLs.units="0-100 %"
      CLs.long_name="Clay content"
      CLs[:,:,:]=CL[:,:,:]
      
      SOMs= soilpara.createVariable("om","f4",dims)
      #SOMs.missing_value="-9999"
      SOMs.units="0-100 %"
      SOMs.long_name="Organic content"
      SOMs[:,:,:]=SOM[:,:,:]

      PSI_Ss= soilpara.createVariable("psi_s","f4",dims)
      #PSI_Ss.missing_value="-9999"
      PSI_Ss.units="m"
      PSI_Ss.long_name="Saturated capillary potential"
      PSI_Ss[:,:,:]=PSI_S[:,:,:]    
      
      LAMBDs= soilpara.createVariable("b","f4",dims)
      #LAMBDs.missing_value="-9999"
      LAMBDs.units="-"
      LAMBDs.long_name="Pore size distribution index"
      LAMBDs[:,:,:]=1.0/LAMBD[:,:,:]

      ZSs= soilpara.createVariable("zs","f4",("depth",))
      #ZSs.missing_value="-9999"
      ZSs.units="m/s"
      ZSs.long_name="Saturate hydraulic conductivity"
      ZSs[:]=Depths[:]                          
      soilpara.close()

