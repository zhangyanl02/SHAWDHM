from netCDF4 import Dataset
import numpy as np
from gisutil import *
import os
from DateUtil import *
from operator import itemgetter, attrgetter 

def SpecificHumidity2RelativeHumidity(SH,T,P):
      RH=0.263*P*SH/(np.exp(17.67*(T-273.15)/(T-29.65)))
      return RH

def ExtVarFromGLDAS(SrcDir,OutDir,year,day,lat,lon,InterpolationType,ncformat,compress,row,col,cellsize,colrowbound,timezone):
    nlatOut=np.shape(lat)[0]
    nlonOut=np.shape(lat)[1]
    solartmp=np.zeros(((24,nlatOut,nlonOut)),dtype=np.float)
    solartmp[:,:,:]=-9999
    if(timezone==0):
        date2=GetDateFromJulianDay(year,day)
        filename2=SrcDir+str(year)+"/"+"/RAD_"+str(year)+str(100+date2.month)[1:]+str(100+date2.day)[1:]+".nc4"
        ncfile2 = Dataset(filename2.strip(), "r") 
        solar_inst=ncfile2.variables["solar"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        solartmp[:,:,:]=solar_inst[:,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
        ncfile2.close()
    elif(timezone>0):
        preday=day-1
        if(preday<=0):
           preyear=year-1
           if(Leap(preyear)):
              preday=366
           else:
              preday=365
        else:
           preyear=year
        date1=GetDateFromJulianDay(preyear,preday)
        date2=GetDateFromJulianDay(year,day)
        filename1=SrcDir+str(year)+"/"+"/RAD_"+str(preyear)+str(100+date1.month)[1:]+str(100+date1.day)[1:]+".nc4"
        filename2=SrcDir+str(year)+"/"+"/RAD_"+str(year)+str(100+date2.month)[1:]+str(100+date2.day)[1:]+".nc4"
        ncfile2 = Dataset(filename2.strip(), "r") 
        solar_inst=ncfile2.variables["solar"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        solartmp[timezone:24,:,:]=solar_inst[0:24-timezone,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
        ncfile2.close()
	if(os.path.exists(filename1)):
            ncfile1 = Dataset(filename1.strip(), "r") 
            solar_inst=ncfile1.variables["solar"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
            solartmp[0:timezone,:,:]=solar_inst[24-timezone:24,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
            ncfile1.close()
        else:
            solartmp[0:timezone,:,:]=-9999


    #solartmp[solartmp==-9999]=np.float("NaN")
    solartmp[solartmp==-9999]=0.0
    #for r in range(nlatOut):
    #    for c in range(nlonOut):
    #       for h in range(5,22):
    #           if(solartmp[h,r,c]==0):
    #               if(solartmp[h-1,r,c]!=0 and solartmp[h+1,r,c]!=0):
    #                  solartmp[h,r,c]=(solartmp[h-1,r,c]+solartmp[h+1,r,c])/2.0
    #               if(solartmp[h-1,r,c]!=0 and solartmp[h+1,r,c]==0 and solartmp[h+2,r,c]!=0):
    #                  solartmp[h,r,c]=solartmp[h-1,r,c]+1/3.0*(solartmp[h+2,r,c]-solartmp[h-1,r,c])
            

    outputncfilename=OutDir+"/"+str(year)+"-"+str(year*1000+day)[-3:]+".nc"
    #print outputncfilename
    nc4file = Dataset(outputncfilename, "w", format=ncformat)
    NY     = nc4file.createDimension("NY", nlatOut)
    NX     = nc4file.createDimension("NX", nlonOut)
    time    = nc4file.createDimension("time",24)
    lats= nc4file.createVariable("lat","f4",("NY","NX",))
    lats.units="degrees_north"
    lats.long_name="Latitude"
    lats[:,:]=lat[:,:]
    
    lons= nc4file.createVariable("lon","f4",("NY","NX",))
    lons.units="degrees_east"
    lons.long_name="Longitude"
    lons[:,:]=lon[:,:]

    solar= nc4file.createVariable("solar","f4",("time","NY","NX",),zlib=compress)
    solar.units="Wm-2"
    solar[:,:,:]=solartmp
    nc4file.close()


def GetRowColInGLDAS2(SrcDir,startYear,lat,lon):
      nlat=np.shape(lat)[0]
      nlon=np.shape(lat)[1]
      ncfilename=SrcDir+str(2007)+"/RAD_20070101.nc4"
      ncfile = Dataset(ncfilename.strip(), "r")
      nlat_in=len(ncfile.dimensions["NY"])
      nlon_in=len(ncfile.dimensions["NC"])
      lat0=ncfile.variables["lat"][0,0]
      latn=ncfile.variables["lat"][nlat_in-1,0]
      lon0=ncfile.variables["lon"][0,0]
      lonn=ncfile.variables["lon"][0,nlon_in-1]
      cellsize=[0.0,0.0]
      cellsize[0]=(latn-lat0)/(nlat_in-1)
      cellsize[1]=(lonn-lon0)/(nlon_in-1)
      row=np.zeros((nlat,nlon),dtype=np.int)
      col=np.zeros((nlat,nlon),dtype=np.int)
      row=((lat[:,:]-lat0+0.5*cellsize[0])//cellsize[0]).astype(int)
      col=((lon[:,:]-lon0+0.5*cellsize[1])//cellsize[1]).astype(int)
      colrowbound=[0,0,0,0]
      colrowbound[0]=np.max([np.min(row)-1,0])
      colrowbound[1]=np.min([np.max(row)+1,nlat_in-1])
      colrowbound[2]=np.max([np.min(col)-1,0])
      colrowbound[3]=np.min([np.max(col)+1,nlon_in-1])
      tempdata=ncfile.variables["solar"][0,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
      latul=ncfile.variables["lon"][colrowbound[0],colrowbound[2]]
      lonul=ncfile.variables["lat"][colrowbound[0],colrowbound[2]]
      ncfile.close()
      return latul,lonul,colrowbound,col,row,cellsize,tempdata


def ExtractSolarFromTang(SrcDir,OutputDir,startYear,EndYear,lat,lon,interpolationType,timezone,compress,ncformat,DT):
      nlat=np.shape(lat)[0]
      nlon=np.shape(lat)[1]
      latul,lonul,colrowbound,col,row,cellsize,tempdata=GetRowColInGLDAS2(SrcDir,startYear,lat,lon)
      for y in range(startYear,EndYear):
            Dir=OutputDir+str(y)
            if not os.path.exists(Dir):
                  os.makedirs(Dir)
            if(Leap(y)):
                ndays=366
            else:
                ndays=365
            for d in range(1,ndays+1):
		print y,d
                ExtVarFromGLDAS(SrcDir,Dir,y,d,lat,lon,interpolationType,ncformat,compress,row,col,cellsize,colrowbound,timezone)






