from netCDF4 import Dataset
import numpy as np
from gisutil import *
import os
from operator import itemgetter, attrgetter 

import matplotlib.pyplot as plt
import matplotlib.cm as cm

def SpecificHumidity2RelativeHumidity(SH,T,P):
      RH=0.263*P*SH/(np.exp(17.67*(T-273.15)/(T-29.65)))
      return RH

def ExtVarFromYangKun(ncfilename,varname,lat,lon,InterpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights):
      ncfile = Dataset(ncfilename.strip(), "r")
      ntime=len(ncfile.dimensions["time"])
      times=ncfile.variables["time"][:]
      dims=np.shape(lon)
      extVar=np.zeros((ntime,dims[0],dims[1]))
      #extVar1=np.zeros((ntime,dims[0],dims[1]))
      tempdata=ncfile.variables[varname][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
      latul=ncfile.variables["lat"][colrowbound[0]]
      lonul=ncfile.variables["lon"][colrowbound[2]]
      ncfile.close()

      #extVar1[:,:,:]=tempdata[:,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
      #extVar1[extVar1<-9999.0]=np.float("nan")
      #np.savetxt("extVar1.txt",extVar1[0,:,:])
      # No interpolation
      if(InterpolationType==0):
          extVar[:,:,:]=tempdata[:,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
          extVar[extVar<-9999.0]=np.float("nan")
      else:
      # Inverse distance interpolation       
          extVar[:,:,:]=IDWWeights[0,:,:]*tempdata[:,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*tempdata[:,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*tempdata[:,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*tempdata[:,NearPointRow[3,:,:],NearPointCol[3,:,:]]
          #for r in range(dims[0]):
          #      for c in range(dims[1]):
          #            # Inverse distance interpolation       
          #            extVar1[:,r,c]=IDWWeights[0,r,c]*tempdata[:,NearPointRow[0,r,c],NearPointCol[0,r,c]]+IDWWeights[1,r,c]*tempdata[:,NearPointRow[1,r,c],NearPointCol[1,r,c]]+IDWWeights[2,r,c]*tempdata[:,NearPointRow[2,r,c],NearPointCol[2,r,c]]+IDWWeights[3,r,c]*tempdata[:,NearPointRow[3,r,c],NearPointCol[3,r,c]]
      #np.savetxt("extVar.txt",extVar[0,:,:])
      #np.savetxt("extVar.txt",extVar[0,:,:])
      #stop
      return times,extVar

def ExtVarFromYangKun2(ncfilename,varname,lat,lon,InterpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights):
      ncfile = Dataset(ncfilename.strip(), "r")
      times=ncfile.variables["time"][0]
      dims=np.shape(lon)
      extVar=np.zeros((dims[0],dims[1]))
      tempdata=ncfile.variables[varname][0,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
      latul=ncfile.variables["lat"][colrowbound[0]]
      lonul=ncfile.variables["lon"][colrowbound[2]]
      ncfile.close()      
      # No interpolation
      if(InterpolationType==0):
          extVar[:,:]=tempdata[row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
          extVar[extVar<-9999.0]=np.float("nan")
      else:
      # Inverse distance interpolation       
          extVar[:,:]=IDWWeights[0,:,:]*tempdata[NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*tempdata[NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*tempdata[NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*tempdata[NearPointRow[3,:,:],NearPointCol[3,:,:]]             
      return times,extVar

def ConvertThreeHourlyVar2Hourly(varsin):
      ntime=np.shape(varsin)[0]
      nlat=np.shape(varsin)[1]
      nlon=np.shape(varsin)[2]
      HourlyVar=np.zeros((24,nlat,nlon))
      for t in range(ntime-1):
            HourlyVar[3*t,:,:]=varsin[t,:,:]
            HourlyVar[3*t+1,:,:]=varsin[t,:,:]*2.0/3.0+varsin[t+1,:,:]*1.0/3.0            
            HourlyVar[3*t+2,:,:]=varsin[t,:,:]*1.0/3.0+varsin[t+1,:,:]*2.0/3.0
      if(ntime==8):
      	for i in range(3):
            HourlyVar[3*(ntime-1)+i,:,:]=varsin[ntime-1,:,:]
      return HourlyVar

def GetHourlyTimes(ThreeHourlyTimes):
      ntime=np.shape(ThreeHourlyTimes)[0]
      timeHourly=np.zeros((24),dtype=np.int)
      for t in range(ntime-1):
            timeHourly[3*t]=ThreeHourlyTimes[t]
            timeHourly[3*t+1]=ThreeHourlyTimes[t]+1                  
            timeHourly[3*t+2]=ThreeHourlyTimes[t]+2   
      if(ntime==8):
        for i in range(3):
            timeHourly[3*(ntime-1)+i]=ThreeHourlyTimes[ntime-1]+i
      return timeHourly  


def GetVariable(SrcDir,VariableName,datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights):
            ncfilename=SrcDir+VariableName+"/"+VariableName+"_ITPCAS-CMFD_V0106_B-01_03hr_010deg_"+datestr+".nc"
            timesout,VarData=ExtVarFromYangKun(ncfilename,VariableName,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
            return VarData,timesout


def GetNextMonthVariable(SrcDir,VariableName,datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights):
            ncfilename=SrcDir+VariableName+"/"+VariableName+"_ITPCAS-CMFD_V0106_B-01_03hr_010deg_"+datestr+".nc"
            timesout,VarData=ExtVarFromYangKun2(ncfilename,VariableName,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
            return VarData,timesout


def GetRowColInYangkun(SrcDir,startYear,lat,lon):
      nlat=np.shape(lat)[0]
      nlon=np.shape(lat)[1]
      ncfilename=SrcDir+"temp"+"/"+"temp"+"_ITPCAS-CMFD_V0106_B-01_03hr_010deg_"+str(startYear*100+1)+".nc"
      ncfile = Dataset(ncfilename.strip(), "r")
      nlat_in=len(ncfile.dimensions["lat"])
      nlon_in=len(ncfile.dimensions["lon"])
      lat0=ncfile.variables["lat"][0]
      latn=ncfile.variables["lat"][nlat_in-1]
      lon0=ncfile.variables["lon"][0]
      lonn=ncfile.variables["lon"][nlon_in-1]      
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

      tempdata=ncfile.variables["temp"][0,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
      latul=ncfile.variables["lat"][colrowbound[0]]
      lonul=ncfile.variables["lon"][colrowbound[2]]
      ncfile.close()
      return latul,lonul,colrowbound,col,row,cellsize,tempdata

def GetNearPoints(lat,lon,latul,lonul,colrowbound,cellsize,tempdata):
      nlat=np.shape(lat)[0]
      nlon=np.shape(lat)[1]
      NearPointRow=np.zeros((4,nlat,nlon),dtype=np.int) #0 for INROW, 1 for INROW2, 2 for INCOL, 3 for INCOL2
      NearPointCol=np.zeros((4,nlat,nlon),dtype=np.int) #0 for INROW, 1 for INROW2, 2 for INCOL, 3 for INCOL2
      NearPointDis=np.zeros((4,nlat,nlon))
      IDWWeights=np.zeros((4,nlat,nlon))
      for r in range(nlat):
          for c in range(nlon):
                        IX=lon[r,c]
                        IY=lat[r,c]
                        NearPointCol[0,r,c]=np.int((IX-lonul+0.5*cellsize[1])/cellsize[1])
                        NearPointRow[0,r,c]=np.int((IY-latul+0.5*cellsize[0])/cellsize[0])
                        distance=np.zeros((9,3))
                        count=0
                        for i in range(-1,2):
                           for j in range(-1,2):
                               rr=NearPointRow[0,r,c]+i
                               cc=NearPointCol[0,r,c]+j
                               if(rr>=0 and rr<=colrowbound[1]-colrowbound[0] and cc>=0 and cc<=colrowbound[3]-colrowbound[2]):
                                   #distance[count,0]=((rr+0.5)*cellsize[0]+latul-IY)**2.0+((cc+0.5)*cellsize[1]+lonul-IX)**2.0
                                   distance[count,0]=((rr)*cellsize[0]+latul-IY)**2.0+((cc)*cellsize[1]+lonul-IX)**2.0
                                   if(np.isnan(np.float(tempdata[rr,cc]))):
                                      distance[count,0]=9999.0
                                   distance[count,1]=rr
                                   distance[count,2]=cc
                                   count=count+1
                        sortedDis=np.asarray(sorted(distance[:,0:count],key=itemgetter(0)))
                        NearPointRow[:,r,c]=sortedDis[0:4,1]
                        NearPointCol[:,r,c]=sortedDis[0:4,2]
                        NearPointDis[:,r,c]=sortedDis[0:4,0]
                        if(sortedDis[0,0]<=0.00001):
                              IDWWeights[0,r,c]=1.0
                              IDWWeights[1:,r,c]=0.0
                        else:
                              IDWWeights[:,r,c]=1.0/NearPointDis[:,r,c]/(np.sum(1.0/NearPointDis[:,r,c]))
      return NearPointRow[0:4,:,:],NearPointCol[0:4,:,:],IDWWeights[0:4,:,:]

def ExtractHourlyDataFromYangKunDataset(SrcDir,ForingOutputDir,startYear,EndYear,lat,lon,interpolationType,timezone,compress):
      nlat=np.shape(lat)[0]
      nlon=np.shape(lat)[1]
      ## calculate the cellsize, column and row numbers, and the min and max in col and row numbers
      latul,lonul,colrowbound,col,row,cellsize,tempdata=GetRowColInYangkun(SrcDir,startYear,lat,lon)
      NearPointRow,NearPointCol,IDWWeights=GetNearPoints(lat,lon,latul,lonul,colrowbound,cellsize,tempdata)

      
      for y in range(startYear,EndYear):
            jd=1
            Dir=ForingOutputDir+str(y)
            if not os.path.exists(Dir):
                  os.makedirs(Dir)
            for m in range(1,13):
                  datestr=str(y*100+m)
                  tempsThreeHour,timeThreeHour=GetVariable(SrcDir,"temp",datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
                  precsThreeHour,timeThreeHour=GetVariable(SrcDir,"prec",datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
                  lradsThreeHour,timeThreeHour=GetVariable(SrcDir,"lrad",datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
                  sradsThreeHour,timeThreeHour=GetVariable(SrcDir,"srad",datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
                  pressThreeHour,timeThreeHour=GetVariable(SrcDir,"pres",datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
                  shumsThreeHour,timeThreeHour=GetVariable(SrcDir,"shum",datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
                  windsThreeHour,timeThreeHour=GetVariable(SrcDir,"wind",datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
                  rhsThreeHour=SpecificHumidity2RelativeHumidity(shumsThreeHour,tempsThreeHour,pressThreeHour)

                  ndays=np.int(np.shape(timeThreeHour)[0]/8)
                  for d in range(ndays):
                        print y,m,d
                        outputncfilename=Dir+"/"+str(y)+"-"+str(y*1000+jd)[-3:]+".nc"
                        nc4file = Dataset(outputncfilename, "w", format="NETCDF3_CLASSIC")
                        NY     = nc4file.createDimension("NY", nlat)
                        NX     = nc4file.createDimension("NX", nlon)
                        time    = nc4file.createDimension("time",24)
      
                        lats= nc4file.createVariable("lats","f4",("NY","NX",))
                        lats.units="degrees_north"
                        lats.long_name="Latitude"
                        lats[:,:]=lat[:,:]
            
                        lons= nc4file.createVariable("lons","f4",("NY","NX",))
                        lons.units="degrees_east"
                        lons.long_name="Longitude"
                        lons[:,:]=lon[:,:]

                        times= nc4file.createVariable("time","i4",("time",))
                        times.units="UTC times since 1900-01-01 00:00:0.0"
                        times.long_name="Time"

                        prec= nc4file.createVariable("prec","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #prec.scale_factor=0.0025
                        #prec.add_offset=50.0
                        prec.missing_value=-32767
                        prec.units="mm hr-1"
                        prec.long_name="Precipitation rate"

                        temp= nc4file.createVariable("temp","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #temp.scale_factor=0.01
                        #temp.add_offset=0.0
                        temp.missing_value=-32767
                        temp.units="K"
                        temp.long_name="Near surface air temperature"

                        lrad= nc4file.createVariable("lrad","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #lrad.scale_factor=0.25
                        #lrad.add_offset=685.0
                        lrad.missing_value=-32767
                        lrad.units="W m-2"
                        lrad.long_name="Surface downward longwave radiation"
                  
                        srad= nc4file.createVariable("srad","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #srad.scale_factor=0.25
                        #srad.add_offset=685.0
                        srad.missing_value=-32767
                        srad.units="W m-2"
                        srad.long_name="Surface downward shortwave radiation"
                  
                        pres= nc4file.createVariable("pres","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #pres.scale_factor=2.0
                        #pres.add_offset=63500.
                        pres.missing_value=-32767
                        pres.units="Pa"
                        pres.long_name="Near surface air pressure"                  

                        shum= nc4file.createVariable("shum","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #shum.scale_factor=0.000001
                        #shum.add_offset=0.025
                        shum.missing_value=-32767
                        shum.units="kg kg-1"
                        shum.long_name="Near surface air specific humidity"
                        
                        wind= nc4file.createVariable("wind","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #wind.scale_factor=0.002
                        #wind.add_offset=60.0
                        wind.missing_value=-32767
                        wind.units="m s-1"
                        wind.long_name="Near surface wind speed"

                        RH= nc4file.createVariable("rhum","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #RH.scale_factor=0.01            
                        #RH.missing_value=-32767
                        RH.units="%"
                        RH.long_name="Near surface Relative Humidity"
                        if(d<ndays-1):
		                times[:]   =GetHourlyTimes(timeThreeHour[d*8:(d+1)*8+1])[:]
		                prec[:,:,:]=ConvertThreeHourlyVar2Hourly(precsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
		                temp[:,:,:]=ConvertThreeHourlyVar2Hourly(tempsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]-273.15
		                lrad[:,:,:]=ConvertThreeHourlyVar2Hourly(lradsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
		                srad[:,:,:]=ConvertThreeHourlyVar2Hourly(sradsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
		                pres[:,:,:]=ConvertThreeHourlyVar2Hourly(pressThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
		                shum[:,:,:]=ConvertThreeHourlyVar2Hourly(shumsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
		                wind[:,:,:]=ConvertThreeHourlyVar2Hourly(windsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
		                RH[:,:,:]  =ConvertThreeHourlyVar2Hourly(rhsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]/100.0 
                        else:
                                m1=m+1
                                y1=y
                                if(m1>12):
                                   y1=y+1
                                   m1=1
                                if(y1==2016):
				        times[:]   =GetHourlyTimes(timeThreeHour[d*8:(d+1)*8+1])[:]
				        prec[:,:,:]=ConvertThreeHourlyVar2Hourly(precsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
				        temp[:,:,:]=ConvertThreeHourlyVar2Hourly(tempsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]-273.15
				        lrad[:,:,:]=ConvertThreeHourlyVar2Hourly(lradsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
				        srad[:,:,:]=ConvertThreeHourlyVar2Hourly(sradsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
				        pres[:,:,:]=ConvertThreeHourlyVar2Hourly(pressThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
				        shum[:,:,:]=ConvertThreeHourlyVar2Hourly(shumsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
				        wind[:,:,:]=ConvertThreeHourlyVar2Hourly(windsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
				        RH[:,:,:]  =ConvertThreeHourlyVar2Hourly(rhsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]/100.0 
                                else:
				        datestr1=str(y1*100+m1)
				        precsThreeHourNextMonth,timeThreeHourNextMonth=GetNextMonthVariable(SrcDir,"prec",datestr1,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
				        tempsThreeHourNextMonth,timeThreeHourNextMonth=GetNextMonthVariable(SrcDir,"temp",datestr1,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
				        lradsThreeHourNextMonth,timeThreeHourNextMonth=GetNextMonthVariable(SrcDir,"lrad",datestr1,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
				        sradsThreeHourNextMonth,timeThreeHourNextMonth=GetNextMonthVariable(SrcDir,"srad",datestr1,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
				        pressThreeHourNextMonth,timeThreeHourNextMonth=GetNextMonthVariable(SrcDir,"pres",datestr1,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
				        shumsThreeHourNextMonth,timeThreeHourNextMonth=GetNextMonthVariable(SrcDir,"shum",datestr1,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
				        windsThreeHourNextMonth,timeThreeHourNextMonth=GetNextMonthVariable(SrcDir,"wind",datestr1,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
				        rhsThreeHourNextMonth=SpecificHumidity2RelativeHumidity(shumsThreeHourNextMonth,tempsThreeHourNextMonth,pressThreeHourNextMonth)

		                        timeThreeHour1=np.zeros((9))                                
		                        precsThreeHour1=np.zeros((9,nlat,nlon))
		                        tempsThreeHour1=np.zeros((9,nlat,nlon))
		                        lradsThreeHour1=np.zeros((9,nlat,nlon))
		                        sradsThreeHour1=np.zeros((9,nlat,nlon))
		                        pressThreeHour1=np.zeros((9,nlat,nlon))
		                        shumsThreeHour1=np.zeros((9,nlat,nlon))
		                        windsThreeHour1=np.zeros((9,nlat,nlon))
		                        rhsThreeHour1  =np.zeros((9,nlat,nlon))

		                        timeThreeHour1[0:8]=timeThreeHour[d*8:(d+1)*8]
		                        timeThreeHour1[8]=timeThreeHourNextMonth
		                        precsThreeHour1[0:8,:,:]=precsThreeHour[d*8:(d+1)*8,:,:]
		                        precsThreeHour1[8,:,:]=precsThreeHourNextMonth[:,:]
		                        tempsThreeHour1[0:8,:,:]=tempsThreeHour[d*8:(d+1)*8,:,:]
		                        tempsThreeHour1[8,:,:]=tempsThreeHourNextMonth[:,:]

		                        lradsThreeHour1[0:8,:,:]=lradsThreeHour[d*8:(d+1)*8,:,:]
		                        lradsThreeHour1[8,:,:]=lradsThreeHourNextMonth[:,:]

		                        sradsThreeHour1[0:8,:,:]=sradsThreeHour[d*8:(d+1)*8,:,:]
		                        sradsThreeHour1[8,:,:]=sradsThreeHourNextMonth[:,:]

		                        pressThreeHour1[0:8,:,:]=pressThreeHour[d*8:(d+1)*8,:,:]
		                        pressThreeHour1[8,:,:]=pressThreeHourNextMonth[:,:]

		                        shumsThreeHour1[0:8,:,:]=shumsThreeHour[d*8:(d+1)*8,:,:]
		                        shumsThreeHour1[8,:,:]=shumsThreeHourNextMonth[:,:]

		                        windsThreeHour1[0:8,:,:]=windsThreeHour[d*8:(d+1)*8,:,:]
		                        windsThreeHour1[8,:,:]=windsThreeHourNextMonth[:,:]

		                        rhsThreeHour1[0:8,:,:]=rhsThreeHour[d*8:(d+1)*8,:,:]
		                        rhsThreeHour1[8,:,:]=rhsThreeHourNextMonth[:,:]

				        times[:]   =GetHourlyTimes(timeThreeHour1)[:]
				        prec[:,:,:]=ConvertThreeHourlyVar2Hourly(precsThreeHour1)[:,:,:]
				        temp[:,:,:]=ConvertThreeHourlyVar2Hourly(tempsThreeHour1)[:,:,:]-273.15
				        lrad[:,:,:]=ConvertThreeHourlyVar2Hourly(lradsThreeHour1)[:,:,:]
				        srad[:,:,:]=ConvertThreeHourlyVar2Hourly(sradsThreeHour1)[:,:,:]
				        pres[:,:,:]=ConvertThreeHourlyVar2Hourly(pressThreeHour1)[:,:,:]
				        shum[:,:,:]=ConvertThreeHourlyVar2Hourly(shumsThreeHour1)[:,:,:]
				        wind[:,:,:]=ConvertThreeHourlyVar2Hourly(windsThreeHour1)[:,:,:]
				        RH[:,:,:]  =ConvertThreeHourlyVar2Hourly(rhsThreeHour1)[:,:,:]/100.0                      
                        nc4file.close()
                        jd=jd+1


def ExtractHourlyDataFromYangKunDataset2(SrcDir,ForingOutputDir,startYear,EndYear,lat,lon,interpolationType,timezone,compress):
      nlat=np.shape(lat)[0]
      nlon=np.shape(lat)[1]
      ## calculate the cellsize, column and row numbers, and the min and max in col and row numbers
      latul,lonul,colrowbound,col,row,cellsize,tempdata=GetRowColInYangkun(SrcDir,startYear,lat,lon)
      NearPointRow,NearPointCol,IDWWeights=GetNearPoints(lat,lon,latul,lonul,colrowbound,cellsize,tempdata)

      
      for y in range(startYear,EndYear):
            jd=1
            Dir=ForingOutputDir+str(y)
            if not os.path.exists(Dir):
                  os.makedirs(Dir)
            for m in range(1,13):
                  datestr=str(y*100+m)
                  lradsThreeHour,timeThreeHour=GetVariable(SrcDir,"lrad",datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
                  sradsThreeHour,timeThreeHour=GetVariable(SrcDir,"srad",datestr,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
                  ndays=np.int(np.shape(timeThreeHour)[0]/8)
                  for d in range(ndays):
                        print y,m,d
                        outputncfilename=Dir+"/"+str(y)+"-"+str(y*1000+jd)[-3:]+".nc"
                        nc4file = Dataset(outputncfilename, "w", format="NETCDF3_CLASSIC")
                        NY     = nc4file.createDimension("NY", nlat)
                        NX     = nc4file.createDimension("NX", nlon)
                        time    = nc4file.createDimension("time",24)
      
                        lats= nc4file.createVariable("lats","f4",("NY","NX",))
                        lats.units="degrees_north"
                        lats.long_name="Latitude"
                        lats[:,:]=lat[:,:]
            
                        lons= nc4file.createVariable("lons","f4",("NY","NX",))
                        lons.units="degrees_east"
                        lons.long_name="Longitude"
                        lons[:,:]=lon[:,:]

                        times= nc4file.createVariable("time","i4",("time",))
                        times.units="UTC times since 1900-01-01 00:00:0.0"
                        times.long_name="Time"



                        lrad= nc4file.createVariable("lrad","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #lrad.scale_factor=0.25
                        #lrad.add_offset=685.0
                        lrad.missing_value=-32767
                        lrad.units="W m-2"
                        lrad.long_name="Surface downward longwave radiation"
                  
                        srad= nc4file.createVariable("srad","f4",("time","NY","NX",),fill_value=-32767,zlib=compress)
                        #srad.scale_factor=0.25
                        #srad.add_offset=685.0
                        srad.missing_value=-32767
                        srad.units="W m-2"
                        srad.long_name="Surface downward shortwave radiation"

                        if(d<ndays-1):
		                times[:]   =GetHourlyTimes(timeThreeHour[d*8:(d+1)*8+1])[:]
		                lrad[:,:,:]=ConvertThreeHourlyVar2Hourly(lradsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
		                srad[:,:,:]=ConvertThreeHourlyVar2Hourly(sradsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
                        else:
                                m1=m+1
                                y1=y
                                if(m1>12):
                                   y1=y+1
                                   m1=1
                                if(y1==2016):
				        times[:]   =GetHourlyTimes(timeThreeHour[d*8:(d+1)*8+1])[:]
				        lrad[:,:,:]=ConvertThreeHourlyVar2Hourly(lradsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
				        srad[:,:,:]=ConvertThreeHourlyVar2Hourly(sradsThreeHour[d*8:(d+1)*8+1,:,:])[:,:,:]
                                else:
				        datestr1=str(y1*100+m1)
				        
				        
				        lradsThreeHourNextMonth,timeThreeHourNextMonth=GetNextMonthVariable(SrcDir,"lrad",datestr1,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)
				        sradsThreeHourNextMonth,timeThreeHourNextMonth=GetNextMonthVariable(SrcDir,"srad",datestr1,lat,lon,interpolationType,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)

		                        timeThreeHour1=np.zeros((9))                                
		                        precsThreeHour1=np.zeros((9,nlat,nlon))
		                        tempsThreeHour1=np.zeros((9,nlat,nlon))
		                        lradsThreeHour1=np.zeros((9,nlat,nlon))
		                        sradsThreeHour1=np.zeros((9,nlat,nlon))
		                        pressThreeHour1=np.zeros((9,nlat,nlon))
		                        shumsThreeHour1=np.zeros((9,nlat,nlon))
		                        windsThreeHour1=np.zeros((9,nlat,nlon))
		                        rhsThreeHour1  =np.zeros((9,nlat,nlon))

		                        timeThreeHour1[0:8]=timeThreeHour[d*8:(d+1)*8]
		                        timeThreeHour1[8]=timeThreeHourNextMonth

		                        lradsThreeHour1[0:8,:,:]=lradsThreeHour[d*8:(d+1)*8,:,:]
		                        lradsThreeHour1[8,:,:]=lradsThreeHourNextMonth[:,:]

		                        sradsThreeHour1[0:8,:,:]=sradsThreeHour[d*8:(d+1)*8,:,:]
		                        sradsThreeHour1[8,:,:]=sradsThreeHourNextMonth[:,:]
				        times[:]   =GetHourlyTimes(timeThreeHour1)[:]
				        lrad[:,:,:]=ConvertThreeHourlyVar2Hourly(lradsThreeHour1)[:,:,:]
				        srad[:,:,:]=ConvertThreeHourlyVar2Hourly(sradsThreeHour1)[:,:,:]
                   
                        nc4file.close()
                        jd=jd+1



