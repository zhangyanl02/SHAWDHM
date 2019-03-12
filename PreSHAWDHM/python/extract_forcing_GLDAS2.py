from netCDF4 import Dataset
import numpy as np
from gisutil import *
import os
from DateUtil import *
from operator import itemgetter, attrgetter 

def SpecificHumidity2RelativeHumidity(SH,T,P):
      RH=0.263*P*SH/(np.exp(17.67*(T-273.15)/(T-29.65)))
      return RH

def ExtVarFromGLDAS(SrcDir,OutDir,year,day,lat,lon,InterpolationType,offset,datastr,ncformat,compress,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights):
    ListLen=len(datastr)
    nlatOut=np.shape(lat)[0]
    nlonOut=np.shape(lat)[1]
    precsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    tempsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    #tempsTmp1=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)####
    lradsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    sradsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    pressTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    shumsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    windsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    rhsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),  dtype=np.float)
    albedotmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),  dtype=np.float)
    precsTmp[:,:,:]=-9999.0
    tempsTmp[:,:,:]=-9999.0
    lradsTmp[:,:,:]=-9999.0
    sradsTmp[:,:,:]=-9999.0
    pressTmp[:,:,:]=-9999.0
    shumsTmp[:,:,:]=-9999.0
    windsTmp[:,:,:]=-9999.0
    rhsTmp[:,:,:]  =-9999.0

    for itime in range(len(datastr)):
        year1=int(datastr[itime][0:4])
        if(year1>=2010):
            ncfilename=SrcDir+str(year1)+"/GLDAS_NOAH025_3H.A"+datastr[itime]+".021.nc4"
        else:
            ncfilename=SrcDir+str(year1)+"/GLDAS_NOAH025_3H.A"+datastr[itime]+".020.nc4"
        ncfile = Dataset(ncfilename.strip(), "r")
        #minrow=  colrowbound[0]
        #maxrow=  colrowbound[1]
        #mincol=  colrowbound[2]
        #maxcol=  colrowbound[3]                                                   
        Albedo_inst=ncfile.variables["Albedo_inst"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        Wind_f_inst=ncfile.variables["Wind_f_inst"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        Tair_f_inst=ncfile.variables["Tair_f_inst"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        Qair_f_inst=ncfile.variables["Qair_f_inst"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        Psurf_f_inst=ncfile.variables["Psurf_f_inst"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        Rainf_f_tavg=ncfile.variables["Rainf_f_tavg"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        SWdown_f_tavg=ncfile.variables["SWdown_f_tavg"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        LWdown_f_tavg=ncfile.variables["LWdown_f_tavg"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        ncfile.close()

        #tempsTmp1[3*itime,:,:] =Tair_f_inst[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
        #np.savetxt("temp_NR.txt",tempsTmp1[3*itime,:,:])
      # No interpolation
        if(InterpolationType==0):
              albedotmp[3*itime,:,:]=Albedo_inst[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
              windsTmp[3*itime,:,:] =Wind_f_inst[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
              tempsTmp[3*itime,:,:] =Tair_f_inst[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
              shumsTmp[3*itime,:,:] =Qair_f_inst[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
              pressTmp[3*itime,:,:] =Psurf_f_inst[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
              precsTmp[3*itime,:,:] =Rainf_f_tavg[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
              sradsTmp[3*itime,:,:] =SWdown_f_tavg[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
              lradsTmp[3*itime,:,:] =LWdown_f_tavg[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
              rhsTmp[3*itime,:,:]   =0.263*pressTmp[3*itime,:,:]*shumsTmp[3*itime,:,:]/(np.exp(17.67*(tempsTmp[3*itime,:,:]-273.15)/( tempsTmp[3*itime,:,:]-29.65)))
        else:
        # Inverse distance interpolation       
              albedotmp[3*itime,:,:]=IDWWeights[0,:,:]*Albedo_inst[0,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*Albedo_inst[0,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*Albedo_inst[0,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*Albedo_inst[0,NearPointRow[3,:,:],NearPointCol[3,:,:]]
              windsTmp[3*itime,:,:]=IDWWeights[0,:,:]*Wind_f_inst[0,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*Wind_f_inst[0,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*Wind_f_inst[0,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*Wind_f_inst[0,NearPointRow[3,:,:],NearPointCol[3,:,:]]
              tempsTmp[3*itime,:,:]=IDWWeights[0,:,:]*Tair_f_inst[0,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*Tair_f_inst[0,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*Tair_f_inst[0,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*Tair_f_inst[0,NearPointRow[3,:,:],NearPointCol[3,:,:]]
              shumsTmp[3*itime,:,:]=IDWWeights[0,:,:]*Qair_f_inst[0,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*Qair_f_inst[0,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*Qair_f_inst[0,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*Qair_f_inst[0,NearPointRow[3,:,:],NearPointCol[3,:,:]]
              pressTmp[3*itime,:,:]=IDWWeights[0,:,:]*Psurf_f_inst[0,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*Psurf_f_inst[0,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*Psurf_f_inst[0,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*Psurf_f_inst[0,NearPointRow[3,:,:],NearPointCol[3,:,:]]
              precsTmp[3*itime,:,:]=IDWWeights[0,:,:]*Rainf_f_tavg[0,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*Rainf_f_tavg[0,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*Rainf_f_tavg[0,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*Rainf_f_tavg[0,NearPointRow[3,:,:],NearPointCol[3,:,:]]
              sradsTmp[3*itime,:,:]=IDWWeights[0,:,:]*SWdown_f_tavg[0,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*SWdown_f_tavg[0,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*SWdown_f_tavg[0,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*SWdown_f_tavg[0,NearPointRow[3,:,:],NearPointCol[3,:,:]]
              lradsTmp[3*itime,:,:]=IDWWeights[0,:,:]*LWdown_f_tavg[0,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*LWdown_f_tavg[0,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*LWdown_f_tavg[0,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*LWdown_f_tavg[0,NearPointRow[3,:,:],NearPointCol[3,:,:]]
              rhsTmp[3*itime,:,:]   =0.263*pressTmp[3*itime,:,:]*shumsTmp[3*itime,:,:]/(np.exp(17.67*(tempsTmp[3*itime,:,:]-273.15)/( tempsTmp[3*itime,:,:]-29.65)))
              #np.savetxt("temp_IWD.txt",tempsTmp[3*itime,:,:])
              #stop
    for t in range(len(datastr)-1):
        albedotmp[3*t+1,:,:]=albedotmp[3*t,:,:]*2.0/3.0+albedotmp[3*(t+1),:,:]*1.0/3.0
        albedotmp[3*t+2,:,:]=albedotmp[3*t,:,:]*1.0/3.0+albedotmp[3*(t+1),:,:]*2.0/3.0

        windsTmp[3*t+1,:,:]=windsTmp[3*t,:,:]*2.0/3.0+windsTmp[3*(t+1),:,:]*1.0/3.0
        windsTmp[3*t+2,:,:]=windsTmp[3*t,:,:]*1.0/3.0+windsTmp[3*(t+1),:,:]*2.0/3.0
        
        tempsTmp[3*t+1,:,:]=tempsTmp[3*t,:,:]*2.0/3.0+tempsTmp[3*(t+1),:,:]*1.0/3.0
        tempsTmp[3*t+2,:,:]=tempsTmp[3*t,:,:]*1.0/3.0+tempsTmp[3*(t+1),:,:]*2.0/3.0
        
        shumsTmp[3*t+1,:,:]=shumsTmp[3*t,:,:]*2.0/3.0+shumsTmp[3*(t+1),:,:]*1.0/3.0
        shumsTmp[3*t+2,:,:]=shumsTmp[3*t,:,:]*1.0/3.0+shumsTmp[3*(t+1),:,:]*2.0/3.0
        
        pressTmp[3*t+1,:,:]=pressTmp[3*t,:,:]*2.0/3.0+pressTmp[3*(t+1),:,:]*1.0/3.0
        pressTmp[3*t+2,:,:]=pressTmp[3*t,:,:]*1.0/3.0+pressTmp[3*(t+1),:,:]*2.0/3.0
        
        precsTmp[3*t+1,:,:]=precsTmp[3*t,:,:]*2.0/3.0+precsTmp[3*(t+1),:,:]*1.0/3.0
        precsTmp[3*t+2,:,:]=precsTmp[3*t,:,:]*1.0/3.0+precsTmp[3*(t+1),:,:]*2.0/3.0
        
        sradsTmp[3*t+1,:,:]=sradsTmp[3*t,:,:]*2.0/3.0+sradsTmp[3*(t+1),:,:]*1.0/3.0
        sradsTmp[3*t+2,:,:]=sradsTmp[3*t,:,:]*1.0/3.0+sradsTmp[3*(t+1),:,:]*2.0/3.0
        
        lradsTmp[3*t+1,:,:]=lradsTmp[3*t,:,:]*2.0/3.0+lradsTmp[3*(t+1),:,:]*1.0/3.0
        lradsTmp[3*t+2,:,:]=lradsTmp[3*t,:,:]*1.0/3.0+lradsTmp[3*(t+1),:,:]*2.0/3.0

        rhsTmp[3*t+1,:,:]=rhsTmp[3*t,:,:]*2.0/3.0+rhsTmp[3*(t+1),:,:]*1.0/3.0
        rhsTmp[3*t+2,:,:]=rhsTmp[3*t,:,:]*1.0/3.0+rhsTmp[3*(t+1),:,:]*2.0/3.0

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

    times= nc4file.createVariable("time","i4",("time",))
    times.units="hours"
    times.long_name="Time"

    prec= nc4file.createVariable("prec","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #prec.scale_factor=0.0025
    #prec.add_offset=50.0
    prec.missing_value=-9999
    prec.units="mm hr-1"  #kg m-2 s-1 in GLDAS
    prec.long_name="Total precipitation rate"

    temp= nc4file.createVariable("temp","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #temp.scale_factor=0.01
    #temp.add_offset=0.0
    temp.missing_value=-9999
    temp.units="K"
    temp.long_name="air_temperature"

    lrad= nc4file.createVariable("lrad","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #lrad.scale_factor=0.25
    #lrad.add_offset=685.0
    lrad.missing_value=-9999
    lrad.units="W m-2"
    lrad.long_name="surface_downwelling_longwave_flux_in_air"
    
    srad= nc4file.createVariable("srad","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #srad.scale_factor=0.25
    #srad.add_offset=685.0
    srad.missing_value=-9999
    srad.units="W m-2"
    srad.long_name="surface_downwelling_shortwave_flux_in_air"

    pres= nc4file.createVariable("pres","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #pres.scale_factor=2.0
    #pres.add_offset=63500.
    pres.missing_value=-9999
    pres.units="Pa"
    pres.long_name="surface_air_pressure"                  

    shum= nc4file.createVariable("shum","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #shum.scale_factor=0.000001
    #shum.add_offset=0.025
    shum.missing_value=-9999
    shum.units="kg kg-1"
    shum.long_name="specific_humidity"
    
    wind= nc4file.createVariable("wind","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #wind.scale_factor=0.002
    #wind.add_offset=60.0
    wind.missing_value=-9999
    wind.units="m s-1"
    wind.long_name="Near surface wind speed"

    RH= nc4file.createVariable("rhum","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #RH.scale_factor=0.01            
    RH.missing_value=-9999
    RH.units="0~1"
    RH.long_name="Near surface Relative Humidity"

    albedo= nc4file.createVariable("albedo","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #RH.scale_factor=0.01            
    albedo.missing_value=-9999
    albedo.units="%"
    albedo.long_name="surface albedo"
    
    prec[:,:,:]=precsTmp[offset:offset+24,:,:]*3600.0
    temp[:,:,:]=tempsTmp[offset:offset+24,:,:]-273.15
    lrad[:,:,:]=lradsTmp[offset:offset+24,:,:]
    srad[:,:,:]=sradsTmp[offset:offset+24,:,:]
    pres[:,:,:]=pressTmp[offset:offset+24,:,:]
    shum[:,:,:]=shumsTmp[offset:offset+24,:,:]
    wind[:,:,:]=windsTmp[offset:offset+24,:,:]
    RH[:,:,:]=rhsTmp[offset:offset+24,:,:]/100.0
    albedo[:,:,:]=albedotmp[offset:offset+24,:,:]
    nc4file.close()


def GetThreeHourlyFileNames(year,day,timezone):
    if(timezone!=0):
	    FirstUCTHour=0-timezone
	    FirstDay=day
	    FirstYear=year
	    if(FirstUCTHour<0):
		FirstDay=day-1
		FirstUCTHour=FirstUCTHour+24
	    if(FirstDay<=0):
		FirstYear=year-1
		if(Leap(FirstYear)):
		    FirstDay=366
		else:
		    FirstDay=365
	    for i in range(8):
		if(FirstUCTHour>=3*i and FirstUCTHour<3*(i+1)):
		    firstHour=3*i
	    offset=FirstUCTHour-firstHour
	    datestrList=[]
	    for i in range(10):
		tmpday=FirstDay
		tmpyear=FirstYear
		hour=firstHour+i*3
		if(hour>=24):
		    hour=hour-24
		    tmpday=FirstDay+1
		    if(Leap(FirstYear)):
		        ndays=366
		    else:
		        ndays=365
		    if(tmpday>ndays):
		        tmpyear=FirstYear+1
		        tmpday=1
		date1=GetDateFromJulianDay(tmpyear,tmpday)
		datestr=str(tmpyear)+str(100+date1.month)[1:]+str(100+date1.day)[1:]+"."+str(10000+hour*100)[1:]
		datestrList.append(datestr)
    else:
	    datestrList=[]
	    offset=0
	    tmpday=day
	    tmpyear=year
	    for i in range(9):
		hour=i*3
		if(hour>=24):
		    hour=hour-24
		    tmpday=day+1
		    if(Leap(year)):
		        ndays=366
		    else:
		        ndays=365
		    if(tmpday>ndays):
		        tmpyear=year+1
		        tmpday=1
		date1=GetDateFromJulianDay(tmpyear,tmpday)
		datestr=str(tmpyear)+str(100+date1.month)[1:]+str(100+date1.day)[1:]+"."+str(10000+hour*100)[1:]
		datestrList.append(datestr)
    return offset,datestrList


def GetRowColInGLDAS2(SrcDir,startYear,lat,lon):
      nlat=np.shape(lat)[0]
      nlon=np.shape(lat)[1]
      if(startYear>=2011):
            ncfilename=SrcDir+str(startYear)+"/GLDAS_NOAH025_3H.A"+str(startYear)+"0101.0000"+".021.nc4"
      else:
            ncfilename=SrcDir+str(startYear)+"/GLDAS_NOAH025_3H.A"+str(startYear)+"0101.0000"+".020.nc4"
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
      tempdata=ncfile.variables["Tair_f_inst"][0,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
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


def ExtractHourlyDataFromGLDAS(SrcDir,ForingOutputDir,startYear,EndYear,lat,lon,interpolationType,timezone,compress,ncformat):
      nlat=np.shape(lat)[0]
      nlon=np.shape(lat)[1]
      latul,lonul,colrowbound,col,row,cellsize,tempdata=GetRowColInGLDAS2(SrcDir,startYear,lat,lon)
      NearPointRow,NearPointCol,IDWWeights=GetNearPoints(lat,lon,latul,lonul,colrowbound,cellsize,tempdata)

      for y in range(startYear,EndYear):
            Dir=ForingOutputDir+str(y)
            if not os.path.exists(Dir):
                  os.makedirs(Dir)
            if(Leap(y)):
                ndays=366
            else:
                ndays=365
            for d in range(1,ndays+1):
		print y,d
                offset,datestrList=GetThreeHourlyFileNames(y,d,timezone)
                ExtVarFromGLDAS(SrcDir,Dir,y,d,lat,lon,interpolationType,offset,datestrList,ncformat,compress,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights)




def ExtVarFromGLDAS2(SrcDir,OutDir,year,day,lat,lon,InterpolationType,offset,datastr,ncformat,compress,row,col,cellsize,colrowbound,NearPointRow,NearPointCol,IDWWeights):
    ListLen=len(datastr)
    nlatOut=np.shape(lat)[0]
    nlonOut=np.shape(lat)[1]
    precsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    tempsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    #tempsTmp1=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)####
    lradsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    sradsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    pressTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    shumsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    windsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),dtype=np.float)
    rhsTmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),  dtype=np.float)
    albedotmp=np.zeros(((ListLen-1)*3+1,nlatOut,nlonOut),  dtype=np.float)
    precsTmp[:,:,:]=-9999.0
    tempsTmp[:,:,:]=-9999.0
    lradsTmp[:,:,:]=-9999.0
    sradsTmp[:,:,:]=-9999.0
    pressTmp[:,:,:]=-9999.0
    shumsTmp[:,:,:]=-9999.0
    windsTmp[:,:,:]=-9999.0
    rhsTmp[:,:,:]  =-9999.0

    for itime in range(len(datastr)):
        year1=int(datastr[itime][0:4])
        if(year1>=2010):
            ncfilename=SrcDir+str(year1)+"/GLDAS_NOAH025_3H.A"+datastr[itime]+".021.nc4"
        else:
            ncfilename=SrcDir+str(year1)+"/GLDAS_NOAH025_3H.A"+datastr[itime]+".020.nc4"
        ncfile = Dataset(ncfilename.strip(), "r")
        #minrow=  colrowbound[0]
        #maxrow=  colrowbound[1]
        #mincol=  colrowbound[2]
        #maxcol=  colrowbound[3]
        SWdown_f_tavg=ncfile.variables["SWdown_f_tavg"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        LWdown_f_tavg=ncfile.variables["LWdown_f_tavg"][:,colrowbound[0]:colrowbound[1]+1,colrowbound[2]:colrowbound[3]+1]
        ncfile.close()

        #tempsTmp1[3*itime,:,:] =Tair_f_inst[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
        #np.savetxt("temp_NR.txt",tempsTmp1[3*itime,:,:])
      # No interpolation
        if(InterpolationType==0):
              sradsTmp[3*itime,:,:] =SWdown_f_tavg[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
              lradsTmp[3*itime,:,:] =LWdown_f_tavg[0,row[:,:]-colrowbound[0],col[:,:]-colrowbound[2]]
        else:
        # Inverse distance interpolation       
              sradsTmp[3*itime,:,:]=IDWWeights[0,:,:]*SWdown_f_tavg[0,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*SWdown_f_tavg[0,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*SWdown_f_tavg[0,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*SWdown_f_tavg[0,NearPointRow[3,:,:],NearPointCol[3,:,:]]
              lradsTmp[3*itime,:,:]=IDWWeights[0,:,:]*LWdown_f_tavg[0,NearPointRow[0,:,:],NearPointCol[0,:,:]]+IDWWeights[1,:,:]*LWdown_f_tavg[0,NearPointRow[1,:,:],NearPointCol[1,:,:]]+IDWWeights[2,:,:]*LWdown_f_tavg[0,NearPointRow[2,:,:],NearPointCol[2,:,:]]+IDWWeights[3,:,:]*LWdown_f_tavg[0,NearPointRow[3,:,:],NearPointCol[3,:,:]]
              #np.savetxt("temp_IWD.txt",tempsTmp[3*itime,:,:])
              #stop
    for t in range(len(datastr)-1):      
        sradsTmp[3*t+1,:,:]=sradsTmp[3*t,:,:]*2.0/3.0+sradsTmp[3*(t+1),:,:]*1.0/3.0
        sradsTmp[3*t+2,:,:]=sradsTmp[3*t,:,:]*1.0/3.0+sradsTmp[3*(t+1),:,:]*2.0/3.0        
        lradsTmp[3*t+1,:,:]=lradsTmp[3*t,:,:]*2.0/3.0+lradsTmp[3*(t+1),:,:]*1.0/3.0
        lradsTmp[3*t+2,:,:]=lradsTmp[3*t,:,:]*1.0/3.0+lradsTmp[3*(t+1),:,:]*2.0/3.0


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

    times= nc4file.createVariable("time","i4",("time",))
    times.units="hours"
    times.long_name="Time"


    lrad= nc4file.createVariable("lrad","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #lrad.scale_factor=0.25
    #lrad.add_offset=685.0
    lrad.missing_value=-9999
    lrad.units="W m-2"
    lrad.long_name="surface_downwelling_longwave_flux_in_air"
    
    srad= nc4file.createVariable("srad","f4",("time","NY","NX",),fill_value=-9999,zlib=compress)
    #srad.scale_factor=0.25
    #srad.add_offset=685.0
    srad.missing_value=-9999
    srad.units="W m-2"
    srad.long_name="surface_downwelling_shortwave_flux_in_air"    


    lrad[:,:,:]=lradsTmp[offset:offset+24,:,:]
    srad[:,:,:]=sradsTmp[offset:offset+24,:,:]

    nc4file.close()


