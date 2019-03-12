from extract_soil_mod import *
from gisutil import *
from extract_lai import *
from extract_landuse import *
from extract_forcing_Yangkun_UTC import *
import extract_Radiation
from tsavg import *
import matplotlib.pyplot as plt
import extract_forcing_GLDAS2
import time

landuse_para=[]
landuse_para.append("iland	Kcanopy(i)	LAImax(i)	kcrop(i)	root(i)	anik(i)	Sstmax(i)	surfn(i)	landusename")
landuse_para.append("1	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("2	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("3	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("4	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("5	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("6	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("7	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("8	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("9	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("10	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("11	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("12	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("13	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("14	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("15	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")
landuse_para.append("16	0.5	4.0	2.0	2.0	1.0	1.0	1.0	grassland")



#  Get the latitude and longitude
#model_para_dir="/disk2/SHAWDHM_model/dongbei/hulunbe/input/"
#InputDemFileName="/disk2/SHAWDHM_model/dongbei/hulunbe/input/hlbe_th_dem_5km.tif"
#wsfile="/disk2/SHAWDHM_model/dongbei/hulunbe/input/mask.asc"


model_para_dir="/disk3/Solar_Radiation/heihe_YK/"
InputDemFileName="/disk3/Solar_Radiation/heihe.tif"
wsfile="/disk3/Solar_Radiation/heihe.asc"
ws=readasc(wsfile)
proj4string=""
latfilename=model_para_dir+"lat.asc"
lonfilename=model_para_dir+"lon.asc"
nc,nr,xllcorner,yllcorner,cellsize,nodata=readaschead(wsfile)
lat,lon=GetGridLatLon2(InputDemFileName,model_para_dir,wsfile,proj4string,latfilename,lonfilename)




soil=False
if(soil==True):
      #  Extract Soil Parameters
      sourceSoilDir="/disk2/data/soil/dai/nc4/"
      OutputSoilFile="/disk2/SHAWDHM_model/dongbei/hulunbe/input/soilpara.nc"
      #SoilNodeDepths=[0.025,0.1,0.2,0.5,0.8,1.4,2.0,3.0,5.0,8.0,12.0,15.0,17.0,18.0]#0, 0.04, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.4, 4, 6, 9, 12, 14, 15, 15.5
      SoilNodeDepths=[0.1,0.2,0.4,0.8,1.2,2.0,3.2,4.8,8.0,12.0,16.0,20.0,24.0,25.0]#0, 0.04, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.4, 4, 6, 9, 12, 14, 15, 15.5
      baseDepth=4.0
      print "extracting soil parameter"
      ExtractSoilDataFromDaiYongJiu(sourceSoilDir,OutputSoilFile,SoilNodeDepths,lat,lon,baseDepth)

ncformat="NETCDF3_CLASSIC"
landuse=False
if (landuse==True):
      #  Extract land use cover
      print "extracting landuse cover"
      SourceLanduse="/disk2/data/landuse/landuse.nc"
      OutputLanduse="/disk2/SHAWDHM_model/dongbei/hulunbe/input/landuse.nc"
      OutputVegPara="/disk2/SHAWDHM_model/dongbei/hulunbe/input/veg_para.dat"
      landuse=Extract_Landuse(SourceLanduse,lat,lon)
      vegParaFile=open(OutputVegPara,'w')
      vegParaFile.write("vege para\n")
      vegParaFile.write("iland, Kcanopy(i), LAImax(i), kcrop(i),root(i),anik(i), Sstmax(i), surfn(i), landusename\n")
      ii=1
      for i in range(1,17):
        if(i in landuse):
            vegParaFile.write(str(ii)+" "+str(landuse_para[i][:])+"\n")
            ii=ii+1
      vegParaFile.close()
      
      writeascfile("/disk2/SHAWDHM_model/dongbei/hulunbe/input/landuse.asc",nr,nc,xllcorner,yllcorner,cellsize,nodata,landuse,"int")
      dims=np.shape(landuse)
      landusefile = Dataset(OutputLanduse,"w",format="NETCDF3_CLASSIC")
      nlat   = landusefile.createDimension("NY", dims[0])
      nlon   = landusefile.createDimension("NX", dims[1])
      nland  = landusefile.createDimension("NV", ii-1)
      print "nland",nland
      landuses= landusefile.createVariable("landuse","i1",("NY","NX",))
      landuse_ratio= landusefile.createVariable("landratio","f4",("NY","NX","NV",))
      landuses.missing_value="-9999"
      landuses.units=""
      landuses.long_name="USGS Land cover"
      landuses[:,:]=landuse[:,:]
      landuse_ratio[:,:,:]=1/ii
      landusefile.close()

StartYear=2003
EndYear=2016
lai=False
ncformat="NETCDF3_CLASSIC"
#ncformat="NETCDF4"
if (lai==True):
      #  Extract lai
      print "extracting lai"
      for y in range(StartYear,EndYear): 
            print y
            SourceLAIname="/disk2/data/LAI/nc4/global_30s_"+str(y)+".nc"
            LAI=Extract_LAI(SourceLAIname,lat,lon)
            plt.imshow(LAI[0,:,:],interpolation='nearest')
            plt.show()
            dims=np.shape(LAI)
            laifile = Dataset("/disk2/SHAWDHM_model/dongbei/hulunbe/input/lai/"+str(y)+".nc","w",format=ncformat)
            nlat   = laifile.createDimension("lat", dims[1])
            nlon   = laifile.createDimension("lon", dims[2])
            ntime  = laifile.createDimension("time", dims[0])
            Lais= laifile.createVariable("LAI","f4",("time","lat","lon",))
            Lais.missing_value="-9999"
            Lais.units="-"
            Lais.long_name="Leaf area index"
            Lais[:,:,:]=LAI[:,:,:]
            laifile.close()

StartYear=2007
EndYear=2015
ncformat="NETCDF4"
forcing=True
if(forcing==True):
      ForingOutputDir="/disk3/Solar_Radiation/heihe_YK/"
      SrcDir="/disk3/data/MeteoForcing/Yangkun/Data_forcing_03hr_010deg_nc4/"
      interpolationType=0
      timezone=8
      compress=True
      time1=time.time()
      ExtractHourlyDataFromYangKunDataset2(SrcDir,ForingOutputDir,StartYear,EndYear,lat,lon,interpolationType,timezone,compress)
      print time.time()-time1


StartYear=2009
EndYear=2015
ncformat="NETCDF4"
forcing=False
if(forcing==True):
      ForingOutputDir="/disk3/Solar_Radiation/QTP_GLDAS/"
      SrcDir="/media/zhangyanlin-t7610/disk41/GLDAS2.0/"
      interpolationType=1
      timezone=0
      compress=True
      if (ncformat=="NETCDF3_CLASSIC"):
          compress=False
      extract_forcing_GLDAS2.ExtractHourlyDataFromGLDAS2(SrcDir,ForingOutputDir,StartYear,EndYear,lat,lon,interpolationType,timezone,compress,ncformat)


StartYear=2007
EndYear=2015
ncformat="NETCDF4"
solar=False
DT="daily" # "hourly"
if(solar==True):
      ForingOutputDir="/disk3/Solar_Radiation/bbh/"
      SrcDir="/disk3/Solar_Radiation/Tang/nc4/"
      interpolationType=0
      timezone=8
      compress=True
      if (ncformat=="NETCDF3_CLASSIC"):
          compress=False
      extract_Radiation.ExtractSolarFromTang(SrcDir,ForingOutputDir,StartYear,EndYear,lat,lon,interpolationType,timezone,compress,ncformat,"daily")



SourceDir="/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/forcing/Yangkun/"
SourceDir2="/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/output/"

tsavg=False
StartYear=2010
EndYear=2015
if(tsavg==True):
      #tsavg=ExtTsavg(SourceDir,StartYear,EndYear,nr,nc)
      tsavg=ExtTsavgFromCalculatedSoil(SourceDir2,StartYear,EndYear,nr,nc)
      writeascfile("/home/zhangyanlin-t7610/model/Horton_SHAWDHM/BBH/input/tsavg5.asc",nr,nc,xllcorner,yllcorner,cellsize,nodata,tsavg,"real")

