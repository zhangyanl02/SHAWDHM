import numpy as np
import pyproj
import subprocess
###import rpy2.robjects.packages as rpackages    linux

def writeascfile(filename,nr,nc,xllcorner,yllcorner,cellsize,nodata,var,int_real):
      ascfile=open(filename,'w')
      ascfile.write("ncols        "+str(nc)+"  \n")
      ascfile.write("nrows        "+str(nr)+"  \n")
      ascfile.write("xllcorner    "+str(xllcorner)+"  \n")
      ascfile.write("yllcorner    "+str(yllcorner)+"  \n")
      ascfile.write("cellsize     "+str(cellsize)+"  \n")
      ascfile.write("nodata       "+str(nodata)+"  \n")
      if(int_real=="int"):
            np.savetxt(ascfile,var,"%5i")
      elif(int_real=="real"):
            np.savetxt(ascfile,var,"%10.3f")
      ascfile.close()

def readaschead(filename):
      assascfile=open(filename,'r')
      temp=assascfile.readline()
      NX=int(temp.split()[1])
      
      temp=assascfile.readline()
      NY=int(temp.split()[1])
      
      temp=assascfile.readline()
      xllcorner=float(temp.split()[1])
      
      temp=assascfile.readline()
      yllcorner=float(temp.split()[1])
      
      temp=assascfile.readline()
      cellsize=float(temp.split()[1])
      
      temp=assascfile.readline()
      nodata=float(temp.split()[1])
      assascfile.close()
      return NX,NY,xllcorner,yllcorner,cellsize,nodata


def readasc(filename):
      assascfile=open(filename,'r')
      for i in range(6):
            temp=assascfile.readline()
      data=np.loadtxt(assascfile)
      assascfile.close()
      return data

def SplitByChar(inputstring,splitchar):
      list1=inputstring.replace('\n','').split(splitchar)
      listlen=len(list1)
      list2=[]
      for i in range(listlen):
            if(len(list1[i])!=0):
                  list2.append(list1[i])            
      return list2      

def GetTiffProj4String_win(filename):
      pipe = subprocess.Popen("gdalinfo -proj4 "+filename+" | grep proj", shell=True,stdout=subprocess.PIPE)
      result = pipe.communicate()
      list1=SplitByChar(result[0],'\'')
      if(len(list1)<=0):
            print "there is no proj4 string, please give the string mannually"
            proj4string=input("Input: ")
      else:
            for i in range(len(list1)):
                  if(len(list1[i])>5):
                        proj4string=list1[i]
      print "projr string : "  ,proj4string
      return proj4string
      
#def GetTiffProj4String2_linux(filename):
#      rgdal = rpackages.importr('rgdal')
#      sp = rpackages.importr('sp')
#      g1=rgdal.readGDAL(filename)
#      prjstring=sp.proj4string(g1)
#      proj4string=prjstring[0]    
#      return proj4string

def GetAscHead(asciifilename):
      asciifile=open(asciifilename,'r')
      temp=asciifile.readline()
      rlist=SplitByChar(temp,' ')
      strncols=rlist[0]
      ncols=np.int(rlist[1])

      temp=asciifile.readline()
      rlist=SplitByChar(temp,' ')
      strnrows=rlist[0]
      nrows=np.int(rlist[1])

      temp=asciifile.readline()
      rlist=SplitByChar(temp,' ')
      strxllcorner=rlist[0]
      xllcorner=np.float(rlist[1])

      temp=asciifile.readline()
      rlist=SplitByChar(temp,' ')
      stryllcorner=rlist[0]
      yllcorner=np.float(rlist[1])

      temp=asciifile.readline()
      rlist=SplitByChar(temp,' ')
      strcellsize=rlist[0]
      cellsize=np.float(rlist[1])

      temp=asciifile.readline()
      rlist=SplitByChar(temp,' ')
      strnodata=rlist[0]
      nodata=np.float(rlist[1])
      asciifile.close()      
      return ncols,nrows,xllcorner,yllcorner,cellsize,nodata

def GetGridLatLon(ncols,nrows,xllcorner,yllcorner,cellsize,proj4string):
      lat=np.zeros((nrows,ncols),dtype=np.float)
      lon=np.zeros((nrows,ncols),dtype=np.float)
      project2latlon = pyproj.Proj(proj4string)
      for r in range(nrows):
            for c in range(ncols):
                  X=xllcorner+cellsize/2.0+c*cellsize
                  Y=yllcorner+cellsize/2.0+(nrows-1-r)*cellsize
                  lon[r,c],lat[r,c]=project2latlon(X, Y, inverse=True) 
      return lat,lon
      

def GetGridLatLon2(InputDemFileName,model_para_dir,wsfile,projstring,latfilename,lonfilename):
      asciifilename=wsfile
      if (len(projstring.strip())==0):      
            proj4string=GetTiffProj4String_win(InputDemFileName)
            #proj4string2=GetTiffProj4String2(InputDemFileName)
      else:
            proj4string=projstring.strip()
      ncols,nrows,xllcorner,yllcorner,cellsize,nodata= GetAscHead(asciifilename)
      lat=np.zeros((nrows,ncols),dtype=np.float)
      lon=np.zeros((nrows,ncols),dtype=np.float)
      project2latlon = pyproj.Proj(proj4string)
      for r in range(nrows):
            for c in range(ncols):
                  X=xllcorner+cellsize/2.0+c*cellsize
                  Y=yllcorner+cellsize/2.0+(nrows-1-r)*cellsize
                  lon[r,c],lat[r,c]=project2latlon(X, Y, inverse=True)
      latfile=open(latfilename,'w')
      latfile.write("ncols        "+str(ncols)+"\n")
      latfile.write("nrows        "+str(nrows)+"\n")
      latfile.write("xllcorner    "+str(xllcorner)+"\n")
      latfile.write("yllcorner    "+str(yllcorner)+"\n")
      latfile.write("cellsize     "+str(cellsize)+"\n")
      latfile.write("nodata       "+str(nodata)+"\n")
      np.savetxt(latfile,lat,fmt="%12.8f")
      latfile.close()
      lonfile=open(lonfilename,'w')
      lonfile.write("ncols        "+str(ncols)+"\n")
      lonfile.write("nrows        "+str(nrows)+"\n")
      lonfile.write("xllcorner    "+str(xllcorner)+"\n")
      lonfile.write("yllcorner    "+str(yllcorner)+"\n")
      lonfile.write("cellsize     "+str(cellsize)+"\n")
      lonfile.write("nodata       "+str(nodata)+"\n")
      np.savetxt(lonfile,lon,fmt="%12.8f")
      lonfile.close()      
      return lat,lon
      
def MakeSrc():
      pipe = subprocess.Popen("make", shell=True,stdout=subprocess.PIPE)
      pipe.wait()
      print pipe.communicate()[0]
