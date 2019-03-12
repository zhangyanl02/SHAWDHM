import pygrib
import glob
import os
SourceDir="/disk2/data/Forcing/gldas/"
StartYear=2010
EndYear=2011
for Y in range(StartYear,EndYear):
      for jd in range(200,201):
            for H in range(1):
                  FileName=SourceDir+str(Y)+"/"+"GLDAS_NOAH025SUBP_3H.A"+str(Y*1000+jd)+"."+str(10000+3*H*100)[1:]+".001.grb"
                  grbs = pygrib.open(FileName)
                  for grb in grbs:
                        print grb

