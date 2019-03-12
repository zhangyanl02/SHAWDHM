from netCDF4 import Dataset
import numpy as np
import os.path
from matplotlib import pyplot as plt
from matplotlib import dates
import datetime as dt
# Arou ROW=30  COL=32

NROW=66
NCOL=94
Julian=366
first_LAI=np.zeros((NROW,NCOL),dtype='f')
second_LAI=np.zeros((NROW,NCOL),dtype='f')

for year in range(2000,2010):
	firstday=1
	secondday=1
	LAI=np.zeros((366,NROW,NCOL),dtype='f')
	day=0
	for j in range(1,47):
		firstday=(j-1)*8+1
		secondday=secondday+8
		first_filename="./ascii/"+str(year*1000+firstday)+".asc"
		if(secondday>366):
			second_filename="./ascii/"+str((year+1)*1000+1)+".asc"
		else:
			second_filename="./ascii/"+str(year*1000+secondday)+".asc"
		first_file=open(first_filename,"r")
		for k in range(6):
			first_file.readline()
		first_LAI=np.loadtxt(first_file)[:,:]
		#print first_LAI
		second_file=open(second_filename,"r")
		for k in range(6):
			second_file.readline()
		second_LAI=np.loadtxt(second_file)[:,:]
		first_file.close()
		second_file.close()
		for i in range(8):
			print "  ",year,j,day
			if(day>365):
				break		
			LAI[day,:,:]=first_LAI[:,:]+(second_LAI[:,:]-first_LAI[:,:])/8.0*i
			print LAI[day,0,0]
			day=day+1
	outfile =Dataset(str(year)+".nc", 'w')
	dimNY="NY"
	dimNX="NX"
	dimDay="Julian"
	outfile.createDimension(dimNY, NROW) 
	outfile.createDimension(dimNX, NCOL) 
	outfile.createDimension(dimDay, Julian) 
	dimLAI = ('Julian','NY','NX', )    
	varLAI = outfile.createVariable('LAI','f', dimLAI )
	LAI[LAI>200]=0.0
	varLAI[:,:,:] = LAI/10.0*4.0/5.0
	outfile.close()
