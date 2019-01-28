#!/bin/bash
HM          	= .
LHM		= $(HM)
BINPATH 	= $(LHM)
NAME		= SHAWDHM
BINS		= $(BINPATH)/$(NAME)

CFLAGS	= -O3   #-mcmodel=medium

ifeq ($(FC),f77)
	FC=gfortran
endif


ifeq ($(FC),gfortran)
	LIB_NETCDF_INC =/usr/include
	LIB_NETCDF_LIB =/usr/lib
	FLAGS          =-fbounds-check -ffpe-trap=invalid,zero,overflow -fbacktrace -fopenmp
	NETCDF_FLAG    =-I$(LIB_NETCDF_INC) -L$(LIB_NETCDF_LIB) -lnetcdf -lnetcdff
endif


SCRSPATH1	= .


OBJ		= 	$(SCRSPATH1)/shr_kind_mod.mod $(SCRSPATH1)/dims_mod.mod\
		  	$(SCRSPATH1)/constvar_mod.mod $(SCRSPATH1)/controlpara_mod.mod\
		  	$(SCRSPATH1)/soilproperty_mod.mod $(SCRSPATH1)/matrix_mod.mod\
		  	$(SCRSPATH1)/rsparm_mod.mod $(SCRSPATH1)/clayrs_mod.mod\
		  	$(SCRSPATH1)/writeit_mod.mod $(SCRSPATH1)/waterbal_mod.mod\
			$(SCRSPATH1)/windv_mod.mod $(SCRSPATH1)/radcan_mod.mod\
			$(SCRSPATH1)/radres_mod.mod $(SCRSPATH1)/radout_mod.mod\
			$(SCRSPATH1)/swrcoe_mod.mod $(SCRSPATH1)/spwatr_mod.mod\
			$(SCRSPATH1)/savedata_mod.mod $(SCRSPATH1)/residu_mod.mod\
			$(SCRSPATH1)/statevar_mod.mod $(SCRSPATH1)/calsun_mod.mod\
			$(SCRSPATH1)/shaw27_mod.mod $(SCRSPATH1)/input_mod.mod\
			$(SCRSPATH1)/goshaw_mod.mod $(SCRSPATH1)/hydro_mod.mod
			
OBJ_O		= 	$(SCRSPATH1)/shr_kind_mod.o $(SCRSPATH1)/dims_mod.o\
		  	$(SCRSPATH1)/constvar_mod.o $(SCRSPATH1)/controlpara_mod.o\
		  	$(SCRSPATH1)/soilproperty_mod.o $(SCRSPATH1)/matrix_mod.o\
		  	$(SCRSPATH1)/rsparm_mod.o $(SCRSPATH1)/clayrs_mod.o\
		  	$(SCRSPATH1)/writeit_mod.o $(SCRSPATH1)/waterbal_mod.o\
			$(SCRSPATH1)/windv_mod.o $(SCRSPATH1)/radcan_mod.o\
			$(SCRSPATH1)/radres_mod.o $(SCRSPATH1)/radout_mod.o\
			$(SCRSPATH1)/swrcoe_mod.o $(SCRSPATH1)/spwatr_mod.o\
			$(SCRSPATH1)/savedata_mod.o $(SCRSPATH1)/residu_mod.o\
			$(SCRSPATH1)/statevar_mod.o $(SCRSPATH1)/calsun_mod.o\
			$(SCRSPATH1)/shaw27_mod.o $(SCRSPATH1)/input_mod.o\
			$(SCRSPATH1)/goshaw_mod.o $(SCRSPATH1)/hydro_mod.o			




DEBUG   = -g -Wall
CC      = gcc $(DEBUG)




%.mod: %.F90 
	$(FC) -c $(CFLAGS) $(FLAGS) $(NETCDF_FLAG) $<
#	$(FC) -c -fbounds-check -g -ffpe-trap=invalid,zero,overflow -fbacktrace $<
%.o: %.F90 
	$(FC) -c $(CFLAGS) $(FLAGS) $(NETCDF_FLAG) $<

	
	
all: SHAWDHM

SHAWDHM: $(OBJ) 
	$(FC) -o ./SHAWDHM.exe $(OBJ_O) ./testPrgram.F90 $(CFLAGS) $(FLAGS) $(NETCDF_FLAG)
	
clean:
	rm -rf *.o *~ $(OBJ)
