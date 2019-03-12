!@  Developed by Dr. Yanlin Zhang,Hunan University of Science and Technology, Xiangtan, Hunan province,411201
!@  Apr 2017

      program PreGBHM
      use gisutil
      use mod_preprocess
      !integer,parameter::NX=538,NY=698
      !real,dimension(NY,NX)::flowlen
      !integer,dimension(NY,NX)::dir
      integer::flag,i,j,openstat,subcount,iostatus
      character(200)::temp,filename,Inputfile,Outputfile,maskfile,outputfileasc
      character(200)::subbasinfile,atmp
      character(5)::postfix
      real::cellsize,nodata,start,finish
      real,allocatable::indata(:,:)
      integer,allocatable::indataint(:,:),subbasin(:,:)
      logical::paradir,file_e
      integer::NY,NX,r,c,tmpi,nsub,fileunit,d1,d2,a1,a2
      real(8)::X0,Y0  
      integer::number
      character*200::setupfile    
      NAMELIST /Dir_and_files/ model_para_dir,DemFile,Outlet,skyviewfile,midfileformat
      NAMELIST /para/ ThresHold,cpus,ModelRes,DemRes,smallestWSshed,dx_max,river_parameter_only
      
      CALL GETARG(1, setupfile)  
      print*,number,setupfile

      cellsize=1000
      nodata=-9999
      flag=1
      !OPEN CONTROL FILE FOR INPUT MODEL CONTROL PARAMETERS
      fileunit=10
      print*,setupfile
      call strlen(setupfile,a1,a2)
      open(unit=fileunit,file=setupfile(a1:a2),status='OLD',err=3,iostat=openstat)
      read(fileunit,Dir_and_files,END=102)
      !read(fileunit,para,END=101)
      !close(fileunit)

      !open(unit=12,file=setupfile,status='OLD',err=3,iostat=openstat)
      read(fileunit,para,END=102)
      close(fileunit)
3     if (openstat.NE.0) then
        print*,openstat
        print *, "Error in opening the setup.txt file"
        STOP
      end if


      print*,river_parameter_only
      call cpu_time(start)

      call strlen(model_para_dir,d1,d2) 
      hydro_para_dir=trim(model_para_dir(d1:d2))//'hydropara'
      inquire(file=hydro_para_dir,exist=paradir)
      if(.not. paradir) call system("mkdir "//hydro_para_dir)
     
      call strlen(hydro_para_dir,d1,d2) 
      inquire(file=trim(hydro_para_dir(d1:d2))//'/riverpara',exist=paradir)


      if(paradir) then
         call system("rm -r "//trim(hydro_para_dir(d1:d2))//'/riverpara/')
		 print*,"rm -r"//trim(hydro_para_dir(d1:d2))//'/riverpara/'
      end if
      call system("mkdir "//trim(hydro_para_dir(d1:d2))//"/riverpara")
	  
      if(river_parameter_only.eq.0) then
        call TauDem()
        call priver_horton()
        call deminput()
        call morph_river()
        call strlen(hydro_para_dir,d1,d2)
        subcatchment2file=trim(hydro_para_dir(d1:d2))//"/subcatchment.txt"
        inquire(file=trim(hydro_para_dir(d1:d2))//"/subcatchment.txt", exist=file_e)
        if(file_e) then
          call river_parameter()
        else
          print*,subcatchment2file//' does not exist, please provide it and run the program again!'
          stop
        end if
      else
        call strlen(hydro_para_dir,d1,d2)
        outputfile=trim(hydro_para_dir(d1:d2))//"/ws.asc"
        call read_asc_head(outputfile,NY,NX,x0,y0,cellsize,nodata)
        allocate(subbasin(NY,NX))
        call read_asc_int(outputfile,NY,NX,subbasin)
        tmpi=0
        nsub=0
        do r=1,NY
          do c = 1,NX
            if( subbasin(r,c) > tmpi) then
              nsub=subbasin(r,c)
              tmpi=subbasin(r,c)
            endif
          enddo
        enddo
        deallocate(subbasin)

        call strlen(hydro_para_dir,l1,l2)
        subbasinfile=trim(hydro_para_dir(l1:l2))//'/subbasin.dat'
        open(14,file = subbasinfile, status='old')
        subcount=0
        iostatus=1
        do while(iostatus.ge.0)
          read(14,*,iostat=iostatus) atmp
          subcount=subcount+1
        end do
        close(14)
        subcount=subcount-1

        if(nsub.ne.subcount) then
           print*, "nsub in ws.asc is not equals to that in subbasin.dat"
           print*,"The program stopped"
           stop
        end if

        NumberOfSubWatersheds=nsub
        print*,'river_parameter_only'
        call strlen(hydro_para_dir,d1,d2)
        subcatchment2file=trim(hydro_para_dir(d1:d2))//"/subcatchment.txt"
        inquire(file=trim(hydro_para_dir(d1:d2))//"/subcatchment.txt", exist=file_e)
        if(file_e) then
          call river_parameter()
        else
          print*,subcatchment2file//' does not exist, please provide it and run the program again!'
          stop
        end if
      end if
      call cpu_time(finish)
      print*,'runing time: ', finish-start
      return
 102  PRINT *,'Error reading NAMELIST file. Program stopped.'
      STOP
      
      end program PreGBHM
