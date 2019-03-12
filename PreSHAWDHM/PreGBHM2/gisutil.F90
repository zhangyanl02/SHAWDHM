!@  Developed by Dr. Yanlin Zhang,Hunan University of Science and Technology, Xiangtan, Hunan province,411201
!@  Apr 2017

   module gisutil
      public
      integer,allocatable,save::dirtemp(:,:)
      real,allocatable,save::flowlentemp(:,:)
        
      contains
        
      subroutine strlen(str,l1,l2)
        implicit none
        character(200)::str
        integer::i,l1,l2,k
        k=0
        do i = 1, 200
          if(k.eq.0 .and. str(i:i).NE.' ') then
            l1=i
            k=1
          elseif(k.eq.1 .and. str(i:i).EQ.' ') then
            l2 = i-1
            exit
          endif
        end do
      end subroutine

      subroutine filepostfix(str,l1,l2)
        implicit none
        character(*)::str
        integer::i,l1,l2,k
        call strlen(str,l1,l2)
        k=0
        do i= l2,1,-1
          if(k.eq.0 .and. str(i:i).NE.'.') then
            l2=i
            k=1
          elseif(k.eq.1 .and. str(i:i).EQ.'.') then
            l1 = i+1
            exit
          endif
        end do
      end subroutine
        
! most time is cost on writing the asc files.
!      subroutine check(status)
!        use netcdf
!        implicit none
!        integer, intent ( in) :: status
!        if(status /= nf90_noerr) then
!          print *, trim(nf90_strerror(status))
!          stop "Stopped"
!        end if
!      end subroutine check
!      subroutine read_nc_head(inputfile,NY,NX,X0,Y0,cellsize,nodata)
!        use netcdf
!        implicit none
!      ! This is the name of the data file we will create.
!        character(*)::inputfile
!        integer::NY,NX
!        real(8)::X0,Y0
!        real::cellsize,nodata
!        integer::ncid,data_id
!        call check(nf90_open(trim(inputfile),nf90_nowrite,ncid))
!        call check(nf90_inq_varid(ncid,"data",data_id))
!        call check( nf90_get_att(ncid, data_id, "ncols", NX))
!        call check( nf90_get_att(ncid, data_id, "nrows", NY))
!        call check( nf90_get_att(ncid, data_id, "xllcorner", x0))
!        call check( nf90_get_att(ncid, data_id, "yllcorner", y0))
!        call check( nf90_get_att(ncid, data_id, "cellsize", cellsize))
!        call check( nf90_get_att(ncid, data_id, "nodata", nodata))
!        call check( nf90_close(ncid) )
!      end subroutine
!
!      subroutine read_nc_int(inputfile,NY,NX,var)
!        use netcdf
!        implicit none
!      ! This is the name of the data file we will create.
!        character(*)::inputfile
!        integer::NY,NX
!        integer,dimension(:,:)::var
!        integer::ncid,var_id
!        integer,dimension(2)::dimids_2d
!        call check(nf90_open(trim(inputfile),nf90_nowrite,ncid))
!        call check(nf90_inq_varid(ncid,"data",var_id))
!        call check(nf90_get_var(ncid,var_id,var))
!        call check( nf90_close(ncid) )
!      end subroutine
!
!      subroutine read_nc_real(inputfile,NY,NX,var)
!        use netcdf
!        implicit none
!      ! This is the name of the data file we will create.
!        character(*)::inputfile
!        integer::NY,NX
!        real,dimension(:,:)::var
!        integer::ncid,var_id
!        integer,dimension(2)::dimids_2d
!        call check(nf90_open(trim(inputfile),nf90_nowrite,ncid))
!        call check(nf90_inq_varid(ncid,"data",var_id))
!        call check(nf90_get_var(ncid,var_id,var))
!        call check( nf90_close(ncid) )
!      end subroutine
!
!      subroutine write_nc_int(inputfile,NY,NX,X0,Y0,cellsize,nodata,var)
!        use netcdf
!        implicit none
!      ! This is the name of the data file we will create.
!        character(*)::inputfile
!        integer::NY,NX
!        real(8),intent(in)::X0,Y0
!        real::cellsize,nodata
!        integer::var(:,:)
!        integer::ncid,x_dimid,y_dimid,var_id
!        integer,dimension(2)::dimids_2d
!
!        call check(nf90_create(trim(inputfile), NF90_CLOBBER, ncid) )
!        ! Define the dimensions. NetCDF will hand back an ID for each.
!        call check(nf90_def_dim(ncid, "NX", NX, x_dimid))
!        call check(nf90_def_dim(ncid, "NY", NY, y_dimid))
!        dimids_2d(1)=y_dimid
!        dimids_2d(2)=x_dimid
!        ! Define the variable.
!        call check(nf90_def_var(ncid, "data", NF90_INT, dimids_2d, var_id))
!        call check(nf90_put_att(ncid,var_id,"ncols",NX))
!        call check(nf90_put_att(ncid,var_id,"nrows",NY))
!        call check(nf90_put_att(ncid,var_id,"xllcorner",x0))
!        call check(nf90_put_att(ncid,var_id,"yllcorner",y0))
!        call check(nf90_put_att(ncid,var_id,"cellsize",cellsize))
!        call check(nf90_put_att(ncid,var_id,"nodata",nodata))
!        call check(nf90_enddef(ncid))
!        call check( nf90_inq_varid(ncid, "data", var_id))
!        call check( nf90_put_var(ncid, var_id, var(1:NY,1:NX)))
!        call check( nf90_close(ncid) )
!      end subroutine
!
!      subroutine write_nc_real(inputfile,NY,NX,X0,Y0,cellsize,nodata,var)
!        use netcdf
!        implicit none
!      ! This is the name of the data file we will create.
!        character(*)::inputfile
!        integer::NY,NX
!        real(8),intent(in)::X0,Y0
!        real::cellsize,nodata
!        real::var(:,:)
!        integer::ncid,x_dimid,y_dimid,var_id
!        integer,dimension(2)::dimids_2d
!
!        call check(nf90_create(trim(inputfile), NF90_CLOBBER, ncid) )
!        ! Define the dimensions. NetCDF will hand back an ID for each.
!        call check(nf90_def_dim(ncid, "NX", NX, x_dimid))
!        call check(nf90_def_dim(ncid, "NY", NY, y_dimid))
!        dimids_2d(1)=y_dimid
!        dimids_2d(2)=x_dimid
!        ! Define the variable.
!        call check(nf90_def_var(ncid, "data", NF90_REAL, dimids_2d, var_id))
!        call check(nf90_put_att(ncid,var_id,"ncols",NX))
!        call check(nf90_put_att(ncid,var_id,"nrows",NY))
!        call check(nf90_put_att(ncid,var_id,"xllcorner",x0))
!        call check(nf90_put_att(ncid,var_id,"yllcorner",y0))
!        call check(nf90_put_att(ncid,var_id,"cellsize",cellsize))
!        call check(nf90_put_att(ncid,var_id,"nodata",nodata))
!        call check(nf90_enddef(ncid))
!        call check( nf90_inq_varid(ncid, "data", var_id))
!        call check( nf90_put_var(ncid, var_id, var))
!        call check( nf90_close(ncid) )
!      end subroutine

      subroutine read_asc_head(inputfile,NY,NX,X0,Y0,cellsize,nodata)
        implicit none
        character(*)::inputfile
        integer,intent(inout)::NY,NX
        real(8),intent(inout)::X0,Y0
        real::cellsize,nodata
        character*5 ncols,nrows
        character*12 xllcorner,yllcorner
        character*8 scellsize
        character*12 snodata
        integer::i,j        
        open(1, file=trim(inputfile), status='old')
        read(1,*) ncols,NX
        read(1,*) nrows,NY
        read(1,*) xllcorner,X0
        read(1,*) yllcorner,Y0
        read(1,*) scellsize,cellsize
        read(1,*) snodata,nodata    
        close (1)
      end subroutine
        
      subroutine read_asc_real(inputfile,NY,NX,var)
        implicit none
        character(*)::inputfile
        integer,intent(in)::NY,NX
        real,dimension(:,:)::var
        character(80)::temp
        integer::i,j        
        open(1, file=trim(inputfile), status='old')
        read(1,*) temp
        read(1,*) temp
        read(1,*) temp
        read(1,*) temp
        read(1,*) temp
        read(1,*) temp
        do i=1,NY        
          read(1,*) (var(i,j), j=1,NX)
        end do           
        close (1)
      end subroutine

      subroutine read_asc_int(inputfile,NY,NX,var)
        implicit none
        character(*)::inputfile
        integer,intent(in)::NY,NX
        integer,dimension(:,:)::var
        character(80)::temp
        integer::i,j        
        open(1, file=trim(inputfile), status='old')
        read(1,*) temp
        read(1,*) temp
        read(1,*) temp
        read(1,*) temp
        read(1,*) temp
        read(1,*) temp
        do i=1,NY        
          read(1,*) (var(i,j), j=1,NX)
        end do           
        close (1)
      end subroutine
        
      subroutine write_asc_int(outputfile,NY,NX,X0,Y0,cellsize,nodata,var)
        implicit none
        character(*)::outputfile
        integer,intent(in)::NY,NX
        real(8),intent(in)::X0,Y0
        real::cellsize,nodata
        integer::var(:,:)
        character*5 ncols,nrows
        character*12 xllcorner,yllcorner
        character*8 scellsize
        character*12 snodata
        integer::i,j             
        ncols='ncols'
        nrows='nrows'
        xllcorner='xllcorner'
        yllcorner='yllcorner'
        scellsize='cellsize'
        snodata='NODATA_value'            
        open(1, file=trim(outputfile), status='unknown')
        write(1,'(A,i16)') ncols,NX
        write(1,'(A,i16)') nrows,NY
        write(1,'(A,f15.5)') xllcorner,X0
        write(1,'(A,f15.5)') yllcorner,Y0
        write(1,'(A,f15.5)') scellsize,cellsize
        write(1,'(A,f15.5)') snodata,nodata
        do i=1,NY        
          write(1,*) (var(i,j), j=1,NX)
        end do           
        close (1)
      end subroutine

      subroutine write_asc_real(outputfile,NY,NX,X0,Y0,cellsize,nodata,var)
        implicit none
        character(*)::outputfile
        integer,intent(in)::NY,NX
        real(8),intent(in)::X0,Y0
        real::cellsize,nodata
        real::var(NY,NX)
        character*5 ncols,nrows
        character*12 xllcorner,yllcorner
        character*8 scellsize
        character*12 snodata
        integer::i,j             
        ncols='ncols'
        nrows='nrows'
        xllcorner='xllcorner'
        yllcorner='yllcorner'
        scellsize='cellsize'
        snodata='NODATA_value'            
        open(1, file=trim(outputfile), status='unknown')
        write(1,'(A,i16)') ncols,NX
        write(1,'(A,i16)') nrows,NY
        write(1,'(A,f15.5)') xllcorner,X0
        write(1,'(A,f15.5)') yllcorner,Y0
        write(1,'(A,f15.5)') scellsize,cellsize
        write(1,'(A,f15.5)') snodata,nodata
        do i=1,NY        
          write(1,*) (var(i,j), j=1,NX)
        end do           
        close (1)
      end subroutine
        


      subroutine Blocksum2(inputfile,outputfile,XWidth,YWidth,Outputtype,wsfinefile,wscoarsefile)
        implicit none
        character(*)::inputfile,outputfile !input and output files of *.asc
        character(*)::wsfinefile,wscoarsefile
        integer,intent(in)::XWidth,YWidth,Outputtype  !the width of windows in cells
        integer::InNX,InNY,OutNX,OutNY,i,j
        real(8)::X0,Y0
        real::cellsize,nodata
        real,dimension(:,:),allocatable::Indata,Outdata,windowdata
        integer,allocatable::OutdataInt(:,:),wsfine(:,:),wscoarse(:,:)
        integer::NBR,NBC,ii,jj,cont,r1,r2,c1,c2,l1,l2
        integer::NBR1,NBC1,InNX1,InNY1
        real(8)::X01,Y01
        real::cellsize1,nodata1
        real::sumdata,start,finish
        character(5)::sub

        sub=''
        call filepostfix(inputfile,l1,l2)
        sub=inputfile(l1:l2)
        if(trim(sub).eq.'asc' .or.trim(sub).eq.'txt' ) then
          call read_asc_head(inputfile,InNY,InNX,X0,Y0,cellsize,nodata)
          allocate(Indata(InNY,InNX))
          call read_asc_real(inputfile,InNY,InNX,Indata)
!        elseif(trim(sub).eq.'nc')then
!          call read_nc_head(inputfile,InNY,InNX,X0,Y0,cellsize,nodata)
!          allocate(Indata(InNY,InNX))
!          call read_nc_real(inputfile,InNY,InNX,Indata)
        else
          print*,"file type is not right in block mean"
        end if

        allocate(windowdata(YWidth,XWidth))
        NBR=floor(InNY/YWidth+0.0)+1
        NBC=floor(InNX/XWidth+0.0)+1
        allocate(Outdata(NBR,NBC))
        Outdata=nodata

        call read_asc_head(wsfinefile,InNY1,InNX1,X01,Y01,cellsize1,nodata1)
        if (InNY.eq.InNY1 .and. InNX1.eq.InNX)then
          allocate(wsfine(InNY,InNX))
          call read_asc_int(wsfinefile,InNY,InNX,wsfine)
        else
          print*,"In block mean,the ncols and nrows dismatch between the input and the watershed mask"
          stop
        end if

        print*,"reading coarse pbasin"
        call read_asc_head(wscoarsefile,NBR1,NBC1,X01,Y01,cellsize1,nodata1)
        if (NBR.eq.NBR1 .and. NBC.eq.NBC1)then
          allocate(wscoarse(NBR,NBC))
          print*,"loading data"
          call read_asc_int(wscoarsefile,NBR1,NBC1,wscoarse)
        else
          print*,"In block mean,the ncols and nrows dismatch between the input and the watershed mask"
          stop
        end if
      !$omp parallel private(i,j,r1,r2,c1,c2,sumdata,cont,ii,jj)
      !$omp do
        do i=1,NBR
          do j=1,NBC
           r1=(i-1)*YWidth+1
           r2=i*YWidth
           c1=(j-1)*YWidth+1
           c2=j*YWidth
           if(r2>InNY)r2=InNY
           if(c2>InNX)c2=InNX
           sumdata=0.0
           cont=0
           do ii=1,r2-r1+1
             do jj=1,c2-c1+1
               if(Indata(r1+ii-1,c1+jj-1).ne.-9999 .and. wscoarse(i,j).eq.wsfine(r1+ii-1,c1+jj-1)) then
                 sumdata=sumdata+Indata(r1+ii-1,c1+jj-1)
                 cont=cont+1
               end if
             end do
           end do
           if(cont>0)then
             Outdata(i,j)=sumdata
           else
             Outdata(i,j)=nodata
           end if
          end do
        end do
      !$omp end do
      !$omp end PARALLEL

        Y0=Y0+InNY*cellsize-NBR*cellsize*YWidth
        sub=''
        call filepostfix(outputfile,l1,l2)
        sub=outputfile(l1:l2)
        print*,"sub",sub
        if(Outputtype.eq.1)      then
          allocate(OutdataInt(InNY,InNX))
          OutdataInt(:,:)=int(Outdata(:,:))
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub).eq.'txt' ) then
             call write_asc_int(outputfile,NBR,NBC,X0,Y0,cellsize*YWidth,nodata,OutdataInt)
!          elseif(trim(sub(1:2)).eq.'nc')then
!             call write_nc_int(outputfile,NBR,NBC,X0,Y0,cellsize*YWidth,nodata,OutdataInt)
          else
            print*,"output file type is not right in block sum"
          end if
        else
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub).eq.'txt' ) then
            call write_asc_real(outputfile,NBR,NBC,X0,Y0,cellsize*YWidth,nodata,Outdata)
!          elseif(trim(sub(1:2)).eq.'nc')then
!            call write_nc_real(outputfile,NBR,NBC,X0,Y0,cellsize*YWidth,nodata,Outdata)
          else
           print*,"output file type is not right in block sum"
          end if
        end if
        if(allocated(OutdataInt))deallocate(OutdataInt)
        if(allocated(Outdata))deallocate(Outdata)
        if(allocated(Indata))deallocate(Indata)
        if(allocated(windowdata))deallocate(windowdata)
        if(allocated(wscoarse))deallocate(wscoarse)
        if(allocated(wsfine))deallocate(wsfine)
      end subroutine


      subroutine BlockSum(inputfile,outputfile,XWidth,YWidth,Outputtype)
        implicit none
        character(*)::inputfile,outputfile !input and output files of *.asc
        integer,intent(in)::XWidth,YWidth  !the width of windows in cells
        integer,intent(in)::Outputtype
        real,allocatable::Indata(:,:)
        real,allocatable::Outdata(:,:)
        real,allocatable::windowdata(:,:)
        integer,allocatable::OutdataInt(:,:)
        integer::NX,NY,i,j,l1,l2
        real(8)::X0,Y0
        real::cellsize,nodata
        integer::NBR,NBC,ii,jj,cont,r1,r2,c1,c2
        real::sumdata
        character(5)::sub            
            
        sub=''            
        call filepostfix(inputfile,l1,l2)
        sub=inputfile(l1:l2)            
        if(trim(sub(1:3)).eq.'asc' .or.trim(sub(1:3)).eq.'txt' ) then
          call read_asc_head(inputfile,NY,NX,X0,Y0,cellsize,nodata)
          allocate(Indata(NY,NX))
          call read_asc_real(inputfile,NY,NX,Indata)            
!        elseif(trim(sub(1:2)).eq.'nc')then
!          call read_nc_head(inputfile,NY,NX,X0,Y0,cellsize,nodata)
!          allocate(Indata(NY,NX))
!          call read_nc_real(inputfile,NY,NX,Indata)
        else
          print*,"Input file type is not right in block sum"
        end if
            
        allocate(windowdata(YWidth,XWidth))
        allocate(Outdata(NY,NX))            
        Outdata=nodata
        NBR=floor(NY/YWidth+0.0)+1
        NBC=floor(NX/XWidth+0.0)+1
      !$omp parallel private(i,j,r1,r2,c1,c2,sumdata,cont,ii,jj)
      !$omp do      
        do i=1,NBR
          do j=1,NBC
           r1=(i-1)*YWidth+1
           r2=i*YWidth
           c1=(j-1)*YWidth+1
           c2=j*YWidth
           if(r2>NY)r2=NY
           if(c2>NX)c2=NX
           windowdata(:,:)=Indata(r1:r2,c1:c2)
           sumdata=0.0
           cont=0
           do ii=1,r2-r1+1
             do jj=1,c2-c1+1
               if(Indata(r1+ii-1,c1+jj-1).ne.-9999) then
                 sumdata=sumdata+Indata(r1+ii-1,c1+jj-1)
                 cont=cont+1
               end if
             end do
           end do
           if(cont>0)then
             Outdata(r1:r2,c1:c2)=sumdata
           else
             Outdata(r1:r2,c1:c2)=nodata
           end if
          end do
        end do
      !$omp end do
      !$omp end PARALLEL
        sub=''
        call filepostfix(outputfile,l1,l2)
        sub=outputfile(l1:l2)
        if(Outputtype.eq.1)      then      
          allocate(OutdataInt(NY,NX))
          OutdataInt(:,:)=int(Outdata(:,:))            
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub).eq.'txt' ) then
             call write_asc_int(outputfile,NY,NX,X0,Y0,cellsize,nodata,OutdataInt)
!          elseif(trim(sub(1:2)).eq.'nc')then
!             call write_nc_int(outputfile,NY,NX,X0,Y0,cellsize,nodata,OutdataInt)
          else
            print*,"output file type is not right in block sum"
          end if
        else
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub).eq.'txt' ) then
            call write_asc_real(outputfile,NY,NX,X0,Y0,cellsize,nodata,Outdata)
!          elseif(trim(sub(1:2)).eq.'nc')then
!            call write_nc_real(outputfile,NY,NX,X0,Y0,cellsize,nodata,Outdata)
          else
           print*,"output file type is not right in block sum"
          end if
        end if            
        if(allocated(OutdataInt))deallocate(OutdataInt)
        if(allocated(Outdata))deallocate(Outdata)
        if(allocated(Indata))deallocate(Indata)
        if(allocated(windowdata))deallocate(windowdata)
      end subroutine

    subroutine sort2(level,n)
      implicit none
      integer::level(n,2)
      integer::n,i,j,tmp
      do i=1,n-1
           do j=i+1,n
            if(level(j,2)>level(i,2)) then
                tmp =level(i,2)
                level(i,2)=level(j,2)
                level(j,2)=tmp

                tmp =level(i,1)
                level(i,1)=level(j,1)
                level(j,1)=tmp
            end if
         end do
      end do
    end subroutine

      subroutine BlockMajor(inputfile,outputfile,XWidth,YWidth)
        implicit none
        character(*)::inputfile,outputfile !input and output files of *.asc
        integer,intent(in)::XWidth,YWidth  !the width of windows in cells
        integer(4),allocatable::Indata(:,:)
        integer(4),allocatable::Outdata(:,:)
        integer(4),allocatable::windowdata(:,:)
        integer::NX,NY,i,j,l1,l2
        real(8)::X0,Y0
        real::cellsize,nodata
        integer::NBR,NBC,ii,jj,cont,r1,r2,c1,c2,cc,cont2,flag
        real::sumdata
        character(5)::sub
        integer,allocatable::gridcount(:,:)

        sub=''
        call filepostfix(inputfile,l1,l2)
        sub=inputfile(l1:l2)
        if(trim(sub(1:3)).eq.'asc' .or.trim(sub(1:3)).eq.'txt' ) then
          call read_asc_head(inputfile,NY,NX,X0,Y0,cellsize,nodata)
          allocate(Indata(NY,NX))
          call read_asc_int(inputfile,NY,NX,Indata)
!        elseif(trim(sub(1:2)).eq.'nc')then
!          call read_nc_head(inputfile,NY,NX,X0,Y0,cellsize,nodata)
!          allocate(Indata(NY,NX))
!          call read_nc_int(inputfile,NY,NX,Indata)
        else
          print*,"Input file type is not right in block major"
        end if


        NBR=floor(NY/YWidth+0.0)+1
        NBC=floor(NX/XWidth+0.0)+1
        allocate(windowdata(YWidth,XWidth))
        allocate(Outdata(NBR,NBC))
        allocate(gridcount(YWidth*XWidth,2))
        Outdata=nodata
      !!!!$omp parallel private(i,j,r1,r2,c1,c2,sumdata,cont,ii,jj)
      !!!!$omp do
        do i=1,NBR
          do j=1,NBC
           r1=(i-1)*YWidth+1
           r2=i*YWidth
           c1=(j-1)*YWidth+1
           c2=j*YWidth
           if(r2>NY)r2=NY
           if(c2>NX)c2=NX
           windowdata(:,:)=Indata(r1:r2,c1:c2)
           sumdata=0.0

           gridcount(:,1)=-9999
           gridcount(:,2)=0
           cont2=1
           do ii=1,r2-r1+1
             do jj=1,c2-c1+1
               flag=0
               do cc=1,cont2
                 if(windowdata(ii,jj) .eq. gridcount(cc,1)) then
                   gridcount(cc,2)=gridcount(cc,2)+1
                   flag=1
                   exit
                 end if
               end do
               if(flag.eq.0)then
                 cont2=cont2+1
                 gridcount(cont2,1)=windowdata(ii,jj)
                 gridcount(cont2,2)=1
               end if
             end do
           end do
           call sort2(gridcount(1:cont2,:),cont2)
           Outdata(i,j)=gridcount(1,1)
          end do
        end do
      !!!!$omp end do
      !!!!$omp end PARALLEL

        Y0=Y0+NY*cellsize-NBR*cellsize*YWidth
        sub=''
        call filepostfix(outputfile,l1,l2)
        sub=outputfile(l1:l2)
        !print*,"NBR,NBC",NBR,NBC,X0,Y0,cellsize*YWidth,nodata
        !print*,outputfile

        if(trim(sub(1:3)).eq.'asc' .or.trim(sub).eq.'txt' ) then
             call write_asc_int(outputfile,NBR,NBC,X0,Y0,cellsize*YWidth,nodata,Outdata)
!        elseif(trim(sub(1:2)).eq.'nc')then
!             call write_nc_int(outputfile,NBR,NBC,X0,Y0,cellsize*YWidth,nodata,Outdata)
        else
            print*,"output file type is not right in block major"
        end if
        if(allocated(Outdata))deallocate(Outdata)
        if(allocated(Indata))deallocate(Indata)
        if(allocated(windowdata))deallocate(windowdata)
      end subroutine

      subroutine Blockmean2(inputfile,outputfile,XWidth,YWidth,Outputtype,wsfinefile,wscoarsefile)
        implicit none
        character(*)::inputfile,outputfile !input and output files of *.asc
        character(*)::wsfinefile,wscoarsefile
        integer,intent(in)::XWidth,YWidth,Outputtype  !the width of windows in cells
        integer::InNX,InNY,OutNX,OutNY,i,j
        real(8)::X0,Y0
        real::cellsize,nodata
        real,dimension(:,:),allocatable::Indata,Outdata,windowdata
        integer,allocatable::OutdataInt(:,:),wsfine(:,:),wscoarse(:,:)
        integer::NBR,NBC,ii,jj,cont,r1,r2,c1,c2,l1,l2
        integer::NBR1,NBC1,InNX1,InNY1
        real(8)::X01,Y01
        real::cellsize1,nodata1
        real::sumdata,start,finish
        character(5)::sub

        call read_asc_head(inputfile,InNY,InNX,X0,Y0,cellsize,nodata)
        allocate(Indata(InNY,InNX))
        call read_asc_real(inputfile,InNY,InNX,Indata)

        call read_asc_head(wsfinefile,InNY1,InNX1,X01,Y01,cellsize1,nodata1)
        if (InNY.eq.InNY1 .and. InNX1.eq.InNX)then
          allocate(wsfine(InNY,InNX))
          call read_asc_int(wsfinefile,InNY,InNX,wsfine)
        else
          print*,"In block mean,the ncols and nrows dismatch between the input and the watershed mask"
          print*,wsfinefile,InNY1,InNX1
          print*,inputfile,InNY,InNX
          stop
        end if

        allocate(windowdata(YWidth,XWidth))
        NBR=floor(InNY/YWidth+0.0)+1
        NBC=floor(InNX/XWidth+0.0)+1
        allocate(Outdata(NBR,NBC))
        Outdata=nodata

        !print*,"reading head"
        call read_asc_head(wscoarsefile,NBR1,NBC1,X01,Y01,cellsize1,nodata1)
        if (NBR.eq.NBR1 .and. NBC.eq.NBC1)then
          allocate(wscoarse(NBR,NBC))
          !print*,"reading data"
          call read_asc_int(wscoarsefile,NBR,NBC,wscoarse)
        else
          print*,"In block mean,the ncols and nrows dismatch between the input and the watershed mask"
          print*,wscoarsefile,NBR1,NBC1
          print*,inputfile,NBR,NBC
          stop
        end if


      !$omp parallel private(i,j,r1,r2,c1,c2,sumdata,cont,ii,jj)
      !$omp do
        do i=1,NBR
          do j=1,NBC
           r1=(i-1)*YWidth+1
           r2=i*YWidth
           c1=(j-1)*YWidth+1
           c2=j*YWidth
           if(r2>InNY)r2=InNY
           if(c2>InNX)c2=InNX
           sumdata=0.0
           cont=0
           do ii=1,r2-r1+1
             do jj=1,c2-c1+1
               if(Indata(r1+ii-1,c1+jj-1).ne.-9999 .and. wscoarse(i,j).eq.wsfine(r1+ii-1,c1+jj-1)) then
                 sumdata=sumdata+Indata(r1+ii-1,c1+jj-1)
                 cont=cont+1
               end if
             end do
           end do
           if(cont>0)then
             Outdata(i,j)=sumdata/cont
           else
             Outdata(i,j)=nodata
           end if
          end do
        end do
      !$omp end do
      !$omp end PARALLEL

        Y0=Y0+InNY*cellsize-NBR*cellsize*YWidth
        sub=''
        call filepostfix(outputfile,l1,l2)
        sub=outputfile(l1:l2)
        !print*,"sub",sub
        if(Outputtype.eq.1)      then
          allocate(OutdataInt(InNY,InNX))
          OutdataInt(:,:)=int(Outdata(:,:))
          call write_asc_int(outputfile,NBR,NBC,X0,Y0,cellsize*YWidth,nodata,OutdataInt)
        else
          call write_asc_real(outputfile,NBR,NBC,X0,Y0,cellsize*YWidth,nodata,Outdata)
        end if
        if(allocated(OutdataInt))deallocate(OutdataInt)
        if(allocated(Outdata))deallocate(Outdata)
        if(allocated(Indata))deallocate(Indata)
        if(allocated(windowdata))deallocate(windowdata)
        if(allocated(wscoarse))deallocate(wscoarse)
        if(allocated(wsfine))deallocate(wsfine)
      end subroutine


      subroutine Blockmean(inputfile,outputfile,XWidth,YWidth,Outputtype)
        implicit none
        character(*)::inputfile,outputfile !input and output files of *.asc
        integer,intent(in)::XWidth,YWidth,Outputtype  !the width of windows in cells
        integer::InNX,InNY,OutNX,OutNY,i,j
        real(8)::X0,Y0
        real::cellsize,nodata
        real,dimension(:,:),allocatable::Indata,Outdata,windowdata
        integer,allocatable::OutdataInt(:,:)
        integer::NBR,NBC,ii,jj,cont,r1,r2,c1,c2,l1,l2
        real::sumdata,start,finish
        character(5)::sub            
              
        sub=''
        call filepostfix(inputfile,l1,l2)
        sub=inputfile(l1:l2)            
        if(trim(sub).eq.'asc' .or.trim(sub).eq.'txt' ) then
          call read_asc_head(inputfile,InNY,InNX,X0,Y0,cellsize,nodata)
          allocate(Indata(InNY,InNX))
          call read_asc_real(inputfile,InNY,InNX,Indata)
!        elseif(trim(sub).eq.'nc')then
!          call read_nc_head(inputfile,InNY,InNX,X0,Y0,cellsize,nodata)
!          allocate(Indata(InNY,InNX))
!          call read_nc_real(inputfile,InNY,InNX,Indata)
        else
          print*,"file type is not right in block mean"
        end if                  
        allocate(windowdata(YWidth,XWidth))            
        allocate(Outdata(InNY,InNX))
        Outdata=nodata
        NBR=floor(InNY/YWidth+0.0)+1
        NBC=floor(InNX/XWidth+0.0)+1
      !$omp parallel private(i,j,r1,r2,c1,c2,sumdata,cont,ii,jj)
      !$omp do      
        do i=1,NBR
          do j=1,NBC
           r1=(i-1)*YWidth+1
           r2=i*YWidth
           c1=(j-1)*YWidth+1
           c2=j*YWidth
           if(r2>InNY)r2=InNY
           if(c2>InNX)c2=InNX
           windowdata(:,:)=Indata(r1:r2,c1:c2)
           sumdata=0.0
           cont=0
           do ii=1,r2-r1+1
             do jj=1,c2-c1+1
               if(Indata(r1+ii-1,c1+jj-1).ne.-9999) then
                 sumdata=sumdata+Indata(r1+ii-1,c1+jj-1)
                 cont=cont+1
               end if
             end do
           end do
           if(cont>0)then
             Outdata(r1:r2,c1:c2)=sumdata/cont
           else
             Outdata(r1:r2,c1:c2)=nodata
           end if
          end do
        end do
      !$omp end do
      !$omp end PARALLEL            
        sub=''
        call filepostfix(outputfile,l1,l2)
        sub=outputfile(l1:l2)
        !print*,"sub",sub
        if(Outputtype.eq.1)      then      
          allocate(OutdataInt(InNY,InNX))
          OutdataInt(:,:)=int(Outdata(:,:))            
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub).eq.'txt' ) then
             call write_asc_int(outputfile,InNY,InNX,X0,Y0,cellsize,nodata,OutdataInt)
!          elseif(trim(sub(1:2)).eq.'nc')then
!             call write_nc_int(outputfile,InNY,InNX,X0,Y0,cellsize,nodata,OutdataInt)
          else
            print*,"output file type is not right in block sum"
          end if
        else
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub).eq.'txt' ) then
            call write_asc_real(outputfile,InNY,InNX,X0,Y0,cellsize,nodata,Outdata)
!          elseif(trim(sub(1:2)).eq.'nc')then
!            call write_nc_real(outputfile,InNY,InNX,X0,Y0,cellsize,nodata,Outdata)
          else
           print*,"output file type is not right in block sum"
          end if
        end if            
        if(allocated(OutdataInt))deallocate(OutdataInt)
        if(allocated(Outdata))deallocate(Outdata)
        if(allocated(Indata))deallocate(Indata)
        if(allocated(windowdata))deallocate(windowdata)                        
      end subroutine
        
      subroutine Blockmean_memory(InNY,InNX,XWidth,YWidth,indata,outdata)
        implicit none
        integer::InNY,InNX
        integer,intent(in)::XWidth,YWidth  !the width of windows in cells
        real,dimension(InNY,InNX)::indata,outdata            
        integer::i,j
        real,dimension(YWidth,XWidth)::windowdata
        integer::NBR,NBC,ii,jj,cont,r1,r2,c1,c2,l1,l2
        real::sumdata,start,finish
        character(5)::sub
        Outdata=-9999
        NBR=floor(InNY/YWidth+0.0)+1
        NBC=floor(InNX/XWidth+0.0)+1
        do i=1,NBR
          do j=1,NBC
           r1=(i-1)*YWidth+1
           r2=i*YWidth
           c1=(j-1)*YWidth+1
           c2=j*YWidth
           if(r2>InNY)r2=InNY
           if(c2>InNX)c2=InNX
           windowdata(:,:)=Indata(r1:r2,c1:c2)
           sumdata=0.0
           cont=0
           do ii=1,YWidth
             do jj=1,XWidth
               if(windowdata(ii,jj).ne.-9999) then
                 sumdata=sumdata+windowdata(ii,jj)
                 cont=cont+1
               end if
             end do
           end do
           if(cont>0)then
             Outdata(r1:r2,c1:c2)=sumdata/cont
           else
             Outdata(r1:r2,c1:c2)=-9999
           end if
          end do
        end do            
      end subroutine
        

      subroutine Blocksum_memory(InNY,InNX,XWidth,YWidth,indata,outdata)
        implicit none
        integer::InNY,InNX
        integer,intent(in)::XWidth,YWidth  !the width of windows in cells
        real,dimension(InNY,InNX)::indata,outdata,outdata2      
        integer::i,j
        real,dimension(YWidth,XWidth)::windowdata
        integer::NBR,NBC,ii,jj,cont,r1,r2,c1,c2,l1,l2
        real::sumdata,start,finish
        character(5)::sub
        Outdata(:,:)=-9999
        NBR=floor(InNY/YWidth+0.0)+1
        NBC=floor(InNX/XWidth+0.0)+1
      !$omp parallel private(i,j,r1,r2,c1,c2,sumdata,cont,ii,jj)
      !$omp do            
        do i=1,NBR
          do j=1,NBC
           r1=(i-1)*YWidth+1
           r2=i*YWidth
           c1=(j-1)*YWidth+1
           c2=j*YWidth
           if(r2>InNY)r2=InNY
           if(c2>InNX)c2=InNX
           sumdata=0.0
           cont=0
           do ii=1,r2-r1+1
             do jj=1,c2-c1+1
               if(Indata(r1+ii-1,c1+jj-1).ne.-9999) then
                 sumdata=sumdata+Indata(r1+ii-1,c1+jj-1)
                 cont=cont+1
               end if
             end do
           end do
           if(cont>0)then
             Outdata(r1:r2,c1:c2)=sumdata
           else
             Outdata(r1:r2,c1:c2)=-9999
           end if
          end do
        end do
      !$omp end do
      !$omp end PARALLEL              
      end subroutine
        
      subroutine Resample(inputfile,outputfile,outcellsize,InterType,Outputtype)
        implicit none
        character(*),intent(in)::inputfile,outputfile
        real,intent(in)::outcellsize
        integer,intent(in)::InterType! 1 for nearest, 2 for binear
        integer,intent(in)::Outputtype   ! 1 for integer variable; and other for real variable
        integer::InNX,InNY,OutNX,OutNY,i,j
        real(8)::X0,Y0
        real::Insize,nodata
        real,allocatable::Indata(:,:)
        real,allocatable::Outdata(:,:)
        integer,allocatable::OutdataInt(:,:)
        integer::InCOL,INROW,InCOL2,InRow2,l1,l2
        real::x1,x2,y1,y2,R1,R2,IX,IY,DeltaX,DeltaY
        real::start,finish
        character(5)::sub
            
        sub=''            
        call filepostfix(inputfile,l1,l2)
        sub=inputfile(l1:l2)      
        if(trim(sub(1:3))=="asc") then 
          call read_asc_head(inputfile,InNY,InNX,X0,Y0,Insize,nodata)
          allocate(Indata(InNY,InNX))
          call read_asc_real(inputfile,InNY,InNX,Indata)            
!        elseif(trim(sub(1:2)).eq.'nc')then
!          call read_nc_head(inputfile,InNY,InNX,X0,Y0,Insize,nodata)
!          allocate(Indata(InNY,InNX))
!          call read_nc_real(inputfile,InNY,InNX,Indata)
        else
          print*,"Input file type is not right in resampling"
        end if
        OutNX=int(ANINT(Insize*InNX/outcellsize))
        OutNY=int(ANINT(Insize*InNY/outcellsize))
        allocate(Outdata(OutNY,OutNX))
        Outdata(:,:)=nodata
      !$omp parallel private(i,j,IX,IY,InCOL,InRow,DeltaX,DeltaY,INCOL2,INROW2,X1,X2,Y1,Y2,R1,R2)
      !$omp do
        DO i=1,OutNY
          DO j=1,OutNX
            IX=x0+(j-0.5)*outcellsize
            IY=y0+Insize*InNY-(i-0.5)*outcellsize!y0+(OutNY-i+0.5)*outcellsize
            InCOL=floor((IX-x0)/Insize)+1
            InRow=InNY-floor((IY-y0)/Insize)
            DeltaX=x0+InCOL*Insize-IX
            DeltaY=y0+(InNY-InRow+1)*Insize-IY
            !print*,DeltaX,DeltaY
            IF(DeltaX<Insize/2.0) then
                     INCOL2=InCOL+1
            ELSE
               INCOL2=INCOL-1
            END IF                  
            IF(DeltaY<Insize/2.0) then
                     INROW2=InROW-1
            ELSE
               INROW2=InROW+1
            END IF                  
            IF(INROW2==0) INROW2=InROW+1
            IF(INCOL2==0) INCOL2=INCOL+1
            IF(INROW2>InNY) INROW2=InROW-1
            IF(INCOL2>InNX) INCOL2=INCOL-1                  
            X1=x0+(InCOL-0.5)*Insize
            X2=x0+(InCOL2-0.5)*Insize
            Y1=y0+(InNY-INROW+0.5)*Insize
            Y2=y0+(InNY-InRow2+0.5)*Insize
            IF(InterType.eq.1) then
              IF(Indata(InRow,InCOL).eq.nodata)then
                Outdata(i,j)=nodata
              else
                Outdata(i,j)=Indata(InRow,InCOL)
              end if
            elseif(InterType.eq.2) then
              IF(Indata(InRow,InCOL).eq.nodata)then
                Outdata(i,j)=nodata
              else
                if(Indata(InRow,INCOL2).eq.nodata .or. &
                 Indata(INROW2,InCOL).eq.nodata.or. Indata(INROW2,INCOL2).eq.nodata)then
                  Outdata(i,j)=Indata(InRow,InCOL)
                else
                  IF(DeltaX.eq.Insize/2.0 .or. DeltaY.eq.Insize/2.0) then
                    IF(DeltaX.eq.Insize/2.0 .and. DeltaY.eq.Insize/2.0)then
                      Outdata(i,j)=Indata(InRow,InCOL)
                    else
                      IF(DeltaX.eq.Insize/2.0) then
                        Outdata(i,j)=((Y2-IY)*Indata(INROW,InCOL)+(IX-Y1)*Indata(INROW2,INCOL))/(Y2-Y1)
                      end if
                      IF(DeltaY.eq.Insize/2.0) then
                        Outdata(i,j)=((X2-IX)*Indata(INROW,InCOL)+(IX-X1)*Indata(INROW,INCOL2))/(x2-x1)
                      end if
                    end if
                  else
                    R1=((X2-IX)*Indata(InRow,InCOL)+(IX-X1)*Indata(InRow,INCOL2))/(x2-x1)
                    R2=((X2-IX)*Indata(INROW2,InCOL)+(IX-X1)*Indata(INROW2,INCOL2))/(x2-x1)
                    Outdata(i,j)=((Y2-IY)*R1+(IY-Y1)*R2)/(Y2-Y1)
                  end if
                end if
              end if
            end if
            if(Outdata(i,j)<0 .and. Outdata(i,j).ne.-9999) then
            print*,"X,Y",IX,IY
            print*,X1,Y1,INROW,INCOL,Indata(InRow,INCOL)
            print*,X2,Y1,INROW,INCOL2,Indata(InRow,INCOL2)
            print*,X1,Y2,INROW2,INCOL,Indata(INROW2,INCOL)
            print*,X2,Y2,INROW2,INCOL2,Indata(INROW2,INCOL2)
            print*,Outdata(i,j),i,j,x0+(InCOL-0.5)*Insize,DeltaX,DeltaY
            end if                  
          END DO
        END DO
      !$omp end do
      !$omp end PARALLEL  
        sub=''
        call filepostfix(outputfile,l1,l2)
        sub=outputfile(l1:l2)
        Y0=Y0+InNY*Insize-OutNY*outcellsize
        if(Outputtype.eq.1) then
          allocate(OutdataInt(OutNY,OutNX))
          OutdataInt(:,:)=int(Outdata(:,:))
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub(1:3)).eq.'txt' ) then
            call write_asc_int(outputfile,OutNY,OutNX,X0,Y0,outcellsize,nodata,OutdataInt)
!          elseif(trim(sub(1:2)).eq.'nc')then
!            call write_nc_int(outputfile,OutNY,OutNX,X0,Y0,outcellsize,nodata,OutdataInt)
          else
            print*,"Output file type is not right in resampling"
          end if
        else
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub(1:3)).eq.'txt' ) then
            call write_asc_real(outputfile,OutNY,OutNX,X0,Y0,outcellsize,nodata,Outdata)
!          elseif(trim(sub(1:2)).eq.'nc')then
!            call write_nc_real(outputfile,OutNY,OutNX,X0,Y0,outcellsize,nodata,Outdata)
          else
            print*,"Output file type is not right in resampling"
          end if
        end if            
        if(allocated(Outdata))deallocate(Outdata)
        if(allocated(OutdataInt))deallocate(OutdataInt)
        if(allocated(Indata))deallocate(Indata)
      end subroutine

        
      subroutine ExtByMask(inputfile,maskfile,outputfile,OutputType)
        implicit none
        character(*)::inputfile,maskfile,outputfile
        integer::OutputType
        integer::InNX1,InNY1,InNX2,InNY2
        real::nodata1,nodata2,Insize1,Insize2
        real(8)::x01,y01,x02,y02
        real,dimension(:,:),allocatable::Indata,Outdata,maskdata
        integer,dimension(:,:),allocatable::OutdataInt
        integer::i,j,INRow,InCOL,l1,l2
        real::xx,yy,start
        character(5)::sub

        sub=''            
        call filepostfix(inputfile,l1,l2)
        sub=inputfile(l1:l2)      
        if(trim(sub(1:3)).eq.'asc' .or.trim(sub).eq.'txt' ) then
          call read_asc_head(inputfile,InNY1,InNX1,x01,y01,Insize1,nodata1)
          allocate(Indata(InNY1,InNX1))
          call read_asc_real(inputfile,InNY1,InNX1,Indata)      
!        elseif(trim(sub(1:2)).eq.'nc')then
!          call read_nc_head(inputfile,InNY1,InNX1,x01,y01,Insize1,nodata1)
!          allocate(Indata(InNY1,InNX1))
!          call read_nc_real(inputfile,InNY1,InNX1,Indata)
        else
          print*,"file type is not right for input file in ExtByMask"
        end if            
        sub=''            
        call filepostfix(maskfile,l1,l2)
        sub=maskfile(l1:l2)      
        if(trim(sub).eq.'asc' .or.trim(sub).eq.'txt' ) then
          call read_asc_head(maskfile,InNY2,InNX2,x02,y02,Insize2,nodata2)
          allocate(maskdata(InNY2,InNX2))
          call read_asc_real(maskfile,InNY2,InNX2,maskdata)
!        elseif(trim(sub).eq.'nc')then
!          call read_nc_head(maskfile,InNY2,InNX2,x02,y02,Insize2,nodata2)
!          allocate(maskdata(InNY2,InNX2))
!          call read_nc_real(maskfile,InNY2,InNX2,maskdata)
        else
          print*,"file type is not right for maskfile in ExtByMask"
        end if

        allocate(Outdata(InNY2,InNX2))      
      !$omp parallel private(i,j,xx,yy,InCOL,InRow)
      !$omp do            
        do i=1,InNY2
          do j=1,InNX2
             if(maskdata(i,j).eq.nodata2) then
               Outdata(i,j)=nodata2
             else
               xx=x02+(j-1+0.5)*Insize2
               yy=y02+(InNY2-i+0.5)*Insize2
               InCOL=floor((xx-x01)/Insize1)+1
               InRow=InNY1-floor((yy-y01)/Insize1)
               Outdata(i,j)=Indata(InRow,InCOL)
               if(Outdata(i,j).eq.nodata2 .and. maskdata(i,j).ne.1) then
                 print*,"warning: there is no data in the watershed at row: ",InRow,"and col:",InCOL
               end if
            end if
          end do
        end do
      !$omp end do
      !$omp end PARALLEL  
        call filepostfix(outputfile,l1,l2)
        sub=outputfile(l1:l2)              
        if(OutputType.eq.1) then
          allocate(OutdataInt(InNY2,InNX2))
          OutdataInt(:,:)=int(Outdata(:,:))
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub(1:3)).eq.'txt' ) then
            call write_asc_int(outputfile,InNY2,InNX2,X02,Y02,Insize2,nodata2,OutdataInt)
!          elseif(trim(sub(1:2)).eq.'nc')then
!            call write_nc_int(outputfile,InNY2,InNX2,X02,Y02,Insize2,nodata2,OutdataInt)
          else
           print*,"file type is not right for output file in ExtByMask"
          end if
        else
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub(1:3)).eq.'txt' ) then
            call write_asc_real(outputfile,InNY2,InNX2,X02,Y02,Insize2,nodata2,Outdata)
!          elseif(trim(sub(1:2)).eq.'nc')then
!            call write_nc_real(outputfile,InNY2,InNX2,X02,Y02,Insize2,nodata2,Outdata)
          else
           print*,"file type is not right for output file in ExtByMask"
          end if
        end if              
        if(allocated(Outdata)) deallocate(Outdata)
        if(allocated(Indata)) deallocate(Indata)      
        if(allocated(maskdata))deallocate(maskdata)
        if(allocated(OutdataInt))deallocate(OutdataInt)
      end subroutine        
        
      subroutine slope(dem,slp,NX,NY,cellsize,nodata)
      implicit none
      integer,intent(in)::NX,NY
      real,intent(in)::cellsize,nodata
      real,dimension(NY,NX),intent(in)::dem
      real,dimension(NY,NX),intent(inout)::slp
      integer::row,col,rr,cc,valid
      real::dz_dx,dz_dy,rise_run
      slp(:,:)=nodata
      !$omp parallel private(row,col,valid,rr,cc,dz_dx,dz_dy,rise_run)
      !$omp do
      do row=2,NY-1
          do col=2,NX-1
             !print*,row,col
             if(dem(row,col).ne.nodata) then
                valid=1
                do rr=-1,1,1
                   do cc=-1,1,1
                     if(dem(row+rr,col+cc).eq.nodata) valid=0
                   end do
                end do
                if(valid.eq.1)then
                 dz_dx =(dem(row-1,col+1)+2*dem(row,col+1)+dem(row+1,col+1)-&
                   (dem(row-1,col-1)+2*dem(row,col-1)+dem(row+1,col-1)))/8.0/cellsize
                 dz_dy =(dem(row+1,col-1)+2*dem(row+1,col)+dem(row+1,col+1)-&
                   (dem(row-1,col-1)+2*dem(row-1,col)+dem(row-1,col+1)))/8.0/cellsize
                 slp(row,col)=sqrt(dz_dx**2.0+dz_dy**2.0)
                 !slp(row,col)=ATAN (rise_run)*57.29578
                end if
             end if
          end do
      end do      
      !$omp end do
      !$omp end PARALLEL  
      end subroutine
        
      subroutine aspect(dem,asp,NX,NY,cellsize,nodata)
      implicit none
      integer,intent(in)::NX,NY
      real,intent(in)::cellsize,nodata
      real,dimension(NY,NX),intent(in)::dem
      real,dimension(NY,NX),intent(inout)::asp
      integer::row,col,rr,cc,valid
      real::dz_dx,dz_dy,cell
        asp=nodata
      do row=2,NY-1
          do col=2,NX-1
             if(dem(row,col).ne.nodata) then
                valid=1
                do rr=-1,1,1
                   do cc=-1,1,1
                     if(dem(row+rr,col+cc).eq.nodata) valid=0
                   end do
                end do
                if(valid.eq.1)then
                 dz_dx =(dem(row-1,col+1)+2*dem(row,col+1)+dem(row+1,col+1)-&
                   (dem(row-1,col-1)+2*dem(row,col-1)+dem(row+1,col-1)))/8.0/cellsize
                 dz_dy =(dem(row+1,col-1)+2*dem(row+1,col)+dem(row+1,col+1)-&
                   (dem(row-1,col-1)+2*dem(row-1,col)+dem(row-1,col+1)))/8.0/cellsize
                 asp(row,col)=57.29578 * atan2 (dz_dy, -dz_dx)
                 if(asp(row,col).lt.0) then
                   cell=90.-asp(row,col)
                 elseif(asp(row,col).gt.90)  then
                   cell = 360.0 - asp(row,col) + 90.0
                 else
                   cell = 90.0 - asp(row,col)
                 end if
                 asp(row,col)=cell
                 if((dz_dx**2.0+dz_dy**2.0).lt. 0.000001) then 
                   asp(row,col)=-1
                 end if
                end if
             end if
          end do
      end do        
      end subroutine
        
      subroutine flowlength(dir,flowlen,NX,NY,cellsize,nodata)
      implicit none
      integer,intent(in)::NX,NY
      integer,intent(in),dimension(NY,NX)::dir
      real,intent(inout),dimension(NY,NX)::flowlen

      real,intent(in)::cellsize,nodata
      integer::row,col
        
      allocate(dirtemp(NY,NX))
      allocate(flowlentemp(NY,NX))
      dirtemp(1:NY,1:NX)=dir(1:NY,1:NX)
      flowlentemp(1:NY,1:NX)=nodata
      !$omp parallel private(row,col)
      !$omp do
      do row=1,NY
          do col=1,NX              
             if(dirtemp(row,col).ne.nodata) then
               call calflowlen(row,col,NX,NY,cellsize,nodata)
             end if
          end do
      end do
      !$omp end do
      !$omp end PARALLEL  
      flowlen(1:NY,1:NX)=flowlentemp(1:NY,1:NX) 
      deallocate(dirtemp,flowlentemp)
      end subroutine

      recursive subroutine calflowlen(row,col,NX,NY,cellsize,nodata)
      implicit none
      integer,intent(in)::row,col,NX,NY
      real,intent(in)::cellsize,nodata
      if(flowlentemp(row,col).ne.nodata) return 
      if(dirtemp(row,col).eq.1) then
        if(col+1>NX) then
          flowlentemp(row,col)=0.0
          return
        else
          if(dirtemp(row,col+1).ne.nodata) then
            if(flowlentemp(row,col+1).eq.nodata) call calflowlen(row,col+1,NX,NY,cellsize,nodata)
            flowlentemp(row,col)=flowlentemp(row,col+1)+cellsize
          else
            flowlentemp(row,col)=0.0
          end if
        end if
      end if

      if(dirtemp(row,col).eq.2) then
        if(col+1>NX .or. row+1>NY) then
          flowlentemp(row,col)=0.0
          return
        else
          if(dirtemp(row+1,col+1).ne.nodata) then
            if(flowlentemp(row+1,col+1).eq.nodata) call calflowlen(row+1,col+1,NX,NY,cellsize,nodata)
            flowlentemp(row,col)=flowlentemp(row+1,col+1)+cellsize*sqrt(2.0)
          else
            flowlentemp(row,col)=0.0
          end if
        end if
      end if        

      if(dirtemp(row,col).eq.4) then
        if(row+1>NY) then
          flowlentemp(row,col)=0.0
          return
        else
          if(dirtemp(row+1,col).ne.nodata) then
            if(flowlentemp(row+1,col).eq.nodata) call calflowlen(row+1,col,NX,NY,cellsize,nodata)
            flowlentemp(row,col)=flowlentemp(row+1,col)+cellsize
          else
            flowlentemp(row,col)=0.0
          end if
        end if
      end if
        
      if(dirtemp(row,col).eq.8) then
        if(col-1<0 .or. row+1>NY) then
          flowlentemp(row,col)=0.0
          return
        else
          if(dirtemp(row+1,col-1).ne.nodata) then
            if(flowlentemp(row+1,col-1).eq.nodata) call calflowlen(row+1,col-1,NX,NY,cellsize,nodata)
            flowlentemp(row,col)=flowlentemp(row+1,col-1)+cellsize*sqrt(2.0)
          else
            flowlentemp(row,col)=0.0
          end if
        end if
      end if

      if(dirtemp(row,col).eq.16) then
        if(col-1<0) then
          flowlentemp(row,col)=0.0
          return
        else
          if(dirtemp(row,col-1).ne.nodata) then
            if(flowlentemp(row,col-1).eq.nodata) call calflowlen(row,col-1,NX,NY,cellsize,nodata)
            flowlentemp(row,col)=flowlentemp(row,col-1)+cellsize
          else
            flowlentemp(row,col)=0.0
          end if
        end if
      end if        

      if(dirtemp(row,col).eq.32) then
        if(col-1<0 .or. row-1<0) then
          flowlentemp(row,col)=0.0
          return
        else
          if(dirtemp(row-1,col-1).ne.nodata) then
            if(flowlentemp(row-1,col-1).eq.nodata) call calflowlen(row-1,col-1,NX,NY,cellsize,nodata)
            flowlentemp(row,col)=flowlentemp(row-1,col-1)+cellsize*sqrt(2.0)
          else
            flowlentemp(row,col)=0.0
          end if
        end if
      end if      

      if(dirtemp(row,col).eq.64) then
        if(row-1<0) then
          flowlentemp(row,col)=0.0
          return
        else
          if(dirtemp(row-1,col).ne.nodata) then
            if(flowlentemp(row-1,col).eq.nodata)call calflowlen(row-1,col,NX,NY,cellsize,nodata)
            flowlentemp(row,col)=flowlentemp(row-1,col)+cellsize
          else
            flowlentemp(row,col)=0.0
          end if
        end if
      end if

      if(dirtemp(row,col).eq.128) then
        if(row-1<0 .or. col+1>NX) then
          flowlentemp(row,col)=0.0
          return
        else
          if(dirtemp(row-1,col+1).ne.nodata) then
            if(flowlentemp(row-1,col+1).eq.nodata)call calflowlen(row-1,col+1,NX,NY,cellsize,nodata)
            flowlentemp(row,col)=flowlentemp(row-1,col+1)+cellsize*sqrt(2.0)
          else
            flowlentemp(row,col)=0.0
          end if
        end if
      end if
      end subroutine
        
      subroutine smallextent(infile,outfile,nodataWidth,OutputType)
        implicit none
        character(*)::infile,outfile !input and output files of *.asc
        integer,intent(in)::nodataWidth  !the width of nodata rows or columns surrounding the active region
        integer,intent(in)::OutputType !OutputType of 1 for integer variable,other for real
        integer::l1,l2
        character(5)::sub
        integer::InNY,InNX,OutNY,OutNX
        real(8)::X01,Y01,X02,Y02
        real::cellsize,nodata
        real,allocatable::Indata(:,:),Outdata(:,:)
        integer,allocatable::OutdataInt(:,:)
        integer::r1,r2,c1,c2,i,j,colmin
            
        sub=''            
        call filepostfix(infile,l1,l2)
        sub=infile(l1:l2)
        if(trim(sub(1:3))=="asc") then 
          call read_asc_head(infile,InNY,InNX,X01,Y01,cellsize,nodata)
          allocate(Indata(InNY,InNX))
          call read_asc_real(infile,InNY,InNX,Indata)            
!        elseif(trim(sub(1:2)).eq.'nc')then
!          call read_nc_head(infile,InNY,InNX,X01,Y01,cellsize,nodata)
!          allocate(Indata(InNY,InNX))
!          call read_nc_real(infile,InNY,InNX,Indata)
        else
          print*,"Input file type is not right in small extent"
        end if
        r1=-999
        r2=-999
        c1=99999
        c2=-999            
        do i=1,InNY
          do j=1,InNX
             if(Indata(i,j).ne.nodata .and. r1.eq.-999) r1=i
             if(Indata(i,j).ne.nodata .and. c1>j) c1=j
             if(Indata(i,j).ne.nodata) r2=i
             if(Indata(i,j).ne.nodata .and. j>c2) c2=j
          end do
        end do            
        OutNY=r2-r1+2*nodataWidth+1
        OutNX=c2-c1+2*nodataWidth+1
        allocate(Outdata(OutNY,OutNX))
        Outdata(:,:)=nodata
        Outdata(nodataWidth+1:OutNY-nodataWidth,nodataWidth+1:OutNX-nodataWidth)=Indata(r1:r2,c1:c2)
        X02=X01+(c1-1-nodataWidth)*cellsize
        Y02=Y01+(InNY-r2-nodataWidth)*cellsize
        sub=''
        call filepostfix(outfile,l1,l2)
        sub=outfile(l1:l2)
        if(OutputType.eq.1) then
          allocate(OutdataInt(OutNY,OutNX))
          OutdataInt(:,:)=int(Outdata(:,:))
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub(1:3)).eq.'txt' ) then
            call write_asc_int(outfile,OutNY,OutNX,X02,Y02,cellsize,nodata,OutdataInt)
!          elseif(trim(sub(1:2)).eq.'nc')then
!            call write_nc_int(outfile,OutNY,OutNX,X02,Y02,cellsize,nodata,OutdataInt)
          else
            print*,"Output file type is not right in resampling"
          end if
        else
          if(trim(sub(1:3)).eq.'asc' .or.trim(sub(1:3)).eq.'txt' ) then
            call write_asc_real(outfile,OutNY,OutNX,X02,Y02,cellsize,nodata,Outdata)
!          elseif(trim(sub(1:2)).eq.'nc')then
!             call write_nc_real(outfile,OutNY,OutNX,X02,Y02,cellsize,nodata,Outdata)
          else
            print*,"Output file type is not right in resampling"
          end if
        end if
        if(allocated(Outdata)) deallocate(Outdata)
        if(allocated(Indata)) deallocate(Indata)      
      end subroutine        
  end module
