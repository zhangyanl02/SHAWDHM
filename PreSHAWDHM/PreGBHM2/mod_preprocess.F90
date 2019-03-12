!@  Developed by Dr. Yanlin Zhang,Hunan University of Science and Technology, Xiangtan, Hunan province,411201
!@  Apr 2017
  module mod_preprocess
        character*(200)::midfileformat,hydro_para_dir,subcatchment2file
        character*(200)::model_para_dir,DEM_fine,DemFile,Outlet,ThresHold,lufine,cpus,skyviewfile
        integer::ModelRes,DemRes,smallestWSshed
        integer::NumberOfSubWatersheds,river_parameter_only
        real::dx_max
        integer,parameter::maxlink=15,maxlevel=10
        !integer::maxlink,maxlevel
        !integer::maxlevel

  contains
    subroutine TauDem()
        use gisutil
        implicit none
        integer::fileunit,openstat
        logical::file_e
        integer::NY,NX,i,j,d1,d2
        real(8)::X0,Y0
        real::cellsize,nodata
        integer,allocatable::IndataInt(:,:),ws(:,:),dir(:,:),ws2(:,:)
        real,allocatable::IndataReal(:,:)
        integer::l1,l2,a1,a2
        character(5)::sub
        character(80)::infile,outfile,maskfile
        integer,allocatable::subrename(:,:),mask(:,:),nextsub2(:,:)
        integer::nsub,tmpi
        integer,allocatable::gridscount(:)
        integer::r,c,tmp,nextsub,eliminatednum,redo
        character(200)::resolution,inputfile,outputfile

        sub=''            
        call strlen(DemFile,l1,l2)
        write(resolution,*) DemRes
        call strlen(resolution,a1,a2)

        call system("gdal_translate -of GTiff -tr "//resolution(a1:a2)//" "//resolution(a1:a2)//" "&
                  //DemFile(l1:l2)//" ResampleDEM.tif")
        inquire(file='Fil_fine.tif', exist=file_e)
        if(file_e) call system("rm Fil_fine.tif")
        call system("mpiexec -n "//cpus//" pitremove -z ResampleDEM.tif -fel Fil_fine.tif")
        inquire(file='Fil_fine.tif', exist=file_e)
        if(file_e) call system("gdal_translate -of AAIGrid Fil_fine.tif  FilDEM_fine.asc")
        inquire(file='flowdir_fine.tif', exist=file_e)
        if(file_e) call system("rm flowdir_fine.tif")
        inquire(file='d8slope_fine.tif', exist=file_e)
        if(file_e) call system("rm d8slope_fine.tif")
        call system("mpiexec -n "//cpus//" d8flowdir -p flowdir_fine.tif -sd8 d8slope_fine.tif -fel Fil_fine.tif")
         
        inquire(file='acc_fine.tif', exist=file_e)
        if(file_e) call system("rm acc_fine.tif")
        call system("mpiexec -n "//cpus//" aread8 -p flowdir_fine.tif -ad8 acc_fine.tif -nc")
         
        inquire(file='gord.tif', exist=file_e)
        if(file_e) call system("rm gord.tif")
        inquire(file='plen_fine.tif', exist=file_e)
        if(file_e) call system("rm plen_fine.tif")
        inquire(file='tlen_fine.tif', exist=file_e)
        if(file_e) call system("rm tlen_fine.tif")
        call system("mpiexec -n "//cpus//" gridnet -p flowdir_fine.tif -gord gord.tif -plen plen_fine.tif -tlen tlen_fine.tif")
         
        inquire(file='stream_fine.tif', exist=file_e)
        if(file_e) call system("rm stream_fine.tif")
        call system("mpiexec -n "//cpus//" threshold -ssa acc_fine.tif -src stream_fine.tif -thresh "//ThresHold)
         
        inquire(file='ADJOutlet.shp', exist=file_e)
        if(file_e) then
          call system("rm ADJOutlet.shp")
          call system("rm ADJOutlet.shx")
          call system("rm ADJOutlet.prj")
          call system("rm ADJOutlet.dbf")
        end if         
        call system("mpiexec -n "//cpus//" moveoutletstostreams -p flowdir_fine.tif -src stream_fine.tif -o "&
            //Outlet//" -om ADJOutlet.shp")

        inquire(file='acc2__fine.tif', exist=file_e)
        if(file_e) call system("rm acc2__fine.tif")
        call system("mpiexec -n "//cpus//" aread8 -p flowdir_fine.tif -o ADJOutlet.shp -ad8 acc2__fine.tif -nc")

        inquire(file='stream2__fine.tif', exist=file_e)
        if(file_e) call system("rm stream2__fine.tif")
        call system("mpiexec -n "//cpus//" threshold -ssa acc2__fine.tif -src stream2__fine.tif -thresh "//ThresHold)

        inquire(file='ord_fine.tif', exist=file_e)
        if(file_e) call system("rm ord_fine.tif")
        inquire(file='tree.txt', exist=file_e)
        if(file_e) call system("rm tree.txt")
        inquire(file='coord.txt', exist=file_e)
        if(file_e) call system("rm coord.txt")      
        inquire(file='ws_fine.tif', exist=file_e)
        if(file_e) call system("rm ws_fine.tif")
        inquire(file='tree.txt', exist=file_e)
        if(file_e) call system("rm tree.txt")               
        inquire(file='net.shp', exist=file_e)         
        if(file_e) then
          call system("rm net.shp")
          call system("rm net.shx")
          call system("rm net.prj")
          call system("rm net.dbf")
        end if                 
        call system("mpiexec -n "//cpus//" streamnet -fel Fil_fine.tif -p flowdir_fine.tif -ad8 acc_fine.tif -src &
           stream2__fine.tif -o ADJOutlet.shp -ord ord_fine.tif -tree tree.txt -coord coord.txt -net net.shp -w ws_fine.tif")
      
        inquire(file='ws_fine.tif', exist=file_e)
        if(file_e) call system("gdal_translate -of AAIGrid -ot Int16 ws_fine.tif  ws_fine.asc")
        inquire(file='acc_fine.tif', exist=file_e)
        if(file_e) call system("gdal_translate -of AAIGrid -ot Int16 acc_fine.tif  acc_fine.asc")
        inquire(file='flowdir_fine.tif', exist=file_e)
        if(file_e) call system("gdal_translate -of AAIGrid -ot Int16 flowdir_fine.tif  dir_fine.asc")
        inquire(file='stream2__fine.tif', exist=file_e)
        if(file_e) call system("gdal_translate -of AAIGrid -ot Int16 stream2__fine.tif  net_fine.asc")


        call read_asc_head("ws_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
        allocate(IndataINT(NY,NX))
        call read_asc_int("ws_fine.asc",NY,NX,IndataINT)
        if(nodata.ne.-9999) then
          do i=1,NY
            do j=1,NX
              if(IndataINT(i,j).eq.nodata) then
                IndataINT(i,j)=-9999
              endif
            end do
          end do
        end if
        tmpi=0
        do i=1,NY
          do j=1,NX
            if(IndataINT(i,j).ne.-9999) then
              if( IndataINT(i,j)>tmpi) then
                nsub=IndataINT(i,j)
                tmpi=IndataINT(i,j)
              endif
            endif
          end do
        end do

        allocate(subrename(nsub,2))
        subrename(:,:)=-9999
        tmpi=1
        do i=1,NY
          do j=1,NX
            if(IndataINT(i,j).ne.-9999 .and. subrename(IndataINT(i,j),1).eq.-9999) then
               subrename(IndataINT(i,j),1)=IndataINT(i,j)
               subrename(IndataINT(i,j),2)=tmpi
               tmpi=tmpi+1 
            endif
          end do
        end do
        do i=1,NY
          do j=1,NX
            if(IndataINT(i,j).ne.-9999 .and. subrename(IndataINT(i,j),2).ne.-9999) then
              IndataINT(i,j)=subrename(IndataINT(i,j),2)
            endif
          end do
        end do

        nodata=-9999
        call write_asc_int("ws_fine.asc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
        call write_asc_int("ws_fine1.asc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
        if(allocated(subrename))deallocate(subrename)
        if(allocated(IndataINT))deallocate(IndataINT)


        call read_asc_head("acc_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
        if(nodata.ne.-9999) then
          allocate(IndataINT(NY,NX))
          call read_asc_int("acc_fine.asc",NY,NX,IndataINT)
          do i=1,NY
            do j=1,NX
               if(IndataINT(i,j).eq.nodata)IndataINT(i,j)=-9999
            end do
          end do
          nodata=-9999
          call write_asc_int("acc_fine.asc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
          deallocate(IndataINT)
        end if

        call read_asc_head("dir_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
        allocate(IndataINT(NY,NX))
        call read_asc_int("dir_fine.asc",NY,NX,IndataINT)
        do i=1,NY
          do j=1,NX
            if(IndataINT(i,j).eq.nodata) then
              IndataINT(i,j)=-9999
            else
              if(IndataINT(i,j).eq.2) IndataINT(i,j)=128
              if(IndataINT(i,j).eq.3) IndataINT(i,j)=64
              if(IndataINT(i,j).eq.4) IndataINT(i,j)=32
              if(IndataINT(i,j).eq.5) IndataINT(i,j)=16
              if(IndataINT(i,j).eq.7) IndataINT(i,j)=4
              if(IndataINT(i,j).eq.8) IndataINT(i,j)=2
              if(IndataINT(i,j).eq.6) IndataINT(i,j)=8
            end if
          end do
        end do
        nodata=-9999
        call write_asc_int("dir_fine.asc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
        if(allocated(IndataINT))deallocate(IndataINT)



!!!!!there exist a big bug here
!!!!!eliminate small watershed here
        call read_asc_head("ws_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
        allocate(ws(NY,NX))
        allocate(dir(NY,NX))
        call read_asc_int("ws_fine.asc",NY,NX,ws)
        call read_asc_int("dir_fine.asc",NY,NX,dir)
        allocate(gridscount(nsub))            
        gridscount(:)=0
        do i=1,NY
          do j=1,NX
            if(ws(i,j).ne.-9999) then
              gridscount(ws(i,j))=gridscount(ws(i,j))+1
            endif
          end do
        end do
!        eliminatednum=0
!        do i=1,nsub
!          if(gridscount(i)<smallestWSshed*1000.0/DemRes*1000.0/DemRes) then
!            eliminatednum=eliminatednum+1
!            nextsub=0
!            do r=1,NY
!               do c=1,NX
!                  if(ws(r,c).eq.i) then
!                     if(dir(r,c).eq.1)   tmp=ws(r,c+1)
!                     if(dir(r,c).eq.2)   tmp=ws(r+1,c+1)
!                     if(dir(r,c).eq.4)   tmp=ws(r+1,c)
!                     if(dir(r,c).eq.8)   tmp=ws(r+1,c-1)
!                     if(dir(r,c).eq.16)  tmp=ws(r,c-1)
!                     if(dir(r,c).eq.32)  tmp=ws(r-1,c-1)
!                     if(dir(r,c).eq.64)  tmp=ws(r-1,c)
!                     if(dir(r,c).eq.128) tmp=ws(r-1,c+1)
!                     if(tmp.ne.ws(r,c)) nextsub=tmp
!                  end if
!               end do
!            end do
!            do r=1,NY
!               do c=1,NX
!                  if(ws(r,c).eq.i .and. nextsub>0) then
!                     ws(r,c)=nextsub
!                     gridscount(nextsub)=gridscount(nextsub)+1
!                  end if
!               end do
!            end do
!          end if
!        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      allocate(nextsub2(nsub,2))
      redo=1
      do while(redo.eq.1)
        nextsub2(:,:)=0
        do i=1,nsub
            do r=1,NY
               do c=1,NX
                  if(ws(r,c).eq.i) then
                    if(dir(r,c).eq.1)   tmp=ws(r,c+1)
                    if(dir(r,c).eq.2)   tmp=ws(r+1,c+1)
                    if(dir(r,c).eq.4)   tmp=ws(r+1,c)
                    if(dir(r,c).eq.8)   tmp=ws(r+1,c-1)
                    if(dir(r,c).eq.16)  tmp=ws(r,c-1)
                    if(dir(r,c).eq.32)  tmp=ws(r-1,c-1)
                    if(dir(r,c).eq.64)  tmp=ws(r-1,c)
                    if(dir(r,c).eq.128) tmp=ws(r-1,c+1)
                    if(tmp.ne.ws(r,c)) nextsub2(i,1)=tmp
                  end if
               end do
            end do
        end do
        do i=1, nsub
          do j=1,nsub
            if(nextsub2(i,1).eq.j) then
              nextsub2(j,2)=i
            end if
          end do
        end do
        eliminatednum=0
        do i=1,nsub
          if(gridscount(i)<smallestWSshed*1000.0/DemRes*1000.0/DemRes) then
            eliminatednum=eliminatednum+1
            if(nextsub2(i,2).gt.0 .and. gridscount(nextsub2(i,2))<gridscount(nextsub2(i,1))) then
              nextsub=nextsub2(i,2)
            else
              if(nextsub2(i,1)>0) then
                nextsub=nextsub2(i,1)
              else
                nextsub=nextsub2(i,2)
              end if
            end if
            !print*,i,nextsub2(i,:),nextsub,gridscount(nextsub2(i,1)),gridscount(nextsub2(i,2))
            do r=1,NY
               do c=1,NX
                  if(ws(r,c).eq.i .and. nextsub>0) then
                     ws(r,c)=nextsub
                     gridscount(nextsub)=gridscount(nextsub)+1
                  end if
               end do
            end do
          end if
        end do
        gridscount(:)=0
        do i=1,NY
          do j=1,NX
            if(ws(i,j).ne.-9999) then
              gridscount(ws(i,j))=gridscount(ws(i,j))+1
            endif
          end do
        end do
        redo = 0
        do i=1,nsub
          if (gridscount(i)>0 .and. gridscount(i)<smallestWSshed*1000.0/DemRes*1000.0/DemRes) then
            redo=1
          end if
        end do
      end do
      if(allocated(nextsub2))deallocate(nextsub2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        allocate(subrename(nsub,2))
        subrename(:,:)=-9999
        tmpi=1
        do i=1,NY
          do j=1,NX
            if(ws(i,j).ne.-9999 .and. subrename(ws(i,j),1).eq.-9999) then
               subrename(ws(i,j),1)=ws(i,j)
               subrename(ws(i,j),2)=tmpi
               tmpi=tmpi+1 
            endif
          end do
        end do
        do i=1,NY
          do j=1,NX
            if(ws(i,j).ne.-9999 .and. subrename(ws(i,j),2).ne.-9999) then
              ws(i,j)=subrename(ws(i,j),2)
            endif
          end do
        end do



        nodata=-9999
        call write_asc_int("ws_fine.asc",NY,NX,X0,Y0,cellsize,nodata,ws)


!!!!!!!!!!!!!!!!!!!!!!!!
        call read_asc_head("net_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
          allocate(IndataINT(NY,NX))
          call read_asc_int("net_fine.asc",NY,NX,IndataINT)
          do i=1,NY
            do j=1,NX
              if(IndataINT(i,j).eq.nodata) then
                IndataINT(i,j)=-9999
              else
                if (IndataINT(i,j)>0)IndataINT(i,j)=ws(i,j)
              end if
            end do
          end do
          nodata=-9999
          call write_asc_int("net_fine.asc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
        if(allocated(IndataINT))deallocate(IndataINT)
        if(allocated(subrename))deallocate(subrename)

        if(allocated(ws))deallocate(ws)
        if(allocated(ws))deallocate(ws2)
        if(allocated(dir))deallocate(dir)
        if(allocated(gridscount)) deallocate(gridscount)

        infile="ws_fine.asc"
        call smallextent(infile,infile,int(ModelRes/DemRes)*3,1)
        call read_asc_head("ws_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
        allocate(ws2(NY,NX))
        do i=1,NY
          do j=1,NX
              ws2(i,j)=1
          end do
        end do
        call write_asc_int("ws_bound.asc",NY,NX,X0,Y0,cellsize,nodata,ws2)

        infile="acc_fine.asc"
        maskfile="ws_fine.asc"
        outfile="acc_fine.asc"
        call ExtByMask(infile,maskfile,outfile,1)

        infile="dir_fine.asc"
        maskfile="ws_fine.asc"
        outfile="dir_fine.asc"
        call ExtByMask(infile,maskfile,outfile,1)

        infile="net_fine.asc"
        maskfile="ws_fine.asc"
        outfile="net_fine.asc"
        call ExtByMask(infile,maskfile,outfile,1)


        !call read_asc_head("ws_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
        !allocate(IndataINT(NY,NX))
        !call read_asc_int("ws_fine.asc",NY,NX,IndataINT)
        call strlen(hydro_para_dir,d1,d2)
        !print*,trim(hydro_para_dir(d1:d2))
        inputfile="ws_fine.asc"
        outputfile=trim(hydro_para_dir(d1:d2))//"/ws.asc"
        call BlockMajor(inputfile,outputfile,int(ModelRes/DemRes),int(ModelRes/DemRes))
        !call write_asc_int(trim(hydro_para_dir(d1:d2))//"/ws.asc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
        !if(allocated(IndataINT))deallocate(IndataINT)


!        sub=''
!        call filepostfix(midfileformat,l1,l2)
!        sub=midfileformat(l1:l2)
!        if (sub(1:2).eq.'nc') then
!          call read_asc_head("ws_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
!          allocate(IndataINT(NY,NX))
!          call read_asc_int("ws_fine.asc",NY,NX,IndataINT)
!          call write_nc_int("ws_fine.nc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
!          if(allocated(IndataINT))deallocate(IndataINT)
!
!          call read_asc_head("acc_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
!          allocate(IndataINT(NY,NX))
!          call read_asc_int("acc_fine.asc",NY,NX,IndataINT)
!          call write_nc_int("acc_fine.nc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
!          if(allocated(IndataINT))deallocate(IndataINT)
!
!          call read_asc_head("dir_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
!          allocate(IndataINT(NY,NX))
!          call read_asc_int("dir_fine.asc",NY,NX,IndataINT)
!          call write_nc_int("dir_fine.nc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
!          if(allocated(IndataINT))deallocate(IndataINT)
!
!          call read_asc_head("net_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
!          allocate(IndataINT(NY,NX))
!          call read_asc_int("net_fine.asc",NY,NX,IndataINT)
!          call write_nc_int("net_fine.nc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
!          if(allocated(IndataINT))deallocate(IndataINT)
!        end if
        return
end subroutine
  
  
    subroutine priver_horton()
      use gisutil
      implicit none
      character(5)::sub
      integer,allocatable::acc(:,:),dir(:,:),net(:,:)
      integer,allocatable::subbasin(:,:),psubbasin(:,:)
      integer,allocatable::nbasin(:),basinlevel(:)
      integer(4),allocatable::pbasinup(:,:)
      integer,allocatable::netbasinstart(:,:),netbasinend(:,:),nextrowcol(:,:)
      integer,allocatable::netbasinstart2(:,:),netbasinend2(:,:)
      integer,allocatable::basinpnext(:),nbasinup(:),basinup(:,:)
      integer(4),allocatable::pbasin(:)
      integer::nr,nc,ia1,ia2,ib1,ib2,outlet,ny,nx,r,c
      real(8)::x0,y0
      real::nodata,cellsize
      integer::s,l1,l2,i,j,tmpi,nsub,k,tmpj,maxtmp,d1,d2
      integer,allocatable::tmpjj(:),tmpii(:),drainage(:)
      integer::level,p
      character*5:: ncols,nrows
      character*12:: xllcorner,yllcorner
      character*8:: scellsize
      character*12:: snodata
      character *12:: nstr
      character(200)::filename,infile,outputfile,wsfinefile,wscoasefile
      logical::file_e,file_e2

      integer(4)::pbasin_tmp(maxlink)
      integer:: upbasin(maxlink)
      integer:: numlev(maxlevel)

      real,allocatable::distance(:,:),flowlen(:,:)

        
      ncols='ncols'
      nrows='nrows'
      xllcorner='xllcorner'
      yllcorner='yllcorner'
      scellsize='cellsize'
      snodata='NODATA_value'            
!      read input file dir and acc file from ArcGIS
      !sub=''
      !call filepostfix(midfileformat,l1,l2)
      !sub=midfileformat(l1:l2)
      !if(sub(1:3).eq.'asc') then
         call read_asc_head('ws_fine.asc',nr,nc,x0,y0,cellsize,nodata)
         ny=nr
         nx=nc
         allocate(subbasin(nr,nc))
         call read_asc_int('ws_fine.asc',nr,nc,subbasin)
!      elseif(sub(1:2).eq.'nc') then
!         call read_nc_head('ws_fine.nc',nr,nc,x0,y0,cellsize,nodata)
!         ny=nr
!         nx=nc
!         allocate(subbasin(nr,nc))
!         call read_nc_int('ws_fine.nc',nr,nc,subbasin)
!      end if


      !if(sub(1:3).eq.'asc') then
         call read_asc_head('acc_fine.asc',nr,nc,x0,y0,cellsize,nodata)
         if(nr.ne.ny .and. nx.ne.nc) then
            print*,"error: the numbers of rows and cols does not match in inputs"
            stop
         end if
         allocate(acc(nr,nc))
         call read_asc_int('acc_fine.asc',nr,nc,acc)
         ny=nr
         nx=nc
      !elseif(sub(1:2).eq.'nc') then
!         call read_nc_head('acc_fine.nc',nr,nc,x0,y0,cellsize,nodata)
!         if(nr.ne.ny .and. nx.ne.nc) then
!            print*,"error: the numbers of rows and cols does not match in inputs"
!            stop
!         end if
!         allocate(acc(nr,nc))
!         call read_nc_int('acc_fine.nc',nr,nc,acc)
!      end if

      !if(sub(1:3).eq.'asc') then
         call read_asc_head('net_fine.asc',nr,nc,x0,y0,cellsize,nodata)
         if(nr.ne.ny .and. nx.ne.nc) then
            print*,"error: the numbers of rows and cols does not match in inputs"
            stop
         end if
         allocate(net(nr,nc))
         call read_asc_int('net_fine.asc',nr,nc,net)
!      elseif(sub(1:2).eq.'nc') then
!         call read_nc_head('net_fine.nc',nr,nc,x0,y0,cellsize,nodata)
!         if(nr.ne.ny .and. nx.ne.nc) then
!            print*,"error: the numbers of rows and cols does not match in inputs"
!            stop
!         end if
!         allocate(net(nr,nc))
!         call read_nc_int('net_fine.nc',nr,nc,net)
!      end if

      !if(sub(1:3).eq.'asc') then
         call read_asc_head('dir_fine.asc',nr,nc,x0,y0,cellsize,nodata)
         if(nr.ne.ny .and. nx.ne.nc) then
            print*,"error: the numbers of rows and cols does not match in inputs"
            stop
         end if
         allocate(dir(nr,nc))
         call read_asc_int('dir_fine.asc',nr,nc,dir)
!      elseif(sub(1:2).eq.'nc') then
!         call read_nc_head('dir_fine.nc',nr,nc,x0,y0,cellsize,nodata)
!         if(nr.ne.ny .and. nx.ne.nc) then
!            print*,"error: the numbers of rows and cols does not match in inputs"
!            stop
!         end if
!         allocate(dir(nr,nc))
!         call read_nc_int('dir_fine.nc',nr,nc,dir)
!      end if
!      rename the net file with subbasin number
      do i =1,nr
        do j =1,nc
          if(net(i,j) > 0) net(i,j)=subbasin(i,j)
        enddo
      enddo
      !if(sub(1:3).eq.'asc') then
        call write_asc_int('net_fine.asc',nr,nc,x0,y0,cellsize,nodata,net)
!      elseif(sub(1:2).eq.'nc') then
!        call write_nc_int('net_fine.nc',nr,nc,x0,y0,cellsize,nodata,net)
!      end if



!!!!!!!calculate the number of subbasin!!!!!!!!!!!!!!
      tmpi=0
      nsub=0
      do i=1,nr
        do j = 1,nc
          if( subbasin(i,j) > tmpi) then
            nsub=subbasin(i,j)
            tmpi=subbasin(i,j)              
          endif
        enddo
      enddo
      NumberOfSubWatersheds=nsub
      print *, 'number of subbasin', nsub


      allocate(nbasin(nsub))
      allocate(netbasinstart(nsub,2))
      allocate(netbasinend(nsub,2))
      allocate(nextrowcol(nsub,2))
      allocate(tmpii(nsub))
      allocate(tmpjj(nsub))
      allocate(drainage(nsub))

      do k =1,nsub
        nbasin(k)=k
      enddo
      tmpjj(:)=-1
      tmpii(:)=nr*nc
!!!!!!! determine the start and end point of network for each subbasin
      do i = 1,nr
        do j = 1,nc
          if(net(i,j)>0) then
            if(acc(i,j) > tmpjj(net(i,j))) then
              netbasinend(net(i,j),1)=i
              netbasinend(net(i,j),2)=j          
              tmpjj(net(i,j))=acc(i,j)
            endif
            if(acc(i,j) < tmpii(net(i,j))) then
              netbasinstart(net(i,j),1)=i
              netbasinstart(net(i,j),2)=j
              tmpii(net(i,j))=acc(i,j)
            endif
          end if
        enddo
      enddo 

!!!!!! Determining the outlet subbasin
      print*,"Determining the outlet subbasin"
      maxtmp=0
      do k =1,nsub
        i=netbasinend(k,1)
        j=netbasinend(k,2)
        drainage(k)=acc(i,j)
        if(acc(i,j)> maxtmp) then
          maxtmp=acc(i,j)
          outlet=k
        endif
      enddo

!cccccccccccc  detemine the link of each subbasin ccccccccccccccc
      print*,"detemine the link of each subbasin"
      allocate(basinpnext(nsub))
      allocate(nbasinup(nsub))
      allocate(basinlevel(nsub))
      allocate(pbasin(nsub))
      allocate(psubbasin(nr,nc))
      allocate(pbasinup(nsub,maxlink))
      allocate(basinup(nsub,maxlink))
      nbasinup(:)=0
      basinup(:,:)=-9999
      basinlevel(:)=-9999
      pbasinup(:,:)=-9999
      nextrowcol(:,:)=0
      call link(nr,nc,nsub,outlet,dir,acc,net,nbasin,netbasinend,basinpnext,nextrowcol)

      do k=1,nsub
        do i =1,nsub
          if(basinpnext(i)==k) then
            nbasinup(k)= nbasinup(k)+1
            basinup(k,nbasinup(k))=i
          endif
        enddo
      enddo


      print*,"Max number of up basin",maxval(nbasinup(1:nsub))
      do k = 1, nsub
        if(nbasinup(k) == 0) then
          basinlevel(k)=1
        endif
      enddo

      call strlen(hydro_para_dir,d1,d2)
      !open(10,file=trim(hydro_para_dir(d1:d2))//'/basinup.txt',status='replace')
      do k =1,nsub
        upbasin(:)=0
        if(nbasinup(k) > 0) then
          do i =1,nbasinup(k)
            upbasin(i)=basinup(k,i)
          enddo
        endif
        !write(10,*)nbasin(k),nbasinup(k),upbasin(:)
      enddo
      !close(10)


!cccccc determine the level ccccccc
      print*,"recursively determine level"
      call leveldetermine2(nsub,outlet,level,nbasinup,basinup,basinpnext,basinlevel,nbasin,drainage,pbasin,numlev,maxlink)

      pbasinup(:,:)=-9999
      do k =1,nsub
        if(nbasinup(k) > 0) then
          do i =1,nbasinup(k)
            j=basinup(k,i)
            pbasinup(k,i)=pbasin(j)
          enddo
        endif
      enddo

      !open(18,file=trim(hydro_para_dir(d1:d2))//'/basinlevel.txt',status='unknown')
      !do k =1,nsub
      !  write(18,*)nbasin(k),basinlevel(k),pbasin(k)
      !enddo
      !close(18)
      psubbasin(:,:)=-9999
      do i =1,nr
        do j = 1,nc
          do k =1,nsub
            if(subbasin(i,j)==nbasin(k)) then
               psubbasin(i,j)=pbasin(k)
            endif
          enddo
        enddo
      enddo

      call strlen(hydro_para_dir,d1,d2)
      open(18, file='./pbasin.asc', status='unknown')
      write(18,'(A, i16)') ncols,nc
      write(18,'(A, i16)') nrows,nr
      write(18,'(A, f24.8)') xllcorner,x0
      write(18,'(A, f24.8)') yllcorner,y0
      write(18,'(A, f24.8)') scellsize, cellsize
      write(18,'(A, f24.8)') snodata,nodata
      do i=1,nr        
        write(18,*) (psubbasin(i,j), j=1,nc)
      end do           
      close (18) 

      infile='./pbasin.asc'
      outputfile=trim(hydro_para_dir(d1:d2))//'/pbasin.asc'
      call BlockMajor(infile,outputfile,int(ModelRes/DemRes),int(ModelRes/DemRes))


      j=0
      filename=trim(hydro_para_dir(d1:d2))//'/subcatchment.txt'
      open(10,file=trim(filename),status='unknown')
      write  (10,*) 'id wsname  ws  width_min  width_max  height_min  height_max  roughness_min  roughness_max'
      do k =1,level
        do i =1,numlev(k)
          j=j+1
          nstr='ws'
          write  (10,'(i4,3x,a2,i5,3x,i5.5,2x,6f8.3)')j,nstr,k*10000+i,k*10000+i,50.0,50.0,5.0,5.0,0.01,0.01
        enddo
      enddo
      close(10)
      j=0
      open(10,file=trim(hydro_para_dir(d1:d2))//'/subbasin.dat',status='unknown')
      do k =1,level
        do i =1,numlev(k)
          j=j+1
          do p = 1,nsub
            if(pbasin(p)== k*10000+i) then
              if(nbasinup(p) .gt. maxlink) then
                print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                print*,"The max upbasin is greater than ", maxlink
                print*,"The program stopped here"
                print*, "In basin:",pbasin(p),",the number of up basin is:",nbasinup(p)
                print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                stop
              end if
              pbasin_tmp(:)=0
              if (nbasinup(p)>0) then
                pbasin_tmp(1:nbasinup(p))=pbasinup(p,1:nbasinup(p))
              end if
              !write  (10,'(i6,2x,i6,2x,15i7)')j,k*10000+i,pbasin_tmp(:)!,nextrowcol(p,1),nextrowcol(p,2),&
              write  (10,*)j,k*10000+i,pbasin_tmp(:)!,nextrowcol(p,1),nextrowcol(p,2),&
            endif
          enddo
        enddo ! 
      enddo
      close(10)

      !open(10,file=trim(hydro_para_dir(d1:d2))//'/drainage.txt',status='replace')
      do k =1,level
        do i =1,numlev(k)
          j=j+1
          do p = 1,nsub
            if(pbasin(p)== k*10000+i) then
              !write(10,*)pbasin(p),drainage(p)
            endif
          enddo
        enddo ! 
      enddo

      call strlen(hydro_para_dir,d1,d2)
      inquire(file=trim(hydro_para_dir(d1:d2))//'/distance.asc', exist=file_e)
      if(.not. file_e) then
        print*,'generating distance.asc'
        call read_asc_head("dir_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
        allocate(distance(NY,NX))
        allocate(flowlen(NY,NX))
        print*,"calculate flow length"
        call flowlength(dir,flowlen,NX,NY,cellsize,nodata)
        !call write_nc_real('flowlen.nc',NY,NX,X0,Y0,cellsize,nodata,flowlen)
        !infile='flowlen.nc'
        !maskfile='ws_fine.asc'
        !outfile=trim(hydro_para_dir(d1:d2))//'/distance.asc'
        !call ExtByMask(infile,maskfile,outfile,2)
        print*,"calculate flow distance"
        call write_asc_real('flowlen.asc',NY,NX,X0,Y0,cellsize,nodata,flowlen)

        distance(:,:)=-9999
        do r=1,NY
          do c=1,NX
            if(psubbasin(r,c).ne.-9999) then
              do i=1,nsub
                if(pbasin(i) .eq. psubbasin(r,c)) then
                  distance(r,c)=flowlen(r,c)-flowlen(netbasinend(i,1),netbasinend(i,2))
                  !print*,pbasin(i),psubbasin(netbasinend(i,1),netbasinend(i,2)),flowlen(netbasinend(i,1),netbasinend(i,2))
                end if
              end do
            end if
          end do
        end do
        call write_asc_real('distance.asc',NY,NX,X0,Y0,cellsize,nodata,distance)
        infile='./distance.asc'
        outputfile=trim(hydro_para_dir(d1:d2))//'/distance.asc'
        wsfinefile="./pbasin.asc"
        wscoasefile=trim(hydro_para_dir(d1:d2))//'/pbasin.asc'
        call Blockmean2(infile,outputfile,int(ModelRes/DemRes),int(ModelRes/DemRes),0,wsfinefile,wscoasefile)
        print*,"getting outflow interval"
        open(10,file='endinglength.txt',status='unknown')
        write(10,*) "basinNum  EndingBasin  EndingRow  EndingCol  NextRow  NextCol  NextDistance"
        do i = 1, nsub
         if(nextrowcol(i,1).gt.0 .and. nextrowcol(i,1).lt.NY .and. nextrowcol(i,2).gt.0 .and. nextrowcol(i,2).lt.NX) then
         write(10,*) psubbasin(netbasinend(i,1),netbasinend(i,2)),psubbasin(nextrowcol(i,1),nextrowcol(i,2)),&
                     netbasinend(i,1),netbasinend(i,2),nextrowcol(i,1),nextrowcol(i,2),distance(nextrowcol(i,1),nextrowcol(i,2))
         else
         write(10,*) psubbasin(netbasinend(i,1),netbasinend(i,2)),-9999,-9999,&
                     netbasinend(i,1),netbasinend(i,2),-9999,-9999,-9999
         end if
        end do
        close(10)
      end if



      if(allocated(distance))deallocate(distance)
      if(allocated(flowlen))deallocate(flowlen)
      if(allocated(nbasin))deallocate(nbasin)
      if(allocated(psubbasin))deallocate(psubbasin)
      if(allocated(pbasinup))deallocate(pbasinup)
      if(allocated(subbasin))deallocate(subbasin)
      if(allocated(acc))deallocate(acc)
      if(allocated(net))deallocate(net)
      if(allocated(dir))deallocate(dir)
      if(allocated(nbasin))deallocate(nbasin)
      if(allocated(netbasinstart))deallocate(netbasinstart)
      if(allocated(netbasinend))deallocate(netbasinend)
      if(allocated(tmpii))deallocate(tmpii)
      if(allocated(tmpjj))deallocate(tmpjj)
      if(allocated(drainage))deallocate(drainage)
      if(allocated(basinpnext))deallocate(basinpnext)
      if(allocated(nbasinup))deallocate(nbasinup)
      if(allocated(basinup))deallocate(basinup)
      if(allocated(basinlevel))deallocate(basinlevel)
      if(allocated(pbasin))deallocate(pbasin)
      if(allocated(pbasinup))deallocate(pbasinup)
      if(allocated(psubbasin))deallocate(psubbasin)
      if(allocated(nextrowcol))deallocate(nextrowcol)
    end subroutine

    subroutine link(nr,nc,nsub,outlet,dir,acc,net,nbasin,netbasinend,basinpnext,nextrowcol)
      integer::nr,nc
      integer::nsub,outlet
      integer,dimension(nr,nc):: acc, dir,net
      integer,dimension(nsub,2)::netbasinend,nextrowcol
      integer,dimension(nsub)::nbasin,basinpnext
      integer::i,j,k,m,n
      integer::nexti,nextj
      basinpnext(:)=-9999
      do k=1,nsub
        if(k.ne. outlet) then
          m=0
          i=netbasinend(k,1)
          j=netbasinend(k,2)
          do while (m == 0)
            if(dir(i,j) == 64) then
              nexti=i-1
              nextj=j
            elseif(dir(i,j) == 128) then
              nexti=i-1
              nextj=j+1
            elseif(dir(i,j) == 32) then
              nexti=i-1
              nextj=j-1
            elseif(dir(i,j) == 1) then
              nexti=i
              nextj=j+1
            elseif(dir(i,j) == 16) then
              nexti=i
              nextj=j-1
            elseif(dir(i,j) == 8) then
              nexti=i+1
              nextj=j-1
            elseif(dir(i,j) == 4) then
              nexti=i+1
              nextj=j
            elseif(dir(i,j) == 2) then
              nexti=i+1
              nextj=j+1
            endif
            i=nexti
            j=nextj
            if(net(nexti,nextj) > 0 .and. net(nexti,nextj) .ne. nbasin(k)) then
              basinpnext(k)=net(nexti,nextj)
              nextrowcol(k,1)=nexti
              nextrowcol(k,2)=nextj
              m=1
          endif    
          enddo
        endif
      enddo
    end subroutine



    subroutine  order(ac,bs,bn)
      implicit none
      integer,parameter :: maxb=500
      integer::ac(maxb)
      integer::bs(maxb)
      integer::bn
      integer::i,j,tmp
      do i = 1,bn-1
        do j = i+1,bn
          if(ac(j) < ac(i)) then
            tmp =bs(i)
            bs(i)=bs(j)
            bs(j)=tmp
            tmp=ac(i)
            ac(i)=ac(j)
            ac(j)=tmp
          endif
        enddo
      enddo
    end subroutine

    subroutine sort(level,n)
      implicit none
      integer::level(n)
      integer::n,i,j,tmp
      do i=1,n-1
           do j=i+1,n
            if(level(j)<level(i)) then
                tmp =level(i)
                level(i)=level(j)
                level(j)=tmp
            end if
         end do
      end do
    end subroutine

    subroutine recursivelevel(nsub,i,level,basinlevel,nbasinup,basinup,nbasin,maxlink)
      implicit none
      integer::i,nsub,maxlink
      integer::level
      integer,dimension(nsub)::nbasinup ! number of upbasin
      integer,dimension(nsub,maxlink)::basinup! up basin number
      integer,dimension(nsub)::basinlevel ! basin level
      integer,dimension(nsub)::nbasin

      integer::uplevel(maxlink)
      integer::level2(maxlink,2),tmplevel
      integer::j,ii,rr

      if (nbasinup(i).eq.0) then
        basinlevel(i)=1
        return
      end if
      uplevel(:)=0
      do j=1,nbasinup(i)
         if(basinlevel(basinup(i,j))>0)then
           uplevel(j)=basinlevel(basinup(i,j))
         else
           call recursivelevel(nsub,basinup(i,j),basinlevel(basinup(i,j)),basinlevel,nbasinup,basinup,nbasin,maxlink)
           uplevel(j)=basinlevel(basinup(i,j))
         end if
      end do
      call sort(uplevel,nbasinup(i))
      level2(:,:)=-9999
      level2(1,1)=uplevel(1)
      level2(1,2)=1
      rr=1
      do ii=2,nbasinup(i)
        if(uplevel(ii).eq.level2(rr,1)) then
          level2(rr,2)=level2(rr,2)+1
        else
          rr=rr+1
          level2(rr,1)=uplevel(ii)
          level2(rr,2)=1
        end if
      end do
      do ii=1,rr
        if(level2(ii,2)>1)level2(ii,1)=level2(ii,1)+1
      end do
      level=maxval(level2(1:rr,1))
    end subroutine recursivelevel


    subroutine leveldetermine2(nsub,outlet,level,nbasinup,basinup,basinpnext,basinlevel,nbasin,drainage,pbasin,numlev,maxlink)
      implicit none
      integer::nsub,outlet,level,maxlink ! total level
      integer,dimension(nsub)::nbasinup ! number of upbasin
      integer,dimension(nsub,maxlink)::basinup! up basin number
      integer,dimension(nsub)::basinpnext ! next basin number
      integer,dimension(nsub)::basinlevel ! basin level
      integer,dimension(nsub)::nbasin,drainage
      integer(4),dimension(nsub)::pbasin
      integer,dimension(maxlevel) :: numlev
      integer::i,j,k,p,q,temp,r
      integer,dimension(:),allocatable::bs,tmp,ac
      integer::bn,bn2
      integer::w,m,bm,rr,nn,ii


      do i =1, nsub
         if(basinlevel(i)==-9999) then
            call recursivelevel(nsub,i,basinlevel(i),basinlevel,nbasinup,basinup,nbasin,maxlink)
         end if
      end do

      print*,"level",level,basinlevel(outlet),outlet
      do i =1,nsub
        print*,"level",i,basinlevel(i)
      end do
      level=basinlevel(outlet)
      if(level.gt.maxlevel) then
        print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print*,"The level of watershed subdivision is greater than maxlevel:",maxlevel
        print*,"The program stopped"
        print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        stop
      end if

      allocate(bs(nsub+10))
      allocate(tmp(nsub+10))
      allocate(ac(nsub+10))

!    determine  pbasin number according to
      do k =1,1
        bn=0
        do i = 1,nsub
          if(basinlevel(i)==k) then
            bn=bn+1
            pbasin(i)=10000*k+bn
          endif
        enddo
        print *, 'basin number of level' ,k ,'  =',bn
        numlev(k)=bn
      enddo

      do k=2,level
        bn=0
        do i = 1,nsub
          if(basinlevel(i)==k) then
            bn=bn+1
            bs(bn)=i
            ac(bn)=drainage(i)
            pbasin(i)=10000*k+bn
          endif
        enddo
        print *, 'basin number of level' ,k ,'  =',bn
        numlev(k)=bn

        do j = 1,bn-1
          do p = j+1,bn
             if ( p .ne. j ) then
               if(ac(j) > ac(p)) then
                 q=bs(j)
                 r=bs(p)
                 temp=pbasin(q)
                 pbasin(q)=pbasin(r)
                 pbasin(r)=temp
                 temp=q
                 bs(j)=r
                 bs(p)=temp
                 temp=ac(j)
                 ac(j)=ac(p)
                 ac(p)=temp
               endif
             endif
          enddo !p
        enddo !j
      enddo !do 2,level

      if(allocated(bs))deallocate(bs)
      if(allocated(ac))deallocate(ac)
      if(allocated(tmp))deallocate(tmp)
    end subroutine



    subroutine deminput()
      use gisutil
      implicit none
      character(5)::sub,strres,sub2
      integer::l1,l2,a1,a2,d1,d2
      logical::file_e,file_e2
      character(200)::ascfile,infile,outfile,maskfile,filldem
      real(8)::X0,Y0
      integer::NY,NX,i,j,r,c
      real::cellsize,nodata
      real,allocatable::DEM(:,:)!DEM with finer resolution
      real,allocatable::slp(:,:),slpsum(:,:) !slope with finer resolution in percentile
      real,allocatable::asp(:,:),aspsum(:,:) !slope with finer resolution in percentile
      real,allocatable::rivlen(:,:),rivlen2(:,:),slplen(:,:)
      integer,allocatable::dir100(:,:),net100(:,:) !flow direction and river net with finer resolution 
      real,allocatable::cell_area(:,:)
      integer,allocatable::zone(:,:)
      real,allocatable::bedslp(:,:),bedslp2(:,:),flowlen(:,:)
      integer,allocatable::ws(:,:),dir_fine(:,:)
      character(200)::resolution,wsfinefile,wscoasefile
      integer,allocatable::IndataInt(:,:)

       call strlen(hydro_para_dir,d1,d2)
       inquire(file=trim(hydro_para_dir(d1:d2))//'/slope.asc', exist=file_e)
       if(.not. file_e) then
          print*,"generating slope.asc"
          !call system("gdal_translate -of AAIGrid Fil_fine.tif  Fil_fine.asc")
          ascfile='FilDEM_fine.asc'
          call strlen(ascfile,l1,l2)
          call read_asc_head(ascfile(l1:l2),NY,NX,X0,Y0,cellsize,nodata)
          if(nodata.ne.-9999) then
            allocate(DEM(NY,NX))
            call read_asc_real(ascfile(l1:l2),NY,NX,DEM)
            do i=1,NY
              do j=1,NX
                if(DEM(i,j).eq.nodata) then
                  DEM(i,j)=-9999
                endif
              end do
            end do
            nodata=-9999
            call write_asc_real(ascfile(l1:l2),NY,NX,X0,Y0,cellsize,nodata,DEM)
            if(allocated(DEM))deallocate(DEM)
          end if
          infile="./FilDEM_fine.asc"
          maskfile="./ws_bound.asc"
          outfile="./FilDEM_fine.asc"
          call ExtByMask(infile,maskfile,outfile,1)
          infile="./FilDEM_fine.asc"
          call strlen(infile,l1,l2)
          call read_asc_head(infile(l1:l2),NY,NX,X0,Y0,cellsize,nodata)
          allocate(DEM(NY,NX))
          call read_asc_real(infile(l1:l2),NY,NX,DEM)
          allocate(slp(NY,NX))
          allocate(asp(NY,NX))
          call slope(dem,slp,NX,NY,cellsize,nodata)
          call aspect(dem,asp,NX,NY,cellsize,nodata)
          if(allocated(DEM)) deallocate(DEM)

          infile='./slope.asc'
          call write_asc_real(infile,NY,NX,X0,Y0,cellsize,nodata,slp)
          if(allocated(slp))deallocate(slp)
          wsfinefile="pbasin.asc"
          wscoasefile=trim(hydro_para_dir(d1:d2))//'/pbasin.asc'
          outfile=trim(hydro_para_dir(d1:d2))//'/slope.asc'
          call Blockmean2(infile,outfile,int(ModelRes/DemRes),int(ModelRes/DemRes),2,wsfinefile,wscoasefile)

          infile='./aspect.asc'
          call write_asc_real(infile,NY,NX,X0,Y0,cellsize,nodata,asp)
          if(allocated(asp))deallocate(asp)
          outfile=trim(hydro_para_dir(d1:d2))//'/aspect.asc'
          call Blockmean2(infile,outfile,int(ModelRes/DemRes),int(ModelRes/DemRes),2,wsfinefile,wscoasefile)
        endif

        inquire(file=trim(hydro_para_dir(d1:d2))//'/slope_length.asc', exist=file_e)
        if(.not. file_e) then
          print*,'generating slope_length.asc'
          inquire(file='stream2__fine.tif', exist=file_e)
          if(file_e) call system("rm stream2__fine1.tif")
          call system("mpiexec -n "//cpus//" threshold -ssa acc2__fine.tif -src stream2__fine1.tif -thresh 10")
          inquire(file='stream2__fine.tif', exist=file_e)
          if(file_e) call system("gdal_translate -of AAIGrid -ot Int16 stream2__fine1.tif  net_fine1.asc")
          call read_asc_head("net_fine1.asc",NY,NX,X0,Y0,cellsize,nodata)
          allocate(IndataINT(NY,NX))
          call read_asc_int("net_fine1.asc",NY,NX,IndataINT)
          do i=1,NY
            do j=1,NX
              if(IndataINT(i,j).eq.nodata) then
                IndataINT(i,j)=-9999
              end if
            end do
          end do
          nodata=-9999
          call write_asc_int("net_fine1.asc",NY,NX,X0,Y0,cellsize,nodata,IndataINT)
          if(allocated(IndataINT))deallocate(IndataINT)
          infile="net_fine1.asc"
          maskfile="ws_fine.asc"
          outfile="net_fine1.asc"
          call ExtByMask(infile,maskfile,outfile,1)
          call read_asc_head("dir_fine.asc",NY,NX,X0,Y0,cellsize,nodata)
          allocate(dir100(NY,NX))
          allocate(rivlen(NY,NX))
          allocate(net100(NY,NX))
          call read_asc_int("dir_fine.asc",NY,NX,dir100)
          call read_asc_int("net_fine1.asc",NY,NX,net100)
          rivlen(:,:)=0.0
          !$omp parallel private(i,j)
          !$omp do
          do i=1,NY
            do j=1,NX
              if(dir100(i,j).eq.nodata) then
                rivlen(i,j)=-9999
              else
                if(dir100(i,j).eq.1) rivlen(i,j)=net100(i,j)
                if(dir100(i,j).eq.128) rivlen(i,j)=1.4142*net100(i,j)
                if(dir100(i,j).eq.64) rivlen(i,j)=net100(i,j)
                if(dir100(i,j).eq.32) rivlen(i,j)=1.4142*net100(i,j)
                if(dir100(i,j).eq.16) rivlen(i,j)=net100(i,j)
                if(dir100(i,j).eq.8) rivlen(i,j)=1.4142*net100(i,j)
                if(dir100(i,j).eq.4) rivlen(i,j)=net100(i,j)
                if(dir100(i,j).eq.2) rivlen(i,j)=1.4142*net100(i,j)
                if(rivlen(i,j)<0.0)then
                   print*,net100(i,j),dir100(i,j)
                endif
              end if
            end do
          end do
          !$omp end do
          !$omp end PARALLEL
          if (allocated(net100)) deallocate(net100)
          if (allocated(dir100)) deallocate(dir100)
          nodata=-9999
          call write_asc_real("rivlen.asc",NY,NX,X0,Y0,cellsize,nodata,rivlen)
          if (allocated(rivlen))deallocate(rivlen)
          infile="rivlen.asc"
          outfile='rivsum.asc'
          wsfinefile="pbasin.asc"
          wscoasefile=trim(hydro_para_dir(d1:d2))//'/pbasin.asc'
          call BlockSum2(infile,outfile,int(ModelRes/DemRes),int(ModelRes/DemRes),2,wsfinefile,wscoasefile)
          if(allocated(rivlen))deallocate(rivlen)
          if(allocated(rivlen2))deallocate(rivlen2)

          call read_asc_head('rivsum.asc',NY,NX,X0,Y0,cellsize,nodata)
          allocate(slplen(NY,NX))
          call read_asc_real('rivsum.asc',NY,NX,slplen)
          !$omp parallel private(r,c)
          !$omp do
          do r=1,NY
            do c=1,NX
              if(slplen(r,c).ne.nodata) then
               if(slplen(r,c)>0)slplen(r,c)=ModelRes*ModelRes/slplen(r,c)/2.0/DemRes
               if(slplen(r,c)>2*ModelRes)slplen(r,c)=2*ModelRes
              end if
            end do
          end do
          !$omp end do
          !$omp end PARALLEL
          call write_asc_real(trim(hydro_para_dir(d1:d2))//'/slope_length.asc',NY,NX,X0,Y0,cellsize,nodata,slplen)
          if(allocated(slplen))deallocate(slplen)
        endif


        inquire(file=trim(hydro_para_dir(d1:d2))//'/elevation.asc', exist=file_e)
        if(.not. file_e) then
          !generate elevation.asc (1k) from ws_fine.asc
          print*,'generating elevation.asc'
          infile='FilDEM_fine.asc'
          outfile=trim(hydro_para_dir(d1:d2))//'/elevation.asc'
          maskfile='ws_bound.asc'
          call ExtByMask(infile,maskfile,outfile,2)
          call blockmean2(infile,outfile,int(ModelRes/DemRes),int(ModelRes/DemRes),2,wsfinefile,wscoasefile)
        end if


! generate zone.asc from gj1k
        inquire(file=trim(hydro_para_dir(d1:d2))//'/zone.asc', exist=file_e)
        if(.not. file_e) then
          print*,'generating zone.asc'
          call read_asc_head(trim(hydro_para_dir(d1:d2))//'/elevation.asc',NY,NX,X0,Y0,cellsize,nodata)
          if(.not. allocated(DEM))allocate(DEM(NY,NX))
          if(.not. allocated(zone))allocate(zone(NY,NX))
          zone(:,:)=-9999
          call read_asc_real('ws_fine.asc',NY,NX,DEM)
         !$omp parallel private(r,c)
         !$omp do
          do r=1,NY
            do c=1,NX
              if(DEM(r,c)>4000.0) zone(r,c)=1
              if(DEM(r,c).le.4000.0) zone(r,c)=2
              if(DEM(r,c).le.3000.0) zone(r,c)=3
              if(DEM(r,c).le.2000.0) zone(r,c)=4
              if(DEM(r,c).le.1000.0) zone(r,c)=5
            end do
         end do
         !$omp end do
         !$omp end PARALLEL
          call write_asc_int(trim(hydro_para_dir(d1:d2))//'/zone.asc',NY,NX,X0,Y0,cellsize,nodata,zone)
          if(allocated(zone)) deallocate(zone)
          if(allocated(DEM)) deallocate(DEM)
        endif

!generate cell_area.asc from ws_fine
        inquire(file=trim(hydro_para_dir(d1:d2))//'/cell_area.asc', exist=file_e)
        if(.not. file_e) then
          print*,'generating cell_area.asc'
          call read_asc_head(trim(hydro_para_dir(d1:d2))//'/pbasin.asc',NY,NX,X0,Y0,cellsize,nodata)
          if(.not. allocated(cell_area))allocate(cell_area(NY,NX))
          call read_asc_real(trim(hydro_para_dir(d1:d2))//'/pbasin.asc',NY,NX,cell_area)
          !$omp parallel private(r,c)
          !$omp do
          do r=1,NY
            do c=1,NX
              if(cell_area(r,c).ne.nodata) cell_area(r,c)=ModelRes*ModelRes
            end do
          end do
          !$omp end do
          !$omp end PARALLEL
          call write_asc_real(trim(hydro_para_dir(d1:d2))//'/cell_area.asc',NY,NX,X0,Y0,cellsize,nodata,cell_area)
          if(allocated(cell_area)) deallocate(cell_area)
        end if

! generate bedslope.asc from wdh_n50
        inquire(file=trim(hydro_para_dir(d1:d2))//'/bedslope.asc', exist=file_e)
        if(.not. file_e) then
          print*,'generating bedslope.asc'
          infile='slope.asc'
          call read_asc_head(infile,NY,NX,X0,Y0,cellsize,nodata)
          allocate(slp(NY,NX))
          allocate(net100(NY,NX))
          allocate(bedslp(NY,NX))
          allocate(bedslp2(NY,NX))
          bedslp(:,:)=-9999
          call read_asc_real(infile,NY,NX,slp)
          call read_asc_int("net_fine1.asc",NY,NX,net100)
         !$omp parallel private(r,c)
         !$omp do
          do r=1,NY
            do c=1,NX
              if(net100(r,c).eq.1) bedslp(r,c)=slp(r,c)
            end do
          end do
         !$omp end do
         !$omp end PARALLEL
         call write_asc_real('bedslp.asc',NY,NX,X0,Y0,cellsize,nodata,bedslp)
         infile='bedslp.asc'
         outfile='bedslpblcmean1.asc'
         wsfinefile="pbasin.asc"
         wscoasefile=trim(hydro_para_dir(d1:d2))//'/pbasin.asc'
         call Blockmean(infile,outfile,int(ModelRes/DemRes),int(ModelRes/DemRes),2)
         infile='bedslp.asc'
         outfile='bedslpblcmean2.asc'
         call Blockmean(infile,outfile,int(ModelRes/DemRes)*10,int(ModelRes/DemRes)*10,2)
         bedslp(:,:)=-9999
         bedslp2(:,:)=-9999
         call read_asc_real('bedslpblcmean1.asc',NY,NX,bedslp)
         call read_asc_real('bedslpblcmean2.asc',NY,NX,bedslp2)
         !$omp parallel private(r,c)
         !$omp do
         do r = 1,NY
           do c = 1,NX
             if(bedslp(r,c).eq.nodata .and. bedslp2(r,c).ne.nodata) bedslp(r,c)=bedslp2(r,c)
           end do
         end do
         !$omp end do
         !$omp end PARALLEL
         call write_asc_real('bedslp_blcmean.asc',NY,NX,X0,Y0,cellsize,nodata,bedslp)
         infile="bedslp_blcmean.asc"
         outfile='bedslp_resample.asc'
         call Resample(infile,outfile,real(ModelRes),1,2)

         infile='bedslp_resample.asc'
         maskfile=trim(hydro_para_dir(d1:d2))//'/pbasin.asc'
         outfile=trim(hydro_para_dir(d1:d2))//'/bedslope.asc'
         call ExtByMask(infile,maskfile,outfile,2)
         if(allocated(bedslp)) deallocate(bedslp)
         if(allocated(bedslp2)) deallocate(bedslp2)
         if(allocated(slp)) deallocate(slp)
         if(allocated(net100)) deallocate(net100)
       end if

       print*,'checking bedslope.asc and slope_length.asc'
       call strlen(hydro_para_dir,d1,d2)
       call read_asc_head(trim(hydro_para_dir(d1:d2))//'/pbasin.asc',NY,NX,X0,Y0,cellsize,nodata)
       allocate(ws(NY,NX))
       call read_asc_int(trim(hydro_para_dir(d1:d2))//'/pbasin.asc',NY,NX,ws)
       allocate(bedslp(NY,NX))
       call read_asc_real(trim(hydro_para_dir(d1:d2))//'/bedslope.asc',NY,NX,bedslp)
       allocate(slplen(NY,NX))
       call read_asc_real(trim(hydro_para_dir(d1:d2))//'/slope_length.asc',NY,NX,slplen)
       !$omp parallel private(r,c)
       !$omp do      
       do r = 1,NY
         do c = 1,NX
           if(ws(r,c).ne.nodata .and. bedslp(r,c).eq.nodata) then
             print*,'there is error in bedslope.asc, nodata is expected'
             stop
           end if
           if(ws(r,c).ne.nodata .and. slplen(r,c).eq.nodata) then
             print*,'there is error in bedslope.asc, nodata is expected'
             stop
           end if
         end do
       end do
       !$omp end do
       !$omp end PARALLEL
       if(allocated(ws)) deallocate(ws)    
       if(allocated(bedslp)) deallocate(bedslp)
       if(allocated(slplen)) deallocate(slplen)


!generate soil_unit.asc from ws_fine

!generate soil_depth.asc from ws_fine

!generate lu.asc from ws_fine
       if (allocated(flowlen)) deallocate(flowlen)
       if (allocated(dir_fine)) deallocate(dir_fine)
    end subroutine


    subroutine morph_river()
      use gisutil
      implicit none
      integer::nsub
      integer,allocatable::pbasin(:,:)
      real,allocatable::distance(:,:),distance2(:)
      real,allocatable::bedslope(:,:)
      integer,allocatable::subbasin(:,:),ending(:,:)
      integer::NY,NX,r,c,d1,d2
      real(8)::X0,Y0
      real::nodata,cellsize
      integer::IOstatus,i,temp1,temp2,temp3,m,j
      character(5)::rivername
      character(80)::outfile
      character(200)::tmp
      logical::paradir
      integer,allocatable::nextrowcol(:,:)


      call strlen(hydro_para_dir,d1,d2)
      call read_asc_head(trim(hydro_para_dir(d1:d2))//'/pbasin.asc',NY,NX,X0,Y0,cellsize,nodata)
      allocate(pbasin(NY,NX))
      call read_asc_int(trim(hydro_para_dir(d1:d2))//'/pbasin.asc',NY,NX,pbasin)


      allocate(distance(NY,NX))
      call read_asc_real(trim(hydro_para_dir(d1:d2))//'/distance.asc',NY,NX,distance)

      allocate(bedslope(NY,NX))
      call read_asc_real(trim(hydro_para_dir(d1:d2))//'/bedslope.asc',NY,NX,bedslope)

!      open(10,file='./endinglength.txt',status='old')
!      nsub=0
!      IOstatus=1
!      do while(IOstatus.ge.0)
!        nsub=nsub+1
!        READ(10,*,IOSTAT=IOstatus) temp1,temp2,temp3
!      end do
!      close(10)


!      if(NumberOfSubWatersheds.ne. nsub) then
!         print*, 'The number of watersheds is wrong in the endinglength.txt file'
!         stop
!      end if



      nsub=NumberOfSubWatersheds
      allocate(ending(nsub,6))
      allocate(distance2(nsub))
      open(10,file='./endinglength.txt',status='old')
      read(10,*) tmp
      do i=1,nsub
        READ(10,*,IOSTAT=IOstatus) (ending(i,j),j=1,6),distance2(i)
      end do
      close(10)

      !distance2(:,:)=-9999
      !do r=1,NY
      !  do c=1,NX
      !    if(pbasin(r,c).ne.-9999) then
      !      do i=1,nsub
      !        if(pbasin(r,c).eq.ending(i,1)) then
      !          distance2(r,c)=distance(r,c)-distance(ending(i,2),ending(i,3))
      !        end if
      !      end do
      !    end if
      !  end do
      !end do

      !open(10,file='./hydropara/subbasin.dat',status='old')
      !allocate(subbasin(nsub,13))
      !do i = 1, nsub
      !   read(10,*) (subbasin(i,j),j=1,13)
      !end do
      !close(10)

      inquire(file='./morph',exist=paradir)
      if(.not. paradir) then
         call system("mkdir ./morph")
      else
         call system("rm ./morph/*")
      end if

      do i=1,nsub
        write(rivername,'(I5)') ending(i,1)
        outfile='./morph/ws'//rivername//'_morph'
        open(3,file=trim(outfile),status='unknown')
        m=0
        do r=1,NY
          do c=1,NX
            if(pbasin(r,c).ne.-9999 .and. (distance(r,c).eq.-9999 .or. bedslope(r,c).eq.-9999)) then
              print *, 'wrong in basin_morph:',pbasin(r,c),r,c,distance(r,c),bedslope(r,c)
              stop
            end if
            if(pbasin(r,c).ne.-9999 .and. pbasin(r,c).eq. ending(i,1) ) then
              write(3,10) distance(r,c),r,c,bedslope(r,c)
            end if
          end do
        enddo
        close(3)
      end do 
10      format(f15.3,2i8,f15.5)
      if(allocated(pbasin))deallocate(pbasin)
      if(allocated(distance))deallocate(distance)
      if(allocated(distance2))deallocate(distance2)
      if(allocated(bedslope))deallocate(bedslope)
      if(allocated(ending))deallocate(ending)
    end subroutine
      
    subroutine river_parameter()
      use gisutil
      implicit none
      integer::nsub,IOstatus,i,ii,ws,j,maxnum,maxnp
      character(7)::wsname
      real::width_min,width_max,height_min,height_max,roughness_min,roughness_max
      integer,allocatable::outinterval(:,:)
      !real,allocatable::outdistance(:,:)
      integer::NY,NX,r,c,d1,d2
      real(8)::X0,Y0
      real::nodata,cellsize

      real,allocatable::distance2(:)
      integer,allocatable::ending(:,:)

      nsub=NumberOfSubWatersheds
      allocate(ending(nsub,8))
      allocate(distance2(nsub))
      open(10,file='./endinglength.txt',status='old')
      read(10,*)
      do i=1,nsub
        READ(10,*,IOSTAT=IOstatus) (ending(i,j),j=1,6),distance2(i)
      end do
      close(10)

      call strlen(hydro_para_dir,d1,d2)
      call read_asc_head(trim(hydro_para_dir(d1:d2))//'/pbasin.asc',NY,NX,X0,Y0,cellsize,nodata)
      subcatchment2file=trim(hydro_para_dir(d1:d2))//"/subcatchment.txt"
      open(10,file=subcatchment2file,status='old')
      read(10,*)
      nsub=0
      IOstatus=1
      do while(IOstatus.ge.0)        
        READ(10,*,IOSTAT=IOstatus) ii,wsname,ws,width_min,width_max,height_min,height_max,roughness_min,roughness_max
        nsub=nsub+1
      end do
      close(10)

      nsub=nsub-1
      if(NumberOfSubWatersheds.ne. nsub) then
         print*,NumberOfSubWatersheds,nsub
         print*, 'The number of watersheds is wrong in the subcatchment2.dat file'
         stop
      end if

      !outinterval=0


      !allocate(outdistance(nsub,4))
      !open(10,file='./outdistance.txt',status='old')
      !do i = 1, nsub
      !   read(10,*) (outdistance(i,j),j=1,4)
      !end do
      !close(10)
        
      open(11,file=subcatchment2file,status='old')
      read(11,*)
      print *,'Calculating sub-basin parameters...'
      maxnum=0
      maxnp=0
      do i=1,nsub
        !outinterval(i,1)=int(outdistance(i,1))
        !outinterval(i,2)=int(outdistance(i,2))
        !outinterval(i,3)=int(outdistance(i,4))
        READ(11,*) ii,wsname,ws,width_min,width_max,height_min,height_max,roughness_min,roughness_max
        call parameter_model(nsub,ii,wsname,ws,width_min,width_max,height_min,height_max,roughness_min,roughness_max,&
                             ending,maxnum,maxnp,NX,NY,distance2)
      end do
      close(11)

      allocate(outinterval(nsub,maxlink+2))
      open(10,file=trim(hydro_para_dir(d1:d2))//"/subbasin.dat",status='old')
      do i=1,nsub
        READ(10,*) (outinterval(i,j),j=1,maxlink+2)
      end do
      close(10)

      do i=1,nsub
         do j=3,maxlink+2
           if(outinterval(i,j)>0) then
             do ii = 1,nsub
               if(outinterval(i,j).eq.ending(ii,1)) outinterval(i,j)=ending(ii,7)
             end do
           end if
         end do
      end do

      open(11,file="./ending.txt",status='unknown')
      do i=1,nsub
        write(11,*) (ending(i,j),j=1,8)
      end do
      close(11)


      open(11,file=trim(hydro_para_dir(d1:d2))//"/outinterval.dat",status='unknown')
      do i=1,nsub
        write(11,'(17i6)') (outinterval(i,j),j=1,maxlink+2)
      end do
      close(11)

      call strlen(hydro_para_dir,d1,d2)
      open(3,file=trim(hydro_para_dir(d1:d2))//'/riverpara/maxnp',status='unknown')
      write(3,*) maxnum,maxnp,maxlink,maxlevel
      close(3)
      
      !open(10,file=trim(hydro_para_dir(d1:d2))//'/outinterval.txt',status='replace')
      !do i = 1, nsub
      !   write(10,*) (outinterval(i,j),j=1,5)
      !end do
      !close(10)
      if(allocated(distance2))deallocate(distance2)
      if(allocated(ending))deallocate(ending)
      if(allocated(outinterval))deallocate(outinterval)
      !if(allocated(outdistance))deallocate(outdistance)
    end subroutine

    subroutine parameter_model(nsub,jj,wsname,ws,width_min,width_max,height_min,height_max,roughness_min,roughness_max,&
                            ending,maxnum,maxnp,NX,NY,distance2)
      use gisutil
      implicit none
      integer,intent(in)::jj,nsub,ws,NX,NY
      character(7),intent(in)::wsname
      real,intent(in)::width_min,width_max,height_min,height_max,roughness_min,roughness_max
      integer,intent(inout)::ending(nsub,8)
      real,intent(in)::distance2(nsub)
      integer,intent(inout)::maxnum,maxnp
      !real,intent(inout)::outdistance(nsub,4)

      integer::l1,l2,n_totalgrid,IOstatus,tmpr,tmpc,ll1,ll2
      character(200)::basinmorph,output
      real::tmpdis,tmpbedslp,dis_max
      real,allocatable::dis(:),bed_angle(:)
      integer,allocatable::ir(:),ic(:)
      real,allocatable::dx(:),s0(:),b(:),DR(:),roughness(:)
      integer,allocatable::np(:),ipr(:,:),ipc(:,:)
      integer::i,num,j,ii
      real(8)::x0,y0
      real::cellsize,nodata,d1,d2,tmp_angle,distance,dtmp
      

      print *,nsub,jj, '   ',wsname
      !print *,width_min,width_max,height_min,height_max,roughness_min,roughness_max
      !print *        

      basinmorph="./morph/"//wsname//"_morph"
      call strlen(basinmorph,l1,l2)
      !print *, trim(basinmorph(l1:l2))
      call strlen(hydro_para_dir,ll1,ll2)
      output=hydro_para_dir(ll1:ll2)//"/riverpara/"//wsname//'_river'
      dis_max=0.0
      open(10, file=trim(basinmorph(l1:l2)),status='old')       
      n_totalgrid=0
      IOstatus=1
      do while(IOstatus.ge.0)        
        n_totalgrid=n_totalgrid+1
        READ(10,*,IOSTAT=IOstatus) tmpdis,tmpr,tmpc,tmpbedslp
      end do
      close(10)
      n_totalgrid=n_totalgrid-1
      allocate(dis(n_totalgrid))
      allocate(bed_angle(n_totalgrid))
      allocate(ir(n_totalgrid))
      allocate(ic(n_totalgrid))
        
      open(10, file=trim(basinmorph(l1:l2)),status='old')
      do i=1,n_totalgrid
        READ(10,*) dis(i),ir(i),ic(i),bed_angle(i)
        if(ir(i).le.0 .or. ir(i).gt.NY .or. ic(i).le.0 .or. ic(i).gt.NX) then
          print *,'wrong grid number:',i,ir(i),ic(i) 
          stop
        endif
        if(dis(i).gt.dis_max) dis_max=dis(i)
      end do        
      close(10)

      num=int(dis_max/dx_max)+1
      d1=0.0
      d2=0.0
      if(num>maxnum) maxnum=num

      do i = 1, nsub
         if(ending(i,2).eq.ws) then
           if(distance2(i) .ge. 0) ending(i,7)=num-int(distance2(i)/dx_max)
           ending(i,8)=num
         end if
      end do
        
      allocate(dx(num))
      allocate(np(num))
      allocate(b(num))
      allocate(S0(num))
      allocate(Dr(num))
      allocate(roughness(num))
      allocate(ipr(num,n_totalgrid))
      allocate(ipc(num,n_totalgrid))

      dx(:)=0.0
      do i=1,num
        ii=num-i+1
        dx(ii)=dx_max
        d2=d1+dx(ii)
        if(i.eq.num) d2=d1+dx(ii)+1.0
        if(i.eq.num) dx(ii)=dis_max-d1
        np(ii)=0
        tmp_angle=0.0
        do j=1,n_totalgrid
          if(dis(j).ge.d1 .and. dis(j).lt.d2) then
            np(ii)=np(ii)+1
            ipr(ii,np(ii))=ir(j)
            ipc(ii,np(ii))=ic(j)
            if(ir(j).le.0 .or. ic(j).le.0) print *,'wrong:',ii,j,ir(j),ic(j)
            tmp_angle=tmp_angle+bed_angle(j)
          endif
        end do
        !iarea=iarea+np(ii)
        if(np(ii).eq.0) then
           print *,"Zero mode grid in interval:",ii," in basin ",wsname
        endif
        if(dx(ii).eq.0) then
          dx(ii)=500
        endif
        if(np(ii).gt.0)then
          tmp_angle=tmp_angle/float(np(ii))
          s0(ii)=tmp_angle/100.0
        else
          s0(ii)=0.0
        end if
        d1=d2
      end do

      distance=0.0
      do i=1,num
        distance=distance+dx(i)
        if(s0(i).lt.0.0) then
          print *,'Error in calculating river bed slope', i,s0(i)
        end if
      end do
      if(abs(distance-dis_max) .gt. 1.0) then
        print *,'Error in calculating main channel distance IN RIVER',distance,dis_max
      end if
      dtmp=0.0
      do i=1,num
        dtmp=dtmp+dx(i)
        b(i)=width_max+float(num-i)*(width_min-width_max)/float(num)
        roughness(i)=roughness_min+float(num-i)*(roughness_max-roughness_min)/float(num)
        Dr(i)=height_max+float(num-i)*(height_min-height_max)/float(num)
        if(b(i).le.0. .or. roughness(i).le.0. .or. Dr(i).le.0.0) then
          print *,'wrong:  ',wsname
        end if
      end do
        
      call strlen(output,l1,l2)
      open(3,file=output(l1:l2),status='unknown')
      write(3,*) num
      do i=1,num
        if(np(i)>maxnp) maxnp=np(i)
        if(np(i).gt.0)then
          write(3,200) np(i),dx(i),s0(i),b(i),roughness(i),Dr(i)
          write(3,*) (ipr(i,j),ipc(i,j), j=1,np(i))
        elseif(np(i).eq.0)then
          if(i>1 .and. s0(i-1).ne.0.0)then
            s0(i)=s0(i-1)
          end if
          write(3,200) np(i),2000.0,s0(i),b(i),roughness(i),Dr(i)
          write(3,*) -1,-1
        end if
      end do
      close(3)
200     format(i5,f12.2,f12.6,f10.2,f10.5,f10.2)

      deallocate(dx,np,b,s0,dr,roughness,ipr,ipc)
      deallocate(dis,bed_angle,ir,ic)
    end subroutine
      
    subroutine landradio()
      use gisutil
      implicit none
      integer::NY,NX,l,i,j,d1,d2
      real(8)::x0,y0
      real:: nodata,cellsize
      integer,allocatable::lu(:,:)
      real,allocatable::alllu(:,:),alllu2(:,:),lutmp(:,:),lutmp2(:,:)
      real,allocatable::ratio(:,:)
      character(1)::slu
      character(200)::infile,outfile,maskfile        
        
        
      call read_asc_head(lufine,NY,NX,x0,y0,cellsize,nodata)        
      allocate(lu(NY,NX))
      call read_asc_int(lufine,NY,NX,lu)
      allocate(alllu(NY,NX))
      allocate(alllu2(NY,NX))

      alllu(:,:)=-9999
      alllu2(:,:)=-9999
      do i=1,NY
        do j=1,NX
          if(lu(i,j).ne.nodata) alllu(i,j)=1
        end do
      end do
      call Blocksum_memory(NY,NX,int(ModelRes/cellsize),int(ModelRes/cellsize),alllu,alllu2)
      do l=1,9
        allocate(lutmp(Ny,NX))
        allocate(lutmp2(Ny,NX))
        lutmp=-9999
        lutmp2=-9999
        do i=1,NY
          do j=1,NX
            if(lu(i,j).eq.l) lutmp(i,j)=lu(i,j)
          end do
        end do
        call Blocksum_memory(NY,NX,int(ModelRes/cellsize),int(ModelRes/cellsize),lutmp,lutmp2)
        allocate(ratio(NY,NX))
        ratio(:,:)=-9999
        do i=1,NY
          do j=1,NX
            if(lutmp2(i,j).ne.nodata .and. alllu2(i,j).ne.-9999) ratio(i,j)=lutmp2(i,j)/alllu2(i,j)
          end do
        end do
        call write_asc_real('ratiotmp.asc',NY,NX,X0,Y0,cellsize,nodata,ratio)
        deallocate(ratio,lutmp2,lutmp)
        write(slu,'(I1)') l
        infile ='ratiotmp.asc'
        outfile='ratiotmprsmp.asc'
        call Resample(infile,outfile,real(ModelRes),1,2)
        infile ='ratiotmprsmp.asc'
        maskfile='./ws_fine.asc'
        call strlen(hydro_para_dir,d1,d2)
        outfile=trim(hydro_para_dir(d1:d2))//'/lu_'//slu//'.asc'
        call ExtByMask(infile,maskfile,outfile,2)
      end do
      deallocate(lu,alllu,alllu2)
    end subroutine

  end module mod_preprocess
