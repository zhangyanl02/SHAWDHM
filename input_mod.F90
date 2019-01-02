module input_mod
  use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
  use dims_mod,only:nx,ny
  use shaw27_mod
  use statevar_mod,only:nc2d,nsp2d,nr2d,ns2d,inbasin2d,landuse2d,gridarea2d,slopelen2d,skyview2d,dchar2d,&
      tccrit2d,rstom02d,rstexp2d,pleaf02d,rleaf02d,rroot02d,canalb2d,xangle2d,clumpng2d,itype2d,wcandt2d,&
      gmcdt2d,rload2d,zrthik2d,cover2d,albres2d,rescof2d,gmcmax,zs2d,tsdt2d,vlcdt2d,vicdt2d,matdt2d,icesdt2d,&
      albdry2d,albexp2d,zmsrf2d,zhsrf2d,zersrf2d,zmsp2d,zhsp2d,alatud2d,lontitud2d,elevation2d,slope2d,aspect2d,&
      tsavg2d,hrnoon2d,maskflag,nplant2d
      implicit none


!     namelist timepara: in controlpara_mod
!     namelist control
      real(r8)::zmspcm,pondmxcm,zmcm
      real(r8)::albdry,albexp
!     namelist layers
      integer(i4)::ns,nc,nr,nsp,nsoiltype,soilfromtable
!     namelist lvlout
      integer(i4)::lvlout1,lvlout2,lvlout3,lvlout4,lvlout5,lvlout6,lvlout7,lvlout8,&
        lvlout9,lvlout10,lvlout11,lvlout12,lvlout13
!     namelist mapfile
      character(200)::watershed_map,gridarea_map,elevation_map,aspect_map,slopelength_map,&
        slopeangle_map,aquiferdepth_map,met_alt_map,landuse_map,soil_map,met_map,lat_map,lon_map,&
        soilpro_nc,soiltype_map,soilpara_table,tsavg_map,skyview_map
!       namelist initialstate
      character*(200)::soil_ini,snow_ini,veg_ini,res_ini
      integer(i4),dimension(nx,ny)::soil_code
!     namelist vegparameter
      character*(200) vegpara_table,lai_table,hight_table,rootdp_table,weight_table,res_table
!     namelist residue
      real(r8)::rescof,gmcdt1
      


      integer(i4),dimension(nx,ny)::landuse,ws
      real(r8),dimension(16,366)::planth,rtdp,plantw,lai_t
      real(r8),dimension(nx,ny,366)::lai
      real(r8),dimension(16,10)::vegtable
      real(r8),dimension(16,4)::restable

      integer(i4)::ifirst,nstep,ltmpdy,ltmphr,ltmpyr,lvlcdy,lvlchr,lvlcyr
      real(r8)::tmplst,tmp,vlclst,vlc
      real(r8),dimension(nsmax)::tmp1,vlc1,flux
      real(r8),parameter::thresh=0.5
      integer(i4)::lai_from_table
      data ifirst/0/

contains

      subroutine input(setupfile)
        use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4 => shr_kind_r4
        use controlpara_mod
        use soilproperty_mod,only:b2d,entry2d, satk2d,rhob2d,sat2d,sand2d,silt2d,clay2d,om2d,vapexp2d,vapcof2d,&
          thfc2d,thr2d,alpha2d,n2d,l2d
        use constvar_mod
        use netcdf
        implicit none
        character(200),intent(in)::setupfile
        character*(200)::strtemp
!       namelist jobname
        namelist /directories/ hydro_para_dir,output_dir,forcing_dir,lai_dir,shadows_dir
!       namelist timepara
        namelist /timepara/ mtstep,nhrpdt,jstart,hrstar,yrstar,jend,hrend,yrend,hrnoon,timezone
!       namelist control
        namelist /control/ toler,height,pondmxcm,inph2o,mwatrxt,mpltgro,mzcinp,canma,canmb,wcandt1,zmspcm,&
          snotmp,ivlcbc,itmpbc,albdry,albexp,zmcm,restart,shadeef,nsalt,nplant,maskflag,maxstep,hydro_module,deflate,maxiter
          
!       namelist layers
        namelist /layers/ns,nc,nr,nsp,nsoiltype,soilfromtable,lai_from_table
!       namelist mapfile
        namelist /mapfile/ watershed_map,gridarea_map,elevation_map,aspect_map, slopelength_map,&
          slopeangle_map,landuse_map,lat_map,lon_map,soilpro_nc,soiltype_map,soilpara_table,tsavg_map,&
          skyview_map
!       namelist vegparameter
        namelist /vegparameter/ vegpara_table,lai_table,hight_table,rootdp_table,weight_table
!       namelist lvlout
        namelist /levelout/ lvlout1,lvlout2,lvlout3,lvlout4,lvlout5,lvlout6,lvlout7,lvlout8,&
          lvlout9,lvlout10,lvlout11,lvlout12,lvlout13
!       namelist residue
        namelist /residue/ rescof,gmcdt1,res_table
!       namelist initialstate
        namelist /initialstate/ soil_ini,snow_ini,veg_ini,res_ini

        integer(i4)::controlfile,infil,allocatestatus
        integer(i4)::i,j,k,col,row,flunit,ireturn
        real(r8)::textot,cf,dummy
        real(r8)::hum
        integer(i4)::l1,l2,r,c
        logical::paradir


!       open control file for input model control parameters
        call strlen(setupfile,l1,l2)
        controlfile=10
        open(unit=controlfile,file=trim(setupfile(l1:l2)),status='old',err=3,iostat=infil)
3       if (infil.ne.0) then
!         cannot open file for list of input/output files
          print *, infil
          print*,setupfile
          stop
        end if

!       read namelist directories
        print*,"reading namelist directories"
        read(controlfile,directories,err=101,iostat=infil)
101     if(infil.ne.0) then
          print *,'error in reading namelist directories. program stopped.'
          stop
        end if
        
        simulation_dir=output_dir//"RIVER_FLOW"
        call strlen(simulation_dir,l1,l2) 
        inquire(file=trim(simulation_dir(l1:l2)),exist=paradir)
        if(.not. paradir) then
          call system("mkdir "//trim(simulation_dir(l1:l2)))
        end if
        simulation_dir=simulation_dir//"/"

        result1_dir=output_dir//"spatial"
        call strlen(result1_dir,l1,l2) 
        inquire(file=trim(result1_dir(l1:l2)),exist=paradir)
        if(.not. paradir) then
          call system("mkdir "//trim(result1_dir(l1:l2)))
        end if
        result1_dir=result1_dir//"/"
      
        result2_dir=output_dir//"savestate"
        call strlen(result2_dir,l1,l2) 
        inquire(file=trim(result2_dir(l1:l2)),exist=paradir)
        if(.not. paradir) then
          call system("mkdir "//trim(result2_dir(l1:l2)))
        end if
        result2_dir=result2_dir//"/"
      
!       read namelist timepara
        print*,"reading namelist timepara"
        read(controlfile,timepara,err=102,iostat=infil)
102     if(infil.ne.0) then
          print *,'error in reading namelist timepara. program stopped.'
          stop
        end if
        if (hrstar .eq. 0) then
          hrstar=24
          jstart=jstart-1
        end if

        if (mod(24,nhrpdt) .ne. 0) then
          write (6,*)
          write (6,*) ' ************************************************'
          write (6,*) ' parameter for the number of hours per time step'
          write (6,*) ' (nhrpdt; line d of input file) is not valid. '
          write (6,*) ' nhrpdt must be evenly divisible into 24 hours.'
          write (6,*) ' ************************************************'
          stop
        end if
!       check if hrstar is compatilbe with nhrpdt
        if (mod(hrstar,nhrpdt) .ne. 0) then
          write (6,*)
          write (6,*)' **************************************************'
          write (6,*)' the parameter for the number of hours per time'
          write (6,*)' step (nhrpdt; line d of input file) is not'
          write (6,*)' compatible with the beginning hour for the'
          write (6,*)' simulation (hrstar; line b of input file). '
          write (6,*)' hrstar must be zero or evenly divisible by nhrpdt.'
          write (6,*)' **************************************************'
          stop
        end if

!       read namelist control
        print*,"reading namelist control"
        read(controlfile,control,err=103,iostat=infil)
103     if(infil.ne.0) then
          print *,'error in reading namelist control. program stopped.'
          stop
        end if
        pondmx=pondmxcm/100.0
        zmsp2d=zmspcm/100.
        zhsp2d=0.2*zmsp2d
        albdry2d=albdry
        albexp2d=albexp
        zmsrf2d=zmcm/100.
        zhsrf2d=0.2*zmsrf2d
        zersrf2d=0.01

!       read namelist layers
        print*,"reading namelist layers"
        read(controlfile,layers, err=104,iostat=infil)
104     if(infil.ne.0) then
          print *,'error in reading namelist layers. program stopped.'
          stop
        end if
        ns2d=ns
        nc2d=nc
        nr2d=nr
        nsp2d=nsp
        
!       read namelist mapfile
        print*,"reading namelist mapfile"
        read(controlfile,mapfile,err=105,iostat=infil)
105     if(infil.ne.0) then
          print *,'error in reading namelist mapfile. program stopped.'
          stop
        end if


!       read mapfiles
        print*,"reading mapfile ",watershed_map
        call readmapi(nx,ny,ws,trim(watershed_map),ireturn)
        if(ireturn.eq.0)then 
          print *, watershed_map//" doese not exist"
          stop
        endif
        do r = 1,ny
          do c=1,nx
            if (ws(c,r)>0)then
              inbasin2d(c,r)=1
            else
              inbasin2d(c,r)=-9999
            end if
          enddo
        enddo

        print*,"reading mapfile ",gridarea_map
        call readmapr(nx,ny,gridarea2d,gridarea_map,ireturn)
        if(ireturn.eq.0) then  ! no data is available for gridarea
          print *, gridarea_map//" does not exist, default value is used"
          gridarea2d=ddx*ddy
        endif

        print*,"reading mapfile ",elevation_map
        call readmapr(nx,ny,elevation2d,elevation_map,ireturn)
        if(ireturn.eq.0)then 
          print *, elevation_map//" does not exist"
          stop
        endif

        print*,"reading mapfile ",aspect_map
        call readmapr(nx,ny,aspect2d,aspect_map,ireturn)
        if(ireturn.eq.0)then 
          print *, aspect_map//" does not exist"
          stop
        else
          do r=1,ny
            do c=1,nx
              if(aspect2d(c,r)<0.0) then
                 aspect2d(c,r)=0.0
              end if
            end do
          end do
          aspect2d=aspect2d*3.14159/180.
        endif

        if(hydro_module.eq.1) then
          print*,"reading mapfile ",slopelength_map
          call readmapr(nx,ny,slopelen2d,slopelength_map,ireturn)
          if(ireturn.eq.0)then 
            print *, slopelength_map//" does not exist"
            stop
          endif
        end if

        print*,"reading mapfile ",slopeangle_map
        call readmapr(nx,ny,slope2d,slopeangle_map,ireturn)
        if(ireturn.eq.0)then  
          print *, slopeangle_map//" does not exist"
          stop
        else
          slope2d=atan(slope2d)
        endif

        print*,"reading mapfile ",landuse_map
        call readmapi(nx,ny,landuse2d,landuse_map,ireturn)
        if(ireturn.eq.0)then  
          print *, landuse_map//" does not exist"
          landuse2d=10
        endif
        do r=1,ny
           do c=1,nx
             if (landuse2d(c,r).eq.0) then
               inbasin2d(c,r)=-9999
               landuse2d(c,r)=10
             end if
           end do
        end do
        landuse=landuse2d

        print*,"reading mapfile ",lat_map
        call readmapr(nx,ny,alatud2d,lat_map,ireturn)
        if(ireturn.eq.0)then  ! no data is available for zref
          print *, lat_map//" does not exist"
          stop
        else
          alatud2d=alatud2d*3.14159/180.
        endif

        print*,"reading mapfile ",lon_map
        call readmapr(nx,ny,lontitud2d,lon_map,ireturn)
        if(ireturn.eq.0)then  ! no data is available for zref
          print *, lon_map//" does not exist"
          stop
        else
          hrnoon2d=12.0-(lontitud2d-timezone*15.0)/15.0
          lontitud2d=lontitud2d*3.14159/180.
        endif

        print*,"reading mapfile ",soiltype_map
        if (soilfromtable.eq.1) then
          call readmapi(nx,ny,soil_code,soiltype_map,ireturn)
          if(ireturn.eq.0)then  
            print *, soiltype_map//" does not exist"
            stop
          endif
        endif

        print*,"reading mapfile ",tsavg_map
        call readmapr(nx,ny,tsavg2d,tsavg_map,ireturn)
        if(ireturn.eq.0)then 
          print *, tsavg_map//" doese not exist"
          tsavg2d(:,:)=0.0
        endif

        do r=1,ny
          do c=1,nx
            if(inbasin2d(c,r).eq.1 .and. tsavg2d(c,r) >200) then
              print*,"The input tsavg is wrong, it should be in celius degree"
              tsavg2d(c,r) =tsavg2d(c,r) -273.15
            end if
          end do
        end do

        if (shadeef.eq.1) then
          print*,"reading mapfile ",skyview_map
          call readmapr(nx,ny,skyview2d,skyview_map,ireturn)
          if(ireturn.eq.0)then 
            print *, skyview_map//" does not exist"
            stop
          endif
        end if

!       read soil parameter from the paratable based on soil type
        if(soilfromtable.eq.0)  then
          call ReadSoilFromNC(zs2d,soilpro_nc,ns)
        else
          call ReadSoilFromTable(soilpara_table,nsoiltype,ns,soil_code)
        end if


!       reading namelist vegparameter
        print *, "reading namelist vegparameter"
        read(controlfile,vegparameter,err=106,iostat=infil)
106     if(infil.ne.0) then
          print *,'error in reading namelist vegparameter. program stopped.'
          stop
        end if
        call ReadVegPara()


        read(controlfile,residue,err=107,iostat=infil)
107     if(infil.ne.0) then
          print *,'error in reading namelist residue. program stopped.'
          stop
        end if
        if (nr.gt.0) then
          call ReadResiduePara()
        end if

        read(controlfile,initialstate,err=108,iostat=infil)
108     if(infil.ne.0) then
          print *,'error in reading namelist initialstate. program stopped.'
          stop
        end if
        call read_initial_state(restart)






    read(controlfile,levelout,end=103)
      lvlout(1)=lvlout1
      lvlout(2)=lvlout2
      lvlout(3)=lvlout3
      lvlout(4)=lvlout4
      lvlout(5)=lvlout5
      lvlout(6)=lvlout6
      lvlout(7)=lvlout7
      lvlout(8)=lvlout8
      lvlout(9)=lvlout9
      lvlout(10)=lvlout10
      lvlout(11)=lvlout11
      lvlout(12)=lvlout12
      lvlout(13)=lvlout13
    print*,"lvlout",lvlout(1:13)


    close(controlfile)
    return

end  subroutine input




subroutine readmapr(nx,ny,var,flnm,ireturn)
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4 => shr_kind_r4
    implicit none
    integer nx
    integer ny
    real(r8),dimension(nx,ny)::    var
    character*(*) flnm
    integer ireturn ! flag for data reading =1, data was read
                      ! =0, failed reading datac
    integer flunit
    integer ldirnam
    integer lenfl
    logical fexist
    character*5 ncols,nrows
    character*9 xllcorner,yllcorner
    character*8 cellsize
    character*12 nodata
    integer nxin,nyin
    real    x0,y0
    real    gridsize, znodata
    character*200 filename
    integer  i, j



    lenfl = len(flnm)
    call strlnth(flnm,lenfl)
    filename = flnm(1:lenfl)
    inquire(file=filename(1:lenfl),exist = fexist )
    if(.not.fexist ) then
      write(6,'(/1x,a,/1x,a/)') 'file '//filename(1:lenfl)//' is not found.', &
       ' default value will be used.'
       ireturn = 0
       return
    end if

    call getunit(flunit)
    open(flunit,file=filename(1:lenfl))
    read (flunit,*,err=998) ncols,nxin
    read (flunit,*,err=998) nrows,nyin
    read (flunit,*,err=998) xllcorner,x0
    read (flunit,*,err=998) yllcorner,y0
    read(flunit,*,err=998) cellsize,gridsize
    read(flunit,*,err=998) nodata,znodata

    if( (nx.ne.nxin) .or. (ny.ne.nyin) ) then
      write(6,'(a/a/a,i5,a,i5/a/a,i5,a,i5/a)') ' array dimension(s) inconsistent with model definitions,',&
        ' dimensions in the input file were',&
        ' nx = ', nxin, ', ny = ', nyin,&
        ' the model definitions were',&
        ' nx = ', nx,   ', ny = ', ny,&
        ' default value will be used.'
      ireturn = 0
      return
    endif


    do j = 1, ny
      read (flunit,*,err=998) (var(i,j), i = 1, nx)
    enddo
    close(flunit)
    call retunit(flunit)
    ireturn = 1
    return

998 write (6,'(/a,i2/a)')&
     '     read error in data file ' &
     //filename(1:lenfl)//' with the i/o unit ',flunit,&
     ' default value will be used.'
    close(flunit)
    call retunit(flunit)
    ireturn = 0
    return
end subroutine readmapr




 subroutine readmapi(nx,ny,var,flnm,ireturn)
    implicit none

    integer(i4),intent(in):: nx,ny     ! number of grid points in the x-dir (east/west) and y-dir(north/south)
    integer(i4),intent(inout)::var(nx,ny)
    character*(*),intent(in)::flnm
    integer(i4),intent(inout):: ireturn ! flag for data reading  :=1, data was read  :=0, failed reading data 
  
    integer(i4):: flunit
    integer(i4):: lenfl,ldirnam
    logical:: fexist
    character*(5):: ncols,nrows
    character*(9):: xllcorner,yllcorner
    character*(8):: cellsize
    character*(12):: nodata
    integer(i4):: nxin,nyin
    real(r8)::    x0,y0
    real(r8)::    gridsize, znodata
    character*(80):: filename
    integer(i4)::  i, j


    lenfl = len(flnm)
    call strlnth(flnm,lenfl) 

    filename = flnm(1:lenfl)
    inquire(file=filename(1:lenfl),exist = fexist )

    if(.not.fexist ) then
      write(6,'(/1x,a,/1x,a/)') 'file '//filename(1:lenfl)//' is not found.', &
        ' default value will be used.'
      ireturn = 0
      return
    end if

    call getunit(flunit) 
    open(flunit,file=filename(1:lenfl),status='old')

    read (flunit,*,err=998) ncols,nxin
    read (flunit,*,err=998) nrows,nyin
    read (flunit,*,err=998) xllcorner,x0
    read (flunit,*,err=998) yllcorner,y0
    read(flunit,*,err=998) cellsize,gridsize
    read(flunit,*,err=998) nodata,znodata

    if( (nx.ne.nxin) .or. (ny.ne.nyin) ) then
      write(6,'(a/a/a,i5,a,i5/a/a,i5,a,i5/a)') ' array dimension(s) inconsistent with model definitions,', &
        ' dimensions in the input file were', &
        ' nx = ', nxin, ', ny = ', nyin,&
        ' the model definitions were',&
        ' nx = ', nx,   ', ny = ', ny,&
        ' default value will be used.'
      ireturn = 0
      return
    endif
    do j = 1, ny
      read (flunit,*,err=998) (var(i,j), i = 1, nx)
    enddo 
    close(flunit)
    call retunit(flunit)
    ireturn = 1        
    return
998 write (6,'(/a,i2/a)') '     read error in data file ' &
    //filename(1:lenfl)//' with the i/o unit ',flunit, &
    ' default value will be used.'
    close(flunit)
    call retunit(flunit)
    ireturn = 0
end subroutine readmapi


subroutine strlnth(string, length )
    implicit none
    character string*(*)
    integer length
    integer i

    do 100 i = length,1,-1
      if(string(i:i) .ne. ' ') goto 200
    100   continue
200 continue
    length = max(1,i)
    return
end subroutine strlnth



subroutine getunit(nunit )
    implicit none
    integer nunit
    logical used
    integer list(10), nfree
    save list, nfree
    data list /40,41,42,43,44,45,46,47,48,49/
    data nfree /10/
5   if( nfree .lt. 1) then
      write(6,'(1x,a,a)') 'no more unit number is available from the list. ',&
        'job stopped in getunit.'
      stop
    endif
    nunit = list(nfree)
    nfree = nfree-1
    inquire( unit=nunit, opened=used)
      if( used ) goto 5
    return


    entry retunit( nunit )
    inquire( unit=nunit, opened=used)
    if(used) return
    nfree = nfree + 1
    if( nfree .le. 10 ) list( nfree ) = nunit
      nfree = min( nfree, 10)
    return
end subroutine getunit

      subroutine strlen(str,l1,l2)
        implicit none
        character(200)::str
        integer::i,l1,l2,k
        k=0
        do i = 1, 200
          if(k.eq.0 .and. str(i:i).ne.' ') then
            l1=i
            k=1
          elseif(k.eq.1 .and. str(i:i).eq.' ') then
            l2 = i-1
            exit
          endif
        end do
      end subroutine


!***********************************************************************
!
subroutine dayinp2 (julian,year,lastyear,mtstep,mpltgro,mwatrxt,sunhor,tmpday,winday,humday,&
    precip,snoden,soilxt,nplant,plthgt,dchar,pltwgt,pltlai,rootdp,presur,shadow)
!
!     this subroutine reads the time varying data required to run the
!     program, and stores it in arrays with values for the end of each
!     hour.  data such as soil temperature and moisture content at
!     depth, which are seldom available at every hour, are
!     interpolated to get hourly values.
!
!***********************************************************************
    use netcdf
    use dims_mod
    use controlpara_mod,only:lai_dir,forcing_dir,shadows_dir
    use controlpara_mod,only:shadeef
    implicit none


!   input argument
    integer(i4),intent(in)::julian,year,lastyear,mtstep,mpltgro,mwatrxt
!   forcing in a day
    real(r8),dimension(nx,ny,24),intent(inout)::sunhor,tmpday,winday,humday,precip,snoden,presur,shadow
    integer(i4),intent(in)::nplant
!   
    real(r8),dimension(nx,ny,nsmax,24),intent(inout)::soilxt
!
    real(r8),dimension(nx,ny,npmax),intent(inout)::plthgt,dchar,pltwgt,pltlai,rootdp


! this is the name of the data file we will create.
    character (len = *), parameter :: asc_file = 'shaw.wea'
    character (len = *), parameter :: ele_file = 'elevation.asc'

    integer, parameter :: ndims=3
    character (len=4)::syear
    character (len=3)::sjulian


    integer::ncid
    integer::temp_id,win_id,hum_id,pre_id,snow_id,pres_id,solar_id,lai_id,shadow_id
    character*(200)::ncfilename

    integer::col,row,i,j,r,c,h

!****************read weather data
!
!     check flag (mtstep) for type of weather data:  0 = hourly; 1 = daily;   2 ==> data matches nhrpdt

!     daily input weather data

    if (mtstep.eq.1) then
      print *, "the model need hourly forcing data"
      print *, "the daily data is not support yet"
      stop
    end if


    if (mtstep .eq. 0 ) then
      write(syear,"(i4.4)") year
      write(sjulian,"(i3.3)") julian
      ncfilename=trim(forcing_dir)//trim("/"//syear//"/"//syear//"-"//sjulian//".nc")
      call check( nf90_open(ncfilename, nf90_nowrite, ncid) )
      call check( nf90_inq_varid(ncid, "temp", temp_id) )
      call check( nf90_inq_varid(ncid, "wind", win_id) )
      call check( nf90_inq_varid(ncid, "rhum", hum_id) )
      call check( nf90_inq_varid(ncid, "prec",  pre_id) )
      call check( nf90_inq_varid(ncid, "srad",  solar_id ) )
      call check( nf90_inq_varid(ncid, "pres",  pres_id ) )

    ! define the dimensions. netcdf will hand back an id for each. 
      call check( nf90_get_var(ncid, temp_id, tmpday) )
      call check( nf90_get_var(ncid, win_id, winday) )
      call check( nf90_get_var(ncid, hum_id, humday) )
      call check( nf90_get_var(ncid, pre_id, precip) )
      call check( nf90_get_var(ncid, solar_id, sunhor) )
      call check( nf90_get_var(ncid, pres_id, presur) )
      snoden=0.0
      !tmpday=tmpday-273.15
      precip=precip/1000.0
      humday=humday/100.0
      do r=1,ny
        do c=1,nx
          if (inbasin2d(c,r).eq.1) then
            do h=1,24
              if(tmpday(c,r,h)< -200.0 .or. tmpday(c,r,h)> 200.0) then
                print*,"The input air temperature is in degree",tmpday(c,r,h)
                inbasin2d(c,r)=-9999
              end if
              if(presur(c,r,h)< 40000 .or. presur(c,r,h)< 0) then
                print*,"The input pressure is wrong"
                stop
              end if
              if(winday(c,r,h)< 0) then
                print*,"The input winday is wrong"
                stop
              end if
              if(humday(c,r,h)< 0.05) humday(c,r,h)= 0.05
              if(humday(c,r,h)< 0 .or. humday(c,r,h)> 2) then
                print*,"The input humday is wrong"
                stop
              end if
              if(precip(c,r,h)< 0) then
                print*,"The input precip is wrong"
                stop
              end if
              if(sunhor(c,r,h)< 0) then
                print*,"The input sunhor is wrong"
                stop
              end if
            end do
          end if
        end do
      end do
      call check( nf90_close(ncid) )
    end if

!   unit convert
    do col=1,nx
      do row=1,ny
        do i=1,24
          if (winday(col,row,i).lt.thresh) winday(col,row,i)=thresh
        enddo
      enddo
    enddo

      write(syear,"(i4.4)") year
      write(sjulian,"(i3.3)") julian
    if(shadeef.eq.1)then
      ncfilename=trim(shadows_dir)//trim("2014"//"-"//sjulian//".nc")
      call check( nf90_open(ncfilename, nf90_nowrite, ncid) )
      call check( nf90_inq_varid(ncid, "shawdow", shadow_id) )
      call check( nf90_get_var(ncid, shadow_id, shadow) )
      call check( nf90_close(ncid) )
    end if
!***  read canopy information
    if (lai_from_table.eq.0 .and. lastyear .ne. year) then
      print*,"reading lai file"
      ncfilename=trim(lai_dir)//trim("/"//syear//".nc")
      call check( nf90_open(ncfilename, nf90_nowrite, ncid) )
      call check( nf90_inq_varid(ncid, "LAI", lai_id) )
    ! define the dimensions. netcdf will hand back an id for each. 
      call check( nf90_get_var(ncid, lai_id, lai) )
      call check( nf90_close(ncid) )
    endif

  200 if (nplant .eq. 0 .or. mpltgro .eq. 0) go to 250
      do 215 j=1,nplant
         do col=1,nx
         do row=1,ny
          plthgt(col,row,j) = planth(landuse(col,row),julian)
          dchar (col,row,j) = dchar (col,row,j)
          pltwgt(col,row,j) = plantw(landuse(col,row),julian)
          if(lai_from_table .eq. 1) then
            pltlai(col,row,j) = lai_t(landuse(col,row),julian)
          else
            pltlai(col,row,j) = lai(col,row,julian)
          endif
          rootdp(col,row,j) = rtdp(landuse(col,row),julian)
         end do
         end do
  215 continue
  250 ifirst=1
    return

 end subroutine dayinp2


!***********************************************************************
!
 subroutine cloudy2 (clouds,alatud,declin,hafday,sunhor,julian,nhrpdt)
!
!     this subroutine totals the solar radiation for the day.  it then
!     estimates the cloud cover, which will be used in the calculation
!     for the emmissivity of the atmosphere.  (declination and half-day
!     length for the current julian day are also calculated.)
!
!***********************************************************************
    use swrcoe_mod
    implicit none
!input
    real(r8),intent(inout)::clouds,hafday,declin
    real(r8),intent(in)::alatud
    integer(i4)::julian,nhrpdt
    real(r8),dimension(24)::sunhor
!temp
    real(r8)::totsun,coshaf,sunmax,ttotal
    integer(i4)::i,j

    totsun=0.0
    do 10, i=nhrpdt,24
        if (sunhor(i) .lt. 0.0) sunhor(i)=0.0
        totsun=totsun + sunhor(i)
    10 continue
    declin=0.409*sin(2*3.14159*(julian-81)/365.)
    coshaf=-tan(alatud)*tan(declin)
    if (abs(coshaf) .ge. 1.0) then
        if (coshaf .ge. 1.0) then
!           sun does not come up on this day (winter in arctic circle)
            hafday=0.0
        else
!           sun does not set on this day (summer in the arctic circle)
            hafday=3.14159
        end if
    else
        hafday=acos(coshaf)
    end if
    sunmax=24.*solcon*(hafday*sin(alatud)*sin(declin) &
        +cos(alatud)*cos(declin)*sin(hafday))/3.14159
    ttotal=totsun/sunmax
!     use cloud equation based flerchinger & yu (2007), ag & forest met
    clouds=1.333 - 1.666*ttotal
    if (clouds .gt. 1.0) clouds=1.0
    if (clouds .lt. 0.0) clouds=0.0
!
    return
end subroutine cloudy2



  subroutine check(status)
     use netcdf
     implicit none   
     integer, intent ( in) :: status
 
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "stopped"
    end if
  end subroutine check 


      subroutine ReadSoilFromNC(zs,soiltexture,ns)
        use netcdf
        use dims_mod
        use soilproperty_mod
        use constvar_mod
        implicit none

        integer(i4),intent(in)::ns
        character(80),intent(in)::soiltexture
        real(r8),dimension(nx,ny,nsmax)::zs

        integer(i4):: ns1
!       integer inbasin(nx,ny)

!       netcdf variables
        real(r8),dimension(:,:,:),allocatable::sand,silt,clay,om,rhob,sat,satk,b,entry,vapexp,vapcof,fc,thr,&
                                           l,n,alpha

        integer(i4):: ncid, varid
        integer(i4):: x, y, retval
        integer(i4):: row,col,i,allocatestatus,j
        real(r8),dimension(:),allocatable::zs1

        real(r8):: cf,textot
        real(r8)::satcon
    
        print *, "reading soil texture files"
        call check(nf90_open(soiltexture,nf90_nowrite,ncid))

        !call check(nf90_inq_varid(ncid,"ns",varid))
        !call check(nf90_get_var(ncid,varid,ns1))
        !if(ns1.ne.ns) then
        !  print*,"the ns in the netcdf file is:",ns1," is not equal to the model require:",ns
        !  stop
        !endif

        print*,"reading zs"
        allocate(zs1(ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"zs",varid))
        call check(nf90_get_var(ncid,varid,zs1))
        do i=1,nx
          do j=1,ny
            zs(i,j,1:ns)=zs1
          enddo
        enddo

        print*,"reading sa"
        allocate(sand(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"sa",varid))
        call check(nf90_get_var(ncid,varid,sand))
        sand2d(:,:,1:ns)=sand

        print*,"reading si"
        allocate(silt(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"si",varid))
        call check(nf90_get_var(ncid,varid,silt))
        silt2d(:,:,1:ns)=silt

        print*,"reading cl"
        allocate(clay(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"cl",varid))
        call check(nf90_get_var(ncid,varid,clay))
        clay2d(:,:,1:ns)=clay

        print*,"reading bd"
        allocate(rhob(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"bd",varid))
        call check(nf90_get_var(ncid,varid,rhob))
        rhob2d(:,:,1:ns)=rhob

        print*,"reading om"
        allocate(om(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"om",varid))
        call check(nf90_get_var(ncid,varid,om))
        om2d(:,:,1:ns)=om

        print*,"reading thetas"
        allocate(sat(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"thetas",varid))
        call check(nf90_get_var(ncid,varid,sat))
        sat2d(:,:,1:ns)=sat

        print*,"reading ks"
        allocate(satk(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"ks",varid))
        call check(nf90_get_var(ncid,varid,satk))
        satk2d(:,:,1:ns)=satk

        print*,"reading b"
        allocate(b(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"b",varid))
        call check(nf90_get_var(ncid,varid,b))
        b2d(:,:,1:ns)=b

        print*,"reading fc"
        allocate(fc(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"fc",varid))
        call check(nf90_get_var(ncid,varid,fc))
        thfc2d(:,:,1:ns)=fc

        print*,"reading thr"
        allocate(thr(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"thr",varid))
        call check(nf90_get_var(ncid,varid,thr))
        thr2d(:,:,1:ns)=thr

        print*,"reading psi_s"
        allocate(entry(nx,ny,ns),stat=allocatestatus)
        if (allocatestatus /= 0) stop "*** not enough memory for zs(nx,ny,ns)***"
        call check(nf90_inq_varid(ncid,"psi_s",varid))
        call check(nf90_get_var(ncid,varid,entry))
        entry2d(:,:,1:ns)=entry 

        do col=1,nx
          do row=1,ny
            do i=1,ns
              om2d(col,row,i)  =om2d(col,row,i)/100.0
              textot=sand2d(col,row,i)+silt2d(col,row,i)+clay2d(col,row,i)
              cf=1.0/textot
              sand2d(col,row,i)=sand2d(col,row,i)*cf
              silt2d(col,row,i)=silt2d(col,row,i)*cf
              clay2d(col,row,i)=clay2d(col,row,i)*cf
              ! (adjust specific density for organic matter)
              if (sat2d(col,row,i) .gt. (1.- rhob2d(col,row,i)*((1.-om2d(col,row,i))/rhom+om2d(col,row,i)/rhoom))) then
                sat2d(col,row,i)=1.- rhob2d(col,row,i)*((1.-om2d(col,row,i))/rhom+om2d(col,row,i)/rhoom)
              end if
            enddo
          enddo
        enddo

        ! define the dimensions. netcdf will hand back an id for each. 
        call check(nf90_close(ncid))
        deallocate(sand,silt,clay,om,rhob,sat,satk,b,entry)
        print *, "read soil texture successfully"
      end subroutine ReadSoilFromNC
 
 
 
      subroutine ReadSoilFromTable(soilpara_table,nsoiltype,ns,soil_code)
        use dims_mod
        use soilproperty_mod
        use constvar_mod
        implicit none
        character(200),intent(in)::soilpara_table
        integer(i4),intent(in)::nsoiltype,ns
        integer(i4),dimension(nx,ny),intent(in)::soil_code
        
        integer(i4)::l1,l2,flunit,infil,i,j,k,col,row
        character(200)::strtemp
        integer(i4),dimension(50,2)::soilindex
        real(r8),dimension(50,nsmax,20)::soil_para
        real(r8)::cf,textot
        
        print*,"Reading soil parameters from table according to the soil code"
        call getunit(flunit)
        call strlen(soilpara_table,l1,l2)
        open(flunit,file = soilpara_table(l1:l2), status='old',err=10000,iostat=infil)
        read(flunit,*) strtemp
        do i =1,nsoiltype
          read(flunit,*) soilindex(i,1)
          soilindex(i,2)=i
          do j =1,ns
            read(flunit,*) (soil_para(i,j,k),k=1,10)
          end do
        end do
        close(flunit)
10000   if (infil.ne.0) then
          print*,"The following input file cannot be opened: ",soilpara_table
          stop
        end if

        do col=1,nx
          do row=1,ny
            j=0
            do i =1,nsoiltype
              if(soil_code(col,row).eq.soilindex(i,1))then
                j=soilindex(i,2)
              end if
            end do
            if (j.eq.0 .or. j.gt.nsoiltype) then
              print*,"the parameter for soil type:",soil_code(col,row),"does not exist"
              stop
            end if
            do i=1,ns
              zs2d(col,row,i)   = soil_para(j,i,1)
              b2d(col,row,i)    = soil_para(j,i,2)
              entry2d(col,row,i)=soil_para(j,i,3)
              satk2d(col,row,i) =soil_para(j,i,4)/360000.
              rhob2d(col,row,i) =soil_para(j,i,5)
              sat2d(col,row,i)  =soil_para(j,i,6)
              sand2d(col,row,i) =soil_para(j,i,7)
              silt2d(col,row,i) =soil_para(j,i,8)
              clay2d(col,row,i) =soil_para(j,i,9)
              om2d(col,row,i)   =soil_para(j,i,10)/100.0
              thfc2d(col,row,i)=0.25
              thr2d(col,row,i)=0.05
              l2d(col,row,i)=0.0
              alpha2d(col,row,i)=0.0
              n2d(col,row,i)=0.0
              textot=sand2d(col,row,i)+silt2d(col,row,i)+clay2d(col,row,i)
              cf=1.0/textot
              sand2d(col,row,i)=sand2d(col,row,i)*cf
              silt2d(col,row,i)=silt2d(col,row,i)*cf
              clay2d(col,row,i)=clay2d(col,row,i)*cf
              
              ! (adjust specific density for organic matter)
              if (sat2d(col,row,i) .gt. (1.- rhob2d(col,row,i)*((1.-om2d(col,row,i))/rhom+om2d(col,row,i)/rhoom)))then
                 sat2d(col,row,i)=1.- rhob2d(col,row,i)*((1.-om2d(col,row,i))/rhom+om2d(col,row,i)/rhoom)
              end if
            enddo
          enddo
        enddo
        vapcof2d=0.66
        vapexp2d=1.0
      end subroutine ReadSoilFromTable
      
      
      subroutine ReadVegPara()
        use controlpara_mod
        use constvar_mod
        implicit none
        
        integer(i4)::flunit
        integer(i4)::l1,l2,infil,i,j,row,col,k
        character(200)::strtemp
        real(r8)::hum,dummy
!       read veg parameter tables
!       1 parameter table
        call getunit(flunit)
        call strlen(vegpara_table,l1,l2)
        open(flunit,file = vegpara_table(l1:l2), status='old',err=10001,iostat=infil)
        read(flunit,*) strtemp
        do i = 1, 16
          read(flunit,*) (vegtable(i,j),j=1,10)
        end do
        close(flunit)
!       can not open some file of veg tables
10001   if (infil.ne.0) then
          print*,vegpara_table
          stop
        end if

!       monthly lai
        print *,"reading lai table"
        call strlen(lai_table,l1,l2)
        open(flunit,file = lai_table(l1:l2), status='old',err=10002,iostat=infil)
        read(flunit,*) strtemp
        print*,strtemp
        do i = 1, 366
          read(flunit,*) (lai_t(j,i),j=1,16)
        end do
        close(flunit)
!       can not open some file of veg tables
10002   if (infil.ne.0) then
          print*,lai_table
          stop
        end if

!       monthly height
        call strlen(hight_table,l1,l2)
        open(flunit,file = hight_table(l1:l2), status='old',err=10003,iostat=infil)
        read(flunit,*) strtemp
        do i = 1, 366
          read(flunit,*) (planth(j,i),j=1,16)
        end do
        close(flunit)
!       can not open some file of veg tables
10003   if (infil.ne.0) then
          print*,hight_table
          stop
        end if

!       monthly root depth
        call strlen(rootdp_table,l1,l2)
        open(flunit,file = rootdp_table(l1:l2), status='old',err=10004,iostat=infil)
        read(flunit,*) strtemp
        do i = 1, 366
          read(flunit,*) (rtdp(j,i),j=1,16)
        end do
        close(flunit)
!       can not open some file of veg tables
10004   if (infil.ne.0) then
          print*,rootdp_table
          stop
        end if

!       monthly weight
        call strlen(weight_table,l1,l2)
        open(flunit,file = weight_table(l1:l2), status='old',err=10005,iostat=infil)
        read(flunit,*) strtemp
        do i = 1, 366
          read(flunit,*) (plantw(j,i),j=1,16)
        end do
        close(flunit)
        call retunit(flunit)
!       can not open some file of veg tables
10005   if (infil.ne.0) then
          print*,weight_table
          stop
        end if

        hum = 0.999
!       configure plant data
        wcandt2d(:,:,1)=wcandt1
        if (nplant .gt. 0) then
          nplant2d=nplant
          call canhum (2,hum,dummy,wcmax,0.0_r8,canma,canmb)
          do row = 1,ny
            do col = 1, nx          
              if (wcandt2d(col,row,1) .gt. wcmax) then
                wcandt2d(col,row,1)=wcmax
              end if
              do k=1,nplant
                dchar2d (col,row,k)=vegtable(landuse2d(col,row),1)
                dchar2d (col,row,k)=dchar2d (col,row,k)/100.0
                canalb2d(col,row,k)=vegtable(landuse2d(col,row),2)
                tccrit2d(col,row,k)=vegtable(landuse2d(col,row),3)
                pleaf02d(col,row,k)=vegtable(landuse2d(col,row),4)
                rstom02d(col,row,k)=vegtable(landuse2d(col,row),5)
                rstexp2d(col,row,k)=vegtable(landuse2d(col,row),6)
                rleaf02d(col,row,k)=vegtable(landuse2d(col,row),7)
                rroot02d(col,row,k)=vegtable(landuse2d(col,row),8)
                rleaf02d(col,row,k)=1.0/rleaf02d(col,row,k)
                rroot02d(col,row,k)=1.0/rroot02d(col,row,k)
                xangle2d(col,row,k)=vegtable(landuse2d(col,row),9)
                clumpng2d(col,row,k)=vegtable(landuse2d(col,row),10)
                itype2d (col,row,k)=1
              end do
              if(landuse2d(col,row).eq.16) nplant2d(col,row)=0
            end do
          end do
        endif
      end subroutine ReadVegPara
      
      subroutine ReadResiduePara()
        implicit none
        real(r8)::hum
        integer(i4)::flunit,l1,l2,i,j,infil,row,col
        character(200)::strtemp
        real(r8)::dummy

        hum=0.999
        call getunit(flunit)
        call strlen(res_table,l1,l2)
        open(flunit,file = res_table(l1:l2), status='old',err=1004,iostat=infil)
1004    if (infil.ne.0) then
          print*, "residue table file does not exist",res_table(l1:l2)
          stop
        end if
        read(flunit,*) strtemp
        do i = 1, 16
          read(flunit,*) (restable(i,j),j=1,4)
        end do
        close(flunit)
        call retunit(flunit)
        call reshum (2,hum,dummy,gmcmax,0.0_r8)
        gmcdt2d=gmcdt1
        do col=1,nx
          do row=1,ny
!           determine maximum water content of residue (at 99.9% rh)
            if (gmcdt2d(col,row,1) .gt. gmcmax) gmcdt2d(col,row,1)=gmcmax
            cover2d(col,row)=restable(landuse2d(col,row),1)
            albres2d(col,row)=restable(landuse2d(col,row),2)
            rload2d(col,row)=restable(landuse2d(col,row),3)
            rload2d(col,row)=rload2d(col,row)/10000.
            zrthik2d(col,row)=restable(landuse2d(col,row),4)
            zrthik2d(col,row)=zrthik2d(col,row)/100.
            rescof2d(col,row)=rescof
          end do
        end do
      end subroutine ReadResiduePara
      
      
      subroutine read_initial_state(restart)
        use netcdf
        use soilproperty_mod,only:sat2d
        implicit none
        integer(i4),intent(inout)::restart
        real(r8),dimension(:,:,:),allocatable::tsdt,vlcdt,matdt,vicdt
        integer(i4),dimension(:,:,:),allocatable::icedt
        integer(i4)::allocatestatus
        integer(i4)::r,c,countn,col,row,i,j,l1,l2
        integer::tsdt_id,vlcdt_id,vicdt_id,icedt_id,matdt_id,ncid,m
        real::tsmean
        
        
        print*,"reading intial state variables"
        if(restart.gt.0) then
          print*,"Reading initial states from the stored files"
          allocate(tsdt(nx,ny,ns),stat=allocatestatus)
          if (allocatestatus /= 0) stop "*** not enough memory for tsdt(nx,ny,ns)***"
          allocate(vlcdt(nx,ny,ns),stat=allocatestatus)
          if (allocatestatus /= 0) stop "*** not enough memory for vlcdt(nx,ny,ns)***"
          allocate(vicdt(nx,ny,ns),stat=allocatestatus)
          if (allocatestatus /= 0) stop "*** not enough memory for vlcdt(nx,ny,ns)***"
          allocate(matdt(nx,ny,ns),stat=allocatestatus)
          if (allocatestatus /= 0) stop "*** not enough memory for vlcdt(nx,ny,ns)***"
          allocate(icedt(nx,ny,ns),stat=allocatestatus)
          if (allocatestatus /= 0) stop "*** not enough memory for vlcdt(nx,ny,ns)***"

          call strlen(soil_ini,l1,l2)
          call check( nf90_open(trim(soil_ini(l1:l2)), nf90_nowrite, ncid) )
          call check( nf90_inq_varid(ncid, "TSDT2d", tsdt_id) )
          call check( nf90_inq_varid(ncid, "VLCDT2d",  vlcdt_id) )
          call check( nf90_inq_varid(ncid, "VICDT2d",  vicdt_id) )
          call check( nf90_get_var(ncid, tsdt_id, tsdt))
          call check( nf90_get_var(ncid, vlcdt_id, vlcdt))
          call check( nf90_get_var(ncid, vicdt_id, vicdt))
          call check( nf90_close(ncid) )
          tsdt2d(:,:,1:ns)=tsdt
          vlcdt2d(:,:,1:ns)=vlcdt
          vicdt2d(:,:,1:ns)=vicdt
          do col=2,nx-1
            do row=2,ny-1
              if(inbasin2d(col,row)==1) then
              do i=1,ns2d(col,row)                
                tsmean=0.0
                m=0
                if(inbasin2d(col-1,row)==1) then
                  tsmean=tsmean+tsdt2d(col-1,row,i)
                  m=m+1
                end if

                if(inbasin2d(col+1,row)==1) then
                  tsmean=tsmean+tsdt2d(col+1,row,i)
                  m=m+1
                end if

                if(inbasin2d(col,row-1)==1) then
                  tsmean=tsmean+tsdt2d(col,row-1,i)
                  m=m+1
                end if

                if(inbasin2d(col,row+1)==1) then
                  tsmean=tsmean+tsdt2d(col,row+1,i)
                  m=m+1
                end if

                if(m>0) then
                tsmean=tsmean/m
                  if(abs(tsdt2d(col,row,i)-tsmean)>8.0) then
                    tsdt2d(col,row,i)=tsmean
                  end if
                end if
                if (vlcdt2d(col,row,i) .gt. sat2d(col,row,i)) vlcdt2d(col,row,i)=sat2d(col,row,i)
                call matvl1 (i,matdt2d(col,row,i),vlcdt2d(col,row,i),col,row)
                if(vicdt2d(col,row,i).gt. 0.0) then
                   icesdt2d(:,:,:)=1
                else
                   icesdt2d(:,:,:)=0
                end if
              enddo
              end if
            end do
          end do
          !tsdt2d(35,34,:)=tsdt2d(35,34,:)+4.0  !!AR
          !tsdt2d(31,41,:)=tsdt2d(31,41,:)+4.0  !!ARNF
          !tsdt2d(40,29,:)=tsdt2d(40,29,:)+4.0  !!ARSF
          !tsdt2d(74,44,:)=tsdt2d(74,44,:)+4.0  !!EB
          !tsdt2d(58,38,:)=tsdt2d(58,38,:)+3.0  !!HCG
          !tsdt2d(91,55,:)=tsdt2d(91,55,:)+4.0  !!JYL
          matdt2d(:,:,1:ns)=matdt
          icesdt2d(:,:,1:ns)=icedt
          !do col=1,nx
          !  do row=1,ny
          !    if(tsdt2d(col,row,1) .eq. -9999 ) then
          !      tsdt2d(col,row,:)=0.0
          !      countn=0
          !      do r=-1,1
          !        do c=-1,1
          !          if(tsdt2d(col+c,row+r,1).ne.-9999) then
          !            countn=countn+1
          !            do i=1,ns
          !              tsdt2d(col,row,i)=tsdt2d(col,row,i)+tsdt2d(col+c,row+r,i)
          !            end do
          !          end if
          !        end do
          !      end do
          !      tsdt2d(col,row,:)=tsdt2d(col,row,:)/countn
          !    endif
          !  end do
          !end do
        else
          print*,"Setting default initial states for the variable"
          tsdt2d(:,:,1:ns)=0.0
          vlcdt2d(:,:,1:ns)=0.2
          vicdt2d(:,:,1:ns)=0.0
          icesdt2d(:,:,:)=0
          do col=1,nx
            do row=1,ny
              do i=1,ns2d(col,row)
                if (vlcdt2d(col,row,i) .gt. sat2d(col,row,i)) vlcdt2d(col,row,i)=sat2d(col,row,i)
                call matvl1 (i,matdt2d(col,row,i),vlcdt2d(col,row,i),col,row)
              enddo
            end do
          end do
        end if
      end subroutine read_initial_state


 subroutine saxtonhydraulic(sand,clay,silt,organic,e,b,ks,theta_sat)
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
    implicit none
    real(r8),intent(inout)::sand,clay,silt,organic,e,b,ks,theta_sat
    real(r8):: theta1500t,theta1500,theta33t,theta33,theta_s_33_t,theta_s_33,lamda
    real(r8):: a

    if (organic>8.0) organic=8.0
    if (clay>0.60)   clay=0.60
    if (clay<0.02)   clay=0.02
    if (sand<0.02)   sand=0.02
    if (sand>0.96)   sand=0.96

    theta1500t = -0.024*sand + 0.487*clay + 0.006*organic + 0.005* &
      (sand * organic)-0.013*(clay*organic)+0.068*(sand*clay)+0.031
      theta1500 =theta1500t + 0.14*theta1500t- 0.02

    theta33t=-0.251*sand+0.195*clay+0.011*organic+0.006* &
      (sand*organic)-0.027*(clay * organic) + 0.452*(sand * clay) + 0.299
    theta33 = theta33t + (1.283*(theta33t**2) - 0.374*theta33t- 0.015)

    theta_s_33_t=0.278*sand+0.034*clay+0.022*organic-0.018* &
      (sand*organic)-0.027*(clay*organic)-0.584*(sand*clay)+0.078

    theta_s_33=theta_s_33_t+0.636*theta_s_33_t-0.107

    theta_sat = theta33 + theta_s_33 - 0.097*sand + 0.043
!      psiet=-21.67*sand -27.93*clay-81.97*theta_s_33+71.12*sand*theta_s_33+8.29*clay*theta_s_33+14.05*sand*clay+27.16
!      psie=psiet+(0.02*psiet**2-0.113*psiet-0.70)
    b=(log(1500.0)-log(33.0))/(log(theta33)-log(theta1500))
    lamda=1.0/b
    a=exp(log(33.0)+b*log(theta33))
    e=-a*theta_sat**(-b)/9.8
!    #calculated aa,b,psie_1986 according to the equations of saxton (1986)
!    #aa=10.19*numpy.exp(-4.396-0.0715*clay*100.0-0.000488*(sand*100.0)**2.0-0.00004285*clay*100.0*(100.0*sand)**2.0)
!    #b=3.14+0.00222*(clay*100.0)**2.0+0.00003484*clay*100.0*(100.0*sand)**2.0
!    #psie_1986=-aa*theta_sat**(-b)
!    #calculated aa,b,psie_1986 according to the equations of saxton (1986)
    ks=1930.0*(theta_sat-theta33)**(3.0-lamda)
  end subroutine saxtonhydraulic
  
end module input_mod
