module hydro_mod
 use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8
 use dims_mod
 use input_mod,only:strlnth,GETUNIT,RETUNIT
 implicit none
!  num_sub:    number of subcatchment
!  num_flow:   number of flow interval for each sub-catchment
!  sub_catch:  name of the subcatchment using Pfafstetter numbering system
!  river_file: name of river parameter file for each sub-catchment
!
!  ngrid:     number of grids in this flow interval
!  grid_row:  row number of grid in this flow interval
!  grid_col:  column number of grid in this flow interval
!  rdx:       length of flow intervals (m)
!  s0:        river-bed slope of the flow interval
!  b:         river width of the flow interval (m)
!  roughness: river roughness (Manning's coefficient) of the flow interval
!  Dr:        river depth of the flow interval (m)
!
!  qin:       lateral inflow of one flow interval (m3/sec/m)
!  hr:        water depth in river (m)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    public
    save
    integer(i4),parameter::NSUBMAX=100  !sub catchment numbers    
    integer(i4),parameter::NF=300
    integer(i4),parameter::NG=1000
    CHARACTER(5),dimension(100)::SUB_CATCH
    integer(i4)::NUM_SUB
    integer(i4),dimension(100)::NUM_FLOW
    integer(i4),dimension(100,300)::ngrid
    integer(i4),dimension(100,300,1000)::grid_row,grid_col
    real(r8),dimension(100,300)::rdx,s0,riverb,Drr,roughness,hhr
    real(r8),dimension(300)::qin
    integer(i4),dimension(9)::kfs1
    integer(i4),dimension(9,9)::kfs2
    integer(i4),dimension(9,9,9)::kfs3
    integer(i4),dimension(9)::p1,p2,p3
!      COMMON /subcatchment/  NUM_SUB,sub_catch
!      COMMON /flow_interval/ num_flow,ngrid,grid_row,grid_col
!      COMMON /river/         rdx,s0,riverb,roughness,Drr
!      COMMON /inflow/        qin                    
!      COMMON /river_water/   hhr
!      COMMON   /pfaf/  p1,p2,p3,kfs1,kfs2,kfs3     ! Pfaftetter basin number  

    integer(i4),parameter::ntl=3 ! Pfafstetter numbering system
    integer(i4),parameter::SCEflag = 0 ! =1 for SCE optimization; = 0 for no SCE optimization. 

    real(r8),dimension(NSUBMAX,300)::q1,q2
    real(r8),dimension(NSUBMAX,300)::qlin1,qlin2
!   common /discharge/q1,q2,qlin1,qlin2


!   DATA BLOCK
    DATA      kfs1/9*0/
    DATA      kfs2/81*0/
    DATA      kfs3/729*0/
!   DATA BLOCK END
  contains

  SUBROUTINE HYDRO_PARA(para_dir)
!    ##################################################################
!    ##################################################################
!    ######                                                      ######
!    ######                 SUBROUTINE hydro_para                ######
!    ######                                                      ######
!    ######                     Developed by                     ######
!    ######     River and Environmental Engineering Laboratory   ######
!    ######                University of Tokyo                   ######
!    ######                                                      ######
!    ##################################################################
!    ##################################################################
!************************************************************************
!                                                                       *
!     read sub-basin parameter                                  *
!                                                                       *
!************************************************************************
!     num_flow:   number of flow interval for each sub-catchment
!     ngrid:     number of grids in this flow interval
!     rdx:       length of flow intervals (m)
!     s0:        river-bed slope of the flow interval
!     b:         river width of the flow interval (m)
!     roughness: river roughness (Manning's coefficient) of the flow interval
!     Dr:        river depth of the flow interval (m)
!     grid_row:  row number of grid in this flow interval
!     grid_col:  column number of grid in this flow interval
!-----------------------------------------------------------
    use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
    use dims_mod
    implicit none
!input
    character(80),intent(in)::para_dir
!temp variable
    character(200)::RIVER_FILE
    integer(i4)::i,j, ii, LFILE, LDIR, FLUNIT, PARUNIT,NLEVEL,LEVEL
    character(200)::ATMP
    integer(i4)::l1,l2,l3,nsub,IOstatus

!------------------------ reading Horton river topology
      call strlen(para_dir,l1,l2)
      open(1,file = para_dir(l1:l2)//'subbasin.dat', status='old')
      nsub=0
      IOstatus=1
      do while(IOstatus.ge.0)	  
        READ(1,*,IOSTAT=IOstatus) ATMP
        nsub=nsub+1
      end do
      close(1)
      nsub=nsub-1
      print*,nsub

      call strlen(para_dir,l1,l2)
      open(1,file = para_dir(l1:l2)//'subbasin.dat', status='old')
      do i = 1,nsub
          read(1,*) num_sub,psubbasin(i),(pbasinup(i,j),j=1,8)
      enddo
      close(1)


!------------------------

!-------------------------------------------------------------------
!           generate subcatchment name from kfs.dat
!-------------------------------------------------------------------
      LDIR = LEN(PARA_DIR)
      CALL STRLNTH(PARA_DIR,LDIR)
      CALL GETUNIT(PARUNIT)
      OPEN(PARUNIT,FILE=PARA_DIR(1:LDIR)//'kfs.dat',STATUS='OLD')
      READ(PARUNIT,*) ATMP, NLEVEL
      READ(PARUNIT,*)
      DO LEVEL = 1, NLEVEL-1
        IF(LEVEL .EQ.1) THEN                 ! Level 1
          READ(PARUNIT,*)
          READ(PARUNIT,*) ATMP, (kfs1(l1), l1=1,9)
        ELSEIF(level .EQ. 2) THEN            ! Level 2
          READ(PARUNIT,*)
          DO l1=1,9
            IF(kfs1(l1).EQ.1) READ(PARUNIT,*)ATMP,(kfs2(l1,l2), l2=1,9)
          END DO
        ELSEIF(level .EQ. 3) THEN            ! Level 3
          READ(PARUNIT,*)
          DO l1=1,9
            IF(kfs1(l1).EQ.1) THEN
              DO l2=1,9
              IF(kfs2(l1,l2) .EQ. 1)  READ(PARUNIT,*) ATMP,(kfs3(l1,l2,l3), l3=1,9)
              END DO
            END IF
          END DO
        END IF
      END DO
      CLOSE(PARUNIT)
      CALL RETUNIT(PARUNIT)
      NUM_SUB=0
      DO l1 = 9, 1, -1                 ! Level 1
        IF(kfs1(l1) .EQ. 0) THEN
          NUM_SUB=NUM_SUB+1
          sub_catch(NUM_SUB)='ws'//char(48+l1)//char(48+0)//char(48+0)
!	  PRINT *,NUM_SUB,'   ',sub_catch(NUM_SUB)
        ELSE
          DO l2 = 9, 1, -1                 ! Level 2
            IF(kfs2(l1,l2) .EQ. 0) THEN
            NUM_SUB=NUM_SUB+1
            sub_catch(NUM_SUB)='ws'//char(48+l1)//char(48+l2)//char(48+0)
!	    PRINT *,NUM_SUB,'   ',sub_catch(NUM_SUB)
            ELSE
            DO l3 = 9, 1, -1                 ! Level 3
            NUM_SUB=NUM_SUB+1
            sub_catch(NUM_SUB)='ws'//char(48+l1)//char(48+l2)//char(48+l3)
!            PRINT *,NUM_SUB,'   ',sub_catch(NUM_SUB)
            END DO
            END IF
          END DO
        END IF
      END DO
      PRINT *,'NUM_SUB:',NUM_SUB
!	  read sub_catchment parameters
!-------------------------------------------------------------------
      PRINT *,'Reading catchment parameters ...'
      DO ii=1,NUM_SUB
        RIVER_FILE = sub_catch(ii)//'_river'
!        PRINT *,ii,sub_catch(ii)
        LFILE = LEN(RIVER_FILE)
        CALL STRLNTH(RIVER_FILE,LFILE)
        LDIR = LEN(PARA_DIR)
        CALL STRLNTH(PARA_DIR,LDIR)
        CALL GETUNIT(FLUNIT)    !add
        OPEN(FLUNIT,FILE=PARA_DIR(1:LDIR)//RIVER_FILE(1:LFILE),STATUS='OLD')
        READ(FLUNIT,*) num_flow(ii)
        DO i=1,num_flow(ii)
          READ(FLUNIT,*) ngrid(ii,i),rdx(ii,i),s0(ii,i),riverb(ii,i),roughness(ii,i),Drr(ii,i)

          IF(ngrid(ii,i) .GT. 1000)PRINT *,'more than 1000 grids:',ngrid(ii,i),ii,i
          READ(FLUNIT,*)(grid_row(ii,i,j),grid_col(ii,i,j), j=1,ngrid(ii,i))
          IF(s0(ii,i).LE.0.001) s0(ii,i)=0.001
        END DO
        CLOSE(FLUNIT)
        CALL RETUNIT(FLUNIT)
      ENDDO
  END SUBROUTINE HYDRO_PARA




  SUBROUTINE LATERAL_ROUTING(JustStart,INITAL,dthydro,inicon,AREA,SLP2D,LEN2D,Dr2d,Drw2d,&
	runoff_inter, runoff_g,waterflx,JULIAN,HOUR,YEAR,YREND,JEND,HourEND,simulation_dir,RIVERNET,&
	RUNOFF,POND,PONDMX,totsurflow,totintflow,totgflow)
!    ##################################################################
!    ##################################################################
!    ######                                                      ######
!    ######                 SUBROUTINE hydro_para                ######
!    ######                                                      ######
!    ######                     Developed by                     ######
!    ######     River and Environmental Engineering Laboratory   ######
!    ######                University of Tokyo                   ######
!    ######                                                      ######
!    ##################################################################
!    ##################################################################
!    ##################################################################
!              (Mudule of Lateral Water Flow in sib-gbhm)
!           (1) Using basin unit (Pfafstetter Basin Scheme)
!           (2) Top-morphology parameterization (River-Hillslope Concept)
!           (3) Kinematic water (surface and river) routing
!    ##################################################################
    use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
    use dims_mod,only:NX,NY
    implicit none
!input
    integer(i4),intent(in)::JustStart,inital,inicon
    real(r8),intent(in)::dthydro
    real(r8),dimension(NX,NY),intent(in)::AREA 		! area of each GRID (m2)
    real(r8),dimension(NX,NY),intent(in)::SLP2D		! slope of ground surface (ND)
    real(r8),dimension(NX,NY),intent(in)::LEN2D		! hillslope length of each GRID (m)
    real(r8),dimension(NX,NY),intent(inout)::Dr2d 	! river depth in each GRID (m)
    real(r8),dimension(NX,NY),intent(inout)::Drw2d 	! water depth in the river of each GRID (m)
    real(r8),dimension(NX,NY),intent(inout)::runoff_inter,runoff_g,waterflx,RUNOFF,POND
    real(r8),intent(in)::PONDMX
    real(r8),intent(inout)::totsurflow,totintflow,totgflow
    integer(i4),intent(in)::JULIAN,HOUR,YEAR,YREND,JEND,HourEND
    character(80),intent(in):: simulation_dir
    integer(i4),dimension(NX,NY),intent(in)::RIVERNET
!temp variable
    integer(i4):: i, ig, ii, l1,l2,l3,j,outflag
    integer(i4):: ix,iy		!Column number and row number for each grid in the domain
!    character(100)::atmp
    integer(i4)::flowunit
    integer(i4),dimension(NX,NY)::flowno !which flow interval of the grid belongs in the subcatchment
    integer(i4):: m, num_grid(NSUBMAX)  !the number of grids in each subcatchment
    integer(i4):: endflow          !for output Murakami-gauge's data 12/25/2006
!    REAL(r8),dimension(NX,NY)::Sst, Sstmax0
    real(r8)::dx,dy,xsw,ysw
    integer(i4),dimension(12)::dayinmonth
    DATA       dayinmonth/31,28,31,30,31,30,31,31,30,31,30,31/
    character(80)::outfilename

    outfilename="./output/flow.asc"
    outflag = 1

      IF(outflag == 1)THEN		
        IF(inital .EQ. 0) THEN
          CALL getunit(flowunit)
          OPEN(flowunit,FILE=outfilename, status='unknown')
          DO iy = 1, NY
            DO ix = 1, NX
	           flowno(ix,iy) = -9999
            END DO
          END DO


          DO ii=1,num_sub
            num_grid(ii) = 0
	    DO i=1,num_flow(ii)
	      DO ig=1,ngrid(ii,i)
		iy=grid_row(ii,i,ig)
		ix=grid_col(ii,i,ig)
	        Dr2d(ix,iy)=Drr(ii,i)
                flowno(ix,iy)= i
                num_grid(ii) = num_grid(ii)+1
	      END DO
            END DO
          END DO
          WRITE(flowunit,*) 'ncols',NX
          WRITE(flowunit,*) 'nrows',NY
          WRITE(flowunit,*) 'xllcorner',xsw
          WRITE(flowunit,*) 'yllcorner',ysw
          WRITE(flowunit,*) 'cellsize',dx
          WRITE(flowunit,*) 'NODATA_value  -9999'
          Do iy = 1, ny
            WRITE(flowunit,*)  (flowno(ix,iy),ix = 1, nx)
          END DO
          CLOSE(flowunit)
          CALL retunit(flowunit)
        END IF
      END IF


      ii=0
      DO l1 = 9, 1, -1                                           ! Level 1
        IF(kfs1(l1) .EQ. 0) THEN
          ii=ii+1
!     Just considering specific subcathments
	  IF(ii>=1 .AND. ii<=num_sub) THEN
           CALL SURFACE_ROUTING(JustStart,ii,AREA,SLP2D,LEN2D,runoff_inter,runoff_g,waterflx,dthydro,&
		RIVERNET,RUNOFF,POND,PONDMX,totsurflow,totintflow,totgflow)
           CALL RIVER_ROUTING(dthydro,INITAL,JustStart,inicon,ii,1,l1,l2,l3,JULIAN,HOUR,YEAR,YREND,JEND,HourEND,&
		simulation_dir,totsurflow,totintflow,totgflow)
          ENDIF
        ELSE
          DO l2 = 9, 1, -1                                       ! Level 2
            IF(kfs2(l1,l2) .EQ. 0) THEN
              ii=ii+1
	        IF(ii>=1 .AND. ii<=num_sub) THEN
                CALL SURFACE_ROUTING(JustStart,ii,AREA,SLP2D,LEN2D,runoff_inter,runoff_g,waterflx,dthydro,&
		RIVERNET,RUNOFF,POND,PONDMX,totsurflow,totintflow,totgflow)
                CALL RIVER_ROUTING(dthydro,INITAL,JustStart,inicon,ii,2,l1,l2,l3,JULIAN,HOUR,YEAR,YREND,JEND,HourEND,&
		simulation_dir,totsurflow,totintflow,totgflow)
                END IF
            ELSE
              DO l3 = 9, 1, -1                                     ! Level 3
                ii=ii+1
	          IF(ii>=1 .AND. ii<=num_sub) THEN
                  CALL SURFACE_ROUTING(JustStart,ii,AREA,SLP2D,LEN2D,runoff_inter,runoff_g,waterflx,dthydro,RIVERNET,&
		RUNOFF,POND,PONDMX,totsurflow,totintflow,totgflow)
                  CALL RIVER_ROUTING(dthydro,INITAL,JustStart,inicon,ii,3,l1,l2,l3,JULIAN,HOUR,YEAR,YREND,JEND,HourEND,&
		simulation_dir,totsurflow,totintflow,totgflow)
                  END IF
		  p3(l3)=ii
              END DO
            END IF
            p2(l2)=ii
          END DO
        ENDIF
        p1(l1)=ii
      END DO

!     After the surface routing and the river routing, Update the river water depth
      DO ii=1,num_sub
        DO i=1,num_flow(ii)
          DO ig=1,ngrid(ii,i)
            iy=grid_row(ii,i,ig)
            ix=grid_col(ii,i,ig)
            Drw2d(ix,iy)=hhr(ii,i)
          END DO
        END DO
      END DO
  END SUBROUTINE LATERAL_ROUTING


  SUBROUTINE SURFACE_ROUTING(JustStart,ii,AREA,SLP2D,LEN2D,runoff_inter,runoff_g,waterflx,&
	dthydro,RIVERNET,RUNOFF,POND,PONDMX,totsurflow,totintflow,totgflow)
!**********************************************************************
!     Input:
!     	ii,nx,ny,nvtyps,inbasin,area,slp2d,len2d,vtypfrct,vegtyp,runoff_inter,runoff_g,waterflx
!     	Sstmax0
!     INOUT:
!     	Sst
!     OUT:
!     	qin(ii)
!**********************************************************************
!     Variables Definition:
!	  ii:	indicator of subcatchment loop
!	  i:    indicator of flow interval loop
!         ig:   indicator of grid loop within a flow interval
!
!
!	  qg_hr:       subsurface runoff of one flow interval (m3/s/m)
!	  q_hr:        surface runoff of one simulation unit (m/s, before flow into river)
!	  q_hillslope: surface runoff of one simulation unit (m3/s/m, flow into river)
!
!	  runoff_sub: subsurface runoff of a grid (m)
!	  runoff_sur: surface runoff of a grid (m)
!	  Sst:     surface storage of each grid (m)
!	  Sst_max: maximum surface detension storage of each grid (m)
!
!	  year:  calendar year
!	  month: month in calendar
!	  day:   day in calendar
!         hour:  hour of the reference day
!         minute:minute of the reference hour
!         second:second of the reference minute
!         jday:  Julian day starting from Jan. 1st of the year
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
    use dims_mod,only:NX,NY
    implicit none
!input
    integer(i4),intent(in)::JustStart,ii
    real(r8),dimension(NX,NY),intent(in)::AREA,slp2d,LEN2D
    real(r8),dimension(NX,NY),intent(in)::runoff_inter!m/s
    real(r8),dimension(NX,NY),intent(in)::runoff_g    !m/s
    real(r8),dimension(NX,NY),intent(in)::waterflx    !m/s
    real(r8),intent(in)::dthydro
    integer(i4),dimension(NX,NY),intent(in)::RIVERNET
    real(r8),dimension(NX,NY),intent(inout)::RUNOFF
    real(r8),dimension(nx,ny),intent(inout)::POND
    real(r8),intent(in)::PONDMX
    real(r8),intent(inout)::  totsurflow,totintflow,totgflow
!temp variable
    integer(i4):: i, ig,ic,ir,ix,iy
    real(r8):: dt
    real(r8):: basin_area
    real(r8),dimension(NSUBMAX)::area_sub ! basin-area, subcatchment area
    real(r8)::  qin_inter            !Subsurface flow In each grid
    real(r8)::  qin_g                !exchanges between groundwater and river In each grid
    real(r8)::  qin_sfc              !surface runoff from a flow interval In each grid
    real(r8)::  qin_tmp              !totoal runoff from a flow interval

    real(r8),dimension(NX,NY)::Sst,Sstmax0
    SAVE  basin_area,area_sub
    real(r8)::power,q_hillslope,q_surf,roughness_n,s0_slope,water_depth


    dt = dthydro ! second

!**********************************************************
!                     model initiation
!***********************************************************
      IF(JustStart == 1) THEN
        IF(ii .EQ. 1) THEN
          basin_area = 0.0
        ENDIF
!     make a sum of the total area of the model domain and subcatments
        area_sub(ii)=0.0
        DO i=1,num_flow(ii)
          DO ig=1,ngrid(ii,i)
             ir=grid_row(ii,i,ig)
             ic=grid_col(ii,i,ig)
             IF(area(ic,ir).EQ.-9999.0) print *,'wrong in grid-area:',ic,ir,area(ic,ir)
             basin_area=basin_area+area(ic,ir)
             area_sub(ii)=area_sub(ii)+area(ic,ir)
          END DO
        END DO
      END IF
      IF(JustStart == 1) RETURN




!****************************************************************
!                       start simulation
!****************************************************************

      DO 999 i=1,num_flow(ii)
        qin(i)	  =  0.0            ! total lateral inflow (m3/s/m)
        qin_tmp   =  0.0            ! lateral inflow from present flow-interval
        DO 600 ig=1,ngrid(ii,i)
          iy=grid_row(ii,i,ig)
          ix=grid_col(ii,i,ig)

          qin_inter =  0.0
          qin_g     =  0.0
          qin_sfc   =  0.0

!**************************************************
!         for water body
!**************************************************

!          RIVERNET(ix,iy)=0
          IF(RIVERNET(ix,iy).EQ. 1 .AND. waterflx(ix,iy) .GT. 0.1E-10) THEN          ! = prcp - evap on water surface (m)
            PRINT*,"RUNOFF RIVER"
            qin_sfc   = waterflx(ix,iy)*area(ix,iy)/dt   ! m3/s  for waterbody
            qin_inter = 0.0
            qin_g     = 0.0
            qin_tmp   = qin_tmp + (qin_sfc + qin_inter + qin_g)		! m3/s

!**************************************************
!			for	other land use
!**************************************************
          ELSE
!	 1. subsurface runoff
            qin_inter =  runoff_inter(ix,iy) 			!m/s
            qin_g     =  runoff_g(ix,iy)		        !m/s

            
!        2. surface runoff

            q_surf=amax1(0.0,  runoff(ix,iy))!in m, surface runoff

!        hillsope routing: steady constant sheet flow
            q_hillslope=0.0
            IF(q_surf .LE. 0.1E-8) GOTO 500	    !note: important 2007.01.07
            water_depth = q_surf				! in meter
            runoff(ix,iy) =  runoff(ix,iy) - water_depth
!need change
            roughness_n = 0.15		!calibrate

            IF( len2d(ix,iy).LE. 0.0)print *, 'wrong len2d:', ix,iy,len2d(ix,iy)
              s0_slope    = TAN(slp2d(ix,iy)) +water_depth/len2d(ix,iy)
            
            IF(s0_slope .LT. 0.1E-8) s0_slope = 0.1E-8
            power=1.6667
!           q_hillslope: m^3/m
            q_hillslope = dt*sqrt(s0_slope)*   water_depth**power/roughness_n

            q_hillslope = amin1(q_hillslope, water_depth*len2d(ix,iy))
            water_depth = water_depth - q_hillslope/len2d(ix,iy)
            runoff(ix,iy) =  runoff(ix,iy) + water_depth		!update surface storage
           water_depth = 0.0

!            IF (water_depth .GT. 0.0)
!     >          POND(ix,iy) =  POND(ix,iy) + water_depth       !update surface storage
!            IF (POND(ix,iy).GT.PONDMX)
!     >          POND(ix,iy)=PONDMX
!            q_hillslope=q_hillslope+(PONDMX-POND(ix,iy))*len2d(ix,iy)
!            water_depth = 0.0
!	surface runoff
           qin_sfc   = q_hillslope/dt		!m^3/m/s


500        qin_tmp = qin_tmp +(qin_sfc+ qin_inter + qin_g)* area(ix,iy)/len2d(ix,iy)		! m3/s单宽流量 
           totsurflow = totsurflow+qin_sfc*area(ix,iy)/len2d(ix,iy)
           totintflow = totintflow+qin_inter*area(ix,iy)/len2d(ix,iy)
           totgflow   = totgflow  +qin_g*area(ix,iy)/len2d(ix,iy)
!---------------------------------------------------------------------------------------------------
!       I think this should be
! 500         qin_tmp = qin_tmp +
!     >          (qin_inter + qin_g)* area(ix,iy)+ qin_sfc*area(ix,iy)/len2d(ix,iy)	! m3/s单宽流量 
!--------------------------------------------------------------------------------------------------
          END IF
600     CONTINUE					!do 600 ig=1,ngrid(ii,i)
!
!	----------------------------------------------------------------------
        qin(i) = qin_tmp / rdx(ii,i)			! m^3/m/s, total lateral inflow
 999  CONTINUE
      RETURN
  END SUBROUTINE SURFACE_ROUTING








  SUBROUTINE RIVER_ROUTING(dthydro,INITAL,JustStart,inicon,ii,level,l1,l2,l3,JULIAN,HOUR,YEAR,YREND,JEND,&
	HourEND,simulation_dir,totsurflow,totintflow,totgflow)
!**********************************************************************
!                                                                     *
!            ***** Kinematic River Routing Model *****                *
!                                                                     *
!**********************************************************************
!     Variables:
!
!     dt_couple: time step for coupling (s), from 'globcst.inc'
!     qin:       lateral inflow of one flow interval (m3/sec/m), from 'hydro.inc'
!
!     beta:      Manning's equation parameter
!     dt:        time step used in routing (s)
!     q1:        discharge of last time step (m^3/s)
!     q2:        discharge of current time step (m^3/s)
!     qlin1:     lateral inflow of last time step (m^3/m/s)
!     qlin2:     lateral inflow of current time step (m^3/m/s)
!     level:     Pfafstetter level
!     l1-3:      Pfafstetter number of local level (1-9)
!     p1:        basin name (number of whole system) of nine basin of level 1
!     p2:        basin name (number of whole system) of nine basin of level 2
!     p3:        basin name (number of whole system) of nine basin of level 3
!
!     Qd:	    daily average discharge (m^3/s)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
    implicit none
!      INCLUDE 'dims.inc'
!      INCLUDE 'hydro.inc'
    real(r8),intent(in)::dthydro
    integer(i4),intent(in)::INITAL,JustStart
    integer(i4),intent(in)::inicon,ii,level,l1,l2,l3,JULIAN,HOUR,YEAR,YREND,JEND,HourEND
    character(80),intent(in):: simulation_dir
    real(r8),intent(in)::totsurflow,totintflow,totgflow
! temp
    integer(i4)::i!indicator for subcatchment and flow interval respectivelysubroutine STRLNTH
    real(r8)::dt
    integer(i4)::ldir_d, ldir_s, ldir_r, lfile
    integer(i4)::smlunit, wsunit, obsunit, relunit
    integer(i4)::iyy,imm,idd,ihh
    integer(i4):: j, itmp, k , mtime,ii1,ii2
    real(r8)::criterion, h1, tmp, f,df,h2,h,q1_upper,q2_upper,rs,value
    real(r8)::q2tmp
    integer(i4)::outflag

    dt =  dthydro          ! in second


!-----initialize model


      IF (JustStart .EQ. 1) THEN
        DO i=1,num_flow(ii)
          qlin1(ii,i) = 0!lateral inflow last time step
        END DO

!      specify initial condition
!      integer inicon   ! = 0, Give arbitrary value  ! = 1, Import from data file

!     Read the initial discharge in the river
        IF(inicon .EQ. 1) THEN
	  ldir_s = LEN(simulation_dir)
	  CALL STRLNTH(simulation_dir, ldir_s)
          IF(ii == 1)print *, 'reading initial condition from simulation_dir'
          CALL getunit(smlunit)
          OPEN(smlunit,file=simulation_dir(1:ldir_s)//sub_catch(ii)//'I_flow',status='old')
          READ (smlunit,*)(itmp,q1(ii,i), i=1,num_flow(ii))
          CLOSE(smlunit)
	  CALL retunit(smlunit)
!     calculate initial river discharge
        ELSE
          q1(ii,1)=0.1  !initial waterflow in river channel of flow intervals
          DO i=2,num_flow(ii)
            q1(ii,i)=q1(ii,i-1)+0.01			!initial flow
          END DO

        END IF

!     calculate initial river water depth from discharge for each flow interval
        DO i=1,num_flow(ii)
          criterion=0.01
          k=1
          h1=q1(ii,i)/riverb(ii,i)
!         b: river width of the flow interval (m)
 5        tmp=roughness(ii,i)*q1(ii,i)/sqrt(s0(ii,i))
          f=riverb(ii,i)*h1-(tmp**0.6)*((riverb(ii,i)+2.0*h1)**0.4)
          IF (k .GT. 1 .AND. abs(f) .LT. criterion) GOTO 8
          df=riverb(ii,i)-0.8*((tmp/(riverb(ii,i)+2.0*h1))**0.6)
          h2=h1-f/df
          h1=h2
          IF(k .GE. 10) GOTO 8
          k=k+1
          GOTO 5
 8        h=h2
          hhr(ii,i)=h      ! depth of river water
        ENDDO
      END IF
      IF(JustStart .EQ. 1) RETURN		!Initialization finished


!#######################################################################
!
!     lateral inflow
!
!#######################################################################
!     Calculate the discharge nees the lateral flow in the last step
      IF (INITAL .EQ. 0 ) THEN
        DO i=1,num_flow(ii)
          qlin1(ii,i)=qin(i)
        END DO
      END IF

      DO i=1,num_flow(ii)
         qlin2(ii,i)=qin(i)
      END DO
!
!#######################################################################
!
!     River routing

      mtime = 1!int(dt_couple/dt)
      DO 555 j = 1, mtime        !time loop
        DO 333 i=1,num_flow(ii)  !river segment loop

!         junction boundary condition: define the river network
          q1_upper=0.
          q2_upper=0.
!         Calculate discharge inflow from the upper interflow
          IF(i .EQ. 1) THEN
            IF(level .EQ. 1) THEN
            !for the first interflow in head water subcatchments
              IF(l1 .EQ. 9 .OR. mod(l1,2).EQ.0) THEN
                q1_upper=0.
                q2_upper=0.
              ELSE
                ii1=p1(l1+1)
                ii2=p1(l1+2)
                q1_upper=q1(ii1,num_flow(ii1))+q1(ii2,num_flow(ii2))
                q2_upper=q2(ii1,num_flow(ii1))+q2(ii2,num_flow(ii2))
              ENDIF
            ELSEIF(level .EQ. 2) THEN
              IF(l2 .EQ. 9 .or. mod(l2,2).EQ.0) THEN
                q1_upper=0.
                q2_upper=0.

                IF(l2.EQ.9 .AND. (l1.NE.9 .AND. mod(l1,2).NE.0)) THEN !taken from upper level
                  ii1=p1(l1+1)
                  ii2=p1(l1+2)
                  q1_upper=q1(ii1,num_flow(ii1))+q1(ii2,num_flow(ii2))
                  q2_upper=q2(ii1,num_flow(ii1))+q2(ii2,num_flow(ii2))
                ENDIF
              ELSE
                ii1=p2(l2+1)
                ii2=p2(l2+2)
                q1_upper=q1(ii1,num_flow(ii1))+q1(ii2,num_flow(ii2))
                q2_upper=q2(ii1,num_flow(ii1))+q2(ii2,num_flow(ii2))
              ENDIF
            ELSEIF(level .EQ. 3) THEN
              IF(l3 .EQ. 9 .OR. mod(l3,2) .EQ. 0) THEN
                q1_upper=0.
                q2_upper=0.

                IF(l3.EQ.9 .AND. (l2.NE.9 .AND. mod(l2,2).NE.0)) THEN !taken from upper level
                  ii1=p2(l2+1)
                  ii2=p2(l2+2)
                  q1_upper=q1(ii1,num_flow(ii1))+q1(ii2,num_flow(ii2))
                  q2_upper=q2(ii1,num_flow(ii1))+q2(ii2,num_flow(ii2))
                ENDIF
              ELSE
                ii1=p3(l3+1)
                ii2=p3(l3+2)
                q1_upper=q1(ii1,num_flow(ii1))+q1(ii2,num_flow(ii2))
                q2_upper=q2(ii1,num_flow(ii1))+q2(ii2,num_flow(ii2))
              ENDIF
            ENDIF
          ELSE
            q1_upper=q1(ii,i-1)
            q2_upper=q2(ii,i-1)
          ENDIF
!-----end of river network definition

!
          rs=roughness(ii,i)
111       call nkws(dt,rdx(ii,i),riverb(ii,i),s0(ii,i),rs,qlin1(ii,i),qlin2(ii,i),q1_upper,q1(ii,i),q2_upper,q2(ii,i),hhr(ii,i))

333     CONTINUE

!     Update states for the next step
        DO i=1,num_flow(ii)
          qlin1(ii,i)=qlin2(ii,i)
          q1(ii,i)=q2(ii,i)
        ENDDO
 555  CONTINUE

!     OUTPUT THE DISCHARGE AT OUTLET EACH HOUR
      IF(INITAL .EQ. 0) THEN
        IF(sub_catch(ii) .EQ.'ws100' ) THEN          !Outlet
          ldir_s = LEN(simulation_dir)
          CALL strlnth(simulation_dir, ldir_s)  
          OPEN(51,file=simulation_dir(1:ldir_s)//sub_catch(ii)//'.discharge', status='unknown')
          WRITE(51,"(3I5,4F15.7)") YEAR,JULIAN,HOUR,q2(ii,num_flow(ii)),totsurflow,totintflow,totgflow
        END IF
      END IF

      IF(sub_catch(ii) .EQ. 'ws100' ) THEN          !Outlet
         WRITE(51,"(3I5,4F15.7)")YEAR,JULIAN,HOUR,q2(ii,num_flow(ii)),totsurflow,totintflow,totgflow
      END IF

      IF(sub_catch(ii) .EQ. 'ws100' ) THEN          !Outlet
        IF(YEAR .EQ. YREND  .AND.  JULIAN .EQ. JEND.AND. HOUR+1 .GT. HourEND)THEN
          CLOSE(51)
        ENDIF
      END IF
  END SUBROUTINE RIVER_ROUTING


!******************************************************************************
!     Sub-program of Nonlinear Kinematic Wave Scheme solving of river routing
!     Calculate water depth using Newton's method
!
!     Definition of Variales:
!     q0 -- discharge of last time step;
!     q -- discharge of current time step;
!     qlin0 -- lateral inflow of last time step;   qlin1
!     qlin -- lateral inflow of current time step; qlin2
!     h -- water depth;
!     p -- wetted perimeter;
!     b -- river width;
!     s0 -- river slope
!     dx -- Flow interval width
!     tempq01 --q1_upper
!     tempq02 --q1
!     tempq1 -- q2_upper
!     qtr -- q2
!
!     Input:
!     	dt,dx,b,s0,rn,qlin1,qlin2,q1,q1_upper,q2_upper
!     Output:
!     	qtr(q2),h
  SUBROUTINE NKWS(dt,dx,b,s0,rn,qlin0,qlin,tmpq01,tmpq02,tmpq1,qtr,h)
    use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
    implicit none
!   input
    real(r8),intent(in)::dt,dx,b,s0,rn,qlin0,qlin,tmpq01,tmpq02,tmpq1
    real(r8),intent(inout)::qtr,h
!   temp
    real(r8),dimension(2)::q0,q
    real(r8)::roughness,beta
    real(r8)::aa,alfa,bb,cc,criterion,ctmp,df,f,h1,h2,p,qq1,qq2,tmp
    integer(i4)::k

      roughness=rn
      beta=0.6
      q0(1)=tmpq01!q1_upper
      q0(2)=tmpq02!q1
      q(1)=tmpq1	!q2_upper
      criterion=0.00001

!     Calculate the water depth from discharge at the beginning of the simulation
      k=1
      h1=q0(2)/b
15    tmp=roughness*q0(2)/sqrt(s0)
      f=b*h1-(tmp**0.6)*((b+2.0*h1)**0.4)
      IF(k .GT. 1 .AND. abs(f) .LT. criterion) GOTO 18
      df=b-0.8*((tmp/(b+2.0*h1))**0.6)
      h2=h1-f/df
      h1=h2
      IF(k .GE. 20) GOTO 18
      k=k+1
      GOTO 15
18    h=h2
      p=b+2.0*h

!     the initial discharge estimated using linear scheme
      alfa=(roughness*p**(2.0/3.0)/sqrt(s0))**0.6
!
      IF((q0(2)+q(1)) .LE. 0.0) THEN
        cc=0.0
        GOTO 22
      ENDIF
      cc=(0.5*(q0(2)+q(1)))**(beta-1.0)
22    aa=dt*q(1)/dx+alfa*beta*q0(2)*cc+0.5*(qlin0+qlin)*dt
      bb=dt/dx+alfa*beta*cc
      qq1=aa/bb
!
      if(qq1 .le. 0.1e-6) then
        qq2=0.0
        goto 30
      endif
!
!     Using Newton's method to calculate discharge
      criterion=0.00005
      ctmp=dt*q(1)/dx+alfa*q0(2)**beta+0.5*dt*(qlin+qlin0)
      k=1
20    f=dt*qq1/dx+alfa*qq1**beta-ctmp
      IF((k .GT. 1 ).AND.(abs(f).LE.criterion)) GOTO 30
      df=dt/dx+alfa*beta*qq1**(beta-1.0)
      qq2=qq1-f/df
!
      IF(qq2 .LE. 0.1e-6) THEN
        qq2=0.0
        GOTO 30
      ENDIF
!
      qq1=qq2
      IF(k.GE.20) GOTO 30
      k=k+1
      GOTO 20
30    q(2)=qq2
      IF(q(2) .LE. 0.1e-6) q(2)=0.0

      qtr=q(2)

!     Calculate the water depth according to the discharge at the present time step
      criterion=0.00001
      k=1
      h1=q(2)/b
45    tmp=roughness*q(2)/sqrt(s0)
      f=b*h1-(tmp**0.6)*((b+2.0*h1)**0.4)
      IF(k.GT.1 .AND. abs(f) .LT. criterion) GOTO 50
      df=b-0.8*((tmp/(b+2.0*h1))**0.6)
      h2=h1-f/df
      h1=h2
      IF(k.GE.20) GOTO 50
      k=k+1
      GOTO 45
50    h=h2     
  END SUBROUTINE NKWS



  SUBROUTINE GWRIV(Dtg,length,slope,Ds,Dg,Dr,Drw,kgs,Q)
!
!     The datum is sitted at the bottom of the unconfined aquifer
!     Varibles:
!     Dtg:      depth to groundwater (m)
!     length:   length of hillslope (m)
!     slope:    slope of hillslope (m)
!     Ds:       depth of top soil (m)
!     Dg:       depth of unconfined groundwater acquifer below topsoil (m)
!     Dr:       depth of river (m)
!     Drw:      depth of river water (m)
!     kgs:      hydraulic conductivity (m/sec)
!     conlen:   contact length between the river and aquifer (m)
!     grad:     gradient of water head
!     Q:        discharge exchanged between aquifer and river (m^3/sec)
!
    use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
    implicit none
    real(r8),intent(in)::Dtg,length,slope,Ds,Dg,Dr,kgs
    real(r8),intent(inout)::Drw,Q
    real(r8):: H1,H2,hs1,hs2,grad, conlen, Hrd

    Drw = amax1(Drw, 0.0)
    Hrd = Ds+Dg-Dr                     ! distance from datum to riverbed (m)
	if(Hrd.ge.0.0) then
        H1 =0.5*length*slope+Ds+Dg-Dtg   ! waterhead of groundwater (m)
        hs1=H1-Hrd                       ! saturated acquifer depth (m)
        H2 =Hrd+Drw                      ! waterhead of river (m)
        hs2=Drw                          ! water depth in river (m)
	  grad   =(H1-H2)/(0.5*length)     ! gradient of waterhead
	  conlen =0.5*(abs(hs1)+hs2)       ! contact length between the river
                                         ! and aquifer	(m)
	else
        H1 =0.5*length*slope+Ds+Dg-Dtg
        hs1=Ds+Dg-Dtg
        H2 =amax1(Hrd+Drw, 0.0)
	  hs2=H2
	  grad   =(H1-H2)/(0.5*length)
	  conlen =0.5*(abs(hs1)+hs2)
	endif

        Q=kgs*grad*conlen			  ! discharge per unit width of a hillslope
                                    ! (m3/sec/m)
	if(Drw.le.5E-3 .and. Q.lt.0.0) Q=0.0		! 2010.4.26 ge -> le
	return
  end SUBROUTINE GWRIV










end module hydro_mod
