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
           CALL lateral_inflow(JustStart,ii,AREA,SLP2D,LEN2D,runoff_inter,runoff_g,waterflx,dthydro,&
		RIVERNET,RUNOFF,POND,PONDMX,totsurflow,totintflow,totgflow)
           CALL RIVER_ROUTING(dthydro,INITAL,JustStart,inicon,ii,1,l1,l2,l3,JULIAN,HOUR,YEAR,YREND,JEND,HourEND,&
		simulation_dir,totsurflow,totintflow,totgflow)
          ENDIF
        ELSE
          DO l2 = 9, 1, -1                                       ! Level 2
            IF(kfs2(l1,l2) .EQ. 0) THEN
              ii=ii+1
	        IF(ii>=1 .AND. ii<=num_sub) THEN
                CALL lateral_inflow(JustStart,ii,AREA,SLP2D,LEN2D,runoff_inter,runoff_g,waterflx,dthydro,&
		RIVERNET,RUNOFF,POND,PONDMX,totsurflow,totintflow,totgflow)
                CALL RIVER_ROUTING(dthydro,INITAL,JustStart,inicon,ii,2,l1,l2,l3,JULIAN,HOUR,YEAR,YREND,JEND,HourEND,&
		simulation_dir,totsurflow,totintflow,totgflow)
                END IF
            ELSE
              DO l3 = 9, 1, -1                                     ! Level 3
                ii=ii+1
	          IF(ii>=1 .AND. ii<=num_sub) THEN
                  CALL lateral_inflow(JustStart,ii,AREA,SLP2D,LEN2D,runoff_inter,runoff_g,waterflx,dthydro,RIVERNET,&
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



























end module hydro_mod
