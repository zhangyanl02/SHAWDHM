module hydro_mod
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8
    use dims_mod
    use input_mod,only:getunit,retunit,strlen
    implicit none

    public
      save
      integer(i4)::nsub                                      ! total number of subcatchment
      integer(i4)::nflowmax                                  ! max number of flow interval for all subcatchment
      integer(i4)::ngridmax                                  ! max bumber of grid for all flow intervals
      character*6,allocatable:: subbasin(:)                  ! name of each sub catchment
      integer(i4),allocatable:: psubbasin(:)                 ! the ids of sub catchments,e.g., 1001,2003
      integer(i4),allocatable:: pbasinup(:,:)                ! the ids of up-basins for each sub catchment
      integer(i4),allocatable:: nbasinup(:,:)                ! the rank of each up-basin of each sub catchment, e.g., 1,2,3,4
      integer(i4),allocatable:: nextbasin(:)                 ! the next subcatchment of each subcatchment
      integer(i4),allocatable:: nextbasinrow(:)              ! the row number of the grid of next basin
      integer(i4),allocatable:: nextbasincol(:)              ! the col number of the grid of next basin
      integer(i4),allocatable::nflow(:)                      ! the flow interval number of each sub-catchment
      integer(i4),allocatable::ngrid(:,:)                    ! the number of grid in each flow interval                    ----ngrid(nsub,nflow)
      integer(i4),allocatable::grid_row(:,:,:)               ! the row number of each grid                                 ----grid_row(nsub,nflow,ngrid)
      integer(i4),allocatable::grid_col(:,:,:)               ! the col number of each grid                                 ----grid_col(nsub,nflow,ngrid)
      real(r8),allocatable::dx(:,:)                          ! (rdx) length of flow intervals (m)                          ----dx(nsub,nflow)
      real(r8),allocatable::s0(:,:)                          ! river-bed slope of the flow interval                        ----s0(nsub,nflow)
      real(r8),allocatable::b(:,:)                           ! (riverb) river width of the flow interval (m)               ----b(nsub,nflow)
      real(r8),allocatable::roughness(:,:)                   ! river roughness (manning's coefficient) of the flow interval----(nsub,nflow)
      real(r8),allocatable::dr(:,:)                          ! (drr) river depth of the flow interval (m)                  ----dr(nsub,nflow)
      real(r8),allocatable::hrr(:,:)                         ! depth of river water                                        ----hrr(nsub,nflow)
      real(r8),allocatable::qin(:)                           ! lateral inflow of one flow interval (m3/sec/m)              ----qin(nflow)
      real(r8),allocatable::q1(:,:)                          ! discharge of last time step (m^3/s)                         ----q1(nsub,nflow)
      real(r8),allocatable::q2(:,:)                          ! discharge of current time step (m^3/s)                      ----q2(nsub,nflow)
      real(r8),allocatable::qlin1(:,:)                       ! lateral inflow of last time step (m^3/m/s)                  ----qlin1(nsub,nflow)
      real(r8),allocatable::qlin2(:,:)                         ! lateral inflow of current time step (m^3/m/s)               ----qlin2(nsub,nflow)
      real(r8),allocatable::hhr(:,:)                         ! depth of river water                                        ----hhr(nsub,nflow)
      integer(i4),parameter::sceflag = 0                     ! =1 for sce optimization; = 0 for no sce optimization. 
      real(r8):: basin_area
      real(r8),allocatable::area_sub(:)                      ! basin-area, subcatchment area
      real(r8)::dthydro                                      ! the length of time interval, in sencond(s),need to be initianized at the beginning

  contains
    subroutine init_hydro()
    
    end subroutine
    
    subroutine read_hydro_para(para_dir,nhrpdt)
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4      
      implicit none
      character*200,intent(in)::para_dir
      integer(i4),intent(in)::nhrpdt
      
      dthydro=nhrpdt*3600.0
      call initial_subcatchment(para_dir)
      call read_river_para(para_dir)
      return 
    end subroutine read_hydro_para


    subroutine initial_subcatchment(para_dir)
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
      implicit none
      character*200,intent(in)::para_dir
      integer::i,j,k,fileunit
      character(200)::atmp
      integer(i4)::l1,l2,subcount,iostatus
      character*4::ch4 
      character*200::subbasinfile


      call strlen(para_dir,l1,l2)
      subbasinfile=trim(para_dir(l1:l2))//'subbasin.dat'
      call getunit(fileunit)
      open(fileunit,file = subbasinfile, status='old')
      subcount=0
      iostatus=1
      do while(iostatus.ge.0)        
        read(fileunit,*,iostat=iostatus) atmp
        subcount=subcount+1
      end do
      close(fileunit)
      subcount=subcount-1
      nsub=subcount

      call strlen(para_dir, l1, l2)
      open(fileunit,file=trim(para_dir(l1:l2))//'riverpara/'//'maxnp',status='old')
      read(fileunit,*) nflowmax,ngridmax
      close(fileunit)

      allocate(subbasin(nsub))
      allocate(psubbasin(nsub))
      allocate(pbasinup(nsub,8))
      allocate(nbasinup(nsub,8))
      !allocate(nextbasinrow(nsub),nextbasincol(nsub),nextbasin(nsub))
      allocate(nflow(nsub))
      allocate(ngrid(nsub,nflowmax))
      allocate(dx(nsub,nflowmax))
      allocate(s0(nsub,nflowmax))
      allocate(b(nsub,nflowmax))
      allocate(roughness(nsub,nflowmax))
      allocate(Dr(nsub,nflowmax))
      allocate(grid_row(nsub,nflowmax,ngridmax))
      allocate(grid_col(nsub,nflowmax,ngridmax))
      allocate(hrr(nsub,nflowmax))
      allocate(q1(nsub,nflowmax))
      allocate(q2(nsub,nflowmax))      
      allocate(qlin1(nsub,nflowmax))
      allocate(qlin2(nsub,nflowmax))
      allocate(hhr(nsub,nflowmax))
      allocate(qin(nsub))
      allocate(area_sub(nsub))
  
      open(fileunit,file = subbasinfile, status='old')
      do i = 1,nsub
        read(fileunit,*)nsub,psubbasin(i),(pbasinup(i,j),j=1,8)!,nextbasinrow(i),nextbasincol(i),nextbasin(i)
      enddo
      close(fileunit)


      do i =1,nsub
        write(ch4,'(i4.4)') psubbasin(i)
        subbasin(i)='ws'//ch4
      enddo
      
      nbasinup=0     
      do i = 1,nsub
        do k =1,8
          if(pbasinup(i,k)>0) then
            do j =1,nsub
              if(psubbasin(j)== pbasinup(i,k)) then
                nbasinup(i,k)=j 
              endif
            enddo
          endif
        enddo
      enddo  
      call retunit(fileunit)
    end subroutine

  
    subroutine read_river_para(para_dir)
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
      implicit none
      character*200,intent(in)::para_dir
      integer(i4)::isub,l1,l2,ll1,ll2,iflow,j,fileunit
      character*200::infile
      call getunit(fileunit)
      do isub = 1, nsub
          infile = subbasin(isub)//'_river'
        !print *,isub,subbasin(isub),"  ",infile
          call strlen(infile, l1, l2)
          call strlen(para_dir, ll1, ll2)
          open(fileunit,file=trim(para_dir(ll1:ll2))//'riverpara/'//infile(l1:l2),status='old')
          read(fileunit,*) nflow(isub)
          do iflow = 1, nflow(isub)
            read(fileunit,*) ngrid(isub,iflow),dx(isub,iflow),s0(isub,iflow),b(isub,iflow),roughness(isub,iflow),Dr(isub,iflow)
            if(dx(isub,iflow) .lt. 0.1) dx(isub,iflow) = 2000.0
            if(s0(isub,iflow) .le. 0 ) then   
              s0(isub,iflow)=0.00001
              print *, 'wrong in s0',isub,iflow,s0(isub,iflow)  
            endif  
            if(s0(isub,iflow).eq.-9999.0) then 
              s0(isub,iflow)=0.00001 
              print *, 'wrong in s0',isub,iflow,s0(isub,iflow)  
            endif
            if(s0(isub,iflow).le.0.1E-5) s0(isub,iflow) = 0.1E-5

            if(ngrid(isub,iflow) .gt. ngridmax) then
              print *,'more than np grids:',ngrid(isub,iflow),isub,iflow,ngridmax
            endif
            read(fileunit,*) (grid_row(isub,iflow,j),grid_col(isub,iflow,j),j = 1, ngrid(isub,iflow))
          end do
          close(fileunit)
          call retunit(fileunit)
      end do
    end subroutine


    subroutine lateral_inflow(nx,ny,juststart,isub,area2d,slp2d,len2d,runoff_inter,runoff_g,waterflx,&
      rivernet,runoff,pond,pondmx,totsurflow,totintflow,totgflow)  !17 arguments
!**********************************************************************
!  Function:
!     Calculating the lateral inflow of each flow interval in a subcatchment based on the runoff_inter,runoff_g,waterflx 
!  Input:
!     runoff_inter,runoff_g,waterflx
!  OUT:
!           qin(isub)
!**********************************************************************
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
      implicit none

      integer(i4),intent(in)::nx,ny
      integer(i4),intent(in)::juststart                            !
      integer(i4),intent(in)::isub                                 !  indicator of subcatchment loop
      real(r8),dimension(nx,ny),intent(in)::area2d                 !  the area of grid cell
      real(r8),dimension(nx,ny),intent(in)::slp2d                  !  the slope angle of hillslope
      real(r8),dimension(nx,ny),intent(inout)::len2d                  !  the slope lenght of hillslope
      real(r8),dimension(nx,ny),intent(in)::runoff_inter           !  the interflow in each grid cell, in m/s
      real(r8),dimension(nx,ny),intent(in)::runoff_g               !  the base flow in each grid cell, in m/s
      real(r8),dimension(nx,ny),intent(in)::waterflx               !  the surface runoff in each grid cell,in m/s
      integer(i4),dimension(nx,ny),intent(in)::rivernet            !  indicator of river cell
      real(r8),dimension(nx,ny),intent(inout)::runoff              !  
      real(r8),dimension(nx,ny),intent(inout)::pond                !  the depth of suface water ponding
      real(r8),intent(in)::pondmx                                  !  max ponding depth of surface water
      real(r8),intent(inout)::totsurflow                           !  total surface runoff
      real(r8),intent(inout)::totintflow                           !  total interflow
      real(r8),intent(inout)::totgflow                             !  total base flow


      integer(i4):: iflow, ig,ic,ir,ix,iy
      real(r8):: dt                   !time step in second (s)
      real(r8):: qin_inter            !Subsurface flow In each grid
      real(r8):: qin_g                !exchanges between groundwater and river In each grid
      real(r8):: qin_sfc              !surface runoff from a flow interval In each grid
      real(r8):: qin_tmp              !totoal runoff from a flow interval
      real(r8):: power,q_hillslope,q_surf,roughness_n,s0_slope,water_depth


      dt = dthydro 

      if(JustStart == 1) then
        if(isub .EQ. 1) then
          basin_area = 0.0
        endif
        area_sub(isub)=0.0
        do iflow=1,nflow(isub)
          do ig=1,ngrid(isub,iflow)
             ir=grid_row(isub,iflow,ig)
             ic=grid_col(isub,iflow,ig)
             if(area2d(ic,ir).EQ.-9999.0) then
               print *,'wrong in grid-area:',ic,ir,area2d(ic,ir)
               stop
             end if
             if(len2d(ic,ir).LE. 0.0) then
               print *,'wrong in len2d:', ic,ir,len2d(ic,ir)
               len2d(ic,ir)=100.0
             end if
             basin_area=basin_area+area2d(ic,ir)
             area_sub(isub)=area_sub(isub)+area2d(ic,ir)
          enddo
        enddo
        return
      endif


      do iflow=1,nflow(isub)        ! flow interval loop
        qin(iflow)=  0.0            ! total lateral inflow (m3/s/m)
        qin_tmp   =  0.0            ! lateral inflow from present flow-interval
        DO ig=1,ngrid(isub,iflow) !grid loop within a flow interval
          iy=grid_row(isub,iflow,ig)
          ix=grid_col(isub,iflow,ig)
          qin_inter =  0.0
          qin_g     =  0.0
          qin_sfc   =  0.0
          
          !for water body
          if (rivernet(ix,iy).EQ. 1 .and. waterflx(ix,iy) .GT. 0.1E-10) then          ! = prcp - evap on water surface (m)
            qin_sfc   = waterflx(ix,iy)*area2d(ix,iy)/dt   ! m3/s  for waterbody
            qin_inter = 0.0
            qin_g     = 0.0
            qin_tmp   = qin_tmp + (qin_sfc + qin_inter + qin_g)            ! m3/s
         !for other land use
          else
            qin_inter =  runoff_inter(ix,iy)                   !in m/s
            qin_g     =  runoff_g(ix,iy)                        !in m/s 
            q_surf=amax1(0.0,  runoff(ix,iy))               !in m, surface runoff

            !hillsope routing: steady constant sheet flow
            q_hillslope=0.0
            if(q_surf .gt. 0.1E-8) then
              water_depth = q_surf                        ! in meter
              runoff(ix,iy) =  runoff(ix,iy) - water_depth
              roughness_n = 0.2            !need change, calibrate
              s0_slope    = TAN(slp2d(ix,iy)) +water_depth/len2d(ix,iy)            
              if(s0_slope .LT. 0.1E-8) s0_slope = 0.1E-8
              power=1.6667
              q_hillslope = dt*sqrt(s0_slope)*water_depth**power/roughness_n  !in m^3/m
              q_hillslope = amin1(q_hillslope, water_depth*len2d(ix,iy))
              water_depth = water_depth - q_hillslope/len2d(ix,iy)
              runoff(ix,iy) =  runoff(ix,iy) + water_depth            !update surface storage
              water_depth = 0.0
              qin_sfc   = q_hillslope/dt            !m^3/m/s
            end if
            qin_tmp = qin_tmp +(qin_sfc+ qin_inter + qin_g)* area2d(ix,iy)/len2d(ix,iy)            ! m3/s单宽流量 
            totsurflow = totsurflow+qin_sfc*area2d(ix,iy)/len2d(ix,iy)
            totintflow = totintflow+qin_inter*area2d(ix,iy)/len2d(ix,iy)
            totgflow   = totgflow  +qin_g*area2d(ix,iy)/len2d(ix,iy)
          end if
        enddo
        qin(iflow) = qin_tmp / dx(isub,iflow)                  ! m^3/m/s, total lateral inflow
      end do
      return
    end subroutine lateral_inflow


    subroutine river_routing(initial,juststart,inicon,isub,julian,hour,year,yrend,jend,&
      hourend,simulation_dir,totsurflow,totintflow,totgflow)  !14 arguments
!**********************************************************************
!            ***** Kinematic River Routing Model *****                *
!**********************************************************************
!     Variables:
!
!     dt_couple: time step for coupling (s), from 'globcst.inc'
!     qin:       lateral inflow of one flow interval (m3/sec/m), from 'hydro.inc'
!     beta:      Manning's equation parameter
!     dt:        time step used in routing (s)
!     q1:        discharge of last time step (m^3/s)
!     q2:        discharge of current time step (m^3/s)
!     qlin1:     lateral inflow of last time step (m^3/m/s)
!     qlin2:     lateral inflow of current time step (m^3/m/s)
!     Qd:          daily average discharge (m^3/s)
!**********************************************************************
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
      implicit none
      integer(i4),intent(in)::initial,JustStart
      integer(i4),intent(in)::inicon,isub,julian,hour,year,yrend,jend,hourend
      character(200),intent(in):: simulation_dir
      real(r8),intent(in)::totsurflow,totintflow,totgflow

      integer(i4)::iflow!indicator for subcatchment and flow interval respectivelysubroutine STRLNTH
      real(r8)::dt
      integer(i4)::smlunit
      integer(i4):: j, itmp, k , mtime,ii1,ii2
      real(r8)::criterion, tmp, f,df,h2,h1,h,rs
      real(r8)::q1_upper     ! the flow-in discharge from the upper flow interval at the last time step
      real(r8)::q2_upper     ! the flow-in discharge from the upper flow interval at the current time step
      real(r8)::b_tmp,rs_tmp
      integer(i4)::l1,l2

      dt =  dthydro          ! in second        
      if (JustStart .EQ. 1) then
        do iflow=1,nflow(isub)          
          qlin1(isub,iflow) = 0       !lateral inflow at the last time step
        enddo
!      specify initial condition
!      integer inicon   inicon = 0 Give arbitrary value,  inicon = 1 Import from data file
!      Read the initial discharge in the river
        if(inicon .EQ. 1) then
          call strlen(simulation_dir,l1,l2)
          if(isub == 1)  print *, 'reading initial condition from simulation_dir'
          call getunit(smlunit)
          open (smlunit,file=trim(simulation_dir(l1:l2))//subbasin(isub)//'I_flow',status='old')
          read (smlunit,*)(itmp,q1(isub,iflow), iflow=1,nflow(isub))
          close(smlunit)
          call retunit(smlunit)
        else
          q1(isub,1)=0.1  !initial waterflow in river channel of flow intervals
          do iflow=2,nflow(isub)
            q1(isub,iflow)=q1(isub,iflow-1)+0.01                  !initial flow
          enddo
        end if

!     calculate initial river water depth from discharge for each flow interval
        do iflow=1,nflow(isub)
          criterion=0.01
          k=1
          h1=q1(isub,iflow)/b(isub,iflow) !b: river width of the flow interval (m)
          do k=1,10
            tmp=roughness(isub,iflow)*q1(isub,iflow)/sqrt(s0(isub,iflow))
            f=b(isub,iflow)*h1-(tmp**0.6)*((b(isub,iflow)+2.0*h1)**0.4)
            if (k .GT. 1 .AND. abs(f) .LT. criterion) then
              exit
            else
              df=b(isub,iflow)-0.8*((tmp/(b(isub,iflow)+2.0*h1))**0.6)
              h2=h1-f/df
              h1=h2
            end if
          end do
          h=h2
          hhr(isub,iflow)=h    ! depth of river water
        end do
      end if
      if(JustStart .EQ. 1) return            !Initialization finished


!#######################################################################
!
!     lateral inflow
!
!#######################################################################
!     Calculate the discharge needs the lateral flow in the last step
      if (initial .EQ. 0 ) then
        do iflow=1,nflow(isub)
          qlin1(isub,iflow)=qin(iflow)
        enddo
      endif

      do iflow=1,nflow(isub)
         qlin2(isub,iflow)=qin(iflow)
      enddo

!     River routing
      mtime = 1   !int(dt_couple/dt)
      do j = 1, mtime        !time loop
        do iflow=1,nflow(isub)  !river segment loop
!         junction boundary condition: define the river network
          q1_upper=0.
          q2_upper=0.
!         Calculate discharge inflow from the upper interflow
          if(iflow .EQ. 1) then
            if(psubbasin(isub) > 1000 .and. psubbasin(isub) < 2000) then
              q1_upper=0.
              q2_upper=0.
            else
              do ii1=1,8
                if(nbasinup(isub,ii1) > 0) then
                  ii2=nbasinup(isub,ii1)
                  q1_upper=q1_upper+q1(ii2,nflow(ii2))
                  q2_upper=q2_upper+q2(ii2,nflow(ii2))
                endif
              end do
            end if
          else
            q1_upper=q1(isub,iflow-1)
            q2_upper=q2(isub,iflow-1)
          end if
          rs=roughness(isub,iflow)
          b_tmp = b(isub,iflow)
          rs_tmp = rs
          if(hhr(isub,iflow) .gt. 0.5*Dr(isub,iflow)) then
            b_tmp=1.1*b(isub,iflow)
            rs_tmp = 1.5*rs
          endif
          if(hhr(isub,iflow) .gt. 1.0*Dr(isub,iflow)) then
            b_tmp=1.5*b(isub,iflow)
            rs_tmp = 3.0*rs
          endif
          call nkws(dt,dx(isub,iflow),b_tmp,s0(isub,iflow),rs_tmp,qlin1(isub,iflow),qlin2(isub,iflow),&
               q1_upper,q1(isub,iflow),q2_upper,q2(isub,iflow),hhr(isub,iflow))
        end do

!     Update states for the next step
        do iflow=1,nflow(isub)
          qlin1(isub,iflow)=qlin2(isub,iflow)
          q1(isub,iflow)=q2(isub,iflow)
        enddo
      end do

!     output the discharge at outlet each hour
      if(initial .eq. 0) then
        if(isub .eq. nsub) then          !outlet
          call strlen(simulation_dir, l1,l2)  
          open(51,file=trim(simulation_dir(l1:l2))//subbasin(isub)//'.discharge', status='unknown')
          write(51,"(3i5,4f15.7)") year,julian,hour,q2(isub,nflow(isub)),totsurflow,totintflow,totgflow
        end if
      end if

      if(isub .eq. nsub ) then          !outlet
         write(51,"(3i5,4f15.7)")year,julian,hour,q2(isub,nflow(isub)),totsurflow,totintflow,totgflow
      end if

      if(isub .eq. nsub ) then          !outlet
        if(year .eq. yrend  .and.  julian .eq. jend.and. hour+1 .gt. hourend)then
          close(51)
        endif
      end if
    end subroutine river_routing
      
      
    subroutine lateral_routing(nx,ny,juststart,inital,inicon,area,slp2d,len2d,dr2d,drw2d,&
        runoff_inter, runoff_g,waterflx,julian,hour,year,yrend,jend,hourend,simulation_dir,rivernet,&
        runoff,pond,pondmx,totsurflow,totintflow,totgflow)
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
      implicit none
      integer(i4),intent(in)::juststart,inital,inicon
      integer,intent(in)::nx,ny
      real(r8),dimension(nx,ny),intent(in)::area             ! area of each grid (m2)
      real(r8),dimension(nx,ny),intent(in)::slp2d            ! slope of ground surface (nd)
      real(r8),dimension(nx,ny),intent(inout)::len2d            ! hillslope length of each grid (m)
      real(r8),dimension(nx,ny),intent(inout)::dr2d       ! river depth in each grid (m)
      real(r8),dimension(nx,ny),intent(inout)::drw2d       ! water depth in the river of each grid (m)
      real(r8),dimension(nx,ny),intent(inout)::runoff_inter,runoff_g,waterflx,runoff,pond
      real(r8),intent(in)::pondmx
      real(r8),intent(inout)::totsurflow,totintflow,totgflow
      integer(i4),intent(in)::julian,hour,year,yrend,jend,hourend
      character(200),intent(in):: simulation_dir
      integer(i4),dimension(nx,ny),intent(in)::rivernet

      integer(i4):: ig, isub,iflow
      integer(i4):: ix,iy            !column number and row number for each grid in the domain
      do isub = 1,nsub
        call lateral_inflow(nx,ny,juststart,isub,area,slp2d,len2d,runoff_inter,runoff_g,waterflx,&
               rivernet,runoff,pond,pondmx,totsurflow,totintflow,totgflow)
        call river_routing(inital,juststart,inicon,isub,julian,hour,year,yrend,jend,&
               hourend,simulation_dir,totsurflow,totintflow,totgflow)
      enddo
!     after the surface routing and the river routing, update the river water depth
      do isub=1,nsub
        do iflow=1,nflow(isub)
          do ig=1,ngrid(isub,iflow)
            iy=grid_row(isub,iflow,ig)
            ix=grid_col(isub,iflow,ig)
            drw2d(ix,iy)=hhr(isub,iflow)
            dr2d(ix,iy)=dr(isub,iflow)
          enddo
        enddo
      enddo
    end subroutine lateral_routing



!******************************************************************************
!     Sub-program of Nonlinear Kinematic Wave Scheme solving of river routing
!     Calculate water depth using Newton's method
!******************************************************************************
    subroutine NKWS(dt,dx,b,s0,rn,qlin1,qlin2,q1_upper,q1,q2_upper,q2,h)
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
      implicit none
      real(r8),intent(in)::dt             !  time interval
      real(r8),intent(in)::b              !  river width
      real(r8),intent(in)::s0             !  river slope   
      real(r8),intent(in)::dx             !  Flow interval width 
      real(r8),intent(in)::rn             !  roughness           
      real(r8),intent(in)::qlin1          !  lateral inflow of last time step;    qlin1
      real(r8),intent(in)::qlin2           !  lateral inflow of current time step; qlin2
      real(r8),intent(in)::q1_upper         !  q1_upper
      real(r8),intent(in)::q1         !  q1
      real(r8),intent(in)::q2_upper          !  q2_upper
      real(r8),intent(inout)::q2         !  q2
      real(r8),intent(inout)::h           !  water depth
      real(r8)::qqq2          !  discharge of the current time step
      real(r8)::p                         !  wetted perimeter
      real(r8)::roughness,beta
      real(r8)::aa,alfa,bb,cc,criterion,ctmp,df,f,h1,h2,qq1,qq2,tmp
      integer(i4)::k

      roughness=rn
      beta=0.6
      criterion=0.00001

!     Calculate the water depth from discharge at the beginning of the simulation
      k=1
      h1=q1/b
15    tmp=roughness*q1/sqrt(s0)
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
      IF((q1+q2_upper) .LE. 0.0) THEN
        cc=0.0
        GOTO 22
      ENDIF
      cc=(0.5*(q1+q2_upper))**(beta-1.0)
22    aa=dt*q2_upper/dx+alfa*beta*q1*cc+0.5*(qlin1+qlin2)*dt
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
      ctmp=dt*q2_upper/dx+alfa*q1**beta+0.5*dt*(qlin2+qlin1)
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
30    qqq2=qq2
      IF(qqq2 .LE. 0.1e-6) qqq2=0.0

      q2=qqq2

!     Calculate the water depth according to the discharge at the present time step
      criterion=0.00001
      k=1
      h1=qqq2/b
45    tmp=roughness*qqq2/sqrt(s0)
      f=b*h1-(tmp**0.6)*((b+2.0*h1)**0.4)
      IF(k.GT.1 .AND. abs(f) .LT. criterion) GOTO 50
      df=b-0.8*((tmp/(b+2.0*h1))**0.6)
      h2=h1-f/df
      h1=h2
      IF(k.GE.20) GOTO 50
      k=k+1
      GOTO 45
50    h=h2     
  end subroutine nkws



    subroutine gwriv(Dtg,length,slope,Ds,Dg,Dr,Drw,kgs,Q)
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
!     Q:        discharge exchanged between aquifer and river (m^3/sec/m)
      use shr_kind_mod, only:r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
      implicit none
      real(r8),intent(in)::Dtg,length,slope,Ds,Dg,Dr,kgs
      real(r8),intent(inout)::Drw,Q
      real(r8):: h1,h2,hs1,hs2,grad, conlen, hrd

      Drw = amax1(Drw, 0.0)
      Hrd = Ds+Dg-Dr                     ! distance from datum to riverbed (m)
      if(Hrd.ge.0.0) then
        H1 =0.5*length*slope+Ds+Dg-Dtg   ! waterhead of groundwater (m)
        hs1=H1-Hrd                       ! saturated acquifer depth (m)
        H2 =Hrd+Drw                      ! waterhead of river (m)
        hs2=Drw                          ! water depth in river (m)
        grad   =(H1-H2)/(0.5*length)     ! gradient of waterhead
        conlen =0.5*(abs(hs1)+hs2)       ! contact length between the river
                                         ! and aquifer      (m)
      else
        H1 =0.5*length*slope+Ds+Dg-Dtg
        hs1=Ds+Dg-Dtg
        H2 =amax1(Hrd+Drw, 0.0)
        hs2=H2
        grad   =(H1-H2)/(0.5*length)
        conlen =0.5*(abs(hs1)+hs2)
      endif

      Q=kgs*grad*conlen                    ! discharge per unit width of a hillslope (m3/sec/m)
      if(Drw.le.5E-3 .and. Q.lt.0.0) Q=0.0            ! 2010.4.26 ge -> le
      return
    end subroutine gwriv



    subroutine runoff_generation(nx,ny,runoffdis,runoffdepth,dtime,ncount,averunoffdepth)
        use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4 => shr_kind_r4
        use statevar_mod,only:inbasin2d,runoff_g,ns2d,hkdt2d,matdt2d,vlcdt2d,vicdt2d,dgl2d,slopelen2d,slope2d,ds2d,&
                              dg2d,dr2d,drw2d,runoff_inter,zs2d,tsdt2d,runoff12d,runoff2d
        use soilproperty_mod,only:thfc2d
        !use hydro_mod,only:nsub,lateral_inflow,river_routing,read_hydro_para,lateral_routing,gwriv
        use shaw27_mod,only:soilhk
        implicit none
        integer(i4)::nx,ny,ncount
        real(r8)::averunoffdepth
        real(r8)::runoffdis(nx,ny),runoffdepth(nx,ny)
        real(r8)::dtime

        integer(i4)::col,row,i
        real(r8)::tmpqsub,avkg,tempcoe,deltzi,fieldc,tmp

        
        runoffdis(:,:)=0.0
        ncount=0!!added 
        averunoffdepth=0
        !$omp parallel private(col,row,tmpqsub,avkg,tmp,deltzi) default(shared)
        !$omp do
        do col=1,nx
          do row = 1, ny
            if(inbasin2d(col,row).eq. 1)then
              runoffdepth(col,row)=0.0!added
!             calculate groundwater rounoff
              runoff_g(col,row)=0.0
              tmpqsub = 0.0
              ncount=ncount+1
              call soilhk (ns2d(col,row),hkdt2d(col,row,:),matdt2d(col,row,:),vlcdt2d(col,row,:),vicdt2d(col,row,:),col,row)
              avkg = 1.3*hkdt2d(col,row,6)
!             avkg = 0.1*satk(col,row,ns)
              if (dgl2d(col,row).lt.0.0) dgl2d(col,row)=0.0
              call gwriv(dgl2d(col,row),slopelen2d(col,row),tan(slope2d(col,row)),ds2d(col,row),dg2d(col,row), &
              dr2d(col,row),drw2d(col,row),avkg,tmpqsub) ! m3/s/m
              runoff_g(col,row) = tmpqsub * dtime / slopelen2d(col,row)                                   ! m
              !if(abs(runoff_g(col,row)) .lt. 0.1e-20) runoff_g(col,row) = 0.0

              !if(runoff_g(col,row) .lt. 0.0) then      ! river infiltrates into aquifer
              !  if(dgl2d(col,row) .gt. ds2d(col,row)) then
              !    dgl2d(col,row)= dgl2d(col,row)+runoff_g(col,row)/0.3
              !  else 
              !    vlcdt2d(col,row,ns2d(col,row))=vlcdt2d(col,row,ns2d(col,row))-runoff_g(col,row)/(zs2d(col,row,ns2d(col,row))&
              !-zs2d(col,row,ns2d(col,row)-1))/2.0
              !    dgl2d(col,row)=ds2d(col,row)
              !  end if
              !else
              !  dgl2d(col,row)= dgl2d(col,row)+runoff_g(col,row)/0.3
              !end if
              if(dgl2d(col,row).gt.(ds2d(col,row)+dg2d(col,row)))then
                if(runoff_g(col,row).gt. 0.0)then
                   runoff_g(col,row)=runoff_g(col,row)+(ds2d(col,row)+dg2d(col,row)-dgl2d(col,row))*0.3
                end if
                dgl2d(col,row) = ds2d(col,row)+dg2d(col,row)
              end if
              runoff_g(col,row) =runoff_g(col,row)+0.001939/6068* dtime / slopelen2d(col,row)
                  dgl2d(col,row)=dgl2d(col,row)+4.0*3600/2336/1000000.0/0.3
!             caculate subsurface rounoff
              runoff_inter(col,row)=0.0
              tmpqsub= 0.0
              tempcoe= 7.0

              call soilhk (ns2d(col,row),hkdt2d(col,row,:),matdt2d(col,row,:),vlcdt2d(col,row,:),vicdt2d(col,row,:),col,row)
              deltzi=(zs2d(col,row,2)-zs2d(col,row,1))/2.0
              if (tsdt2d(col,row,1)<0.0) then
                  fieldc=thfc2d(col,row,1)*vlcdt2d(col,row,1)/(vlcdt2d(col,row,1)+vicdt2d(col,row,1))
              else
                  fieldc=thfc2d(col,row,1)
              end if
              if((vlcdt2d(col,row,1)-fieldc) .gt. 0.0) then
                if(tsdt2d(col,row,1)<0.0) then
                  tempcoe=7.0
                else
                  tempcoe=7.0
                end if
                tmpqsub = tempcoe*hkdt2d(col,row,1)*sin(slope2d(col,row))*dtime*deltzi
                if(tmpqsub .lt. 0.1e-20) tmpqsub = 0.0
                tmp= (vlcdt2d(col,row,1)-fieldc)*deltzi*slopelen2d(col,row)
                tmpqsub = amin1(tmpqsub, tmp)
                tmpqsub = amax1(tmpqsub, 0.0)
                vlcdt2d(col,row,1)= vlcdt2d(col,row,1)-tmpqsub/deltzi/slopelen2d(col,row)
                runoff_inter(col,row) = runoff_inter(col,row) + tmpqsub
                runoffdepth(col,row)=runoffdepth(col,row)+tmpqsub*zs2d(col,row,1)
              end if

              do i = 2, ns2d(col,row)-1
                tmpqsub=0.0
                deltzi=(zs2d(col,row,i+1)-zs2d(col,row,i-1))/2.0
                if (tsdt2d(col,row,i)<0.0) then
                  fieldc=thfc2d(col,row,i)*vlcdt2d(col,row,1)/(vlcdt2d(col,row,1)+vicdt2d(col,row,1))
                  tempcoe=7.0
                else
                  fieldc=thfc2d(col,row,i)
                  tempcoe=7.0
                end if
                if((vlcdt2d(col,row,i)-fieldc) .gt. 0.0) then
                  call soilhk (ns2d(col,row),hkdt2d(col,row,:),matdt2d(col,row,:),vlcdt2d(col,row,:),vicdt2d(col,row,:),col,row)
                  tmpqsub = tempcoe*hkdt2d(col,row,i)*sin(slope2d(col,row))*dtime*deltzi
                  if(tmpqsub .lt. 0.1e-20) tmpqsub = 0.0
                  tmp     = (vlcdt2d(col,row,i)-fieldc)*deltzi*slopelen2d(col,row)
                  tmpqsub = amin1(tmpqsub, tmp)
                  tmpqsub = amax1(tmpqsub, 0.0)
                  vlcdt2d(col,row,i)= vlcdt2d(col,row,i)-tmpqsub/deltzi/ slopelen2d(col,row)
                  runoff_inter(col,row) = runoff_inter(col,row) + tmpqsub
                end if
                runoffdepth(col,row)=runoffdepth(col,row)+zs2d(col,row,i)*tmpqsub
              enddo
!              runoff_inter(col,row) = runoff_inter(col,row) + hkdt2d(col,row,ns2d(col,row))*dtime/&
!              (gridarea2d(col,row)/slopelen2d(col,row))
              deltzi=(zs2d(col,row,ns2d(col,row))-zs2d(col,row,ns2d(col,row)-1))/2.0
              !runoff_inter(col,row) = runoff_inter(col,row) + bot2d(col,row)*deltzi*slopelen2d(col,row)


              if((vlcdt2d(col,row,ns2d(col,row))-thfc2d(col,row,ns2d(col,row))) .gt. 0.0) then
                tmpqsub = tempcoe*hkdt2d(col,row,ns2d(col,row))*sin(slope2d(col,row))*dtime*deltzi
                if(tmpqsub .lt. 0.1e-20) tmpqsub = 0.0
                tmp= (vlcdt2d(col,row,ns2d(col,row))-thfc2d(col,row,ns2d(col,row)))*deltzi*slopelen2d(col,row)
                tmpqsub = amin1(tmpqsub, tmp)
                tmpqsub = amax1(tmpqsub, 0.0)
                vlcdt2d(col,row,ns2d(col,row))= vlcdt2d(col,row,ns2d(col,row))-tmpqsub/deltzi/slopelen2d(col,row)
                runoff_inter(col,row) = runoff_inter(col,row) + tmpqsub
                runoffdepth(col,row)=runoffdepth(col,row)+tmpqsub*zs2d(col,row,ns2d(col,row))
              end if
              !runoffdepth(col,row)=runoffdepth(col,row)/runoff_inter(col,row)
              !averunoffdepth=averunoffdepth+runoffdepth(col,row)
              runoffdis(col,row)=runoff2d(col,row)*1000.0+runoff_g(col,row)*1000.0+runoff_inter(col,row)/slopelen2d(col,row)*1000.0
            
              runoff_inter(col,row)=runoff_inter(col,row)/dtime
              runoff_g(col,row)=runoff_g(col,row)*slopelen2d(col,row)/dtime
!             roffg     =  roffg * length/dtlsm            ! m --> m3/sec/m
!             roffinter =  roffinter * length/dtlsm      ! m --> m3/sec/m
              runoff12d(col,row)=runoff12d(col,row)+runoff2d(col,row)
            end if
          end do
        end do
        !$omp end do
        !$omp end parallel
   end subroutine


end module hydro_mod
