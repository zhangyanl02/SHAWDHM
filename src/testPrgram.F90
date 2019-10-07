
program testprogram
    use dims_mod
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4 => shr_kind_r4
    use soilproperty_mod
    use controlpara_mod
    use constvar_mod,only:rhoi,rhol,rhom,rhoom
    use goshaw_mod
    use input_mod
    use savedata_mod,only:swe2d,fdepth2d,tdepth2d
    use statevar_mod,only:nc2d,nsp2d,nr2d,ns2d,inbasin2d,gridarea2d,slopelen2d,plthgt2d,pltwgt2d,pltlai2d,rootdp2d,&
      dchar2d,tccrit2d,rstom02d,rstexp2d,pleaf02d,rleaf02d,rroot02d,canalb2d,xangle2d,clumpng2d,pcandt2d,itype2d,&
      zc2d,tcdt2d,vapcdt2d,tlcdt2d,wcandt2d,zsp2d,dzsp2d,rhosp2d,tspdt2d,dlwdt2d,icespt2d,wlag2d,store2d,snowex2d,runoff2d,&
      thflux2d,vflux2d,zr2d,rhor2d,trdt2d,vaprdt2d,gmcdt2d,rload2d,zrthik2d,cover2d,albres2d,rescof2d,gmcmax,zs2d,&
      tsdt2d,vlcdt2d,vicdt2d,matdt2d,concdt2d,saltdt2d,icesdt2d,albdry2d,albexp2d,dgrade2d,sltdif2d,asalt2d,disper2d,&
      zmsrf2d,zhsrf2d,zersrf2d,zmsp2d,zhsp2d,dirres2d,pond2d,alatud2d,lontitud2d,elevation2d,slope2d,aspect2d,clouds2d,&
      tsavg2d,hrnoon2d,maskflag,sunhor2d,tmpday2d,winday2d,precip2d,snoden2d,humday2d,soitmp2d,vlcday2d,soilxt2d,&
      presur2d,etsum2d,evap12d,totflo2d,lwcan2d,swcan2d,lwsnow2d,lwsoil2d,swsoil2d,lwres2d,swres2d,swsnow2d,julian,&
      hour,year,maskgrid,rivernet,runoff_inter,runoff_g,waterflx,hkdt2d,dr2d,drw2d,runoff12d,dgl2d,ds2d,dg2d,bot2d,nplant2d,&
      skyview2d,shadow2d,AbsorbedSW2d,AbsorbedLW2d,snowmelt2d,InDirect2d,InDiffuse2d,TSDT2d_day,VLCDT2d_day,&
      VICDT2d_day,TRDT2d_day,TSPDT2d_day,zsp2d_day,EVAP12d_day,ETSUM2d_day,THFLUX2d_day,NSP2d_day,RUNOFF2d_day,DGL2D_day,&
      RUNOFFDIS_day,infill_day,snowmelt_day
    use hydro_mod,only:nsub,lateral_inflow,river_routing,read_hydro_para,lateral_routing,gwriv,runoff_generation
    use shaw27_mod,only:soilhk
    use soilproperty_mod,only:vapcof2d,vapexp2d
    implicit none

    integer(i4)::col,row,isub,inicon,juststart,hydro,i,k,j,nsp1
    real(r8),dimension(nx,ny)::declin,hafday
    integer(i4)::maxjul
    real(r8)::hum
    integer(i4)::idebug
    integer(i4)::lastyear
    real(r8)::wwdt,fieldc,deltzi

    integer(i4)::timer2,timer1,clock_rate,clock_max
    real(r8)::etime
    integer(i4)::etime_m,etime_h,etime_s

    real(r8),dimension(nx,ny)::runoffdepth,infil2d  !added
    real(r8)::averunoffdepth
    integer(i4)::ncount,l1,l2,year1
    real(r8)::totsurflow,totintflow,totgflow,tmpqsub,avkg,tempcoe,totflowbot,totflowtop
    real(r8)::swe,fdepth,tdepth
    real(r8)::runoffdis(nx,ny)
    character(200)::setupfile
    integer(i4)::jd_local,hour_local,year_local
    real(r8),dimension(nx,ny)::temp_dif
    character*(80):: tmp2  !88888888





!8888888888
    open(unit=12,file="../input/temp_dif.asc",status='old')
    read (12,*) tmp2
    read (12,*) tmp2
    read (12,*) tmp2
    read (12,*) tmp2
    read (12,*) tmp2
    read (12,*) tmp2

    do j = 1, ny
      read (12,*) (temp_dif(i,j), i = 1, nx)
    enddo 
    close(12)
!8888888889
    call system_clock(timer1,clock_rate,clock_max)
    call GETARG(1, setupfile)
    zsp2d(:,:,:)=0.0

    IDEBUG=1
    vapcof2d(:,:,:)=0.66 
    vapexp2d(:,:,:)=1.0
    !call omp_set_num_threads(32)
    call input(setupfile)
    call strlen(output_dir,l1,l2)
    if(hydro_module.eq.1) open(180,File=output_dir(l1:l2)//'RIVER_FLOW/runoff.txt',status='replace')
    satk2d(:,:,1:2)=satk2d(:,:,1:2)*5
    do i = 1,366
        do j = 1,16
           planth(j,i)=planth(j,i)*0.5
           if (planth(j,i)<0.01) planth(j,i)=0.01
        end do
    end do
    ds2d=4.0
    dg2d=6.0-ds2d
    dgl2d=1.6






    simulation_dir=output_dir
    hum=0.999
    lastyear=1
    maskgrid=0      !need to be restored and used to modify the inbasin2d array when restart
    julian=jstart
    hour=hrstar
    year=yrstar
    inital = 0
    if (hrstar .eq. 0) then
      hrstar=24
        jstart=jstart-1
    end if

    if (nhrpdt.eq.1) then
        wwdt=0.6
    else
        wwdt=1.0
    end if


    juststart = 1
    inicon=0
    rivernet(:,:)=0
    pond2d(:,:)=0.0
    runoff_inter(:,:)=0.0
    runoff_g(:,:)=0.0
    runoff12d(:,:)=0.0
    waterflx(:,:)=0.0
    totsurflow=0.0
    totintflow=0.0
    totgflow=0.0
!-------------------------------------
!   for hydro submodel
    if(hydro_module.eq.1)then
      call read_hydro_para(hydro_para_dir,nhrpdt)    
      call lateral_routing(nx,ny,juststart,inital,inicon,gridarea2d,slope2d,slopelen2d,dr2d,drw2d,&
         runoff_inter,runoff_g,waterflx,julian,hour,year,yrend,jend,hrend,simulation_dir,rivernet,&
         runoff12d,pond2d,pondmx,totsurflow,totintflow,totgflow)         
      juststart=0
    end if
!   end hydro submodel
!-------------------------------------

    maxjul=365
    if (mod(year,4) .eq. 0)  maxjul=365

!**** define initial boundary conditions for the day
    if (inph2o.ne.1) then
!        input soil moisture is water content
        vlcday2d(1,1,hour)=vlcdt2d(1,1,ns2d(1,1)) + vicdt2d(1,1,ns2d(1,1))*rhoi/rhol
    else
!        input soil moisture is matric potential
        vlcday2d(1,1,hour)=matdt2d(1,1,ns2d(1,1))
    end if
    do col=1,nx
      do row=1,ny
      vlcday2d(col,row,hour)=vlcday2d(1,1,hour)
      soitmp2d(col,row,hour)=soitmp2d(1,1,hour)
      end do
    enddo
!

    print*,"reading forcing data"
    call dayinp2 (julian,year,lastyear,mtstep,mpltgro,mwatrxt,sunhor2d,tmpday2d,winday2d,humday2d,&
    precip2d,snoden2d,soilxt2d,nplant,plthgt2d,dchar2d,pltwgt2d,pltlai2d,rootdp2d,presur2d,shadow2d)
    precip2d(:,:,:)=precip2d(:,:,:)*nhrpdt
    print*,"ending reading forcing data"
    do col=1,nx
       do row=1,ny
         call cloudy2 (clouds2d(col,row),alatud2d(col,row),declin(col,row),hafday(col,row),sunhor2d(col,row,:),julian,nhrpdt)
       end do
    end do

!**** begin simulation
      if (lvlout(13) .ne. 0) write (6,525)
!
!
!**** start of a new hour - update values assumed constant over the hour
    runoffdis=0.0
  100 dtime=nhrpdt*3600.
    if (hour.eq.3) print*,year,julian
    hour=hour+nhrpdt
    if (hour .gt. 24) then
!        start a new day
        hour=hour-24
        julian=julian + 1
        if (julian .gt. maxjul) then
            julian=julian-maxjul
            year=year + 1
            maxjul=365
            if (mod(year,4) .eq. 0)  maxjul=365
        end if
!need change----->adding hrend 
        if (year .eq. yrend  .and.  julian .gt. jend)   goto 504
        if (year .gt. yrend) goto 504
      call dayinp2 (julian,year,lastyear,mtstep,mpltgro,mwatrxt,sunhor2d,tmpday2d,winday2d,humday2d,&
      precip2d,snoden2d,soilxt2d,nplant,plthgt2d,dchar2d,pltwgt2d,pltlai2d,rootdp2d,presur2d,shadow2d)
      precip2d(:,:,:)=precip2d(:,:,:)*nhrpdt
      do col=1,nx
        do row=1,ny
          call cloudy2 (clouds2d(col,row),alatud2d(col,row),declin(col,row),hafday(col,row),sunhor2d(col,row,:),julian,nhrpdt)
        end do
      end do
    end if

    totflowbot=0.0
    totflowtop=0.0
    swe=0.0
    fdepth=0.0
    tdepth=0.0
    
    open(56,File=output_dir(l1:l2)//'temp.txt',status='unknown')
    open(57,File=output_dir(l1:l2)//'VLCDT.txt',status='unknown')
    tmpday2d(:,:,hour)=tmpday2d(:,:,hour)+temp_dif(:,:)

    !$omp parallel private(col,row) default(shared)
    !$omp do
    do col=1,NX!31,31      
        do row=1,NY!41,41
        if(inbasin2d(col,row).eq.1) then
        call goshaw(col,row,julian,hour,year,nsalt,nc2d(col,row),nsp2d(col,row),nr2d(col,row),ns2d(col,row),mzcinp,inph2o,&
      mwatrxt,ivlcbc,itmpbc,nplant2d(col,row),plthgt2d(col,row,:),pltwgt2d(col,row,:),pltlai2d(col,row,:),rootdp2d(col,row,:),&
      dchar2d(col,row,:),tccrit2d(col,row,:),rstom02d(col,row,:),rstexp2d(col,row,:),pleaf02d(col,row,:),rleaf02d(col,row,:),&
      rroot02d(col,row,:),pcandt2d(col,row,:),canalb2d(col,row,:),xangle2d(col,row,:),clumpng2d(col,row,:),itype2d(col,row,:),&
      zc2d(col,row,:),tcdt2d(col,row,:),tlcdt2d(col,row,:,:),vapcdt2d(col,row,:),wcandt2d(col,row,:),zsp2d(col,row,:),&
      dzsp2d(col,row,:),rhosp2d(col,row,:),tspdt2d(col,row,:),dlwdt2d(col,row,:),icespt2d(col,row,:),wlag2d(col,row,:),&
      store2d(col,row),snowex2d(col,row),zr2d(col,row,:),rhor2d(col,row,:),trdt2d(col,row,:),vaprdt2d(col,row,:),&
      gmcdt2d(col,row,:),gmcmax,rload2d(col,row),zrthik2d(col,row),cover2d(col,row),albres2d(col,row),rescof2d(col,row),&
      dirres2d(col,row),zs2d(col,row,:),tsdt2d(col,row,:),vlcdt2d(col,row,:),vicdt2d(col,row,:),matdt2d(col,row,:),&
      concdt2d(col,row,:,:),icesdt2d(col,row,:),saltdt2d(col,row,:,:),albdry2d(col,row),albexp2d(col,row),&
      dgrade2d(col,row,:),sltdif2d(col,row,:),asalt2d(col,row,:),disper2d(col,row,:),zmsrf2d(col,row),zhsrf2d(col,row),&
      zersrf2d(col,row),zmsp2d(col,row),zhsp2d(col,row),pond2d(col,row),pondmx,alatud2d(col,row),lontitud2d(col,row),&
      elevation2d(col,row),slope2d(col,row),aspect2d(col,row),hrnoon2d(col,row),clouds2d(col,row),declin(col,row),hafday(col,row),&
      sunhor2d(col,row,hour),tmpday2d(col,row,hour),winday2d(col,row,hour),humday2d(col,row,hour),presur2d(col,row,hour),&
      shadow2d(col,row,hour),skyview2d(col,row),precip2d(col,row,hour),snoden2d(col,row,hour),soitmp2d(col,row,hour),&
      vlcday2d(col,row,hour),soilxt2d(col,row,hour,:),inbasin2d(col,row),thflux2d(col,row),vflux2d(col,row),runoff2d(col,row),&
      tsavg2d(col,row),evap12d(col,row),etsum2d(col,row),totflo2d(col,row,:),lwcan2d(col,row,:,:),swcan2d(col,row,:,:),&
      lwsnow2d(col,row),lwsoil2d(col,row),lwres2d(col,row,:),swres2d(col,row,:),swsnow2d(col,row,:),swsoil2d(col,row),&
      wwdt,maskflag,maskgrid(col,row),bot2d(col,row),dgl2d(col,row),infil2d(col,row),AbsorbedSW2d(col,row),AbsorbedLW2d(col,row),&
      snowmelt2d(col,row),InDirect2d(col,row),InDiffuse2d(col,row))
      totflowbot=totflowbot+totflo2d(col,row,ns2d(col,row)-1)*1000000.0
      totflowtop=totflowtop+totflo2d(col,row,1)*1000000.0
      swe=swe+swe2d(col,row)*1000.0
      fdepth=fdepth+fdepth2d(col,row)
      tdepth=tdepth+tdepth2d(col,row)
      if(maskgrid(col,row).eq.1) then
      end if
        if (idebug.eq. 1) then
            if(col.eq.31 .and. row.eq.41) then
                write(56,56) hour,julian,year,tsdt2d(col,row,1:12)
                write(57,57) hour,julian,year,vlcdt2d(col,row,1:12)
                !(58,58),hour,julian,year,matdt2d(col,row,1:12)
                !write(59,59),hour,julian,year,lwcan2d(col,row,1,1)
                !write(62,62),hour,julian,year,swcan2d(col,row,1,1)
                !write(60,60),hour,julian,year,lwres2d(col,row,1)
                !write(63,63),hour,julian,year,swres2d(col,row,1)
                !write(61,61),hour,julian,year,lwsoil2d(col,row),lwsnow2d(col,row)
                !write(64,64),hour,julian,year,swsoil2d(col,row)
                !write(65,65),hour,julian,year,trdt2d(col,row,1:nr2d(col,row))
                !write(67,67),hour,julian,year,tspdt2d(col,row,1:10)
                !write(68,68),hour,julian,year,zsp2d(col,row,nsp2d(col,row)+1),nsp2d(col,row)
                !write(85,85),hour,julian,year,icesdt2d(col,row,1:12)
            endif
        endif
        endif
      end do
    enddo
    !$omp end do
    !$omp end parallel 

    lastyear=year
    inicon=0
    if(hydro_module.eq.1) call runoff_generation(nx,ny,runoffdis,runoffdepth,dtime,ncount,averunoffdepth)
    write(180,*) averunoffdepth,totflowbot,totflowtop,swe,fdepth,tdepth

    if(hydro_module.eq.1) then
      totsurflow = 0.0
      totintflow = 0.0
      totgflow   = 0.0
      averunoffdepth=averunoffdepth/ncount
      swe=swe/ncount
      fdepth=fdepth/ncount
      tdepth=tdepth/ncount        
      call lateral_routing(nx,ny,juststart,inital,inicon,gridarea2d,slope2d,slopelen2d,dr2d,drw2d,&
        runoff_inter,runoff_g,waterflx,julian,hour,year,yrend,jend,hrend,simulation_dir,rivernet,&
        runoff12d,pond2d,pondmx,totsurflow,totintflow,totgflow)
    end if
!   end hydro submodel
!-------------------------------------


    trdt2d_day(:,:,:,int(hour/nhrpdt))=trdt2d(:,:,:)
    tsdt2d_day(:,:,:,int(hour/nhrpdt))=tsdt2d(:,:,1:NS)
    vlcdt2d_day(:,:,:,int(hour/nhrpdt))=vlcdt2d(:,:,1:NS)
    vicdt2d_day(:,:,:,int(hour/nhrpdt))=vicdt2d(:,:,1:NS)
    nsp2d_day(:,:,int(hour/nhrpdt))=nsp2d(:,:)
    tspdt2d_day(:,:,int(hour/nhrpdt))=tspdt2d(:,:,1)
    zsp2d_day(:,:,int(hour/nhrpdt))=zsp2d(:,:,NSPMAX)
    evap12d_day(:,:,int(hour/nhrpdt))=evap12d(:,:)
    etsum2d_day(:,:,int(hour/nhrpdt))=etsum2d(:,:)
    thflux2d_day(:,:,int(hour/nhrpdt))=thflux2d(:,:)
    runoff2d_day(:,:,int(hour/nhrpdt))=runoff2d(:,:)
    dgl2d_day(:,:,int(hour/nhrpdt))=dgl2d(:,:)
    runoffdis_day(:,:,int(hour/nhrpdt))=runoffdis(:,:)
    infill_day(:,:,int(hour/nhrpdt))=infil2d(:,:)
    snowmelt_day(:,:,int(hour/nhrpdt))=snowmelt2d(:,:)

    if (hour .eq. 24) then
      if(daily_output .eq. 0)then
        call output_hourly(trdt2d_day,tsdt2d_day,vlcdt2d_day,vicdt2d_day,NR2d,NS,tspdt2d_day,zsp2d_day,evap12d_day,&
          etsum2d_day,thflux2d_day,year,julian,result1_dir,maskgrid,runoff2d_day,dgl2d_day,runoffdis_day,infill_day,&
          snowmelt_day,deflate,ntimes)
      else 
        call output_daily(trdt2d_day,tsdt2d_day,vlcdt2d_day,vicdt2d_day,NR2d,NS,tspdt2d_day,zsp2d_day,evap12d_day,&
          etsum2d_day,thflux2d_day,year,julian,result1_dir,maskgrid,runoff2d_day,dgl2d_day,runoffdis_day,infill_day,&
          snowmelt_day,deflate,ntimes)
      end if
    end if
    inital=1
    go to 100
 
    
 504 call system_clock(timer2,clock_rate,clock_max)
     etime=(timer2-timer1)/clock_rate
     etime_h = etime/3600
     etime_m = (etime-etime_h*3600)/60.0
     etime_s = (etime-etime_h*3600-etime_m*60)/10*10
!     close(12)
     print*,"elapsed time in seconds: ",etime
     print*,"elapsed time in seconds: ",etime_h,"hour",etime_m,"min",etime_s,"sec"

!    output the soil state of last step
     call output_endstate(tsdt2d,vlcdt2d,vicdt2d,matdt2d,icesdt2d,nx,ny,ns2d(1,1),year,julian,hour,result2_dir)
     do row=1,ny
      write(100,*) (dgl2d(col,row),col=1,nx)
     enddo
     close(56)
     stop



 56 FORMAT (I4,I5,I6,12F13.4)
 57 FORMAT (I4,I5,I6,12F13.4)
 58 FORMAT (I4,I5,I6,12F13.4)
 59 FORMAT (I4,I5,I6,10F13.4)
 60 FORMAT (I4,I5,I6,10F13.4)
 61 FORMAT (I4,I5,I6,2F13.4)
 62 FORMAT (I4,I5,I6,10F13.4)
 63 FORMAT (I4,I5,I6,10F13.4)
 64 FORMAT (I4,I5,I6,2F13.4)
 65 FORMAT (I4,I5,I6,10F13.4)
 67 FORMAT (I4,I5,I6,15F9.4)
 68 FORMAT (I4,I5,I6,F9.4,I5)
 85 FORMAT (I4,I5,I6,12I3)

 525 FORMAT(/20X,'Day   Hour   Year   Min Steps  Max Steps'/)

   90 WRITE (6,*) 'CANNOT FIND SOIL TEMPERATURE DATA FOR INITIAL HOUR'
      WRITE (21,*) 'CANNOT FIND SOIL TEMPERATURE DATA FOR INITIAL HOUR'
      STOP
   95 WRITE (6,*) 'CANNOT FIND SOIL MOISTURE DATA FOR INITIAL HOUR'
      WRITE (21,*) 'CANNOT FIND SOIL MOISTURE DATA FOR INITIAL HOUR'
      STOP

  102 FORMAT (' ',A80)
  105 FORMAT (/5X,'SIMULATION BEGINS ON DAY',I4,', HOUR',I3,', OF ',I4, &
                /5X,'SIMULATION ENDS   ON DAY',I4,', HOUR 24  OF ',I4)
  110 FORMAT (//,' GENERAL SITE DESCRIPTION',//5X,'LATITUDE  :',F5.1, &
        ' DEGREES ',F5.2,' MINUTES',/5X,'SLOPE     : ',F5.2,' %', &
        /5X,'ASPECT    : ',F5.1,' DEGREES FROM NORTH',/5X,'SOLAR NOON: ',&
        F5.2,' HOURS',/5X,'ELEVATION : ',F5.0,' M'/5X,'PRESSURE  : ', &
        F5.1,' KPA')
  115 FORMAT (//,' WIND PROFILE AND SURFACE PARAMETERS',//5X, &
    'ZM =',F5.2,' CM (WHEN SNOW AND CANOPY ARE NOT PRESENT)',/5X,&
    'ZH =',F5.2,' CM (WHEN SNOW AND CANOPY ARE NOT PRESENT)',/5X,&
    'ZERO PLANE DISPLACEMENT =',F5.2,' M',/5X,&
    'HEIGHT OF INSTRUMENTATION =',F5.1,' M',/5X,&
    'MAXIMUM DEPTH OF PONDING =',F4.1,' CM',///1X,&
    'ATMOSPHERIC CONDITIONS',//5X,'MAXIMUM CLEAR-SKY TRANSMISSIVITY ='&
    ,F5.2,/5X,'CLEAR-SKY LONG-WAVE EMISSIVITY PARAMETERS :',&
    F6.3,E15.3)
  158 FORMAT (5X,F5.2,F8.2,F10.2,F9.3,F8.0,2X,4F5.1,F6.1,F8.3)
  159 FORMAT ('JDAY,JH,HYR,TMP,WIND,HUM,SUNHOR:',I4,I3,I5,2X,F8.2,F8.2,F8.2,F8.2,F8.2)
  160 FORMAT (//,' SIMULATION BEGINS',/)
  170 FORMAT (//,' CANOPY PARAMETERS',//5X,&
    'EMISSIVITY OF PLANT MATERIAL :',F5.2,/5X,&
    'MOISTURE PARAMETERS FOR ANY DEAD PLANT MATERIAL :',F7.2,F5.2,&
    /5X,'INITIAL MOISTURE CONTENT FOR DEAD PLANT MATERIAL :',F5.2,&
    ' KG/KG',&
    /5X,'MAXIMUM MOISTURE CONTENT FOR DEAD PLANT MATERIAL :',F5.2,&
    ' KG/KG',//40X,8(A,I1))
end program


subroutine output_endstate(TSDT2d,VLCDT2d,VICDT2d,MATDT2d,ICESDT2d,NX,NY,NS,YEAR,JULIAN,HOUR,output_dir)
    use input_mod,only:check
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4 => shr_kind_r4
    use netcdf
    implicit none
    ! This is the name of the data file we will create.
    integer(i4),intent(in)::NX,NY,NS,YEAR,JULIAN,HOUR
    real(r8),dimension(NX,NY,NS),intent(in)::TSDT2d,VLCDT2d,VICDT2d,MATDT2d
    integer(i4),dimension(NX,NY,NS),intent(in)::ICESDT2d
    character(80),intent(in)::output_dir
    character (len=4)::SYEAR
    character (len=3)::SJULIAN
    character (len=3)::SHOUR
    integer::RETVAL,I,FLAG,TMP
    integer::ncid,xid,yid,zid
    integer,dimension(3)::dimids
    integer::TSDT_ID,VLCDT_ID,VICDT_ID,MATDT_ID,ICESDT_ID,TSW_ID
    integer::x_dimid,y_dimid,z_dimid
    character*(80)::ncfilename


    WRITE(SYEAR,"(I4.4)") YEAR
    WRITE(SJULIAN,"(I3.3)") JULIAN
    WRITE(SHOUR,"(I3.3)") 100+HOUR

    ncfilename=trim(output_dir)//trim(SYEAR)//"-"//trim(SJULIAN)//"-"//trim(SHOUR(2:3))//".nc" 
    call check( nf90_create(ncfilename, NF90_CLOBBER, ncid) )

    ! Define the dimensions. NetCDF will hand back an ID for each. 
    call check( nf90_def_dim(ncid, "NX", NX, x_dimid) )
    call check( nf90_def_dim(ncid, "NY", NY, y_dimid) )
    call check( nf90_def_dim(ncid, "NS", NS, z_dimid) )

    dimids(1) = x_dimid
    dimids(2) = y_dimid
    dimids(3) = z_dimid

    ! Define the variable.
    call check( nf90_def_var(ncid, "TSDT2d", NF90_REAL, dimids, TSDT_ID) )
    call check( nf90_def_var(ncid, "VLCDT2d" , NF90_REAL, dimids, VLCDT_ID ) )
    call check( nf90_def_var(ncid, "VICDT2d" , NF90_REAL, dimids, VICDT_ID ) )
    call check( nf90_def_var(ncid, "MATDT2d" , NF90_REAL, dimids, MATDT_ID ) )
    call check( nf90_def_var(ncid, "ICESDT2d" , NF90_INT, dimids, ICESDT_ID ) )
    call check( nf90_def_var(ncid, "TSW" , NF90_REAL, dimids, TSW_ID ) )
    call check( nf90_enddef(ncid) )
!   Store variable values
    call check( nf90_put_var(ncid, TSDT_ID, TSDT2d) )
    call check( nf90_put_var(ncid, VLCDT_ID, VLCDT2d))
    call check( nf90_put_var(ncid, VICDT_ID, VICDT2d) )
    call check( nf90_put_var(ncid, MATDT_ID, MATDT2d))
    call check( nf90_put_var(ncid, ICESDT_ID, ICESDT2d) )
    call check( nf90_put_var(ncid, TSW_ID, VLCDT2d+VICDT2d*0.9) )
    call check( nf90_close(ncid) )
 end subroutine

subroutine output_hourly(trdt2d,TSDT2d,VLCDT2d,VICDT2d,NR2d,NS,TSPDT2d,ZSPDT2d,EVAP12d,ETSUM2d,THFLUX2d,&
      YEAR,JULIAN,output_dir,maskgrid,RUNOFF2d,DGL2D,RUNOFFDIS,infil2d,snowmelt2d,deflate,ntimes)
    use dims_mod
    use controlpara_mod,only:T_soil_out,VLC_soil_out,VIC_soil_out,T_Residue_out,TSP_out,ZSP_out,DGL_out,EVAP_out,&
                     ET_out,Runoff_out,surfacerunoff_out,SWE_out,THFLUX_out,infill_out,snowmelt_out
    use input_mod,only:check
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4 => shr_kind_r4
    use netcdf
    use savedata_mod,only:SWE2d,FDEPTH2d,TDEPTH2d
    use radout_mod,only:sw_on_soil2d,sw_on_snow2d
    implicit none
    ! This is the name of the data file we will create.
    integer(i4),intent(in)::NS
    integer(i4),intent(in)::ntimes
    integer(i4),dimension(NX,NY),intent(in)::NR2d
    real(r8),dimension(NX,NY,NS,ntimes),intent(inout)::TSDT2d,VLCDT2d,VICDT2d
    real(r8),dimension(NX,NY,NRMAX,ntimes),intent(inout)::TRDT2d
    real(r8),dimension(NX,NY,NSPMAX,ntimes),intent(inout)::TSPDT2d
    real(r8),dimension(NX,NY,ntimes),intent(inout)::ZSPDT2d
    real(r8),dimension(NX,NY,ntimes),intent(in)::EVAP12d,ETSUM2d,THFLUX2d
    integer(i4),intent(in)::YEAR,JULIAN,deflate
    character(80),intent(in)::output_dir
    real(r8),dimension(NX,NY),intent(in)::maskgrid
    real(r8),dimension(NX,NY,ntimes),intent(in)::RUNOFF2d,DGL2D,RUNOFFDIS,infil2d,snowmelt2d
    character (len=4)::SYEAR
    character (len=3)::SJULIAN

    integer::I,FLAG,TMP,COL,ROW,NSP
    integer::ncid,xid,yid,zid
    integer,dimension(4)::dimids_4d
    integer,dimension(3)::dimids_3d
    integer,dimension(2)::dimids_2d
    real(r8),dimension(NX,NY,ntimes)::ZSP



    integer::TRDT_ID,TSDT_ID,VLCDT_ID,TSPDT_ID,ZSPDT_ID,invalidp_ID,evap_id,etsum_id,THFLUX_ID,&
             SWE_ID,FDEPTH_ID,TDEPTH_ID,runoff_ID,VICDT_ID,dgl_id,RUNOFFDIS_id,infil2d_ID,snowmelt_ID,absorbedLW_ID,&
             absorbedSW_ID,InDirect2d_ID,InDiffuse2d_ID,TSP_ID
    integer::x_dimid,y_dimid,z_dimid,h_dimid,nsp_dimid
    character*(80)::ncfilename
    logical::alive
    logical::output_snow


      WRITE(SYEAR,"(I4.4)") YEAR
      WRITE(SJULIAN,"(I3.3)") JULIAN
      output_snow=.TRUE.

      ncfilename=trim(output_dir)//trim(SYEAR)//"-"//trim(SJULIAN)//".nc"
      call check( nf90_create(ncfilename, NF90_NETCDF4, ncid) )
      !call check( nf90_create(ncfilename, NF90_CLOBBER, ncid) )
      ! Define the dimensions. NetCDF will hand back an ID for each. 
      call check( nf90_def_dim(ncid, "NX", NX, x_dimid) )
      call check( nf90_def_dim(ncid, "NY", NY, y_dimid) )
      call check( nf90_def_dim(ncid, "NS", NS, z_dimid) )
      call check (nf90_def_dim(ncid, "HOUR",ntimes,h_dimid))

      dimids_4d(1) = x_dimid
      dimids_4d(2) = y_dimid
      dimids_4d(3) = z_dimid
      dimids_4d(4) = h_dimid

      dimids_2d(1)=x_dimid
      dimids_2d(2)=y_dimid

      dimids_3d(1)=x_dimid
      dimids_3d(2)=y_dimid
      dimids_3d(3)=h_dimid


      ! Define the variable.
      call check( nf90_def_var(ncid, "InvalidPoints", NF90_REAL, dimids_2d, invalidp_ID) )
      if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,invalidp_ID,0,4,4))

      if(T_soil_out .gt. 0)then
        call check( nf90_def_var(ncid, "TSDT2d", NF90_REAL, dimids_4d, TSDT_ID) )
        if(deflate.eq.1) call check( nf90_def_var_deflate(ncid,TSDT_ID,0,4,4))
      end if

      if(T_Residue_out.gt.0)then
        call check( nf90_def_var(ncid, "TRDT2d", NF90_REAL, dimids_4d, TRDT_ID) )
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,TRDT_ID,0,4,4))
      end if

      if(TSP_out .gt. 0)then
        call check( nf90_def_var(ncid, "TSP2d", NF90_REAL, dimids_3d, TSP_ID) )
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,TSP_ID,0,4,4))
      end if

      if(DGL_out.gt.0)then
        call check( nf90_def_var(ncid, "DGL2d", NF90_REAL, dimids_2d, dgl_ID) )
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,dgl_ID,0,4,4))
      end if

      if(VLC_soil_out.gt.0)then
        call check( nf90_def_var(ncid, "VLCDT2d" , NF90_REAL, dimids_4d, VLCDT_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,VLCDT_ID,0,4,4))
      end if

      if(VIC_soil_out .gt. 0) then
        call check( nf90_def_var(ncid, "VICDT2d" , NF90_REAL, dimids_4d, VICDT_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,VICDT_ID,0,4,4))
      end if

      if(EVAP_out .gt. 0) then
        call check( nf90_def_var(ncid, "EVAP" , NF90_REAL, dimids_3d, evap_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,evap_ID,0,4,4))
      end if

      if(ET_out .gt. 0)then
        call check( nf90_def_var(ncid, "ET" , NF90_REAL, dimids_3d, etsum_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,etsum_ID,0,4,4))
      end if

      if(THFLUX_out .gt. 0) then
        call check( nf90_def_var(ncid, "THFLUX" , NF90_REAL, dimids_3d, THFLUX_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,THFLUX_ID,0,4,4))
      end if

      if(SWE_out .gt. 0) then
        call check( nf90_def_var(ncid, "SWE" , NF90_REAL, dimids_3d, SWE_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,SWE_ID,0,4,4))
      end if

      if(surfacerunoff_out .gt. 0) then
        call check( nf90_def_var(ncid, "surfacerunoff" , NF90_REAL, dimids_3d, runoff_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,runoff_ID,0,4,4))
      end if


      if(Runoff_out .gt. 0)then
        call check( nf90_def_var(ncid, "RUNOFFDIS" , NF90_REAL, dimids_3d, RUNOFFDIS_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,RUNOFFDIS_ID,0,4,4))
      end if

      if(infill_out .gt. 0) then
        call check( nf90_def_var(ncid, "infil2d" , NF90_REAL, dimids_3d, infil2d_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,infil2d_ID,0,4,4))
      end if

      if(snowmelt_out .gt. 0) then
        call check( nf90_def_var(ncid, "snowmelt" , NF90_REAL, dimids_3d, snowmelt_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,snowmelt_ID,0,4,4))
      end if

      if(ZSP_out .gt. 0)then
        call check( nf90_def_var(ncid, "ZSP2d", NF90_REAL, dimids_3d, ZSPDT_ID))
      endif
      call check( nf90_enddef(ncid) )

!   Store variable values
      call check( nf90_inq_varid(ncid, "InvalidPoints", invalidp_ID) )
      call check( nf90_put_var(ncid, invalidp_ID, maskgrid))

      if(T_soil_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "TSDT2d", TSDT_ID) )
        call check( nf90_put_var(ncid, TSDT_ID, TSDT2d))
      end if

      if(T_Residue_out.gt.0)then
        call check( nf90_inq_varid(ncid, "TRDT2d", TRDT_ID) )
        call check( nf90_put_var(ncid, TRDT_ID, TRDT2d,start =(/1,1,1,1/),count = (/NX,NY,NR2d(1,1),ntimes/)))
      end if

      if(VLC_soil_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "VLCDT2d", VLCDT_ID) )
        call check( nf90_put_var(ncid, VLCDT_ID, VLCDT2d))
      end if

      if(VIC_soil_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "VICDT2d", VICDT_ID) )
        call check( nf90_put_var(ncid, VICDT_ID, VICDT2d))
      end if

      if(runoff_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "RUNOFFDIS", RUNOFFDIS_ID) )
        call check( nf90_put_var(ncid, RUNOFFDIS_ID, RUNOFFDIS,start =(/1,1,1/),count = (/NX,NY,ntimes/)))
      end if

      if(surfacerunoff_out .gt. 0) then
        call check( nf90_inq_varid(ncid, "surfacerunoff", runoff_ID) )
        call check( nf90_put_var(ncid, runoff_ID, RUNOFF2d,start =(/1,1,1/),count = (/NX,NY,ntimes/)))
      end if

      if (EVAP_out .gt. 0) then
         call check( nf90_inq_varid(ncid, "EVAP", evap_ID) )
         call check( nf90_put_var(ncid, evap_ID, EVAP12d))
      end if

      if (ET_out .gt. 0) then
        call check( nf90_inq_varid(ncid, "ET", etsum_ID) )
        call check( nf90_put_var(ncid, etsum_ID,ETSUM2d))
      end if 

      if(DGL_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "DGL2d", dgl_ID) )
        call check( nf90_put_var(ncid, dgl_ID, DGL2d))      
      end if

      if(THFLUX_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "THFLUX", THFLUX_ID) )
        call check( nf90_put_var(ncid, THFLUX_ID,THFLUX2d))
      end if

      if(SWE_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "SWE", SWE_ID) )
        call check( nf90_put_var(ncid, SWE_ID,SWE2d))
      end if

      if(ZSP_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "ZSP2d", ZSPDT_ID) )
        call check( nf90_put_var(ncid, ZSPDT_ID, ZSPDT2d))
      endif

      if(TSP_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "ZSP2d", ZSPDT_ID) )
        call check( nf90_put_var(ncid, ZSPDT_ID, ZSPDT2d))
      endif

      if(infill_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "infil2d", infil2d_ID) )
        call check( nf90_put_var(ncid, infil2d_ID, infil2d,start =(/1,1,1/),count = (/NX,NY,ntimes/)))
      endif

      if(snowmelt_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "snowmelt", snowmelt_ID))
        call check( nf90_put_var(ncid, snowmelt_ID, snowmelt2d,start =(/1,1,1/),count = (/NX,NY,ntimes/)))
      endif
      call check( nf90_close(ncid) )    

    !call check( nf90_inq_varid(ncid, "absorbedSW", absorbedSW_ID) )
    !call check( nf90_inq_varid(ncid, "absorbedLW", absorbedLW_ID) )
    !call check( nf90_inq_varid(ncid, "InDirect", InDirect2d_ID) )
    !call check( nf90_inq_varid(ncid, "InDiffuse", InDiffuse2d_ID) )
    !call check( nf90_put_var(ncid, absorbedSW_ID,AbsorbedSW2d,start =(/1,1,1/),count = (/NX,NY,ntimes/)))
    !call check( nf90_put_var(ncid, absorbedLW_ID, AbsorbedLW2d,start =(/1,1,1/),count = (/NX,NY,ntimes/)))
    !call check( nf90_put_var(ncid, InDirect2d_ID, InDirect2d,start =(/1,1,1/),count = (/NX,NY,ntimes/)))
    !call check( nf90_put_var(ncid, InDiffuse2d_ID, InDiffuse2d,start =(/1,1,1/),count = (/NX,NY,ntimes/)))
    end subroutine



subroutine output_daily(trdt2d,TSDT2d,VLCDT2d,VICDT2d,NR2d,NS,TSPDT2d,ZSPDT2d,EVAP12d,ETSUM2d,THFLUX2d,&
      YEAR,JULIAN,output_dir,maskgrid,RUNOFF2d,DGL2D,RUNOFFDIS,infil2d,snowmelt2d,deflate,ntimes)


    use dims_mod
    use controlpara_mod,only:T_soil_out,VLC_soil_out,VIC_soil_out,T_Residue_out,TSP_out,ZSP_out,DGL_out,EVAP_out,&
                     ET_out,Runoff_out,surfacerunoff_out,SWE_out,THFLUX_out,infill_out,snowmelt_out
    use input_mod,only:check
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4 => shr_kind_r4
    use netcdf
    use savedata_mod,only:SWE2d,FDEPTH2d,TDEPTH2d
    use radout_mod,only:sw_on_soil2d,sw_on_snow2d
    implicit none
    ! This is the name of the data file we will create.
    integer(i4),intent(in)::NS,ntimes
    integer(i4),dimension(NX,NY),intent(in)::NR2d
    real(r8),dimension(NX,NY,NS,ntimes),intent(inout)::TSDT2d,VLCDT2d,VICDT2d
    real(r8),dimension(NX,NY,NRMAX,ntimes),intent(inout)::TRDT2d
    real(r8),dimension(NX,NY,ntimes),intent(inout)::TSPDT2d,ZSPDT2d
    real(r8),dimension(NX,NY,ntimes),intent(in)::EVAP12d,ETSUM2d,THFLUX2d

    integer(i4),intent(in)::YEAR,JULIAN,deflate
    character(80),intent(in)::output_dir
    real(r8),dimension(NX,NY),intent(in)::maskgrid
    real(r8),dimension(NX,NY,ntimes),intent(in)::RUNOFF2d,DGL2D,RUNOFFDIS,infil2d,snowmelt2d
    character (len=4)::SYEAR
    character (len=3)::SJULIAN

    integer::I,FLAG,TMP,COL,ROW,NSP
    integer::ncid,xid,yid,zid
    integer,dimension(4)::dimids_4d
    integer,dimension(3)::dimids_3d
    integer,dimension(2)::dimids_2d
    real(r8),dimension(NX,NY)::ZSP

    real(r8),dimension(nx,ny,ns)::TSDT2d_daily,VLCDT2d_daily,VICDT2d_daily
    real(r8),dimension(NX,NY)::EVAP12d_daily,ETSUM2d_daily,THFLUX2d_daily,DGL2D_daily
    real(r8),dimension(nx,ny,NRMAX)::TRDT2d_daily
    real(r8),dimension(nx,ny)::ZSP_daily,TSP_daily,RUNOFFDIS_daily,surfacerunoff_daily,infill_daily,snowmelt_daily





    integer::TRDT_ID,TSDT_ID,VLCDT_ID,TSPDT_ID,ZSPDT_ID,invalidp_ID,evap_id,etsum_id,THFLUX_ID,&
             SWE_ID,FDEPTH_ID,TDEPTH_ID,runoff_ID,VICDT_ID,dgl_id,RUNOFFDIS_id,infil2d_ID,snowmelt_ID,absorbedLW_ID,&
             absorbedSW_ID,InDirect2d_ID,InDiffuse2d_ID,h,TSP_ID
    integer::x_dimid,y_dimid,z_dimid,h_dimid,nsp_dimid
    character*(80)::ncfilename
    logical::alive


    WRITE(SYEAR,"(I4.4)") YEAR
    WRITE(SJULIAN,"(I3.3)") JULIAN


      if(T_soil_out .gt. 0) TSDT2d_daily(:,:,:)=0.0
      if(VLC_soil_out .gt. 0) VLCDT2d_daily(:,:,:)=0.0
      if(VIC_soil_out .gt. 0) VICDT2d_daily(:,:,:)=0.0
      if(EVAP_out .gt. 0) EVAP12d_daily(:,:)=0.0
      if(ET_out .gt. 0)ETSUM2d_daily(:,:)=0.0
      if(THFLUX_out .gt. 0)THFLUX2d_daily(:,:)=0.0
      if(DGL_out .gt. 0)DGL2D_daily(:,:)=0.0
      if(T_Residue_out .gt. 0)TRDT2d_daily(:,:,:)=0.0
      if(ZSP_out .gt. 0)  ZSP_daily(:,:)=0.0
      if(TSP_out .gt. 0)   TSP_daily(:,:)=0.0
      if(Runoff_out .gt. 0)RUNOFFDIS_daily(:,:)=0.0
      if(surfacerunoff_out .gt. 0) surfacerunoff_daily(:,:)=0.0
      if(infill_out .gt. 0)  infill_daily(:,:)=0.0
      if(snowmelt_out .gt. 0)  snowmelt_daily(:,:)=0.0
      do h =1,ntimes
        if(T_soil_out .gt. 0)  TSDT2d_daily(:,:,:)=TSDT2d_daily(:,:,:)+TSDT2d(:,:,:,h)/ntimes
        if(VLC_soil_out .gt. 0)VLCDT2d_daily(:,:,:)=VLCDT2d_daily(:,:,:)+VLCDT2d(:,:,:,h)/ntimes
        if(VIC_soil_out .gt. 0)VICDT2d_daily(:,:,:)=VICDT2d_daily(:,:,:)+VICDT2d(:,:,:,h)/ntimes
        if(EVAP_out .gt. 0)EVAP12d_daily(:,:)=EVAP12d_daily(:,:)+EVAP12d(:,:,h)
        if(ET_out .gt. 0)ETSUM2d_daily(:,:)=ETSUM2d_daily(:,:)+ETSUM2d(:,:,h)
        if(THFLUX_out .gt. 0)THFLUX2d_daily(:,:)=THFLUX2d_daily(:,:)+THFLUX2d(:,:,h)/ntimes
        if(DGL_out .gt. 0)DGL2D_daily(:,:)=DGL2D_daily(:,:)+DGL2D(:,:,h)/ntimes
        if(T_Residue_out .gt. 0) TRDT2d_daily(:,:,:)=TRDT2d_daily(:,:,:)+TRDT2d(:,:,:,h)/ntimes
        if(ZSP_out .gt. 0)   ZSP_daily(:,:)=ZSP_daily(:,:)+ZSPDT2d(:,:,h)/ntimes
        if(TSP_out .gt. 0)   TSP_daily(:,:)=TSP_daily(:,:)+TSPDT2d(:,:,h)/ntimes
        if(Runoff_out .gt. 0)RUNOFFDIS_daily(:,:)=RUNOFFDIS_daily(:,:)+RUNOFFDIS(:,:,h)
        if(surfacerunoff_out .gt. 0) surfacerunoff_daily(:,:)=surfacerunoff_daily(:,:)+RUNOFF2d(:,:,h)
        if(infill_out .gt. 0) infill_daily(:,:)=infill_daily(:,:)+infil2d(:,:,h)
        if(snowmelt_out .gt. 0) snowmelt_daily(:,:)=snowmelt_daily(:,:)+snowmelt2d(:,:,h)
      end do

      ncfilename=trim(output_dir)//trim(SYEAR)//"-"//trim(SJULIAN)//".nc"
      call check( nf90_create(ncfilename, NF90_NETCDF4, ncid) )

      call check( nf90_def_dim(ncid, "NX", NX, x_dimid) )
      call check( nf90_def_dim(ncid, "NY", NY, y_dimid) )
      call check( nf90_def_dim(ncid, "NS", NS, z_dimid) )
      dimids_3d(1) = x_dimid
      dimids_3d(2) = y_dimid
      dimids_3d(3) = z_dimid

      dimids_2d(1)=x_dimid
      dimids_2d(2)=y_dimid

      ! Define the variable.
      call check( nf90_def_var(ncid, "InvalidPoints", NF90_REAL, dimids_2d, invalidp_ID) )
      if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,invalidp_ID,0,4,4))

      if(T_soil_out .gt. 0)then
        call check( nf90_def_var(ncid, "TSDT2d", NF90_REAL, dimids_3d, TSDT_ID) )
        if(deflate.eq.1) call check( nf90_def_var_deflate(ncid,TSDT_ID,0,4,4))
      end if

      if(T_Residue_out.gt.0)then
        call check( nf90_def_var(ncid, "TRDT2d", NF90_REAL, dimids_3d, TRDT_ID) )
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,TRDT_ID,0,4,4))
      end if

      if(TSP_out .gt. 0)then
        call check( nf90_def_var(ncid, "TSP2d", NF90_REAL, dimids_2d, TSP_ID) )
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,TSP_ID,0,4,4))
      end if

      if(DGL_out.gt.0)then
        call check( nf90_def_var(ncid, "DGL2d", NF90_REAL, dimids_2d, dgl_ID) )
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,dgl_ID,0,4,4))
      end if

      if(VLC_soil_out.gt.0)then
        call check( nf90_def_var(ncid, "VLCDT2d" , NF90_REAL, dimids_3d, VLCDT_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,VLCDT_ID,0,4,4))
      end if

      if(VIC_soil_out .gt. 0) then
        call check( nf90_def_var(ncid, "VICDT2d" , NF90_REAL, dimids_3d, VICDT_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,VICDT_ID,0,4,4))
      end if

      if(EVAP_out .gt. 0) then
        call check( nf90_def_var(ncid, "EVAP" , NF90_REAL, dimids_2d, evap_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,evap_ID,0,4,4))
      end if

      if(ET_out .gt. 0)then
        call check( nf90_def_var(ncid, "ET" , NF90_REAL, dimids_2d, etsum_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,etsum_ID,0,4,4))
      end if


      if(THFLUX_out .gt. 0) then
        call check( nf90_def_var(ncid, "THFLUX" , NF90_REAL, dimids_2d, THFLUX_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,THFLUX_ID,0,4,4))
      end if

      if(SWE_out .gt. 0) then
        !call check( nf90_def_var(ncid, "SWE" , NF90_REAL, dimids_2d, SWE_ID))
        !if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,SWE_ID,0,4,4))
      end if

      if(surfacerunoff_out .gt. 0) then
        call check( nf90_def_var(ncid, "surfacerunoff" , NF90_REAL, dimids_2d, runoff_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,runoff_ID,0,4,4))
      end if

      if(Runoff_out .gt. 0)then
        call check( nf90_def_var(ncid, "RUNOFFDIS" , NF90_REAL, dimids_2d, RUNOFFDIS_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,RUNOFFDIS_ID,0,4,4))
      end if

      if(infill_out .gt. 0) then
        call check( nf90_def_var(ncid, "infill" , NF90_REAL, dimids_2d, infil2d_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,infil2d_ID,0,4,4))
      end if

      if(snowmelt_out .gt. 0) then
        call check( nf90_def_var(ncid, "snowmelt" , NF90_REAL, dimids_2d, snowmelt_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,snowmelt_ID,0,4,4))
      end if

      if(ZSP_out .gt. 0)then
        call check( nf90_def_var(ncid, "ZSP2d", NF90_REAL, dimids_2d, ZSPDT_ID))
        if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,ZSPDT_ID,0,4,4))
      endif

      call check( nf90_enddef(ncid) )

      !call check( nf90_def_var(ncid, "FDEPTH" , NF90_REAL, dimids_3d, FDEPTH_ID))
      !if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,FDEPTH_ID,0,4,4))
      !call check( nf90_def_var(ncid, "TDEPTH" , NF90_REAL, dimids_3d, TDEPTH_ID))
      !if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,TDEPTH_ID,0,4,4))
      !call check( nf90_def_var(ncid, "absorbedSW" , NF90_REAL, dimids_3d, absorbedSW_ID))
      !if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,absorbedSW_ID,0,4,4))
      !call check( nf90_def_var(ncid, "absorbedLW" , NF90_REAL, dimids_3d, absorbedLW_ID))
      !if(deflate.eq.1)call check( nf90_def_var_deflate(ncid,absorbedLW_ID,0,4,4))
      !call check( nf90_def_var(ncid, "InDirect" , NF90_REAL, dimids_3d, InDirect2d_ID))
      !if(deflate.eq.1) call check( nf90_def_var_deflate(ncid,InDirect2d_ID,0,4,4))
      !call check( nf90_def_var(ncid, "InDiffuse" , NF90_REAL, dimids_3d, InDiffuse2d_ID))
      !if(deflate.eq.1) call check( nf90_def_var_deflate(ncid,InDiffuse2d_ID,0,4,4))




!   Store variable values
      call check( nf90_inq_varid(ncid, "InvalidPoints", invalidp_ID) )
      call check( nf90_put_var(ncid, invalidp_ID, maskgrid))

      if(T_soil_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "TSDT2d", TSDT_ID) )
        call check( nf90_put_var(ncid, TSDT_ID, TSDT2d_daily))
      end if

      if(T_Residue_out.gt.0)then
        call check( nf90_inq_varid(ncid, "TRDT2d", TRDT_ID))
        call check( nf90_put_var(ncid, TRDT_ID, TRDT2d_daily))
      end if

      if(VLC_soil_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "VLCDT2d", VLCDT_ID) )
        call check( nf90_put_var(ncid, VLCDT_ID, VLCDT2d_daily))
      end if

      if(VIC_soil_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "VICDT2d", VICDT_ID) )
        call check( nf90_put_var(ncid, VICDT_ID, VICDT2d_daily))
      end if

      if(runoff_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "RUNOFFDIS", RUNOFFDIS_ID) )
        call check( nf90_put_var(ncid, RUNOFFDIS_ID, RUNOFFDIS_daily))
      end if

      if(surfacerunoff_out .gt. 0) then
        call check( nf90_inq_varid(ncid, "surfacerunoff", runoff_ID) )
        call check( nf90_put_var(ncid, runoff_ID, surfacerunoff_daily))
      end if

      if (EVAP_out .gt. 0) then
         call check( nf90_inq_varid(ncid, "EVAP", evap_ID) )
         call check( nf90_put_var(ncid, evap_ID, EVAP12d_daily))
      end if

      if (ET_out .gt. 0) then
        call check( nf90_inq_varid(ncid, "ET", etsum_ID) )
        call check( nf90_put_var(ncid, etsum_ID,ETSUM2d_daily))
      end if 


      if(DGL_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "DGL2d", dgl_ID) )
        call check( nf90_put_var(ncid, dgl_ID, DGL2d_daily))      
      end if

      if(THFLUX_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "THFLUX", THFLUX_ID) )
        call check( nf90_put_var(ncid, THFLUX_ID,THFLUX2d_daily))
      end if

      if(SWE_out .gt. 0)then
        !call check( nf90_inq_varid(ncid, "SWE", SWE_ID) )
        !call check( nf90_put_var(ncid, SWE_ID,SWE2d_daily))
      end if

      if(ZSP_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "ZSP2d", ZSPDT_ID) )
        call check( nf90_put_var(ncid, ZSPDT_ID, ZSP_daily))
      endif

      if(TSP_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "TSP2d", TSP_ID) )
        call check( nf90_put_var(ncid, TSP_ID, TSP_daily))
      endif

      if(infill_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "infill", infil2d_ID) )
        call check( nf90_put_var(ncid, infil2d_ID, infill_daily))
      endif

      if(snowmelt_out .gt. 0)then
        call check( nf90_inq_varid(ncid, "snowmelt", snowmelt_ID))
        call check( nf90_put_var(ncid, snowmelt_ID, snowmelt_daily))
      endif

    call check( nf90_close(ncid) )
 end subroutine




