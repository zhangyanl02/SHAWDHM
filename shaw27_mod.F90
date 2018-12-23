module shaw27_mod
  use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
  use dims_mod
  implicit none
contains

!***********************************************************************
!
 subroutine rainsl (ns,train,rain,zs,tsdt,vlcdt,vicdt,icesdt,&
    matdt,totflo,saltdt,concdt,pond,pondmx,nsalt,col,row,infiltration)
!
!     this subroutine calculates the infiltration into the soil as well
!     as the final soil temperature and solute concentration resulting
!     from the heat and salts carried by the water.
!     ---- currently the subroutine will not allow the water content to
!          exceed 90% of the saturated value (psat=0.9)
!***********************************************************************
    use statevar_mod,only:maskgrid
    use controlpara_mod,only: wt2d,dt2d,wdt2d
    use dims_mod, only:nsmax,nsaltmax
    use constvar_mod
    use soilproperty_mod, only:rhob2d,saltkq2d,entry2d,sat2d,satk2d
    implicit none

!input
    integer(i4),intent(in)::ns,col,row,nsalt
    integer(i4),dimension(nsmax),intent(in):: icesdt
    real(r8),intent(inout)::train,rain,pond,infiltration
    real(r8),intent(in)::pondmx
    real(r8),dimension(nsmax),intent(in)::zs
    real(r8),dimension(nsmax),intent(inout)::vlcdt,vicdt,matdt,totflo,tsdt
    real(r8),dimension(nsaltmax,nsmax),intent(inout)::saltdt,concdt


!temp variable
    real(r8),dimension(nsmax)::satkmx,matric,cs,psat,vlc9,mat9,satk9
    real(r8)::infmax,infil,vlcmax,dtime,zstar,sumz,sumzk,dtstar,potinf,diminf,sumz2,dz,avheat, &
            cleach,oldvlc,rainht,rqheat,satinf,sheat,skqvlc,sloss,sltnew,time,tlconc,tmpfrz,totpot, &
            tsoil,vicmax,oldvic
    real(r8),parameter::psatk=0.9
    integer(i4)::i,j

!pointer
    real(r8),dimension(:),pointer::rhob,entry,sat,satk
    real(r8),dimension(:,:),pointer::saltkq
    real(r8),pointer::dt,wdt,wt
    real(r8)::temp

!   added for test
    real(r8),dimension(nsmax)::vlctmp,victmp
    real(r8)::poros,fract


    data psat/nsmax*0.95/


    rhob  =>rhob2d(col,row,:)
    entry => entry2d(col,row,:)
    saltkq=>saltkq2d(col,row,:,:)
    sat   =>sat2d(col,row,:)
    satk  =>satk2d(col,row,:)
    dt    =>dt2d(col,row)
    wdt   =>wdt2d(col,row)
    wt    =>wt2d(col,row)


!     determine the maximum conductivity for each node -- this will be
!     defined by the air entry potential unless the soil contains ice
      do 10 i=1,ns
         if (icesdt(i) .ge. 1) then
!           saturated conductivity is limited by the presence of ice
            vlcmax= sat(i) - vicdt(i)
            if (vlcdt(i) .gt. vlcmax) vlcmax=vlcdt(i)
            call matvl1 (i,matric(i),vlcmax,col,row)
            vlc9(i)=psatk*vlcmax
           else
!           set matric potential for computing saturated conductivity 
!           and initialize 90% (or psatk) of maximum water conent
            matric(i)=entry(i)
            vlc9(i)=psatk*sat(i)
         end if
         call matvl1 (i,mat9(i),vlc9(i),col,row)
   10 continue

      call soilhk (ns,satkmx,matric,vlcdt,vicdt,col,row)


!     initialize infiltration parameters
      dtime=dt
      infil=0.0
      zstar=0.0
      sumz=0.0
      sumzk=0.0
      psat(1)=psatk
      i = 1
!     check if conductivity of surface is zero (possibly due to ice)
      if (satkmx(1) .le. 0) go to 40
!     check if wetting front has reached saturated conditions
      !if (matdt(1) .ge. 0.0) go to 40
!     calculate the maximum infiltration that the first layer can hold
      infmax=(psat(1)*sat(1) - vlcdt(1) - vicdt(1))*(zs(2) - zs(1))/2.
!     if current layer is already saturated -- go to next layer
      if (infmax .le. 0.0) go to 30
!     calculate the potential infiltration for the layer
   20 if (zstar .le. 1.0) then
!        assumption of saturated flow behind wetting front is valid
         if(abs(psat(i)*sat(i) - vlcdt(i) - vicdt(i)).lt. 1.0e-10 .or. abs(-matdt(i)+sumz) .lt. 1.0e-10)then
            print*,"error in rainsl, division by zero",psat(i)*sat(i) - vlcdt(i) - vicdt(i),-matdt(i)+sumz
            print*,col,row
            potinf=0.0
         else
           dtstar= satkmx(i)*dtime/(psat(i)*sat(i) - vlcdt(i) - vicdt(i))/(-matdt(i)+sumz)
           temp=(dtstar-2*zstar)**2 +8*dtstar
           if (temp.le.0) then
             temp=0.000001
             maskgrid(col,row)=1
           end if
           if(temp.lt.0.0) then
             print*,"in subroutine rainsl, the temp is less or equal 0.0"
             print*,col,row,temp
             stop
           end if
           potinf=(dtstar-2*zstar + sqrt(temp))/2
           potinf=potinf*(psat(i)*sat(i)-vlcdt(i)-vicdt(i))*(-matdt(i)+sumz)
         end if
       else
!        infiltration flow has become unsaturated -- use flow through 
!        saturated medium 
         potinf=satinf*dtime
      end if
      if (infmax .ge. (rain-infil)) then
!        layer can hold remainder of water - check if all can infiltrate
!
         if (potinf .ge. (rain-infil)) then
!           layer can infiltrate remainder of rain
            infil=rain
           else
!           profile cannot infiltrate all of the rain
            infil=infil+potinf
         end if
!        infiltration calculation is complete - go to energy calculation
         go to 40
      end if
!     layer cannot hold remainder of rain - check if it can hold potinf
      if (potinf .lt. infmax) then
!        profile cannot infiltrate remainder of the rain
         infil=infil+potinf
!        infiltration calculation complete - go to energy calculation
         go to 40
      end if
!     potential infiltration and rain are sufficient to fill layer
!      --- calculate time required to fill the layer and adjust time
      if (zstar .le. 1.0) then
!        assumption of saturated flow behind wetting front is valid     
         if(abs(psat(i)*sat(i) - vlcdt(i) - vicdt(i)).lt. 1.0e-10 .or. abs(-matdt(i)+sumz) .lt. 1.0e-10)then
            print*,"error in rainsl, division by zero",psat(i)*sat(i) - vlcdt(i) - vicdt(i),-matdt(i)+sumz
            print*,col,row
            stop
         end if
         diminf=infmax/(psat(i)*sat(i)-vlcdt(i)-vicdt(i))/(-matdt(i)+sumz)
         if (diminf .gt. 0.01) then
            dtstar= (zstar-1)*log(1. + diminf) + diminf
           else
!           if diminf is small, use alternate method of calculating log
!           term because serious errors may result when taking the log 
!           of a number very close to 1.0, i.e. (1. + 10^-5)
            dtstar= (zstar-1)*2*diminf/(2+diminf) + diminf
         end if
         time= dtstar*(psat(i)*sat(i)-vlcdt(i)-vicdt(i))*(-matdt(i)+sumz)/satkmx(i)
       else
!        infiltration flow has become unsaturated -- use flow through saturated medium
         if(abs(satinf).lt. 1.0E-15)then
             print*,"error in rainsl, division by zero",satinf
             print*,col,row
             stop
         end if
         time=infmax/satinf
      end if
      dtime= dtime-time
      infil=infil+infmax
      if (dtime .le. 0.0) go to 40
!     update parameters for next layer
   30 i = i + 1
!     check if hydraulic conductivity is zero (possibly due to ice)
      if (satkmx(i) .le. 0) go to 40
!     check if wetting front has reached saturated conditions
      if (matdt(i) .gt. sumz) go to 40
      sumz2= (zs(i)+zs(i-1))/2.
      dz= sumz2-sumz
      sumz= sumz2
      sumzk= sumzk + dz/satkmx(i-1)
      if (i.lt.ns) then
         if (zstar .le. 1.0) then
!           saturated flow still exists -- setup zstar for next layer
            zstar=satkmx(i)*sumzk/(-matdt(i) + sumz)
            psat(i)=psatk
            if (zstar .gt. 1.0) then
!              hydraulic conductivity of this layer is sufficiently high
!              that saturated conditions no longer exist -- calculate 
!              infiltration using saturated relationships behind
!              current layer assuming gravity flow.
               satinf=(-matdt(i) + sumz)/sumzk
!              compute conductivity of all layers at psatk (currently 
!              set at 0.90) of maximum water conent for interpolation
               call soilhk (ns,satk9,mat9,vlc9,vicdt,col,row)
!              determine water content of layer such that conductivity 
!              matches flow into layer - use logarithmic interpolation
!              between known points of vlc and conductivity
               psat(i)=10.**(log10(psatk)*log10(satinf/satkmx(i))/log10(satk9(i)/satkmx(i)))
               if (psat(i).gt. psatk) psat(i)=psatk
!              define maximum saturated conductivity for current layer
               satkmx(i)=sumz/sumzk
            end if 
           else
!           saturated conditions no longer exist at wetting front;
!           compare flux to saturated conductivity of current layer
            satinf=(-matdt(i) + sumz)/sumzk
            if (satkmx(i) .lt. satinf) satinf=satkmx(i)
!           determine water content of layer such that conductivity 
!           matches flow into layer - use logarithmic interpolation
!           between known points of vlc and conductivity
            psat(i)=10.**(log10(psatk)*log10(satinf/satkmx(i))/log10(satk9(i)/satkmx(i)))
            if (psat(i).gt. psatk) psat(i)=psatk
!           define maximum saturated conductivity for current layer
            satkmx(i)=sumz/sumzk
         end if
!        calculate the maximum infiltration that layer can hold
         infmax=(psat(i)*sat(i)-vlcdt(i)-vicdt(i))*(zs(i+1)-zs(i-1))/2.
!        if layer cannot hold additional water, go to next layer
         if (infmax .le. 0.0) go to 30
         go to 20
        else
!        infiltration has reached bottom of profile
         if (zstar .le. 1.0) then
!           compute flow through saturated profile      
            infil=infil + dtime*sumz/sumzk
           else
!           compare flux to saturated conducivity of current layer
            if (satkmx(i) .lt. satinf) satinf=satkmx(i)
            infil=infil + dtime*satinf
         end if
         if (infil .gt. rain) infil=rain
      end if
!-----------------------------------------------------------------------
!     infiltration calculation complete -- calculate water to be ponded,
!     then determine temperature, moisture contents and solutes.

   40 pond= rain-infil
      infiltration=infil
      rain= 0.0
      if (pond .gt. pondmx) then
!        set rain equal to the excess rain, i.e. the runoff
         rain= pond-pondmx
         pond= pondmx
      end if
      if (infil.le.0.0) then
!        update of water contents, temperatures, and solutes complete
         infmax=0.0
         go to 110
      end if
!     calculate the specific heat of the layers prior to adding water
      call soilht (ns,cs,vlcdt,vicdt,tsdt,matdt,concdt,nsalt,col,row)
!     define conditions for surface node
!     sloss = salts lost or leached from nodes above
      i = 1
      sloss = 0.0
      dz=(zs(2)-zs(1))/2.
!     start energy and solute calculations for each node
!     calculate the water infiltrated into the current node
   50 infmax=(psat(i)*sat(i) - vlcdt(i) - vicdt(i))*dz
      if (infmax .lt. 0.0) infmax=0.0
      if (infmax .ge. infil) infmax=infil
!     adjust the solute concentration in node to account for leaching
      do 60 j=1,nsalt
!        calculate solute concentration of water coming into node
!        and adjust salt concentration in node for the water
!        required to get to maximum water content
         if (infil .le. 0.0) then
            cleach = 0.0
          else
            cleach = sloss/rhol/infil
         end if
         saltdt(j,i) = saltdt(j,i) + cleach*rhol*infmax/rhob(i)
!        calculate salt present in node after leaching
         skqvlc= saltkq(j,i) + (vlcdt(i)+infmax/dz)*rhol/rhob(i)
         sltnew= (cleach*rhol*(infil-infmax)/dz +saltdt(j,i)*(rhob(i)-wt*rhol*(infil-infmax)/dz/skqvlc))&
              /(rhob(i)+wdt*rhol*(infil-infmax)/dz/skqvlc)
         if (sltnew .lt. 0.0) sltnew=0.0
         sloss= (saltdt(j,i) - sltnew)*rhob(i)*dz
         saltdt(j,i)= sltnew
   60 continue
!     calculate the heat capacity of the rain and of the soil
!     --rain heat capacity is based on total water absorbed and passing
!     through current node (temp water leaving node is soil temp)
      sheat= cs(i)*dz
      rainht= infil*rhol*cl
      if (icesdt(i) .le. 0) then
!        layer is not frozen -- calculate the temperature directly and
!        adjust moisture content, matric potential, and solutes
         tsdt(i)= tsdt(i) + rainht/(rainht+sheat)*(train-tsdt(i))
         vlcdt(i)= vlcdt(i) + infmax/dz
!        do not allow adjustment of saturated potentials
         if (matdt(i).lt.entry(i)) &
            call matvl1 (i,matdt(i),vlcdt(i),col,row)
         do 70 j=1,nsalt
            concdt(j,i)=saltdt(j,i)/(saltkq(j,i)+vlcdt(i)*rhol/rhob(i))
   70    continue
       else
!        layer contains ice -- must account for latent heat of fusion.
!        calculate temperature at which layer will be completely thawed
         oldvlc= vlcdt(i)
         vlcdt(i)= vlcdt(i) + vicdt(i)*rhoi/rhol + infmax/dz
         call matvl1 (i,matdt(i),vlcdt(i),col,row)
         tlconc=0.0
         do 80 j=1,nsalt
            concdt(j,i)= saltdt(j,i)/(saltkq(j,i)+vlcdt(i)*rhol/rhob(i))
            tlconc=tlconc + concdt(j,i)
   80    continue
         totpot= matdt(i) - tlconc*ugas*273.16/g
         tmpfrz= 273.16*totpot/(lf/g-totpot)
!        check if heat of rain is sufficient to melt all of the ice
         avheat= rainht*(train-tmpfrz)
         rqheat= rhoi*lf*vicdt(i)*dz + sheat*(tmpfrz-tsdt(i))
         if (avheat .ge. rqheat) then
!           available heat is more than required to melt the ice.
!           calculate temperature of soil
            tsdt(i)= tmpfrz + (avheat-rqheat)/(sheat+rainht)
            vicdt(i)= 0.0
!           done with this layer -- go to next layer
            go to 110
         end if
!        layer remains frozen - assume temperature is at freezing point
!        and calculate the ice content from an energy balance:
!        (latent heat)*(delta ice) =
!                        rainht*(train-t(new)) + soilht(t(new) - tsdt)
         tsoil= tsdt(i)
         oldvic= vicdt(i)
         vicmax= vicdt(i) + infmax/dz*(rhol/rhoi)
         vicdt(i)= oldvic - (rainht*(train-tmpfrz)-sheat*(tmpfrz-tsoil))/(rhoi*lf*dz)
         if (vicdt(i) .gt. vicmax) vicdt(i) = vicmax
         vlcdt(i)= oldvlc + (vicmax-vicdt(i))*rhoi/rhol
!        determine matric, solutes and temp for this water content
         call matvl1 (i,matdt(i),vlcdt(i),col,row)
         tlconc=0.0
         do 90 j=1,nsalt
            concdt(j,i)= saltdt(j,i)/(saltkq(j,i)+vlcdt(i)*rhol/rhob(i))
            tlconc=tlconc + concdt(j,i)
   90    continue
         totpot= matdt(i) - tlconc*ugas*273.16/g
         tsdt(i)= 273.16*totpot/(lf/g-totpot)
!        calculate ice content using updated temperature
         vicdt(i)= oldvic-(rainht*(train-tsdt(i))-sheat*(tsdt(i)-tsoil))/(rhoi*lf*dz)
         if (vicdt(i) .gt. vicmax) vicdt(i) = vicmax
         vlcdt(i)= oldvlc + (vicmax-vicdt(i))*rhoi/rhol
!        determine matric, solutes and temp for this water content
         call matvl1 (i,matdt(i),vlcdt(i),col,row)
         tlconc=0.0
         do 100 j=1,nsalt
            concdt(j,i)= saltdt(j,i)/(saltkq(j,i)+vlcdt(i)*rhol/rhob(i))
            tlconc=tlconc + concdt(j,i)
  100    continue
         totpot= matdt(i) - tlconc*ugas*273.16/g
         tsdt(i)= 273.16*(lf/g/(lf/g-totpot) - 1.0)
      end if
!     calculate conditions of next layer
  110 i = i + 1
      infil=infil-infmax
!     add infiltration between nodes to total flow between nodes
      totflo(i-1)=totflo(i-1)+infil
      if (i .lt. ns) then
         dz=(zs(i+1) - zs(i-1))/2.
         train = tsdt(i-1)
         if (infil .gt. 0.0) go to 50
      end if
      return
  end subroutine rainsl







!***********************************************************************
!
 subroutine rainrs (nr,train,rain,zr,trdt,gmcdt,rhor, &
 gmcmax,dirres,slope)
!
!     this subroutine calculates the rainfall depth intercepted by the
!     residue layers, and adjusts the temperature to account for the
!     heat (or lack thereof) introduced by the rainfall
!***********************************************************************
    use constvar_mod, only:lf,ls,lv,cl,rhol
    use dims_mod, only:nrmax
    implicit none
! input
    real(r8),dimension(nrmax),intent(in):: zr,rhor
    real(r8),intent(in)::gmcmax,dirres,slope,train
    integer(i4),intent(in)::nr
    real(r8),intent(inout)::rain
    real(r8),dimension(nrmax),intent(inout)::gmcdt,trdt

!temp
    real(r8),dimension(nrmax)::tdirec,tdiffu,cres
    real(r8)::dz,angle,excess,rainht,resiht,rinter
    integer(i4)::i

!
!     use the transmittance to direct radiation to calculate the
!     fraction of rain intercepted by the residue
      angle=1.5708-slope
      call transr (nr,tdirec,tdiffu,zr,rhor,angle,dirres)
!
!     calculate the specific heat of the residue before rain is added
      call resht (nr,cres,gmcdt,rhor)
!
      do i=1,nr
!
         if (i .eq. 1) then
            dz=zr(2)-zr(1)
            if (nr .ne. 1) dz=dz/2.
           else
            if (i .eq. nr) then
               dz=zr(i+1)-zr(i)+(zr(i)-zr(i-1))/2.
              else
               dz=(zr(i+1)-zr(i-1))/2.
            end if
         end if

!****    calculate rainfall intercepted by the node
         rinter=rain*(1.-tdirec(i))
         gmcdt(i)=gmcdt(i) + rinter*rhol/(dz*rhor(1))
         rain=rain-rinter
!        check if the node can hold this much water
         if (gmcdt(i) .gt. gmcmax) then
            excess= rhor(i)*(gmcdt(i)-gmcmax)*dz/rhol
            rinter= rinter-excess
            rain= rain+excess
            gmcdt(i)= gmcmax
         end if

!****    calculate the temperature of the node
         rainht=rinter*rhol*cl
         resiht=cres(i)*dz
         trdt(i)=trdt(i) + rainht/(rainht+resiht)*(train-trdt(i))
      end do
    end subroutine rainrs



!***********************************************************************
!
SUBROUTINE META (NSP,TSPDT,RHOSP,DZSP,DLWDT,col,row)
!
!     COMPUTES THE CHANGE IN DENSITY OF THE SNOW COVER CAUSED BY
!     DESTRUCTIVE (EQUI-TEMPERATURE) METAMORPHISM, COMPACTION, AND THE
!     PRESENCE OF LIQUID-WATER.
!***********************************************************************
    use controlpara_mod, only:dt2d
    use constvar_mod
    use spwatr_mod
    implicit none

!input
    integer(i4),intent(in)::NSP,col,row
    real(r8),dimension(NSPMAX),intent(in)::TSPDT,DLWDT
    real(r8),dimension(NSPMAX),intent(inout)::RHOSP,DZSP

!temp
    integer(i4)::I
    real(r8)::WEIGHT,T1,T2,TERM,WEL
!pointer
    real(r8),pointer::DT
    DT => DT2d(col,row)

!
    IF (CMET1 .LE. 0.0  .OR.  CMET3 .LE. 0.0) THEN
         IF (CMET5 .LE. 0.0) RETURN
    END IF
!
!     WEIGHT IS THE WATER-EQUIVALENT (CM) ABOVE THE LAYER.
      WEIGHT=0.0
!
      DO 10 I=1,NSP
!        IF DENSITY IS THAT OF ICE, DO NOT INCREASE IT ANY MORE
         IF (RHOSP(I) .GE. RHOI) GO TO 10
         WEL= RHOSP(I)*DZSP(I)/RHOL
!
!****    DESTRUCTIVE METAMORPHISM TERM.
!
         TERM=EXP(CMET4*TSPDT(I))
         IF (RHOSP(I) .LE. SNOMAX) THEN
            T1=TERM*CMET3
           ELSE
            T1=TERM*CMET3*EXP(-46.0*(RHOSP(I)-SNOMAX)/RHOL)
         END IF
!
!****    COMPACTION TERM.
!
         TERM=EXP(0.08*TSPDT(I))
         T2=WEIGHT*CMET1*EXP(-CMET2*RHOSP(I)/RHOL)*TERM
!
!***     LIQUID-WATER TERM.
!
         IF (DLWDT(I) .GT. 0.0) T1=CMET5*T1
!
!****    DENSIFICATION OF THE LAYER.  (WATER-EQUIVALENT STAYS THE SAME.)
!
         RHOSP(I)=RHOSP(I)*(1.0 + DT*(T1+T2)/3600.)
         DZSP(I)=RHOL*WEL/RHOSP(I)
         WEIGHT=WEIGHT + (WEL+DLWDT(I))*100.
   10 CONTINUE
      RETURN
END SUBROUTINE 






!***********************************************************************
!
SUBROUTINE SNOMLT (NSP,ICESPT,DZSP,RHOSP,TSPDT,DLWDT,WLAG,&
                        STORE,SCOUT,col,row)
!
!     THIS SUBROUTINE DETERMINES THE SNOW COVER OUTFLOW DURING EACH
!     PERIOD BASED ON THE CHANGE IN LIQUID-WATER DURING THE PERIOD AND
!     THE PHYSICAL CHARTERISTICS OF THE SNOW COVER.
!***********************************************************************
    use controlpara_mod, only:dt2d
    use constvar_mod
    use spwatr_mod
    use savedata_mod,only:IFIRST_snowm2d,NLAG2d
    implicit none

!input
    integer(i4),intent(in)::NSP,col,row
    integer(i4),dimension(NSPMAX),intent(inout)::ICESPT

    real(r8),dimension(NSPMAX),intent(in)::DZSP
    real(r8),dimension(NSPMAX),intent(inout)::RHOSP,TSPDT,DLWDT
    real(r8),dimension(11),intent(inout)::WLAG
    real(r8),intent(inout)::STORE,SCOUT

!temp
    integer(i4)::I,J,K,IDT,NI
    real(r8)::EXCESS,BOTTOM,DENSE,EXCSHR,FJ,FLAG,FLMAX,FN,FREEZE,HOURS,JDT,OUTHR,PLW, &
            POR,R,TDEPTH,W,WE,WEL,WINC,WMAX
    integer(i4)::IHR
! pointer
    integer(i4),pointer::IFIRST,NLAG
    real(r8),pointer::DT

    IFIRST=>IFIRST_snowm2d(col,row)
    NLAG => NLAG2d(col,row)
    DT => DT2d(col,row)




!
      IF(IFIRST .EQ. 0) THEN
!        INITIALIZE VARIABLES.
!        1. STORE IS THE AMOUNT OF LIQUID-WATER(THAT HAS ALREADY BEEN
!           LAGGED) THAT IS IN STORAGE IN THE SNOW-COVER.
!        2. WLAG() IS THE AMOUNT OF LIQUID-WATER IN THE PROCESS OF
!           BEING LAGGED.  NLAG IS THE NUMBER OF ARRAY ELEMENTS USED.
         JDT=CLAG1+0.01
         NLAG=JDT+2
         IF (NLAG.GT.11) NLAG=11
         IFIRST=1
      END IF
!
      IF (NSP .LE. 0) THEN
!        SNOW COVER HAS JUST DISAPPEARED
         DO 10 J=1,NLAG
            SCOUT=SCOUT+WLAG(J)
            WLAG(J)=0.0
   10    CONTINUE
         SCOUT=SCOUT+STORE
         STORE=0.0
         RETURN
      END IF
!
!     DETERMINE THE EXCESS LIQUID-WATER (IN EXCESS OF LIQUID-WATER
!     HOLDING CAPACITY) GENERATED DURING THIS TIME PERIOD.
      EXCESS=0.0
!
!     EXCESS WATER IN BOTTOM LAYER IS NOT TO BE ROUTED.
      IF (NSP .EQ. 1) GO TO 12
      BOTTOM=0.0
      IF (ICESPT(NSP) .EQ. 1) THEN
         WEL=RHOSP(NSP)*DZSP(NSP)/RHOL
         PLW=(PLWMAX-PLWHC)*(PLWDEN-RHOSP(NSP))/PLWDEN+PLWHC
         IF (RHOSP(NSP) .GE. PLWDEN) PLW=PLWHC
         WMAX=PLW*WEL
         IF (DLWDT(NSP) .GT. WMAX) BOTTOM= DLWDT(NSP) - WMAX
      END IF
!
   12 WE=0.0
      TDEPTH=0.0
      DO 20 I=1,NSP
         WEL=RHOSP(I)*DZSP(I)/RHOL
         IF(ICESPT(I) .EQ. 0) THEN
!           LAYER IS BELOW ZERO DEGREES CELSIUS.
            IF(EXCESS .EQ. 0.0) GO TO 15
!           FREEZE SOME OF THE LIQUID-WATER.
            FREEZE=-TSPDT(I)*CI*WEL/LF
            IF (EXCESS .LE. FREEZE) THEN
!              EXCESS IS ALL FROZEN IN THIS LAYER.
               TSPDT(I)=TSPDT(I)+(EXCESS*LF*RHOL)/(CI*WEL*RHOL)
               WEL=WEL+EXCESS
               EXCESS=0.0
               RHOSP(I)=RHOL*WEL/DZSP(I)
               GO TO 15
              ELSE
!              EXCESS EXCEEDS REFREEZE.
               TSPDT(I)= 0.0
               WEL=WEL+FREEZE
               RHOSP(I)=RHOL*WEL/DZSP(I)
               EXCESS=EXCESS-FREEZE
               ICESPT(I)=1
            END IF
         END IF
         PLW=(PLWMAX-PLWHC)*(PLWDEN-RHOSP(I))/PLWDEN+PLWHC
         IF (RHOSP(I) .GE. PLWDEN) PLW=PLWHC
         WMAX= PLW*WEL
         W= DLWDT(I) + EXCESS
         IF (W .LE. WMAX) THEN
!           LIQUID-WATER HOLDING CAPACITY IS NOT SATISFIED.
            DLWDT(I)=W
            EXCESS=0.0
           ELSE
!           LIQUID-WATER HOLDING CAPACITY IS EXCEEDED.
            DLWDT(I)=WMAX
            EXCESS=W-WMAX
         END IF
   15    WE=WE+WEL
         TDEPTH=TDEPTH+DZSP(I)
   20 CONTINUE
!
      IF (NSP .EQ. 1) THEN
!        WATER NOT LAGGED IF ONLY ONE NODE IN THE SNOWPACK
         SCOUT = EXCESS
         DO 50 J=1,NLAG
            SCOUT=SCOUT+WLAG(J)
            WLAG(J)=0.0
   50    CONTINUE
         SCOUT=SCOUT+STORE
         STORE=0.0
         RETURN
      END IF
!
      EXCESS=EXCESS-BOTTOM
      IF (ICESPT(NSP) .EQ. 0) THEN
!        DO NOT ALLOW ANY WATER LAGGED OR STORED TO LEAVE SNOWPACK
!        (EXCESS MUST BE EQUAL 0.0)
         OUTHR=0.0
         SCOUT=0.0
         RETURN
      END IF
!
!     ROUTE EXCESS WATER THROUGH THE SNOW COVER.
!     EMPIRICAL LAG AND ATTENUATION EQUATIONS - ONE HOUR TIME STEP USED.
      EXCSHR = EXCESS*3600./DT
      IDT=DT/3600.
      IF (MOD(DT,3600.).GT.0.00001) IDT=IDT+1
      IF (IDT .EQ. 0) IDT=1
      HOURS=1.0
!
      DENSE = WE/TDEPTH
      SCOUT=0.0
      DO 35 IHR=1,IDT
!        IF THIS IS THE LAST TIME INCREMENT, ACCOUNT FOR ANY REMAINING
!        FRACTION OF AN HOUR IN THE TIME STEP
         IF (IHR.EQ.IDT) HOURS=DT/3600-(IDT-1)
         EXCESS=EXCSHR*HOURS
!
         OUTHR=0.0
!        LAG-FUNCTION OF DEPTH,DENSITY,AND EXCESS WATER.
         IF (EXCSHR .GE. 0.00001) THEN
            NI= ((EXCSHR*10000.0)**0.3)+0.5
            IF (NI .LT. 1) NI=1
            FN=NI
            FLMAX=CLAG1*(1.0-EXP(-0.25*TDEPTH/DENSE))
            DO 25 J=1,NI
               FJ=J
               FLAG=FLMAX/(CLAG2*EXCSHR*100.*(FJ-0.5)/FN+1.0)
               K= FLAG+1.0
               POR=K-FLAG
               WINC=1.0/FN
               WLAG(K)=WLAG(K)+EXCESS*WINC*POR
               WLAG(K+1)=WLAG(K+1)+EXCESS*WINC*(1.0-POR)
   25       CONTINUE
           ELSE
            WLAG(1)=WLAG(1)+EXCESS
         END IF
!
!        ATTENUATION-FUNCTION OF DENSITY AND PREVIOUS OUTFLOW.
         IF ((STORE+WLAG(1)) .NE. 0.0) THEN
            R=1.0/(CLAG3*EXP(-CLAG4*WLAG(1)*DENSE/TDEPTH)+1.0)
            OUTHR=(STORE+WLAG(1))*R
            STORE=STORE+(WLAG(1)-OUTHR)*HOURS
            SCOUT=SCOUT+OUTHR*HOURS
            IF(STORE .LE. 0.00001) THEN
               OUTHR=OUTHR+STORE
               SCOUT=SCOUT+STORE
               STORE=0.0
            END IF
         END IF
         NI=NLAG-1
         DO 30 J=1,NI
            WLAG(J)=WLAG(J)*(1.-HOURS) + WLAG(J+1)*HOURS
            IF (WLAG(J+1).LE. 0.000001) THEN
!              TIME STEPS OF LESS THAN 1 HOUR RESULT IN CONTINUALLY
!              TAKING FRACTIONS OF THE LAGGED DEPTH; DEPTH IS TOO
!              SMALL -- PASS THE ENTIRE DEPTH ONTO THE NEXT LAG
               WLAG(J)=WLAG(J)+WLAG(J+1) + WLAG(J+1)*(1.-HOURS)
               WLAG(J+1)=0.0
            END IF
   30    CONTINUE
         WLAG(NLAG)=0.0
   35 CONTINUE
!
      SCOUT=SCOUT+BOTTOM
      RETURN 
  END SUBROUTINE SNOMLT


!***********************************************************************
SUBROUTINE WBSNOW (NSP,ICESPT,ZSP,DZSP,RHOSP,TSPDT,DLW,DLWDT,RAIN, &
                 TRAIN,TOPSNO,EVAP1,TQVSP,WLAG,STORE,SCOUT,col,row)
!
!     THIS SUBROUTINE PERFORMS A WATER BALANCE OF THE SNOWPACK BY
!     ADSUSTING THE DENSITY FOR VAPOR FLUX AND ANY MELT WHICH OCCURRED
!     OVER THE TIME STEP.  CHECKS ARE MADE TO SEE IF ANY LAYERS HAVE
!     DISAPPEARED DUE TO MELT, OR IF ANY LAYERS ARE OUTSIDE THE
!     ACCEPTABLE THICKNESS RANGE.
!
!***********************************************************************
    use constvar_mod,only:LF,LS,LV,CI,CL,RHOL
    use spwatr_mod
    implicit none


!input
    integer(i4),intent(in)::col,row
    integer(i4),intent(inout)::NSP
    integer(i4),dimension(NSPMAX),intent(inout)::ICESPT
    real(r8),dimension(NSPMAX),intent(in)::DLW,TQVSP
    real(r8),dimension(NSPMAX),intent(inout)::DZSP,RHOSP,TSPDT,DLWDT,ZSP

    real(r8),dimension(11),intent(inout)::WLAG
    real(r8),intent(in)::RAIN,TRAIN
    real(r8),intent(inout)::SCOUT,TOPSNO,EVAP1,STORE

!temp
    real(r8)::CHANGE,ABOVE,EXCESS,WEL,FREEZE,DZ1,TDEPTH,WE,WENL,ZZ,WEI
    integer(i4)::I,NL,J,L,LL,NEXT




!
      SCOUT = 0.0
!
!     ADJUST THE DENSITY FOR VAPOR FLUX
!     (THICKNESS OF LAYER IS ADJUSTED FOR SURFACE LAYER)
!
      CHANGE = (TOPSNO - TQVSP(1))/RHOSP(1)
      DZSP(1) = DZSP(1) + CHANGE
      IF (DZSP(1).LE.0.0 .AND. NSP.EQ.1) THEN
!        SNOWPACK IS LOST TO SUBLIMATION -- ADJUST TOPSNO AND EVAP1 FOR
!        THE WATER BALANCE TO CLOSE AND RETURN TO CALLING SUBROUTINE
         TOPSNO=TOPSNO - DZSP(1)*RHOSP(1)
         EVAP1=EVAP1 - DZSP(1)*RHOSP(1)/RHOL
         NSP=0
         RETURN
      END IF
      DO 10 I=2,NSP
         CHANGE = (TQVSP(I-1) - TQVSP(I))/DZSP(I)
         RHOSP(I) = RHOSP(I) + CHANGE
   10 CONTINUE
!
!     ADJUST LAYERS FOR MELT AND RAINFALL
!        ABOVE = THE RAIN OR THE MELTWATER FROM ABOVE LAYERS
!        EXCESS= THIS IS THE MELTWATER THAWBSNOWT WOULD HAVE BEEN PRODUCED
!                BY ABOVE LAYERS HAD THERE BEEN ENOUGH ICE IN THE LAYER
!                FOR ALL THE ENERGY ABSORBED.  (FOR THE FIRST LAYER,
!                THIS TERM INCLUDES THE ENERGY TRANSFERRED BY RAIN.)
      NL=NSP
      ABOVE= RAIN
      EXCESS= RAIN*(CL*TRAIN+LF)/LF
!
      DO 20 I=1,NSP
      IF (ICESPT(I) .EQ. 0) THEN
!
!        TEMPERATURE IS UNKNOWN--SOME EXCESS LIQUID-WATER FROM
!        THE ABOVE LAYER MUST BE FROZEN.
!
         IF (EXCESS .NE. 0.0) THEN
!           COMPUTE AMOUNT TO BE FROZE IN ORDER TO RAISE THE TEMPERATURE
!           TO ZERO DEGREES CELSIUS. (USE SPECIFIC HEAT OF ICE, CI)
            WEL=RHOSP(I)*DZSP(I)/RHOL
            FREEZE=-TSPDT(I)*WEL*CI/LF
            IF (EXCESS .GT. FREEZE) THEN
!              EXCESS EXCEEDS REFREEZE.
               TSPDT(I)=0.0
               WEL=WEL+FREEZE
               RHOSP(I)=RHOL*WEL/DZSP(I)
               EXCESS=EXCESS-FREEZE
               ABOVE=ABOVE-FREEZE
               ICESPT(I) = 1
             ELSE
!              EXCESS IS ALL FROZEN IN THIS LAYER.
               TSPDT(I)=TSPDT(I) +(EXCESS*LF*RHOL)/(CI*RHOSP(I)*DZSP(I))
               WEL=WEL+EXCESS
               RHOSP(I)=RHOL*WEL/DZSP(I)
               ABOVE = ABOVE - EXCESS
               EXCESS=0.0
            END IF
         END IF
      END IF
!
      DLWDT(I)= DLWDT(I) + EXCESS
      CHANGE= DLWDT(I) - DLW(I) - ABOVE
!
      IF (ABS(CHANGE).LT.0.0000001) THEN
!        THIS IS BEYOND PRECISION OF COMPUTER -- ASSUME ZERO
         ABOVE=0.0
         EXCESS=0.0
       ELSE
         IF (CHANGE .LE. 0.0) THEN
!           SOME LIQUID-WATER IS FROZEN. ADD TO ICE CONTENT OF LAYER.
            WEL=RHOSP(I)*DZSP(I)/RHOL
            WEL=WEL-CHANGE
            RHOSP(I)=RHOL*WEL/DZSP(I)
            ABOVE=0.0
            EXCESS=0.0
          ELSE
!           MELT HAS OCCURRED,SUBTRACT FROM ICE CONTENT
            WEL=RHOSP(I)*DZSP(I)/RHOL
!           IF MELT EXCEEDS ICE CONTENT OF LAYER--LAYER IS GONE.
            IF (CHANGE .GE. WEL) THEN
!              LAYER IS GONE
!              LIQUID-WATER IS ADDED TO THE LAYER BELOW (EXCESS).  THE
!              ICE CONTENT PLUS DLW OF THE LAYER CANNOT BE TAKEN FROM THE
!              NEXT LAYER, BUT IS STILL ADDED TO THE LIQUID-WATER
!              CONTENT OF THE NEXT LAYER (ABOVE).
               NL=NL-1
               ABOVE= ABOVE + WEL + DLW(I)
               EXCESS=DLWDT(I)
               DZSP(I)=0.0
             ELSE
!              LAYER REMAINS
               WEL=WEL-CHANGE
               DZSP(I)=RHOL*WEL/RHOSP(I)
               ABOVE=0.0
               EXCESS=0.0
            END IF
         END IF
      END IF
!
   20 CONTINUE
!
!***********************************************************************
!     CHECK TO SEE IF ENTIRE SNOW COVER IS GONE.
      IF (NL .LE. 0) THEN
         NSP = 0
         WE = 0.0
         TDEPTH = 0.0
         SCOUT = ABOVE
         CALL SNOMLT (NSP,ICESPT,DZSP,RHOSP,TSPDT,DLWDT,WLAG, &
                      STORE,SCOUT,col,row)
         RETURN
      END IF
!
!     ELIMINATE LAYERS WHICH ARE GONE.
      IF (NL .NE. NSP) THEN
         DO 30 I=1,NSP
   23       IF (DZSP(I) .GT. 0.0) GO TO 30
!           LAYER GONE. MOVE OTHER LAYERS UP.
            NEXT=I+1
            DO 25 J=NEXT,NSP
               L=J-1
               DZSP(L)=DZSP(J)
               RHOSP(L)=RHOSP(J)
               TSPDT(L)=TSPDT(J)
               DLWDT(L)=DLWDT(J)
               ICESPT(L)=ICESPT(J)
   25       CONTINUE
            NSP = NSP - 1
            IF (NSP .EQ. NL) GO TO 35
            GO TO 23
   30    CONTINUE
   35    DLWDT(NSP)= DLWDT(NSP) + ABOVE
      END IF
!
!***********************************************************************
!     CHECK THICKNESS OF EACH LAYER. IF NOT WITHIN SPECIFIED LIMITS,
!     DIVIDE (IF TOO LARGE) OR ADD TO AN ADJACENT LAYER (IF TOO SMALL).
      TDEPTH=0.0
      I=1
      IF (NSP .EQ. 1) GO TO 80
!
   40 ZZ=TDEPTH+DZSP(I)*0.5
!     COMPUTED DESIRED THICKNESS FOR THE LAYER.
      DZ1=THICK
      IF (ZZ .GT. 0.30)  DZ1=CTHICK*(ZZ-0.30)+THICK
!     CHECK ACTUAL THICKNESS AGAINST DESIRED THICKNESS.
      IF (DZSP(I) .GT. 1.55*DZ1) GO TO 50
      IF (DZSP(I) .LT. 0.55*DZ1 .AND. NSP.GT.1) GO TO 60
!     THICKNESS IS WITHIN SPECIFIED LIMITS.
      TDEPTH=TDEPTH+DZSP(I)
      GO TO 70
!***********************************************************************
!     THICKNESS IS GREATER THAN SPECIFIED LIMTS.
!          SUB-DIVIDE LAYER,THUS CREATING A NEW LAYER.
!          PROPERTIES ARE THE SAME FOR BOTH LAYERS.
!          MOVE OTHER LAYERS DOWN.
   50 NEXT=I+1
      IF (NEXT .LE. NSP) THEN
         DO 55 J=NEXT,NSP
            L=NSP-J+NEXT
            LL=L+1
            DZSP(LL)=DZSP(L)
            RHOSP(LL)=RHOSP(L)
            TSPDT(LL)=TSPDT(L)
            DLWDT(LL)=DLWDT(L)
            ICESPT(LL)=ICESPT(L)
   55    CONTINUE
      END IF
      TDEPTH=TDEPTH+DZSP(I)
      DZSP(NEXT)=DZSP(I)*0.5
      RHOSP(NEXT)=RHOSP(I)
      TSPDT(NEXT)=TSPDT(I)
      DLWDT(NEXT)=DLWDT(I)*0.5
      ICESPT(NEXT)=ICESPT(I)
      DZSP(I)=DZSP(I)*0.5
      DLWDT(I)=DLWDT(I)*0.5
      I=I+1
      NSP=NSP+1
      GO TO 70
!***********************************************************************
!     THICKNESS IS SMALLER THAN SPECIFIED LIMITS.
!          ADD TO THE SMALLEST ADJACENT LAYER, THUS
!          LOSING A LAYER. PROPERIES OF THE NEW LAYER
!          ARE THE WEIGHTED AVERAGE OF THE TWO FORMER
!          LAYERS. MOVE OTHER LAYERS UP.
   60 IF (I .EQ. 1  .OR.  I .EQ. NSP) THEN
         IF (I .EQ. 1) THEN
            NL=2
           ELSE
            NL=NSP-1
         END IF
       ELSE
         NL=I-1
         IF (DZSP(I+1) .LT. DZSP(I-1)) NL=I+1
      END IF
!
!     NL IS THE SMALLEST ADJACENT LAYER - ADD LAYER I TO LAYER NL
      WEI=RHOSP(I)*DZSP(I)/RHOL
      WENL=RHOSP(NL)*DZSP(NL)/RHOL
      WEL=WEI+WENL
      DZSP(NL)=DZSP(NL)+DZSP(I)
      RHOSP(NL)=RHOL*WEL/DZSP(NL)
      DLWDT(NL)= DLWDT(NL) + DLWDT(I)
      TSPDT(NL)=(WENL*TSPDT(NL)+WEI*TSPDT(I))/WEL
!
      IF (ICESPT(I) .NE. ICESPT(NL)) THEN
!        UNKNOWNS ARE DIFFERENT. COMPUTE THE UNKNOWN FOR THE NEW LAYER.
         FREEZE=-(TSPDT(NL)*CI*RHOSP(NL)*DZSP(NL))/(LF*RHOL)
         IF (DLWDT(NL) .LE. FREEZE) THEN
!           TEMPERATURE IS UNKNOWN
            TSPDT(NL)=TSPDT(NL) + (DLWDT(NL)*LF*RHOL) &
                                  /(CI*RHOSP(NL)*DZSP(NL))
            WEL=WEL+DLWDT(NL)
            RHOSP(NL)=RHOL*WEL/DZSP(NL)
            ICESPT(NL)=0
          ELSE
!           LIQUID-WATER IS UNKNOWN.
            DLWDT(NL)=DLWDT(NL)-FREEZE
            WEL=WEL+FREEZE
            RHOSP(NL)=RHOL*WEL/DZSP(NL)
            ICESPT(NL)=1
         END IF
      END IF
!
      IF (ICESPT(NL) .EQ. 1) TSPDT(NL)=0.0
      IF (ICESPT(NL) .EQ. 0) DLWDT(NL)=0.0
!     MOVE OTHER LAYERS UP.
      IF (NL .LT. I) TDEPTH=TDEPTH+DZSP(I)
      NEXT=I+1
!
      IF (NEXT .LE. NSP) THEN
         DO 65 J=NEXT,NSP
            L=J-1
            DZSP(L)=DZSP(J)
            RHOSP(L)=RHOSP(J)
            TSPDT(L)=TSPDT(J)
            DLWDT(L)=DLWDT(J)
            ICESPT(L)=ICESPT(J)
   65    CONTINUE
      END IF
!
      I = I-1
      NSP = NSP-1
!
!***********************************************************************
!     CHECK FOR LAST LAYER.
   70 IF (I .LT. NSP) THEN
!        START NEXT LAYER
         I=I+1
         GO TO 40
      END IF
!
!***********************************************************************
!     ADJUST THE DENSITY FOR METAMORPHISM AND COMPUTE ANY SNOMELT
   80 CALL META (NSP,TSPDT,RHOSP,DZSP,DLWDT,col,row)
      CALL SNOMLT (NSP,ICESPT,DZSP,RHOSP,TSPDT,DLWDT,WLAG,STORE,SCOUT,col,row)
!
!***********************************************************************
!     CALCULATE THE DEPTH FOR EACH NODE
      ZSP(1) = 0.0
      TDEPTH = DZSP(1)
      WE = DZSP(1)*RHOSP(1)/RHOL + DLWDT(1)
      DO 85 I=2,NSP
         ZSP(I)= TDEPTH + DZSP(I)/2.
         TDEPTH= TDEPTH + DZSP(I)
         WE = WE + DZSP(I)*RHOSP(I)/RHOL + DLWDT(I)
   85 CONTINUE
      ZSP(NSP+1) = TDEPTH
!
      RETURN
      END SUBROUTINE







!***********************************************************************
!
SUBROUTINE NWSNOW (NSP,NC,NR,NS,NPLANT,ICESPT,ZSP,DZSP,RHOSP, &
        TSPDT,DLW,DLWDT,TWBULB,SNOW,SNODEN,TQVSP,MELT,ZC,PCANDT,ZR,TRDT, &
        GMCDT,RHOR,ZS,VLCDT,VICDT,TSDT,MATDT,CONCDT,nsalt,col,row)
!
!     THIS SUBROUTINES ADDS SNOW LAYERS FOR NEWLY FALLEN SNOW
!
!***********************************************************************
    use constvar_mod,only:LF,LS,LV,CI,CL,RHOL
    use spwatr_mod
    implicit none



!input
    integer(i4),intent(inout)::NSP
    integer(i4),intent(in)::NC,NR,NS,NPLANT,col,row,nsalt
    real(r8),dimension(NSPMAX),intent(inout)::DZSP,RHOSP,TSPDT,DLWDT,ZSP,DLW,TQVSP

    real(r8),dimension(NCMAX),intent(in)::ZC
    real(r8),dimension(NRMAX),intent(in)::ZR,GMCDT,RHOR
    real(r8),dimension(NRMAX),intent(inout)::TRDT
    real(r8),dimension(NPMAX),intent(inout)::PCANDT
    real(r8),dimension(NSMAX),intent(in)::ZS,VLCDT,VICDT,MATDT
    real(r8),dimension(NSMAX),intent(inout)::TSDT
    real(r8),dimension(NSALTMAX,NSMAX),intent(in)::CONCDT
    integer(i4),dimension(NSPMAX),intent(inout)::ICESPT
    real(r8),intent(in)::TWBULB,snoden
    real(r8),intent(inout)::MELT,snow





!temp
    real(r8)::DENSTY,DEPTH,TOTAL,CHANGE,dns,down,EXTRA,FREEZE,TDEPTH,WEI,WEL,WENS
    integer(i4)::I,J,NWLAYR
    real(r8),dimension(NRMAX)::CRES
    real(r8),dimension(NSMAX)::CS


!
!     PRECIPITATION IS SNOW
      IF (SNODEN .GT. 0.0) THEN
!        DENSITY OF NEW SNOW IS KNOWN -- CONVERT TO KG/M3
         DENSTY = SNODEN*RHOL
       ELSE
!        COMPUTE DENSITY OF THE NEW SNOW BASED ON WETBULB -3.61199235E-03TEMPERATURE
         IF (TWBULB .LE. -15.) THEN
            DENSTY = 50.
          ELSE
            DENSTY = 50. + 1.7*((TWBULB + 15.)**1.5)
         END IF
      END IF
      DEPTH = SNOW*RHOL/DENSTY
!
      IF (NC .GT. 0) THEN
!        ADJUST CANOPY INTERCEPTION FOR DEPTH OF CANOPY COVERED BY SNOW
!        (NO CONSIDERATION ABOUT WHERE LAI IS WITHIN THE CANOPY)
         TOTAL=0.0
         DO 3 J=1,NPLANT
            CHANGE=PCANDT(J)*DEPTH/ZC(NC+1)
            IF (CHANGE .GT. PCANDT(J)) CHANGE=PCANDT(J)
            PCANDT(J)=PCANDT(J)-CHANGE
            TOTAL=TOTAL+CHANGE
    3    CONTINUE
         SNOW=SNOW+TOTAL
         DENSTY=SNOW*RHOL/DEPTH
      END IF
!
      MELT = 0.0
      IF (NSP .EQ. 0) THEN
!        SNOWCOVER DOES NOT EXIST -- MELT SOME SNOW TO GET UNDERLYING
!        MATERIAL TO 0 DEGREES
         IF (NR .GT. 0) THEN
!           SNOW IS FALLING ON RESIDUE
            CALL RESHT (NR,CRES,GMCDT,RHOR)
            IF (TRDT(1) .LE. 0.0) THEN
!              NO MELTING OF SNOW IS NEEDED TO GET RESIDUE TO 0 DEGREES
               MELT=0.0
              ELSE
               MELT = CRES(1)*ZR(2)/2.*TRDT(1)/(RHOL*LF-RHOL*CI*TWBULB)
               SNOW = SNOW - MELT
               DEPTH = max(DEPTH - MELT*RHOL/DENSTY,0.0)
               IF (SNOW .LE. 0.0) THEN
!                 ALL SNOW THAT FELL IS MELTED - ADJUST RESIDUE FOR
!                 ENERGY REQUIRED TO MELT SNOW (MELTWATER ENERGY WILL
!                 BE CONSIDERED LATER IN SUBROUTINE RAINRS)
                  MELT = MELT + SNOW
                  TRDT(1)= TRDT(1) - (RHOL*LF-RHOL*CI*TWBULB)*MELT &
                                     /(CRES(1)*ZR(2)/2.)
                  RETURN
               END IF
               TRDT(1) = 0.0
            END IF
!
          ELSE
!           SNOW IS FALLING ON BARE SOIL
            IF (TSDT(1) .LE. 0.0) THEN
!              NO MELTING OF SNOW IS NEEDED TO GET SOIL TO 0 DEGREES
               MELT=0.0
              ELSE
               CALL SOILHT(NS,CS,VLCDT,VICDT,TSDT,MATDT,CONCDT,nsalt,col,row)
               MELT = CS(1)*ZS(2)/2.*TSDT(1)/(RHOL*LF-RHOL*CI*TWBULB)
               SNOW = SNOW - MELT
               DEPTH = max(DEPTH - MELT*RHOL/DENSTY,0.0)
               IF (SNOW .LE. 0.0) THEN
!                 ALL SNOW THAT FELL IS MELTED - ADJUST SOIL TEMP FOR
!                 ENERGY REQUIRED TO MELT SNOW (MELTWATER ENERGY WILL
!                 BE CONSIDERED LATER IN SUBROUTINE RAINSL)
                  MELT = MELT + SNOW
                  TSDT(1)= TSDT(1) - (RHOL*LF-RHOL*CI*TWBULB)*MELT &
                                     /(CS(1)*ZS(2)/2.)
                  RETURN
               END IF
               TSDT(1) = 0.0
            END IF
         END IF
!
       ELSE
!        NEW SNOW IS FALLING ON OLD SNOW - CHECK IF NEW SNOW CAUSES THE
!        THE SURFACE LAYER TO EXCEED THE SPECIFIED DEPTH LIMIT.
         IF((DZSP(1) + DEPTH) .LE. (1.55*THICK)) THEN
!           INCORPORATE ALL OF NEW SNOW INTO OLD SURFACE LAYER.
            DNS = DEPTH
            DEPTH = 0.0
            NWLAYR = 0
           ELSE
!           DETERMINE HOW MUCH SHOULD BE COMBINED WITH OLD SURFACE LAYER
            IF (DZSP(1) .GT. THICK) THEN
!              LAYER IS ALREADY GREATER THAN DESIRED THICKNESS - DO NOT
!              ADD NEW SNOW TO OLD SURFACE LAYER
               DNS = 0.0
!              GO SET UP NEW LAYERS FOR THE NEW SNOW
               GO TO 5
            END IF
            DNS = THICK - DZSP(1)
            DEPTH = max(DEPTH - DNS,0.0)
         END IF
!
!        INCORPORATE NEW SNOW INTO OLD SURFACE LAYER - DEFINE PROPERTIES
         WEI=RHOSP(1)*DZSP(1)/RHOL
         WENS=DENSTY*DNS/RHOL
         WEL=WEI+WENS
         DZSP(1) = DZSP(1) + DNS
         RHOSP(1)=RHOL*WEL/DZSP(1)
         TSPDT(1)=(WENS*TWBULB + WEI*TSPDT(1))/WEL
         IF (ICESPT(1) .EQ. 1)  THEN
!           LIQUID-WATER AND SUBFREEZING TEMP MAY EXIST IN NEW LAYER
            FREEZE=-(TSPDT(1)*CI*WEL)/(LF*RHOL)
            IF (DLWDT(1) .LE. FREEZE) THEN
               TSPDT(1)=TSPDT(1)+(DLWDT(1)*LF*RHOL)/(CI*WEL)
               WEL=WEL+DLWDT(1)
               RHOSP(1)=RHOL*WEL/DZSP(1)
               DLWDT(1)=0.0
               ICESPT(1)=0
              ELSE
               DLWDT(1)=DLWDT(1)-FREEZE
               WEL=WEL+FREEZE
               RHOSP(1)=RHOL*WEL/DZSP(1)
               TSPDT(1)= 0.0
            END IF
         END IF
      END IF
!
!
!     COMPUTE THE NUMBER OF LAYERS WITH ENTIRELY NEW SNOW 
    5 NWLAYR = DEPTH/THICK
      if (NWLAYR.lt.0) print*,NWLAYR,DEPTH,THICK,col,row
      EXTRA = DEPTH - NWLAYR*THICK
      IF (EXTRA .GT. 0.55*THICK) NWLAYR = NWLAYR + 1
      IF (EXTRA .GT. 0.0  .AND.  NWLAYR .EQ. 0) NWLAYR = 1
!
      IF (NSP .GT. 0) THEN
         DOWN = DEPTH + DNS
         IF (NWLAYR .EQ. 0) THEN
!           ADJUST DEPTHS OF THE LAYERS AND RETURN
            ZSP(1) = 0.0
            DO 10 I=2,NSP+1
               ZSP(I) = ZSP(I) + DOWN
   10       CONTINUE
            RETURN
          ELSE
!           MOVE LAYERS DOWN TO MAKE ROOM FOR NEW SNOW LAYERS
            ZSP(NSP+NWLAYR+1) = ZSP(NSP+1) + DOWN
            TSPDT(NSP+NWLAYR+1) = TSPDT(NSP+1)
            DO 15 I=NSP,1,-1
               RHOSP(NWLAYR + I)= RHOSP(I)
               TSPDT(NWLAYR + I)= TSPDT(I)
               DLW(NWLAYR + I)= DLW(I)
               DLWDT(NWLAYR + I)= DLWDT(I)
               ICESPT(NWLAYR + I)= ICESPT(I)
               ZSP(NWLAYR + I) = ZSP(I) + DOWN
               DZSP(NWLAYR + I) = DZSP(I)
               TQVSP(NWLAYR + I) = TQVSP(I)
   15       CONTINUE
!           ADJUST DEPTH OF OLD SURFACE LAYER TO BE IN MIDDLE OF LAYER
            ZSP(NWLAYR+1)=ZSP(NWLAYR+1)+DZSP(NWLAYR+1)/2.
         END IF
      END IF
!
!     LEFT-OVER SNOW GOES IN TOP LAYER
      EXTRA=DEPTH - (NWLAYR-1)*THICK
      ZSP(1)=0.0
      DZSP(1)=EXTRA
      RHOSP(1)=DENSTY
      TSPDT(1)=TWBULB
      DLW(1)=0.0
      DLWDT(1)=0.0
      ICESPT(1)= 0
      TQVSP(1)=0.0
      TDEPTH=DZSP(1)
!
!     FULL SIZE NEW SNOW LAYERS.
      DO 20 J=2,NWLAYR
         DZSP(J)=THICK
         ZSP(J)= TDEPTH + DZSP(J)/2.
         TDEPTH=TDEPTH + DZSP(J)
         RHOSP(J)=DENSTY
         TSPDT(J)=TWBULB
         DLW(J)=0.0
         DLWDT(J)=0.0
         ICESPT(J)= 0
         TQVSP(J)=0.0
   20 CONTINUE
!
      IF (NSP .EQ. 0) THEN
!        DEFINE DEPTH OF BOTTOM BOUNDARY OF SNOWPACK
         ZSP(NWLAYR+1) = TDEPTH
!        NEW SNOW COVER - DEFINE TEMP OF BOTTOM BOUNDARY OF SNOWPACK
!        THIS NODE IS SHARED BY THE MATERIAL BELOW
         IF (NR .GT. 0) THEN
!           UNDERLYING MATERIAL IS RESIDUE
            TSPDT(NWLAYR+1) = TRDT(1)
           ELSE
!           UNDERLYING MATERIAL IS SOIL
            TSPDT(NWLAYR+1) = TSDT(1)
         END IF
      END IF
!
      NSP = NSP + NWLAYR
      RETURN
      END SUBROUTINE


!***********************************************************************
!
SUBROUTINE WTBULB (TW,TA,HUM,PRESUR)
!
!     THIS SUBROUTINE CALCULATES THE WETBULB TEMPERATURE (TW IN CELCIUS)
!     FROM THE AIR TEMPERATURE AND RELATIVE HUMIDITY
!     ---- PROCEDURE FOR FINDING WETBULB TEMPERATURE WAS ADOPTED FROM
!          THE ANDERSON SNOWMELT MODEL.
!
!***********************************************************************
    implicit none
    real(r8),intent(in)::PRESUR,TA,HUM
    real(r8),intent(inout)::TW

    real(r8)::PA,EAS,EA,DELTA,TAV,EAV
    integer(i4)::I



!
!     CONVERT ATMOSPHERIC PRESSURE FROM PASCALS TO MBARS
      PA=0.01*PRESUR
!
!     CALCULATE VAPOR PRESSURE AND SATURATED VAPOR PRESSURE IN MBARS
      EAS= 2.7489E8*EXP(-4278.63/(TA+242.792))
      EA=HUM*EAS
!
      DELTA= (4278.63/((TA+242.792)*(TA+242.792)))*EAS
      DO 10 I=1,3
         TW= DELTA*TA+6.6E-4*PA*TA+7.59E-7*PA*TA*TA+EA-EAS
         TW= TW/(DELTA+6.6E-4*PA+7.59E-7*PA*TA)
         TAV= (TA+TW)*0.5
         EAV= 2.7489E8*EXP(-4278.63/(TAV+242.792))
         DELTA= (4278.63/((TAV+242.792)*(TAV+242.792)))*EAV
   10 CONTINUE
!
      RETURN
 END SUBROUTINE WTBULB



!***********************************************************************
!
SUBROUTINE RAINC (NPLANT,NC,XANGLE,CLUMPNG,ITYPE,RAIN,PCANDT, &
     WCANDT,WCMAX,SLOPE,col,row)
!
!     THIS SUBROUTINE CALCULATES THE RAINFALL DEPTH INTERCEPTED BY THE
!     CANOPY LAYERS AND ADJUST THE WATER CONTENT OF THE DEAD PLANT
!
!***********************************************************************
    use constvar_mod
    use clayrs_mod,only:TOTLAI2d,CANLAI2d,DRYCAN2d
    implicit none

! input
    integer(i4),intent(in)::NPLANT,NC,col,row
    real(r8),dimension(NPMAX),intent(in)::XANGLE,CLUMPNG
    integer(i4),dimension(NPMAX),intent(in)::ITYPE
    real(r8),intent(in)::SLOPE,WCMAX

    real(r8),intent(inout)::RAIN
    real(r8),dimension(NPMAX),intent(inout)::PCANDT,WCANDT


!temp
    real(r8)::ANGLE,MAXINT,TRANS,EXCESS,RINTER,RINT,RINTLY,SLAI,SUMLAI,TOTAL
    real(r8),dimension(NCMAX-1)::TDIRCC,TDIFFC
    real(r8),dimension(NPMAX+1,NCMAX-1)::DIRKL,DIFKL
    real(r8),dimension(NPMAX)::FBDU,FDDU
    integer(i4)::I,J

!pointer
    real(r8),dimension(:),pointer::TOTLAI
    real(r8),dimension(:,:),pointer::CANLAI,DRYCAN


    TOTLAI=>TOTLAI2d(col,row,:)
    CANLAI=>CANLAI2d(col,row,:,:)
    DRYCAN=>DRYCAN2d(col,row,:,:)


!
!     USE THE TRANSMITTANCE TO DIRECT RADIATION TO CALCULATE THE
!     FRACTION OF RAIN INTERCEPTED BY THE CANOPY
      ANGLE=1.5708-SLOPE
      CALL TRANSC (NPLANT,NC,XANGLE,CLUMPNG,TDIRCC,TDIFFC,DIFKL,DIRKL, &
                   FBDU,FDDU,ANGLE,1.57_r8,col,row)
!
!**** CALCULATE RAINFALL INTERCEPTED BY THE ENTIRE CANOPY
      TRANS=1.0
      DO 5 I=1,NC
         TRANS=TRANS*TDIRCC(I)
    5 CONTINUE
      RINTER=RAIN*(1.-TRANS)
!
      SUMLAI=0.0
      DO 10 J=1,NPLANT
         SUMLAI=SUMLAI+TOTLAI(J)
   10 CONTINUE
!
      TOTAL = 0.0
!**** DETERMINE HOW MUCH RAINFALL INTERCEPTED BY EACH PLANT SPECIES
      DO 40 J=1,NPLANT
         IF (ITYPE(J) .NE. 0) THEN
!****       PRECIP INTERCEPTED BY LIVING PLANTS
            SLAI=0.0
            DO 20 I=1,NC
               SLAI=SLAI+CANLAI(J,I)
   20       CONTINUE
!
!           CALCULATE PRECIP INTERCEPTED BY ALL LAYERS
            RINT=RINTER*SLAI/SUMLAI
!           ADD TO PRECIP REMAINING ON PLANT LEAVES FROM BEFORE
            PCANDT(J)=PCANDT(J)+RINT
!           ESTIMATE HOW MUCH PRECIP CAN BE HELD ON PLANTS
!           -- ASSUME LIVING PLANTS CAN INTERCEPT AND HOLD 1MM OF WATER
!           FOR EACH UNIT OF LEAF AREA INDEX
            MAXINT=0.001*SLAI
!           CHECK IF PRECIP INTERCEPTED EXCEEDS AMOUNT THAT CAN BE HELD
            IF (PCANDT(J) .GT. MAXINT) THEN
               RINT = MAXINT - (PCANDT(J) - RINT)
               PCANDT(J) = MAXINT
            END IF
            TOTAL = TOTAL + RINT
          ELSE
!****       DEAD PLANT MATERIAL
            EXCESS = 0.0
            DO 30 I=1,NC
!           CHECK IF PLANT EXISTS IN LAYER
            IF (CANLAI(J,I).GT.0.0) THEN
!              PRECIP INTERCEPTED BY LAYER
               RINTLY=RINTER*CANLAI(J,I)/SUMLAI + EXCESS
               WCANDT(I) = WCANDT(I) + RINTLY*RHOL/DRYCAN(J,I)
               EXCESS = 0.0
!              CHECK IF THE LAYER CAN HOLD THIS MUCH WATER
               IF (WCANDT(I) .GT. WCMAX) THEN
                  EXCESS = (WCANDT(I)-WCMAX)*DRYCAN(J,I)/RHOL
                  RINTLY = RINTLY-EXCESS
                  WCANDT(I)=WCMAX
               END IF
               TOTAL=TOTAL+RINTLY
            END IF
   30       CONTINUE
         END IF
   40 CONTINUE
!
!     CALCULATE THE AMOUNT OF PRECIP THAT MAKES IT THROUGH CANOPY
      RAIN = RAIN - TOTAL
!
      RETURN
END SUBROUTINE






!***********************************************************************
!
SUBROUTINE PRECP (NPLANT,NC,NSP,NR,NS,TA,TADT,HUM,HUMDT,&
    ZC,XANGLE,CLUMPNG,ITYPE,WCANDT,WCMAX,PCANDT,&
    ZSP,DZSP,TSPDT,DLW,DLWDT,RHOSP,TQVSP,TOPSNO,WLAG,STORE,SNOWEX,&
    ZR,TRDT,GMCDT,GMCMAX,RHOR,&
    ZS,TSDT,VLCDT,VICDT,MATDT,TOTFLO,SALTDT,CONCDT,ICESPT,ICESDT,&
    RAIN,POND,RUNOFF,EVAP1,PONDMX,SNOTMP,SNODEN,DIRRES,SLOPE,PRESUR,nsalt,col,row,infiltration,snowmelt)
!
!     THIS SUBROUTINE DETERMINES WHERE THE PRECIPITATION AND SNOWMELT
!     SHOULD GO
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d
    use constvar_mod,only:LF,LS,LV,RHOL
    implicit none

!input
    integer(i4),intent(in)::NPLANT,NC,NR,NS,nsalt,col,row
    integer(i4),intent(inout)::NSP
    real(r8),intent(in)::TA,TADT,HUM,HUMDT
    real(r8),dimension(NCMAX),intent(in)::ZC
    real(r8),dimension(NCMAX-1),intent(inout)::WCANDT
    real(r8),dimension(NPMAX),intent(in)::XANGLE,CLUMPNG
    real(r8),dimension(NPMAX),intent(inout)::PCANDT

    integer(i4),dimension(NPMAX),intent(in)::ITYPE
    real(r8),dimension(NSPMAX),intent(inout)::TQVSP
    real(r8),dimension(NSPMAX),intent(inout)::ZSP,DZSP,RHOSP,TSPDT,DLW,DLWDT
    real(r8),dimension(11),intent(inout)::WLAG
    real(r8),dimension(NRMAX),intent(in)::ZR,RHOR
    real(r8),dimension(NRMAX),intent(inout)::TRDT,GMCDT
    real(r8),dimension(NSMAX),intent(in)::ZS
    real(r8),dimension(NSMAX),intent(inout)::TSDT,VLCDT,VICDT,MATDT,TOTFLO
    real(r8),dimension(NSALTMAX,NSMAX),intent(inout)::SALTDT,CONCDT

    integer(i4),dimension(NSPMAX),intent(inout)::ICESPT
    integer(i4),dimension(NSMAX),intent(in)::ICESDT

    real(r8),intent(inout)::RAIN,POND,SNOWEX,RUNOFF,STORE,TOPSNO,EVAP1,infiltration,snowmelt
    real(r8),intent(in)::GMCMAX,PONDMX,SNOTMP,SNODEN,DIRRES,SLOPE,WCMAX
    real(r8),intent(in)::PRESUR


!temp
    real(r8)::MELT,AVGTA,AVGHUM,TRAIN,SNOW,SCOUT
    integer(i4)::I,J
    real(r8)::WE

!pointer
    real(r8),pointer::WT,WDT

    WT  =>  WT2d(col,row)
    WDT =>  WDT2d(col,row)


!
      AVGTA= WT*TA + WDT*TADT
      AVGHUM= WT*HUM + WDT*HUMDT
      TRAIN=AVGTA
      snowmelt=0.0
!
!     RAINFALL INTERCEPTION BY THE CANOPY
      IF (RAIN.GT.0.0 .AND. NC.GT.0) &
         CALL RAINC (NPLANT,NC,XANGLE,CLUMPNG,ITYPE,RAIN,PCANDT, &
                     WCANDT,WCMAX,SLOPE,col,row)
!
      MELT = 0.0
      IF (RAIN .GT. 0.0) THEN
!        RAIN TEMPERATURE IS ASSUMED TO BE EQUAL TO WET BULB TEMPERATURE
         CALL WTBULB (TRAIN,AVGTA,AVGHUM,PRESUR)
         IF (AVGTA .LE. SNOTMP  .OR.  SNODEN .GT. 0.0) THEN
!           PRECIPITATION IS ASSUMED TO BE SNOW
!           (ADD ANY SNOW EXCESS REMAINING FROM ELIMINATING A SHALLOW
!           SNOWCOVER)
            SNOW = RAIN + SNOWEX
            POND = POND - SNOWEX
            SNOWEX = 0.0
            RAIN = 0.0

            IF (TRAIN .GT. 0.0) TRAIN = 0.0
            CALL NWSNOW (NSP,NC,NR,NS,NPLANT,ICESPT,ZSP,DZSP,RHOSP,&
              TSPDT,DLW,DLWDT,TRAIN,SNOW,SNODEN,TQVSP,MELT,ZC,PCANDT,&
              ZR,TRDT,GMCDT,RHOR,ZS,VLCDT,VICDT,TSDT,MATDT,CONCDT,nsalt,col,row)
            snowmelt=MELT+snowmelt
            IF (NSP .LE. 0) THEN
!              SNOW DID NOT STICK (MELT IS NOT INTERCEPTED BY CANOPY)
               TRAIN = 0.0
            END IF
           ELSE
!           PRECIP IS ASSUMED TO BE RAIN - DO NOT ALLOW RAIN TEMP < 0.0
            IF (TRAIN .LT. 0.0) TRAIN = 0.0
         END IF
      END IF
!
!     RAINFALL ABSORPTION BY SNOWPACK AND CALCULATION OF SNOWMELT
      IF (NSP .GT. 0) THEN
         CALL WBSNOW (NSP,ICESPT,ZSP,DZSP,RHOSP,TSPDT,DLW,DLWDT,RAIN, &
                     TRAIN,TOPSNO,EVAP1,TQVSP,WLAG,STORE,SCOUT,col,row)
         snowmelt=SCOUT+snowmelt
         RAIN = SCOUT
         TRAIN = 0.0
      END IF
!
!     DO NOT INFILTRATE PONDED WATER OR SNOW EXCESS (SNOWEX) IF
!     TEMPERATURE IS BELOW FREEZING
      IF ((RAIN + POND + MELT).GT.0.0 .AND. TRAIN.GE.0.0) THEN
!
!        ADD MELT (FROM SNOW NOT STICKING) BACK INTO RAIN
         RAIN=RAIN+MELT
!
!        RAINFALL INTERCEPTION BY RESIDUE LAYERS
         IF(NR.GT.0) CALL RAINRS (NR,TRAIN,RAIN,ZR,TRDT,GMCDT,RHOR, &
                                  GMCMAX,DIRRES,SLOPE)
!
         IF (RAIN .GT. 0.0) THEN
!           ADD ANY EXCESS FROM ELIMINATING SHALLOW SNOWPACK BACK
!           INTO RAIN
            RAIN=RAIN+SNOWEX
            POND=POND-SNOWEX
            SNOWEX=0.0
         END IF
         RAIN = RAIN + POND
         POND = 0.0
!
!        RAINFALL INFILTRATION INTO THE SOIL
         CALL RAINSL (NS,TRAIN,RAIN,ZS,TSDT,VLCDT,VICDT,ICESDT, &
              MATDT,TOTFLO,SALTDT,CONCDT,POND,PONDMX,nsalt,col,row,infiltration)
!
         IF (SNOWEX .GT. 0.0) THEN
!           DO NOT ALLOW EXCESS WATER FROM ELIMINATING SHALLOW SNOWPACK
!           TO RUNOFF - PUT THIS WATER BACK INTO PONDING AND SNOW EXCESS
!           (THIS MAY ALLOW PONDING DEPTH TO BECOME GREATER THAN PONDMX)
            POND = POND + RAIN
            RAIN = 0.0
            IF (POND .LT. SNOWEX) SNOWEX=POND
         END IF
      END IF
!
      IF (NSP .GT. 0) THEN
!
!     CALCULATE WATER EQUIVALENT - IF SUFFICIENTLY SMALL, ASSUME ENTIRE
!     SNOWPACK IS GONE
      WE = 0.0
      DO 20 I=1,NSP
         WE = WE + DZSP(I)*RHOSP(I)/RHOL + DLWDT(I)
   20 CONTINUE
      IF (WE.LT.0.0) THEN
         print*,"WE",WE,col,row
         WE=0.0
      ENDIF
      IF (WE .LT. 0.0005) THEN
!        WATER EQUIVALENT IS SUFFICIENTLY SMALL-ASSUME SNOWPACK IS GONE
         SCOUT = WE
         NSP = 0
!        ADD WATER CURRENTLY BEING LAGGED TO SNOW COVER OUTFLOW
         CALL SNOMLT (NSP,ICESPT,DZSP,RHOSP,TSPDT,DLWDT,WLAG,&
                      STORE,SCOUT,col,row)
         snowmelt=SCOUT+snowmelt
         WRITE (21,*)
         WRITE (21,*)'SNOWCOVER IS ASSUMED GONE: WATER-EQUIVALENT ', &
                     'OF ICE AND TOTAL WATER =',WE,SCOUT,' METERS'
!
!        SET SNOW WATER EXCESS (SNOWEX) TO THE TOTAL WATER EQUIV.
!        REMAINING IN THE SNOW COVER - THIS WILL BE ADDED TO
!        PRECIPITATION OR ALLOWED TO INFITRATE AT A LATER TIME STEP
         SNOWEX=SNOWEX + SCOUT
!        INCLUDE SNOW EXCESS IN AMOUNT PONDED SO IT CAN BE INCLUDED IN
!        WATER BALANCE
         POND=POND + SCOUT
!
        ELSE
!        SNOWPACK STILL PRESENT -- DEFINE TEMPERATURE FOR BOTTOM BOUNDARY
         IF (NR .GT. 0) THEN
!           BOTTOM TEMPERATURE BOUNDARY IS TOP OF RESIDUE
            TSPDT(NSP+1) = TRDT(1)
          ELSE
!           BOTTOM TEMPERATURE BOUNDARY IS SOIL SURFACE
            TSPDT(NSP+1) = TSDT(1)
         END IF
      END IF
!
      END IF
!
      RUNOFF=RAIN
!
      RETURN
END SUBROUTINE PRECP






























!***********************************************************************
!
SUBROUTINE SUMDT (NC,NPLANT,NSP,NR,NS,HFLUX,VFLUX,GFLUX,LWCAN, &
     LWSNOW,LWRES,LWSOIL,SWCAN,SWSNOW,SWRES,SWSOIL,QVC,QVR,QVSP,QSL, &
     QSV,TRNSP,XTRACT,SEEP,TSWSNO,TLWSNO,TSWCAN,TLWCAN,TSWRES,TLWRES, &
     TSWSOI,TLWSOI,THFLUX,TGFLUX,EVAP1,ETSUM,TSEEP,ROOTXT,TOTFLO, &
     NTIMES,TOPSNO,TQVSP,col,row,AbsorbedSW,AbsorbedLW)
!
!     THIS SUBROUTINE SUMS THE NECESSARY FLUXES EACH TIME STEP FOR USE
!     IN THE WATER AND ENERGY BALANCE SUMMARIES AND IN THE SNOWPACK
!     ADJUSTMENTS DUE FOR VAPOR TRANSFER
!
!***********************************************************************
    use controlpara_mod,only:DT2d
    use constvar_mod,only:LF,LS,LV,RHOL
    use savedata_mod,only:TRANSP2d
    implicit none
! input
    integer(i4),intent(in)::NC,NPLANT,NSP,NR,NS,col,row
    real(r8),intent(in)::HFLUX,VFLUX,GFLUX
    real(r8),dimension(NPMAX+1,NCMAX-1),intent(in)::LWCAN,SWCAN
    real(r8),intent(in)::LWSNOW,LWSOIL,SWSOIL
    real(r8),dimension(NRMAX),intent(in)::LWRES,SWRES,QVR
    real(r8),dimension(NSPMAX),intent(in)::SWSNOW,QVSP
    real(r8),dimension(NSPMAX),intent(inout)::TQVSP
    real(r8),dimension(NCMAX-1),intent(in)::QVC
    real(r8),dimension(NSMAX),intent(in)::QSL,QSV,XTRACT
    real(r8),dimension(NSMAX),intent(out)::TOTFLO,ROOTXT
    real(r8),dimension(NPMAX+1),intent(in)::TRNSP
    real(r8),intent(inout)::AbsorbedSW,AbsorbedLW
    real(r8),intent(inout)::SEEP,ETSUM,TSEEP,TOPSNO
    real(r8),intent(inout)::TSWCAN,TLWCAN,TSWSNO,TLWSNO,TSWRES,TLWRES,TSWSOI,TLWSOI,THFLUX,TGFLUX,EVAP1
    integer(i4),intent(in)::NTIMES
!temp
    integer(i4)::I,J
!pointer
    real(r8),dimension(:),pointer::TRANSP
    real(r8),pointer::DT

    TRANSP   =>   TRANSP2d(col,row,:)
    Dt       =>   DT2d(col,row)


!
      IF (NTIMES .EQ. 1) THEN
         TSWCAN=0.0
         TLWCAN=0.0
         TSWSNO=0.0
         TLWSNO=0.0
         TSWRES=0.0
         TLWRES=0.0
         TSWSOI=0.0
         TLWSOI=0.0
         THFLUX=0.0
         TGFLUX=0.0
         EVAP1=0.0
         ETSUM=0.0
         TSEEP=0.0
         TOPSNO=0.0
         AbsorbedSW=0.0
         AbsorbedLW=0.0
         DO 6 J=1,NPLANT
            TRANSP(J)=0.0
   6     CONTINUE
         DO 8 I=1,NSP
            TQVSP(I)=0.0
   8    CONTINUE
        DO 10 I=1,NS
            TOTFLO(I)=0.0
            ROOTXT(I)=0.0
  10    CONTINUE
      END IF
!
!     SUM RADIATION ABSORBED BY THE CANOPY (IN KILOJOULES)
      DO 20 I=1,NC
         TSWCAN = TSWCAN + SWCAN(NPLANT+1,I)*DT/1000.
         TLWCAN = TLWCAN + LWCAN(NPLANT+1,I)*DT/1000.
   20 CONTINUE
      AbsorbedSW=AbsorbedSW+TSWCAN/DT*1000.0
      AbsorbedLW=AbsorbedLW+TLWCAN/DT*1000.0
!
!     SUM RADIATION ABSORBED BY THE SNOWPACK (IN KILOJOULES)
      DO 30 I=1,NSP
         TSWSNO = TSWSNO + SWSNOW(I)*DT/1000.
   30 CONTINUE
      IF (NSP .GT. 0) TLWSNO = TLWSNO + LWSNOW*DT/1000.
      AbsorbedSW=AbsorbedSW+TSWSNO/DT*1000.0
      AbsorbedLW=AbsorbedLW+TLWSNO/DT*1000.0
!     SUM RADIATION ABSORBED BY THE RESIDUE (IN KILOJOULES)
      DO 40 I=1,NR
         TSWRES = TSWRES + SWRES(I)*DT/1000.
         TLWRES = TLWRES + LWRES(I)*DT/1000.
   40 CONTINUE
      AbsorbedSW=AbsorbedSW+TSWRES/DT*1000.0
      AbsorbedLW=AbsorbedLW+TLWRES/DT*1000.0
!     SUM RADIATION ABSORBED BY THE SOIL (IN KILOJOULES)
      TSWSOI = TSWSOI + SWSOIL*DT/1000.
      TLWSOI = TLWSOI + LWSOIL*DT/1000.
      AbsorbedSW=AbsorbedSW+TSWSOI/DT*1000.0
      AbsorbedLW=AbsorbedLW+TLWSOI/DT*1000.0
      
!     SUM SENSIBLE HEAT FLUX AND SOIL HEAT FLUX
      THFLUX = THFLUX + HFLUX*DT/1000.
      TGFLUX = TGFLUX + GFLUX*DT/1000.
!
!     SUM THE EVAPORATION, TRANSPIRATION FROM EACH PLANT SPECIES
      EVAP1 = EVAP1 + VFLUX/RHOL*DT
      DO 50 J=1,NPLANT
         TRANSP(J)=TRANSP(J) + TRNSP(J)*DT
   50 CONTINUE
      ETSUM = ETSUM + TRNSP(NPLANT+1)/RHOL*DT
!
!     SUM THE VAPOR FLUX OCCURRING IN, ABOVE AND BELOW THE SNOWPACK
      IF (NSP .GT. 0) THEN
         IF (NC .GT. 0) THEN
!           VAPOR FLUX ABOVE SNOWPACK IS FROM BOTTOM OF CANOPY
            TOPSNO = TOPSNO + QVC(NC)*DT
         ELSE
!           VAPOR FLUX ABOVE SNOWPACK IS FROM ATMOSPHERE
            TOPSNO = TOPSNO + VFLUX*DT
         END IF
         DO 60 I=1,NSP
            TQVSP(I) = TQVSP(I) + QVSP(I)*DT
   60    CONTINUE
!        TQVSP(NSP) IS THE VAPOR LEAVING BOTTOM OF SNOWPACK
      END IF
!
!     SUM THE ROOT EXTRACTION FROM EACH SOIL LAYER
      IF (NPLANT .GT. 0) THEN
         DO 70 I=1,NS
            ROOTXT(I)=ROOTXT(I)+XTRACT(I)*DT
   70    CONTINUE
      END IF
!
!     SUM THE TOTAL WATER FLUX BETWEEN SOIL NODES
      DO 80 I=1,NS-1
         TOTFLO(I)=TOTFLO(I) + (QSL(I) + QSV(I)/RHOL)*DT
   80 CONTINUE
!
!     SUM THE SEEPAGE FROM THE SURFACE NODE
      TSEEP=SEEP*DT
!
      RETURN
 END SUBROUTINE SUMDT


!***********************************************************************
!
 SUBROUTINE WBSOIL (N,NS,ZS,TS,TSDT,MAT,MATDT,VLC,VLCDT,VIC,VICDT,&
    CONC,CONCDT,ICESDT,QSL,QSV,XTRACT,SEEP,US,ITER,PRESUR,nsalt,col,row)
!
!     THIS SUBROUTINE CALCULATES THE COEFFICIENTS FOR THE SOIL IN THE
!     NEWTON-RAPHSON PROCEDURE OF THE ENERGY BALANCE
!
!***********************************************************************
    use controlpara_mod,only:WDT2d,DT2d
    use constvar_mod
    use soilproperty_mod,only:SAT2d,ENTRY2d
    use matrix_mod,only:A22d,B22d,C22d,D22d
    use savedata_mod,only:QSLT_WB2d,QSVT_WB2d
    implicit none                        
!input
    integer(i4),intent(in)::NS,col,row,ITER,nsalt
    integer(i4),intent(inout)::N
    real(r8),dimension(NSMAX),intent(in)::ZS,TS,TSDT,MAT,MATDT,VLC,VLCDT,VIC,VICDT,XTRACT,US
    real(r8),dimension(NSMAX),intent(inout)::QSL,QSV
    real(r8),dimension(NSALTMAX,NSMAX),intent(in)::CONC,CONCDT
    real(r8),intent(inout)::SEEP
    real(r8),intent(in)::PRESUR
    integer(i4),dimension(NSMAX),intent(in)::ICESDT
!temp
    real(r8),dimension(NSMAX)::QSLDT,HKT,HKDT,CON,CONT,QSVDT,DCONVP
    integer(i4),dimension(NSMAX)::NODE
    integer(i4)::I,J,IFLAG,ISAT,ITOP,JTOP
    real(r8)::QMAX,CLOSE,CLOSVLC,DLDM

!pointer
    real(r8),dimension(:),pointer::QSLT,QSVT
    real(r8),dimension(:),pointer::A2,B2,C2,D2
    real(r8),pointer::WDT,DT
    real(r8),dimension(:),pointer::SAT,ENTRY


    QSLT   =>   QSLT_WB2d(col,row,:)
    QSVT   =>   QSVT_WB2d(col,row,:)
    A2     =>   A22d(col,row,:)
    B2     =>   B22d(col,row,:)
    C2     =>   C22d(col,row,:)
    D2     =>   D22d(col,row,:)
    WDT    =>   WDT2d(col,row)
    DT     =>   DT2d(col,row)
    SAT    =>   SAT2d(col,row,:)
    ENTRY  =>   ENTRY2d(col,row,:)






!     INITIALIZE A2(1) IF SOIL IS FIRST MATERIAL
      IF (N .EQ. 1) A2(N)=0.0
!
!**** DETERMINE THE AVERAGE VAPOR FLUX BETWEEN SOIL NODES
      IF (ITER.EQ.1) CALL QVSOIL (NS,QSVT,ZS,TS,MAT,VLC,VIC,CONC,DCONVP,presur,nsalt,col,row)
      CALL QVSOIL (NS,QSVDT,ZS,TSDT,MATDT,VLCDT,VICDT,CONCDT,DCONVP,presur,nsalt,col,row)
      CALL WEIGHT (NS-1,QSV,QSVT,QSVDT,col,row)
!
!**** DETERMINE THE LIQUID MOISTURE FLUX BETWEEN NODES
      IF (ITER .EQ. 1) THEN
!****    CALCULATE FLUX AT BEGINNING OF TIME STEP
!****    DETERMINE THE HYDRAULIC CONDUCTIVITY OF EACH SOIL NODE
         CALL SOILHK (NS,HKT,MAT,VLC,VIC,col,row)
!****    DETERMINE THE CONDUCTANCE TERM BETWEEN SOIL NODES
         CALL CONDUC (NS,ZS,HKT,CONT,"WBSOIL HKT")
         CALL QLSOIL (NS,ZS,QSLT,CONT,MAT)
      END IF
      CALL SOILHK (NS,HKDT,MATDT,VLCDT,VICDT,col,row)
      CALL CONDUC (NS,ZS,HKDT,CON,"WBSOIL HKDT")
      CALL QLSOIL (NS,ZS,QSLDT,CON,MATDT)
!
!**** OBTAIN THE AVERAGE MOISTURE FLUX OVER THE TIME STEP
      CALL WEIGHT (NS-1,QSL,QSLT,QSLDT,col,row)
!
!     LIMIT NET WATER FLUX INTO NODES IF FLUX FILLS NODE BEYOND 
!     AVAILABLE POROSITY AS A RESULT DOWNWARD FLUX INTO LAYERS WITH
!     LIMITED PERMEABILITY (PERHAPS DUE TO ICE)
!     (START AT BOTTOM OF PROFILE AND WORK UP SO THAT ANY CHANGES IN THE
!     FLUXES CAN WORK THEIR UP THROUGH THE PROFILE)     
      IFLAG=0
      DO 5 I=NS-1,2,-1
         IF (VICDT(I).GT.0.0) THEN
            QMAX=(SAT(I) - VLC(I) - VIC(I))*(ZS(I+1) - ZS(I-1))/2./DT
            IF (QMAX.LT.0.0) QMAX=0.0
            IF (QSL(I-1)-QSL(I) .GT. QMAX) THEN
               IFLAG=1
               IF (QSL(I-1).GT.0.0) THEN
                  CON(I-1)=0.0
                  QSL(I-1)=QSL(I)+QMAX
                  IF (QSL(I-1).LT.0.0) THEN 
                     QSL(I-1)=0.0
                     QSL(I)=QSL(I-1)-QMAX
                     IF (QSL(I).GT.0.0) QSL(I)=0.0
                  END IF
                 ELSE
                  CON(I)=0.0
                  QSL(I)=QSL(I-1)-QMAX
                  IF (QSL(I).GT.0.0) QSL(I)=0.0
               END IF
            END IF
         END IF
    5 CONTINUE
      IF (VICDT(1).GT.0.0) THEN
         QMAX=(SAT(I) - VLC(1) - VIC(1))*(ZS(2) - ZS(1))/2./DT
         IF (QMAX.LT.0.0) QMAX=0.0
         IF (-QSL(1) .GT. QMAX) THEN
            IFLAG=1
            CON(1)=0.0
            QSL(1)=-QMAX
         END IF
      END IF
!
      IF (IFLAG.GT.0) THEN
      DO 10 I=2,NS-1
         IF (VICDT(I).GT.0.0) THEN
            QMAX=(SAT(I) - VLC(I) - VIC(I))*(ZS(I+1) - ZS(I-1))/2./DT
            IF (QMAX.LT.0.0) QMAX=0.0
            IF (QSL(I-1)-QSL(I) .GT. QMAX) THEN
               IF (QSL(I).LT.0.0) THEN
                  CON(I)=0.0
                  QSL(I)=QSL(I-1)-QMAX
                  IF (QSL(I).GT.0.0) THEN 
                     QSL(I)=0.0
                     QSL(I-1)=QSL(I)+QMAX
                     IF (QSL(I-1).LT.0.0) QSL(I-1)=0.0
                  END IF
                 ELSE
                  CON(I-1)=0.0
                  QSL(I-1)=QSL(I)+QMAX
                  IF (QSL(I-1).LT.0.0) QSL(I-1)=0.0
               END IF
            END IF
         END IF
   10 CONTINUE
      END IF
!
!     COUNT NUMBER OF SATURATED NODS
      ISAT=0
      !QSL(NS-1)=QSL(NS-1)*0.2
!
!**** DETERMINE THE COEFFICIENTS FOR THE SURFACE LAYER
      D2(N)=D2(N) - QSL(1) - QSV(1)/RHOL + US(1) - XTRACT(1)/RHOL - &
       (ZS(2)-ZS(1))*(VLCDT(1)-VLC(1)+RHOI/RHOL*(VICDT(1)-VIC(1)))/2./DT
      IF (ICESDT(1) .EQ. 1) THEN
!       ICE IS PRESENT -- CALCULATE COEFFICIENTS FOR ICE CONTENT
         B2(N)=-(ZS(2)-ZS(1))/(2.*DT)*RHOI/RHOL
         A2(N+1)=0.0
         SEEP=0.0
       ELSE
!        NO ICE IS PRESENT -- CALCULATE COEFFICIENTS FOR LIQUID BALANCE
         IF (MATDT(1).LT.0.0 .OR. D2(N).LT.0.0) THEN
!           NO SEEPAGE
            IF (MATDT(1) .GT. ENTRY(1)) THEN
!              NODE IS SATURATED - CHECK IF BOUNDARY IS POORLY DEFINED 
               IF (-B2(N).LT.CON(1)*1.E-6 .OR. CON(1).LE.0.0) THEN
!                 SOIL WATER FLOW OVERWHELMS SURFACE BOUNDARY CONDITION
                  ISAT=ISAT+1
                  NODE(ISAT)=1
               END IF
            END IF
            CALL MATVL3 (1,MATDT(1),VLCDT(1),DLDM,col,row)
            B2(N)=B2(N) - WDT*(CON(1)+DCONVP(1)/RHOL) &
                  - (ZS(2)-ZS(1))/(2.*DT)*DLDM
            A2(N+1)= WDT*(CON(1)+DCONVP(1)/RHOL) 
            SEEP=0.0
           ELSE 
!           ALLOW FOR SEEPAGE FROM SURFACE NODE - BOUNDARY IS PROPERLY 
!           DEFINED (DO NOT INCLUDE AS A SATURATED NODE)
            B2(N)=1.0
            A2(N+1)= WDT*CON(1)
!           NET FLOW INTO NODE IS EQUAL TO SEEPAGE FROM SURFACE
            SEEP=D2(N)
            D2(N)=0.0
         END IF
      END IF
!
!**** DETERMINE THE COEFFICIENTS FOR THE REST OF THE PROFILE
      DO 20 I=N+1,N+NS-2
         J=I-N+1
         D2(I)= QSL(J-1) - QSL(J) + (QSV(J-1)-QSV(J))/RHOL + US(J) &
                - XTRACT(J)/RHOL - (ZS(J+1)-ZS(J-1))&
                /(2.*DT)*(VLCDT(J)-VLC(J) + RHOI/RHOL*(VICDT(J)-VIC(J)))
         IF (ICESDT(J) .EQ. 1) THEN
!           ICE IS PRESENT -- CALCULATE COEFFICIENTS FOR ICE CONTENT
            B2(I)=-(ZS(J+1)-ZS(J-1))/(2.*DT)*RHOI/RHOL
            C2(I-1)=0.0
            A2(I+1)=0.0
          ELSE
!           NO ICE IS PRESENT -- CALCULATE COEFF. FOR LIQUID BALANCE
            IF (MATDT(J) .GT. ENTRY(J)) THEN
!              NODE IS SATURATED
               ISAT=ISAT+1
               NODE(ISAT)=J
            END IF                              
            CALL MATVL3 (J,MATDT(J),VLCDT(J),DLDM,col,row)
            B2(I)= -WDT*(CON(J-1)+CON(J)+(DCONVP(J-1)+DCONVP(J))/RHOL) &
                   -(ZS(J+1)-ZS(J-1))/(2.*DT)*DLDM
            C2(I-1)=WDT*(CON(J-1)+DCONVP(J-1)/RHOL)
            A2(I+1)=WDT*(CON(J)+DCONVP(J)/RHOL)
         END IF
   20 CONTINUE
!
!     ADJUST MATRIX IF SEEPAGE OCCURS
      IF (SEEP.GT.0.0) THEN
         A2(N)=0.0
         C2(N)=0.0
      END IF
!
      IF (ISAT .GT. 0) THEN
!        SATURATED NODES ARE PRESENT - CHECK FOR MATRIX SINGULARITY
         JTOP=0
         NODE(ISAT+1)=NS
         DO 30 J=1,ISAT
            I=NODE(J)+N-1
!           CHECK FOR MATRIX DISCONTINUITY AND SAVE TOP NODE
            IF (NODE(J).EQ.1) THEN
               JTOP=NODE(J)
              ELSE
!xxxx          IF (CON(NODE(J)-1).LE.0.0) JTOP=NODE(J)
!              CHECK IF CONDUCTIVITY OF NODE ABOVE IS SUFFICIENTLY SMALL
!              TO BE CONSIDERED ZERO
               IF (CON(NODE(J)-1).LE.-B2(I)*1.E-7) JTOP=NODE(J)
!              THIS IS TO CHECK IF THE SPECIFIC STORGAGE OF THE 
!              UNSATURATED NODE ABOVE THE SATURATED ZONE IS NEGLIBLE
               IF (-CON(NODE(J)-1)-B2(I-1).LE.-B2(I)*1.E-7 &
                    .AND. JTOP.EQ.0) JTOP=NODE(J)
            END IF
            IF (JTOP.NE.0) THEN
!xxxx          IF (CON(NODE(J)).LE.0.0) THEN
               IF (CON(NODE(J)).LE.-B2(I)*1.E-7) THEN
!                SATURATED NODES SURROUNDED BY ZERO CONDUCTIVITY; 
!                SOLUTION IS UNDEFINED - RESET B2(I) OF TOP SATURATED
!                LAYER AS IF IT IS JUST UNDER AIR ENTRY POTENTIAL
                ITOP=JTOP+N-1
                 CLOSE = 1.0001*ENTRY(JTOP)
                 CALL MATVL3 (JTOP,CLOSE,VLCDT(JTOP),DLDM,col,row)
                 IF (JTOP.EQ.1) THEN
                  B2(ITOP)=B2(ITOP)-(ZS(JTOP+1)-ZS(JTOP))/(2.*DT)*DLDM
                 ELSE
                  B2(ITOP)=B2(ITOP)-(ZS(JTOP+1)-ZS(JTOP-1))/(2.*DT)*DLDM
                 END IF
                 JTOP=0
               END IF
            END IF
!           IF NEXT NODE IS NOT SATURATED, RESET JTOP
            IF (NODE(J)+1.NE.NODE(J+1)) JTOP=0
   30    CONTINUE
      END IF
!
      N=N+NS-2
      RETURN
 END SUBROUTINE WBSOIL

















!***********************************************************************
!
 SUBROUTINE WBRES (N,NR,NSP,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,RHOR,&
        QVR,RHOSP,ICESDT,ITER,PRESUR,col,row)
!
!     THIS SUBOUTINE CALCULATES THE NEWTON-RAPHSON COEFFICIENTS FOR THE
!     WATER BALANCE OF THE RESIDUE LAYERS.
!
!***********************************************************************
    use controlpara_mod,only:DT2d,WDT2d
    use constvar_mod,only:G,UGAS,RHOL
    use matrix_mod,only:A22d,B22d,C22d,D22d
    use residu_mod,only:EVAP2d,EVAPK2d
    !COMMON /RESIDU/ EVAP(10),EVAPK(10)
    implicit none

!input
    integer(i4),intent(in)::NR,NSP,col,row,ITER,ICESDT
    integer(i4),intent(inout)::N
    real(r8),intent(in)::PRESUR
    real(r8),dimension(NRMAX),intent(in)::ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,RHOR
    real(r8),dimension(NRMAX),intent(inout)::QVR
    real(r8),dimension(NSPMAX),intent(in)::RHOSP

!temp
    real(r8),dimension(NRMAX)::CONV,VAPCON,CONVEC,TKRES
    integer(i4)::NR1,I,J

!pointer
    real(r8),dimension(:),pointer::A2,B2,C2,D2
    real(r8),dimension(:),pointer::EVAP,EVAPK
    real(r8),pointer::DT,WDT


    A2   =>  A22d(col,row,:)
    B2   =>  B22d(col,row,:)
    C2   =>  C22d(col,row,:)
    D2   =>  D22d(col,row,:)
    EVAP =>  EVAP2d(col,row,:)
    EVAPK=>  EVAPK2d(col,row,:)
    DT   =>  DT2d(col,row)
    WDT  =>  WDT2d(col,row)


!**** DETERMINE THE HEAT TRANSFER COEFFICIENT FOR EACH NODE
      CALL RESTK (NR,NSP,TKRES,CONVEC,TR,TRDT,GMC,GMCDT,RHOR,RHOSP,col,row)
!
!**** DETERMINE THE VAPOR TRANSPORT FROM THE THERMAL CONVECTION
      CALL RESVK (NR,NSP,TR,TRDT,CONVEC,VAPCON,presur,col,row)
      NR1=NR+1
      VAPCON(NR1)=VAPCON(NR)
!
!**** DETERMINE THE CONDUCTANCE TERM FOR CONVECTIVE VAPOR TRANSPORT
      CALL CONDUC (NR1,ZR,VAPCON,CONV,"WBRES VAPCON")
!
!**** DETERMINE THE VAPOR FLUX BETWEEN RESIDUE NODES
      CALL QVRES (NR,QVR,CONV,VAPR,VAPRDT,col,row)
!
!**** DETERMINE THE COEFFICIENTS FOR THE SURFACE LAYER
      C2(N)=WDT*CONV(1)
      IF (NR.GT.1) THEN
         B2(N)=B2(N) - WDT*(CONV(1) + EVAPK(1)) - (ZR(2)-ZR(1))/(2*DT)
         D2(N)= D2(N) - QVR(1)- EVAP(1) &
                      - (ZR(2)-ZR(1))/(2*DT)*(VAPRDT(1)-VAPR(1))
        ELSE
         B2(N)=B2(N) - WDT*(CONV(1) + EVAPK(1)) - (ZR(2)-ZR(1))/DT
         D2(N)= D2(N) - QVR(1)- EVAP(1) &
                      - (ZR(2)-ZR(1))/DT*(VAPRDT(1)-VAPR(1))
      END IF
!
!**** DETERMINE THE COEFFICIENTS FOR THE REST OF THE RESIDUE LAYERS
      DO 10 I=N+1,N+NR-1
         J=I-N+1
         A2(I)=WDT*CONV(J-1)
         C2(I)=WDT*CONV(J)
         IF (J .NE. NR) THEN
            B2(I)=-WDT*(CONV(J-1)+CONV(J)+EVAPK(J)) &
                  - (ZR(J+1)-ZR(J-1))/2./DT
            D2(I)= QVR(J-1) - QVR(J) - EVAP(J) &
                         - (ZR(J+1)-ZR(J-1))/(2*DT)*(VAPRDT(J)-VAPR(J))
          ELSE
            B2(I)=-WDT*(CONV(J-1)+CONV(J)+EVAPK(J)) &
                  - (ZR(J+1)-ZR(J)+(ZR(J)-ZR(J-1))/2)/DT
            D2(I)= QVR(J-1) - QVR(J) - EVAP(J) &
             - (ZR(J+1)-ZR(J)+(ZR(J)-ZR(J-1))/2)/DT*(VAPRDT(J)-VAPR(J))
         END IF 
   10 CONTINUE
      N=N+NR
!
!**** DETERMINE THE COEFFICIENTS SOIL SURFACE
      A2(N)=WDT*CONV(NR)/RHOL
      IF (ICESDT .EQ. 0) THEN
!        NO ICE IS PRESENT AT SOIL SURFACE
         B2(N)=-WDT*CONV(NR)*VAPRDT(NR+1) &
              *0.018*G/UGAS/(TRDT(NR+1)+273.16)/RHOL
         C2(N-1)=C2(N-1)*VAPRDT(NR+1)*.018*G/UGAS/(TRDT(NR+1)+273.16)
        ELSE
!        ICE IS PRESENT AT SOIL SURFACE - WATER BALANCE FOR ICE CONTENT
         B2(N)=0.0
         C2(N-1)=0.0
      END IF
      D2(N) = QVR(NR)/RHOL
      RETURN
      END subroutine


!***********************************************************************
!
 SUBROUTINE SNOWBC (N,NSP,NR,ZSP,QVSP,VAPSPT,TSDT,TSPDT,ICESDT,PRESUR,col,row)
!
!     THIS SUBROUTINE SETS UP THE UPPER BOUNDARY CONDITION FOR THE WATER
!     BALANCE OF THE RESIDUE-SOIL SYSTEM WHEN THERE IS A SNOWPACK
!***********************************************************************
    use controlpara_mod,only:WDT2d
    use constvar_mod,only:G,P0,RHOL,UGAS
    use matrix_mod,only:B22d,D22d
    use spwatr_mod,only:VAPSPX,VDIFSP
    implicit none

!input
    integer(i4),intent(in)::NSP,NR,col,row,ICESDT
    integer(i4),intent(inout)::N
    real(r8),dimension(NSPMAX),intent(in)::ZSP,QVSP,TSPDT
    real(r8),dimension(NSMAX),intent(in)::TSDT
    real(r8),intent(in)::VAPSPT,PRESUR

!temp
    real(r8)::CON,VDIF

!pointer
    real(r8),dimension(:),pointer::B2,D2
    real(r8),pointer::WDT


    B2=>B22d(col,row,:)
    D2=>D22d(col,row,:)
    WDT=>WDT2d(col,row)

!
!**** DETERMINE THE VAPOR DIFFUSIVITY OF BOTTOM SNOWPACK NODE
      VDIF= VDIFSP*(P0/PRESUR)*(1.+TSPDT(NSP)/273.16)**VAPSPX
      CON = VDIF/(ZSP(NSP+1)-ZSP(NSP))
!
!**** DEFINE COEFFICIENTS
      IF (NR .GT. 0) THEN
!        CALCULATE BOUNDARY COEFFICIENTS FOR SURFACE RESIDUE NODE
         B2(N) = -WDT*CON
         D2(N) = QVSP(NSP)
       ELSE
!        CALCULATE BOUNDARY COEFFICIENTS FOR SURFACE SOIL NODE
         IF (ICESDT .EQ. 0) THEN
!           NO ICE PRESENT AT SOIL SURFACE - MATRIC POTENTIAL UNKNOWN
            B2(N) = -WDT*CON*VAPSPT*0.018*G/(UGAS*(TSDT(1)+273.16))/RHOL
           ELSE
!           ICE PRESENT AT SOIL SURFACE - ICE CONTENT UNKNOWN
            B2(N) = 0.0
         END IF
         D2(N) = QVSP(NSP)/RHOL
      END IF
      RETURN
    END SUBROUTINE SNOWBC









!***********************************************************************
!
 SUBROUTINE WBCAN (N,NPLANT,NC,NSP,NR,NS,ZC,TC,TCDT,TLC,TLCDT,&
    VAPC,VAPCDT,WCAN,WCANDT,PCAN,PCANDT,QVC,TRNSP,MAT,MATDT,XTRACT,&
    SWCAN,LWCAN,CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,ITYPE,&
    ICESDT,ITER,col,row,julian,errorflag)
!
!     THIS SUBROUTINE CALCULATES THE JACOBIAN MATRIX COEFFICIENTS FOR
!     THE CANOPY PORTION OF THE NEWTON-RAPHSON SOLUTION OF THE ENERGY
!     BALANCE
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d,DT2d
    use constvar_mod,only:CA,G,LV,RHOA,RHOL,UGAS
    use matrix_mod,only:B22d,D22d,A22d,C22d
    use writeit_mod,only:xlenc2d
    use savedata_mod,only:CONT_WB2d
    implicit none
!input
    integer(i4),intent(inout)::N
    integer(i4),intent(in)::NPLANT,NC,NSP,NR,NS,col,row,ITER,julian
    real(r8),dimension(NCMAX),intent(in)::ZC,TC,TCDT,VAPC,VAPCDT
    real(r8),dimension(NPMAX,NCMAX-1),intent(in)::TLC
    real(r8),dimension(NPMAX,NCMAX-1),intent(inout)::TLCDT
    real(r8),dimension(NCMAX-1),intent(inout)::WCAN,WCANDT
    real(r8),dimension(NCMAX-1),intent(inout)::QVC
    real(r8),dimension(NPMAX),intent(in)::PCAN
    real(r8),dimension(NPMAX),intent(inout)::PCANDT
    real(r8),dimension(NPMAX+1),intent(inout)::TRNSP
    real(r8),dimension(NSMAX),intent(in)::MAT,MATDT
    real(r8),dimension(NSMAX),intent(inout)::XTRACT
    real(r8),dimension(NPMAX+1,NCMAX-1),intent(in)::SWCAN,LWCAN
    real(r8),intent(in)::CANMA,CANMB
    real(r8),dimension(NPMAX),intent(in)::DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0
    integer(i4),dimension(NPMAX),intent(in)::ITYPE
    integer(i4),intent(in)::ICESDT
    integer(i4),intent(inout)::errorflag

!temp
    real(r8),dimension(NCMAX-1)::CON,CONDT,HEATC,ETLYR,DETLYR,DHEATC,DTLDTC


    real(r8),dimension(:),pointer::CONT
    real(r8),pointer::xlenc
    real(r8),dimension(:),pointer::B2,D2,A2,C2
    real(r8),pointer::DT,WT,WDT

    integer(i4)::I,J

    CONT  =>CONT_WB2d(col,row,:)
    xlenc =>xlenc2d(col,row)
    B2    =>B22d(col,row,:)
    A2    =>A22d(col,row,:)
    C2    =>C22d(col,row,:)
    D2    =>D22d(col,row,:)
    DT    =>DT2d(col,row)
    WDT   =>WDT2d(col,row)
    WT    =>WT2d(col,row)



!
!**** CALCULATE TRANSPIRATION FROM CANOPY

      CALL LEAFT2 (NPLANT,NC,NS,ITER,ITYPE,TC,TCDT,TLC,TLCDT,VAPC,VAPCDT, &
       WCAN,WCANDT,PCAN,PCANDT,MAT,MATDT,TRNSP,XTRACT,SWCAN,LWCAN, &
       HEATC,ETLYR,DETLYR,CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,&
       DHEATC,DTLDTC,col,row,julian,errorflag)
!
!**** DETERMINE THE EDDY CONDUCTANCE TERM BETWEEN NODES
      IF (ITER .EQ. 1) CALL CANTK (NC,CONT,TC,ZC,col,row)
      CALL CANTK (NC,CONDT,TCDT,ZC,col,row)
!
!**** CALCULATE THE AVERAGE CONDUCTANCE TERM OVER THE TIME STEP
      CALL WEIGHT (NC,CON,CONT,CONDT,col,row)
!
!**** CONVERT THERMAL EDDY CONDUCTANCE TO VAPOR CONDUCTANCE
!     AND CALCULATE LIQUID FLUX BETWEEN LAYERS
      DO 5 I=1,NC
         CON(I)=CON(I)/RHOA/CA
         QVC(I)=CON(I)*(WDT*(VAPCDT(I)-VAPCDT(I+1)) &
                       + WT*(VAPC(I)-VAPC(I+1)))
    5 CONTINUE
!
!**** DETERMINE THE MATRIX COEFFICIENTS FOR THE TOP LAYER
!
      A2(N+1)=WDT*CON(1)
      B2(N)=B2(N) - WDT*(CON(1) + DETLYR(1)) - (ZC(2)-ZC(1))/2/DT
      C2(N)=WDT*CON(1)
      D2(N)=D2(N)-QVC(1)+ETLYR(1)-(ZC(2)-ZC(1))*(VAPCDT(1)-VAPC(1))/2/DT
!
!**** DETERMINE THE COEFFICIENTS FOR THE REST OF THE LAYERS
      DO 10 I=N+1,N+NC-1
         J=I-N+1
         A2(I+1)=WDT*CON(J)
         B2(I)=-WDT*(CON(J-1) + CON(J) + DETLYR(J)) &
               - (ZC(J+1)-ZC(J-1))/2/DT
         C2(I)=WDT*CON(J)
         D2(I)=QVC(J-1) - QVC(J) + ETLYR(J) &
              - (ZC(J+1)-ZC(J-1))*(VAPCDT(J)-VAPC(J))/2/DT
   10 CONTINUE
      N=N+NC
!
!**** DETERMINE THE COEFFICIENTS FOR THE TOP LAYER OF THE NEXT MATERIAL
!
      IF (NSP.GT.0) THEN
!****    NEXT MATERIAL IS SNOW -- NO WATER BALANCE MATRIX FOR SNOW
         A2(N) = 0.0
         B2(N) = 0.0
         C2(N-1)=0.0
         D2(N) = 0.0
        ELSE
!
      IF (NR.GT.0) THEN
!****    NEXT MATERIAL IS RESIDUE
         B2(N) = -WDT*CON(NC)
         D2(N) = QVC(NC)
        ELSE
!
!**** NEXT MATERIAL IS SOIL -- WATER BALANCE BASED ON WATER CONTENT AND
!     MATRIC POTENTIAL
      A2(N) = A2(N)/RHOL
      IF (ICESDT .EQ. 0) THEN
!       NO ICE PRESENT AT SOIL SURFACE - WATER BALANCE FOR MATRIC POT.
        B2(N) = -WDT*CON(NC)*VAPCDT(NC+1)*0.018*G &
                 /(UGAS*(TCDT(NC+1)+273.16)*RHOL)
        C2(N-1)=C2(N-1)*VAPCDT(NC+1)*0.018*G/(UGAS*(TCDT(NC+1)+273.16))
       ELSE
!       ICE IS PRESENT AT SOIL SURFACE - WATER BALANCE FOR ICE CONTENT
        B2(N)=0.0
        C2(N-1)=0.0
      END IF
      D2(N) = QVC(NC)/RHOL
!
      END IF
      END IF
!XXXX
      xlenc = LV*qvc(nc)
!
      RETURN
 END SUBROUTINE WBCAN



!
!***********************************************************************
 SUBROUTINE ADJUST(I,VICDT,VLCDT,MATDT,CONCDT,TSDT,SALTDT,ICESDT,nsalt,col,row)
!     THIS SUBROUTINE ADJUSTS THE SOIL LAYER TO ACCOUNT FOR THE LATENT
!     HEAT OF FUSION ON THE FIRST ITERATION THAT THE LAYER CROSSES THE
!     FREEZING POINT.
!***********************************************************************
    !use controlpara_mod,only:WT2d,WDT2d,DT2d
    use constvar_mod,only:G,LF,RHOI,RHOL,UGAS
    use savedata_mod,only:CS2d
    use soilproperty_mod,only:RHOB2d,SALTKQ2d
    implicit none
!input
    integer(i4),intent(in)::I,col,row,nsalt
    integer(i4),dimension(NSMAX),intent(inout)::ICESDT
    real(r8),dimension(NSMAX),intent(inout)::VICDT,VLCDT,MATDT,TSDT
    real(r8),dimension(NSALTMAX,NSMAX),intent(inout)::SALTDT
    real(r8),dimension(NSALTMAX,NSMAX),intent(inout)::CONCDT
!temp
    real(r8)::TLCONC,TOTPOT,TMPFRZ,TMPF,ENERGY,EFFCS,DLDT,DUMMY
    integer(i4)::J

    real(r8),dimension(:),pointer::CS,rhob
    real(r8),dimension(:,:),pointer::SALTKQ
    !real(r8),pointer::WT,WDT,DT



    CS=>CS2d(col,row,:)
    rhob=>rhob2d(col,row,:)
    SALTKQ=>SALTKQ2d(col,row,:,:)
    !WT=>WT2d(col,row)
    !WDT=>WDT2d(col,row)
    !DT=>DT2d(col,row)



!**** INITIALLY ASSUME THAT ALL WATER IS LIQUID
      VLCDT(I)=VLCDT(I) + (RHOI/RHOL)*VICDT(I)
      CALL MATVL1 (I,MATDT(I),VLCDT(I),col,row)
      DO 5 J=1,NSALT
         CONCDT(J,I)=SALTDT(J,I)/(SALTKQ(J,I) + VLCDT(I)*RHOL/RHOB(I))
    5 CONTINUE
      VICDT(I)=0.0
!
!**** CALCULATE THE FREEZING POINT TEMPERATURE ACCORDING TO TOTAL WATER
!     POTENTIAL AT THE END OF THE TIME STEP
      TLCONC=0.0
      DO 10 J=1,NSALT
         TLCONC= TLCONC+CONCDT(J,I)
   10 CONTINUE
      TOTPOT=MATDT(I) - TLCONC*UGAS*273.16/G
      TMPFRZ=273.16*TOTPOT/(LF/G-TOTPOT)
      TMPF=TMPFRZ+273.16
!
!**** CALCULATE THE ENERGY AVAILABLE TO FREEZE WATER
      ENERGY=CS(I)*(TMPFRZ-TSDT(I))
!
!**** CALCULATE THE EFFECTIVE HEAT CAPACITY INCLUDING LATENT HEAT TERM
      CALL FSLOPE (I,DLDT,DUMMY,TMPF,MATDT,CONCDT,VLCDT,nsalt,col,row)
      EFFCS=CS(I)+RHOL*LF*DLDT
!
!**** CALCULATE THE TEMPERATURE AT THE END OF THE TIME STEP
      TSDT(I)=TMPFRZ - ENERGY/EFFCS
!
!**** CALL SUBROUTINE FROZEN TO DETERMINE THE LIQUID, ICE, AND SOLUTES
      CALL FROZEN (I,VLCDT,VICDT,MATDT,CONCDT,TSDT,SALTDT,ICESDT,nsalt,col,row)
      RETURN
 END SUBROUTINE ADJUST


!***********************************************************************
!
 SUBROUTINE BACKUP (NS,NR,NSP,NC,NPLANT,NSALT,ICES,ICESDT,TS,TSDT, &
        MAT,MATDT,CONC,CONCDT,VLC,VLCDT,VIC,VICDT,SALT,SALTDT, &
        TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,TC,TCDT,TLC,TLCDT,VAPC,VAPCDT, &
        WCAN,WCANDT,PCAN,PCANDT,TSP,TSPDT,DLW,DLWDT,ICESP,ICESPT)
!
!     THIS SUBROUTINE SETS END OF TIME STEP VALUES BACK TO THOSE FOR
!     THE BEGINNING OF THE TIME STEP.
!
!***********************************************************************
    implicit none
!input
    integer(i4),intent(in)::NS,NR,NSP,NC,NPLANT,NSALT
    integer(i4),dimension(NSMAX),intent(in)::ICES
    integer(i4),dimension(NSMAX),intent(inout)::ICESDT
    integer(i4),dimension(NSPMAX),intent(in)::ICESP
    integer(i4),dimension(NSPMAX),intent(inout)::ICESPT
    real(r8),dimension(NSALTMAX,NSMAX),intent(in)::CONC,SALT
    real(r8),dimension(NSALTMAX,NSMAX),intent(inout)::CONCDT,SALTDT
    real(r8),dimension(NSMAX),intent(in)::TS,MAT,VLC,VIC
    real(r8),dimension(NSMAX),intent(inout)::TSDT,MATDT,VLCDT,VICDT
    real(r8),dimension(NRMAX),intent(in)::TR,VAPR,GMC
    real(r8),dimension(NRMAX),intent(inout)::TRDT,VAPRDT,GMCDT
    real(r8),dimension(NCMAX),intent(in)::TC,VAPC
    real(r8),dimension(NCMAX),intent(inout)::TCDT,VAPCDT
    real(r8),dimension(NPMAX,NCMAX-1),intent(in)::TLC
    real(r8),dimension(NPMAX,NCMAX-1),intent(inout)::TLCDT
    real(r8),dimension(NCMAX-1),intent(in)::WCAN
    real(r8),dimension(NCMAX-1),intent(inout)::WCANDT
    real(r8),dimension(NPMAX),intent(in)::PCAN
    real(r8),dimension(NPMAX),intent(inout)::PCANDT
    real(r8),dimension(NSPMAX),intent(in)::TSP,DLW
    real(r8),dimension(NSPMAX),intent(inout)::TSPDT,DLWDT

!temp variable
    integer(i4)::I,J

!     SOIL PROPERTIES
      DO 15 I=1,NS
         TSDT(I)=TS(I)
         MATDT(I)=MAT(I)
         VLCDT(I)=VLC(I)
         VICDT(I)=VIC(I)
         ICESDT(I)=ICES(I)
         DO 10 J=1,NSALT
            SALTDT(J,I)=SALT(J,I)
            CONCDT(J,I)=CONC(J,I)
   10    CONTINUE
   15 CONTINUE
!
!     RESIDUE PROPERTIES
      DO 20 I=1,NR
         TRDT(I)=TR(I)
         VAPRDT(I)=VAPR(I)
         GMCDT(I)=GMC(I)
   20 CONTINUE
!
!     SNOWPACK PROPERTIES
      DO 30 I=1,NSP
         ICESPT(I) = ICESP(I)
         TSPDT(I) = TSP(I)
         DLWDT(I) = DLW(I)
   30 CONTINUE
!
!     CANOPY PROPERTIES
      DO 40 I=1,NC
         TCDT(I)=TC(I)
         VAPCDT(I)=VAPC(I)
         WCANDT(I)=WCAN(I)
         DO 35 J=1,NPLANT
            TLCDT(J,I)=TLC(J,I)
   35    CONTINUE
   40 CONTINUE
      IF (NC .GT. 0) THEN
         DO 45 J=1,NPLANT
            PCANDT(J)=PCAN(J)
   45    CONTINUE
      END IF
!
      RETURN
      END SUBROUTINE BACKUP


! no problem checked changed
!***********************************************************************
 SUBROUTINE TDMA (N,A,B,C,D,X,mark)
!
!     THIS SUBROUTINE SOLVES A TRI-DIAGONAL MATRIX OF THE FORM :
!
!              | B1  C1   0   0  .  .   . | |X1|   |D1|
!              | A2  B2  C2   0  .  .   . | |X2|   |D2|
!              |  0  A3  B3  C3  .  .   . | |X3| = |D3|
!              |  .   .   .   .     .   . | | .|   | .|
!              |  0   0   0   0  . AN  BN | |XN|   |DN|
!
!***********************************************************************
    implicit none
    integer(i4),intent(in)::N
    integer(i4),intent(inout)::mark
    real(r8),dimension(N),intent(inout)::A,B,C,D,X
    integer(i4)::I
      DO 10 I=2,N
         C(I-1)=C(I-1)/B(I-1)
         D(I-1)=D(I-1)/B(I-1)
         B(I)=B(I)-A(I)*C(I-1)
         D(I)=D(I)-A(I)*D(I-1)
   10 CONTINUE
      X(N)=D(N)/B(N)
      DO 20 I=N-1,1,-1
         X(I)=D(I)-C(I)*X(I+1)
   20 CONTINUE
      RETURN
 END SUBROUTINE TDMA
 
 
 
    subroutine TDMA2(N,DL,DM,DU,RS,X,mark)
      implicit none
      integer(i4),intent(in)::N
      real(r8),intent(inout)::DL(N)
      real(r8),intent(inout)::DM(N)
      real(r8),intent(inout)::DU(N)
      real(r8),intent(inout)::RS(N)
      real(r8),intent(inout)::X(N)
      integer(i4),intent(inout)::mark
      integer(i4)::markmark,i
      real(r8)::row,d
        MARK = 0
        markmark=0
        row  = abs(DM(1)) + abs(DU(1))
        if (N.Ge.3 .and. row .ne. 0.0)then
           D=1.0E0/ROW
           if(abs(DM(1))*D .gt. 1.0E-20) then
!          checking for strong nonsingularity with N=1
!          factoring A while checking for strong nonsingularity
             DL(1) = 0.0E0
             DU(N) = 0.0E0
             DU(1) = DU(1)/DM(1)
             do I=2,N,1
               row = abs(DL(I)) + abs(DM(I)) + abs(DU(I))
               if (row .eq. 0.0E0) then
                 markmark=1
                 exit
               end if
               D = 1.0E0/ROW
               DM(I) = DM(I) - DL(I) * DU(I-1)
               if (ABS(DM(I))*D .LE. 1.0E-20) then
                 markmark=1
                 exit
               end if
               if (I .LT. N) then
                 DU(I) = DU(I)/DM(I)
               endif
             enddo
           if (markmark.ne.1) then
             MARK=1
           else
             mark=0
           endif
        end if
      end if
!     If MARK = 1, update the right hand side and solve via backsubstitution
      if (MARK .eq. 1) then
        RS(1) = RS(1)/DM(1)
        do I=2,N,1
          RS(I) = (RS(I) - DL(I) * RS(I-1)) / DM(I)
        enddo
!       backsubstitution
        X(N) = RS(N)
        do I=N-1,1,-1
          X(I) = RS(I) - DU(I) * X(I+1)
        enddo
      end if
      return
    end subroutine TDMA2



! no problem checked changed
!***********************************************************************
!
 SUBROUTINE EBSOIL (N,NS,ZS,TS,TSDT,MAT,MATDT,CONC,CONCDT,&
    VLC,VLCDT,VIC,VICDT,ICES,ICESDT,QSL,QSV,SS,GFLUX,ITER,PRESUR,nsalt,col,row)
!     THIS SUBROUTINE CALCULATES THE COEFFICIENTS FOR THE SOIL IN THE
!     NEWTON-RAPHSON PROCEDURE OF THE ENERGY BALANCE
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d,DT2d,itmpbc
    use constvar_mod
    use matrix_mod,only:A12d,B12d,C12d,D12d
    use soilproperty_mod,only:RHOB2d,SALTKQ2d,SAT2d
    use savedata_mod,only:QSLT2d,TKT2d,CST2d,QSVT2d,CS2d
    implicit none
!input
    integer(i4),intent(inout)::N
    integer(i4),intent(in)::NS,col,row,ITER,nsalt
    real(r8),intent(in)::PRESUR
    real(r8),dimension(NSMAX),intent(in)::ZS
    real(r8),dimension(NSMAX),intent(in)::TS,TSDT,MAT,MATDT,VLC,VLCDT,VIC,VICDT
    real(r8),dimension(NSALTMAX,NSMAX),intent(in)::CONC,CONCDT
    integer(i4),dimension(NSMAX),intent(inout)::ICES,ICESDT
    real(r8),dimension(NSMAX),intent(inout)::QSL
    real(r8),dimension(NSMAX),intent(inout)::QSV
    real(r8),dimension(NSMAX),intent(in)::SS
    real(r8),intent(inout)::GFLUX
!temp
    real(r8),dimension(NSMAX)::QSLDT,TK,TKDT,CON,CSDT,HKT,HKDT,CONH,CONHT,QSVDT,DCONVP,DELTA,EPSLON
    integer(i4)::IFLAG,I,J,K
    real(r8)::QMAX,DLDM,DLDT,DLDTDT,DUMMY,T,TLCONC,TDT,TMPFRZ
!2d variable
    real(r8),dimension(:),pointer::A1,B1,C1,D1
    real(r8),dimension(:),pointer::CS
    real(r8),dimension(:),pointer::QSLT,TKT,CST,QSVT
    real(r8),pointer::DT,WT,WDT
    real(r8),dimension(:),pointer::SAT
    real(r8)::tmp1,tmp2

    A1   =>A12d(col,row,:)
    B1   =>B12d(col,row,:)
    C1   =>C12d(col,row,:)
    D1   =>D12d(col,row,:)
    CS   =>CS2d(col,row,:)
    QSLT =>QSLT2d(col,row,:)
    TKT  =>TKT2d(col,row,:)
    CST  =>CST2d(col,row,:)
    QSVT =>QSVT2d(col,row,:)
    DT   =>DT2d(col,row)
    WT   =>WT2d(col,row)
    WDT  =>WDT2d(col,row)
    SAT  =>SAT2d(col,row,:)




!**** DETERMINE THE AVERAGE VAPOR FLUX BETWEEN SOIL NODES
      IF (ITER.EQ.1) CALL QVSOIL (NS,QSVT,ZS,TS,MAT,VLC,VIC,CONC,DCONVP,presur,nsalt,col,row)
      CALL QVSOIL (NS,QSVDT,ZS,TSDT,MATDT,VLCDT,VICDT,CONCDT,DCONVP,presur,nsalt,col,row)
      CALL WEIGHT (NS-1,QSV,QSVT,QSVDT,col,row)
!
!**** DETERMINE THE LIQUID MOISTURE FLUX BETWEEN NODES
      IF (ITER .EQ. 1) THEN
!****    CALCULATE FLUX AT BEGINNING OF TIME STEP
!****    DETERMINE THE HYDRAULIC CONDUCTIVITY OF EACH SOIL NODE
         CALL SOILHK (NS,HKT,MAT,VLC,VIC,col,row)
!****    DETERMINE THE CONDUCTANCE TERM BETWEEN SOIL NODES
         CALL CONDUC (NS,ZS,HKT,CONHT,"EBSOIL HKT")
         CALL QLSOIL (NS,ZS,QSLT,CONHT,MAT)
      END IF
      CALL SOILHK (NS,HKDT,MATDT,VLCDT,VICDT,col,row)
      CALL CONDUC (NS,ZS,HKDT,CONH,"EBSOIL HKDT")
      CALL QLSOIL (NS,ZS,QSLDT,CONH,MATDT)
!
!**** OBTAIN THE AVERAGE MOISTURE FLUX OVER THE TIME STEP
      CALL WEIGHT (NS-1,QSL,QSLT,QSLDT,col,row)
!
!     LIMIT NET WATER FLUX INTO NODES IF FLUX FILLS NODE BEYOND
!     AVAILABLE POROSITY AS A RESULT DOWNWARD FLUX INTO LAYERS WITH
!     LIMITED PERMEABILITY (PERHAPS DUE TO ICE)
!     (START AT BOTTOM OF PROFILE AND WORK UP SO THAT ANY CHANGES IN THE
!     FLUXES CAN WORK THEIR UP THROUGH THE PROFILE)
      IFLAG=0
      DO 5 I=NS-1,2,-1
         IF (VICDT(I).GT.0.0) THEN
            QMAX=(SAT(I) - VLC(I) - VIC(I))*(ZS(I+1) - ZS(I-1))/2./DT
            IF (QMAX.LT.0.0) QMAX=0.0
            IF (QSL(I-1)-QSL(I) .GT. QMAX) THEN
               IFLAG=1
               IF (QSL(I-1).GT.0.0) THEN
                  CONH(I-1)=0.0
                  QSL(I-1)=QSL(I)+QMAX
                  IF (QSL(I-1).LT.0.0) THEN
                     QSL(I-1)=0.0
                     QSL(I)=QSL(I-1)-QMAX
                     IF (QSL(I).GT.0.0) QSL(I)=0.0
                  END IF
                 ELSE
                  CONH(I)=0.0
                  QSL(I)=QSL(I-1)-QMAX
                  IF (QSL(I).GT.0.0) QSL(I)=0.0
               END IF
            END IF
         END IF
!        SET UP LIQUID FLUX COEFFICIENTS SO THAT HEAT IS CARRIED IN
!        DIRECTION OF MOISTURE MOVEMENT.
         IF (QSL(I) .GT. 0.0) THEN
            DELTA(I)=1.0
            EPSLON(I)=0.0
          ELSE
            DELTA(I)=0.0
            EPSLON(I)=-1.0
         END IF
    5 CONTINUE
!
      IF (QSL(1) .GT. 0.0) THEN
         DELTA(1)=1.0
         EPSLON(1)=0.0
       ELSE
         DELTA(1)=0.0
         EPSLON(1)=-1.0
      END IF
!
      IF (VICDT(1).GT.0.0) THEN
         QMAX=(SAT(I) - VLC(1) - VIC(1))*(ZS(2) - ZS(1))/2./DT
         IF (QMAX.LT.0.0) QMAX=0.0
         IF (-QSL(1) .GT. QMAX) THEN
            IFLAG=1
            CONH(1)=0.0
            QSL(1)=-QMAX
         END IF
      END IF
!
      IF (IFLAG.GT.0) THEN
      DO 10 I=2,NS-1
         IF (VICDT(I).GT.0.0) THEN
            QMAX=(SAT(I) - VLC(I) - VIC(I))*(ZS(I+1) - ZS(I-1))/2./DT
            IF (QMAX.LT.0.0) QMAX=0.0
            IF (QSL(I-1)-QSL(I) .GT. QMAX) THEN
               IF (QSL(I).LT.0.0) THEN
                  CONH(I)=0.0
                  QSL(I)=QSL(I-1)-QMAX
                  IF (QSL(I).GT.0.0) THEN
                     QSL(I)=0.0
                     QSL(I-1)=QSL(I)+QMAX
                     IF (QSL(I-1).LT.0.0) QSL(I-1)=0.0
                  END IF
                 ELSE
                  CONH(I-1)=0.0
                  QSL(I-1)=QSL(I)+QMAX
                  IF (QSL(I-1).LT.0.0) QSL(I-1)=0.0
               END IF
            END IF
         END IF
   10 CONTINUE
      END IF
!
!**** CALCULATE THERMAL CONDUCTIVITY
      IF (ITER .EQ. 1) CALL SOILTK (NS,TKT,VLC,VIC,col,row)
      CALL SOILTK (NS,TKDT,VLCDT,VICDT,col,row)
!
!**** OBTAIN AVERAGE CONDUCTIVITY OVER THE TIME STEP
      CALL WEIGHT (NS,TK,TKT,TKDT,col,row)
!
!**** CALCULATE CONDUCTANCE TERMS BETWEEN NODES
      CALL CONDUC (NS,ZS,TK,CON,"EBSOIL TK")
!
!**** CALCULATE THE EFFECTIVE SPECIFIC HEAT AT EACH NODE
      IF (ITER .EQ. 1) CALL SOILHT (NS,CST,VLC,VIC,TS,MAT,CONC,nsalt,col,row)
      CALL SOILHT (NS,CSDT,VLCDT,VICDT,TSDT,MATDT,CONCDT,nsalt,col,row)
!
!**** OBTAIN AVERAGE SPECIFIC HEAT OVER THE TIME STEP
      CALL WEIGHT (NS,CS,CST,CSDT,col,row)
!
!
!**** DETERMINE THE COEFFICIENTS FOR THE SURFACE LAYER
!
      A1(N+1) = WDT*(CON(1) + DELTA(1)*(RHOL*CL*QSL(1)+CV*QSV(1)))
      B1(N) = B1(N)- WDT*(CON(1) + EPSLON(1)*(RHOL*CL*QSL(1)+CV*QSV(1))) &
        - (ZS(2)-ZS(1))*CSDT(1)/(2.*DT)
      C1(N) = WDT*(CON(1) + EPSLON(1)*(RHOL*CL*QSL(1)+CV*QSV(1)))
      D1(N) = D1(N)-(CON(1) + EPSLON(1)*(RHOL*CL*QSL(1)+CV*QSV(1))) &
            *(WT*(TS(1)-TS(2))+ WDT*(TSDT(1)-TSDT(2))) &
            - LV*QSV(1) + SS(1) - (ZS(2)-ZS(1))/(2.*DT) &
            *(CS(1)*(TSDT(1)-TS(1)) - RHOI*LF*(VICDT(1)-VIC(1)))
!
!     COMPUTE SURFACE GROUND HEAT FLUX FOR ENERGY OUTPUT
      GFLUX = (CON(1) + EPSLON(1)*(RHOL*CL*QSL(1)+CV*QSV(1))) &
        *(WT*(TS(1)-TS(2))+ WDT*(TSDT(1)-TSDT(2))) &
        + LV*QSV(1) + (ZS(2)-ZS(1))/(2.*DT) &
        *(CS(1)*(TSDT(1)-TS(1)) - RHOI*LF*(VICDT(1)-VIC(1)))
!
!**** CHECK IF ICE IS PRESENT AT THE END OF THE TIME STEP
      IF (ICESDT(1) .NE. 1) GO TO 20
!
!**** ICE IS PRESENT IN LAYER - - ADJUST COEFFICIENTS FOR LATENT HEAT
!**** TRANSFER AND THE SLOPE OF THE WATER CONTENT-TEMPERATURE CURVE
      IF (ICES(1) .EQ. 1) THEN
!        ICE IS PRESENT FOR THE ENTIRE TIME STEP
         T=273.16+TS(1)
        ELSE
!        ICE IS PRESENT ONLY AT THE END OF THE TIME STEP - DETERMINE
!        THE TEMPERATURE AT WHICH THE SOIL WILL BEGIN TO FREEZE
         TLCONC=0.0
         DO 15 K=1,NSALT
            TLCONC=TLCONC + CONC(K,1)
   15    CONTINUE
         TMPFRZ=273.16*LF/G/(LF/G-MAT(1)+TLCONC*UGAS*(TS(1)+273.17)/G)
         T=TMPFRZ
      END IF
      TDT=273.16+TSDT(1)
!     DETERMINE SLOPE OF LIQUID CONTENT-TEMPERATURE CURVE
      CALL FSLOPE (1,DLDTDT,DUMMY,TDT,MATDT,CONCDT,VLCDT,nsalt,col,row)
      CALL FSLOPE (1,DLDT,DUMMY,T,MAT,CONC,VLC,nsalt,col,row)
!     ENTER MATVLC FOR SLOPE OF LIQUID-MATRIC POTENTIAL CURVE
      CALL MATVL3 (1,MATDT(1),VLCDT(1),DLDM,col,row)
      B1(N)=B1(N) - (ZS(2)-ZS(1))/(2.*DT)* 0.5*RHOL*LF*(DLDT+DLDTDT) &
            - WDT*RHOL*LF*CONH(1)*DLDTDT/DLDM
!
!**** DETERMINE THE MATRIX COEFFICIENTS FOR THE REST OF THE PROFILE
   20 DO 30 I=N+1,N+NS-3
         J=I-N+1
         A1(I+1)=WDT*(CON(J) + DELTA(J)*(RHOL*CL*QSL(J)+CV*QSV(J)))
         B1(I)= -WDT*(CON(J-1)+CON(J) &
            +(DELTA(J-1)*(RHOL*CL*QSL(J-1)+CV*(QSV(J-1))) &
            + EPSLON(J)*(RHOL*CL*QSL(J)+CV*QSV(J)))) &
            - (ZS(J+1)-ZS(J-1))*CSDT(J)/(2.*DT)
         C1(I)=WDT*(CON(J) + EPSLON(J)*(RHOL*CL*QSL(J)+CV*QSV(J)))
         D1(I)=(CON(J-1) + DELTA(J-1)*(RHOL*CL*QSL(J-1) + CV*QSV(J-1))) &
            *(WT*(TS(J-1)-TS(J)) + WDT*(TSDT(J-1)-TSDT(J))) &
            -(CON(J) + EPSLON(J)*(RHOL*CL*QSL(J) + CV*QSV(J))) &
            *(WT*(TS(J)-TS(J+1)) + WDT*(TSDT(J)-TSDT(J+1))) &
            -LV*(QSV(J)-QSV(J-1)) + SS(J) &
            -(ZS(J+1)-ZS(J-1))/(2.*DT) &
            *(CS(J)*(TSDT(J)-TS(J)) - RHOI*LF*(VICDT(J)-VIC(J)))
!
!****    CHECK IF ICE IS PRESENT AT THE END OF THE TIME STEP
         IF (ICESDT(J) .NE. 1) GO TO 30
!
!****    ICE IS PRESENT IN LAYER - - ADJUST COEFFICIENTS FOR LATENT HEAT
!****    TRANSFER AND THE SLOPE OF THE WATER CONTENT-TEMPERATURE CURVE
         IF (ICES(J) .EQ. 1) THEN
!           ICE IS PRESENT FOR THE ENTIRE TIME STEP
            T=273.16+TS(J)
          ELSE
!           ICE IS PRESENT ONLY AT THE END OF THE TIME STEP - DETERMINE
!           THE TEMPERATURE AT WHICH THE SOIL WILL BEGIN TO FREEZE
            TLCONC=0.0
            DO 25 K=1,NSALT
               TLCONC=TLCONC + CONC(K,J)
   25       CONTINUE
            TMPFRZ=273.16*LF/G &
        /(LF/G-MAT(J)+TLCONC*UGAS*(TS(J)+273.16)/G)
            T=TMPFRZ
         END IF
         TDT=273.16+TSDT(J)
!        DETERMINE SLOPE OF LIQUID CONTENT-TEMPERATURE CURVE
         CALL FSLOPE (J,DLDTDT,DUMMY,TDT,MATDT,CONCDT,VLCDT,nsalt,col,row)
         CALL FSLOPE (J,DLDT,DUMMY,T,MAT,CONC,VLC,nsalt,col,row)
!        ENTER MATVLC FOR SLOPE OF LIQUID-MATRIC POTENTIAL CURVE
         CALL MATVL3 (J,MATDT(J),VLCDT(J),DLDM,col,row)
         B1(I)=B1(I) - (ZS(J+1)-ZS(J-1))/(2.*DT) &
            *0.5*RHOL*LF*(DLDT+DLDTDT) &
            - WDT*RHOL*LF*(CONH(J-1)+CONH(J))*DLDTDT/DLDM
   30 CONTINUE

         I=N+NS-2
         J=I-N+1
         if(itmpbc.eq.1)then
           tmp1=0
           tmp2=0
         elseif(itmpbc.eq.0)then
           tmp1=TSDT(J)-TSDT(J+1)
           tmp2=TS(J)-TS(J+1)
         end if
         A1(I+1)=WDT*(CON(J) + DELTA(J)*(RHOL*CL*QSL(J)+CV*QSV(J)))
         B1(I)= -WDT*(CON(J-1)+CON(J) &
            +(DELTA(J-1)*(RHOL*CL*QSL(J-1)+CV*(QSV(J-1))) &
            + EPSLON(J)*(RHOL*CL*QSL(J)+CV*QSV(J)))) &
            - (ZS(J+1)-ZS(J-1))*CSDT(J)/(2.*DT)
         C1(I)=WDT*(CON(J) + EPSLON(J)*(RHOL*CL*QSL(J)+CV*QSV(J)))
         D1(I)=(CON(J-1) + DELTA(J-1)*(RHOL*CL*QSL(J-1) + CV*QSV(J-1))) &
            *(WT*(TS(J-1)-TS(J)) + WDT*(TSDT(J-1)-TSDT(J))) &
            -(CON(J) + EPSLON(J)*(RHOL*CL*QSL(J) + CV*QSV(J))) &
            *(WT*tmp2 + WDT*tmp1) &
            -LV*(QSV(J)-QSV(J-1)) + SS(J) &
            -(ZS(J+1)-ZS(J-1))/(2.*DT) &
            *(CS(J)*(TSDT(J)-TS(J)) - RHOI*LF*(VICDT(J)-VIC(J)))
!
!****    CHECK IF ICE IS PRESENT AT THE END OF THE TIME STEP
         IF (ICESDT(J) .NE. 1) then
           N=N+NS-2
           return
         end if
!
!****    ICE IS PRESENT IN LAYER - - ADJUST COEFFICIENTS FOR LATENT HEAT
!****    TRANSFER AND THE SLOPE OF THE WATER CONTENT-TEMPERATURE CURVE
         IF (ICES(J) .EQ. 1) THEN
!           ICE IS PRESENT FOR THE ENTIRE TIME STEP
            T=273.16+TS(J)
          ELSE
!           ICE IS PRESENT ONLY AT THE END OF THE TIME STEP - DETERMINE
!           THE TEMPERATURE AT WHICH THE SOIL WILL BEGIN TO FREEZE
            TLCONC=0.0
            DO K=1,NSALT
               TLCONC=TLCONC + CONC(K,J)
            enddo
            TMPFRZ=273.16*LF/G &
        /(LF/G-MAT(J)+TLCONC*UGAS*(TS(J)+273.16)/G)
            T=TMPFRZ
         END IF
         TDT=273.16+TSDT(J)
!        DETERMINE SLOPE OF LIQUID CONTENT-TEMPERATURE CURVE
         CALL FSLOPE (J,DLDTDT,DUMMY,TDT,MATDT,CONCDT,VLCDT,nsalt,col,row)
         CALL FSLOPE (J,DLDT,DUMMY,T,MAT,CONC,VLC,nsalt,col,row)
!        ENTER MATVLC FOR SLOPE OF LIQUID-MATRIC POTENTIAL CURVE
         CALL MATVL3 (J,MATDT(J),VLCDT(J),DLDM,col,row)
         B1(I)=B1(I) - (ZS(J+1)-ZS(J-1))/(2.*DT) &
            *0.5*RHOL*LF*(DLDT+DLDTDT) &
            - WDT*RHOL*LF*(CONH(J-1)+CONH(J))*DLDTDT/DLDM
      N=N+NS-2
      RETURN
 END SUBROUTINE EBSOIL



! no problem checked changed
!***********************************************************************
 SUBROUTINE FSLOPE(I,DLDT,PTDLDT,TMP,MAT,CONC,VLC,nsalt,col,row)
!     THIS SUBROUTINE DETERMINES THE SLOPE OF THE LIQUID CONTENT -
!     TEMPERATURE CURVE (DLDT) AS WELL AS THE PARTIAL OF T*DLDT.
!     (OSMOTIC EFFECTS ARE INCLUDED.)
!     ---- TMP MUST BE IN DEGREES KELVIN
!***********************************************************************
    use constvar_mod
    use soilproperty_mod,only:RHOB2d,SALTKQ2d
    implicit none

!input
    integer(i4),intent(in)::I,col,row,nsalt
    real(r8),intent(inout)::DLDT
    real(r8),intent(in)::PTDLDT,TMP
    real(r8),dimension(NSMAX),intent(in)::MAT,VLC
    real(r8),dimension(NSALTMAX,NSMAX),intent(in)::CONC

!temp
    real(r8)::CNCGAS,CRT,CRKQ,DKQ,DLDM
    integer(i4)::J
    real(r8),dimension(:),pointer::RHOB
    real(r8),dimension(:,:),pointer::SALTKQ

    RHOB=>RHOB2d(col,row,:)
    SALTKQ=>SALTKQ2d(col,row,:,:)


      CNCGAS=0.0
      CRT=0.0
      DO 10 J=1,NSALT
         DKQ=SALTKQ(J,I) + VLC(I)*RHOL/RHOB(I)
         CNCGAS=CNCGAS + CONC(J,I)*UGAS
         CRT=CRT + CONC(J,I)*UGAS*TMP*RHOL/RHOB(I)/DKQ
   10 CONTINUE
      CALL MATVL3 (I,MAT(I),VLC(I),DLDM,col,row)
      IF (DLDM .GT. 0.0) THEN
         DLDT=(LF/TMP + CNCGAS)/(G/DLDM + CRT)
        ELSE
!        NO CHANGE IN WATER CONTENT WITH MATRIC POTENTIAL OR TEMPERATURE
!        LAYER IS PROBABLY SATURATED
         DLDT=0.0
      END IF
!
!**** THIS PART OF THE SUBROUTINE CALCULATES THE PARTIAL OF T*DLDT
!
      CRKQ=0.0
      DO 20 J=1,NSALT
         DKQ=SALTKQ(J,I) + VLC(I)*RHOL/RHOB(I)
         CRKQ=CRKQ + &
            CONC(J,I)*UGAS/DKQ*(1.-2.*TMP*(RHOL/RHOB(I))*DLDT/DKQ)
   20 CONTINUE
!XXXX PTDLDT=-LF/ ((-B(I)*G*MAT(I)/VLC(I) + CRT)**2)
!XXXX>      * (B(I)*G*MAT(I)*(B(I)+1)*DLDT/(VLC(I)**2) + CRKQ)
      RETURN
 END SUBROUTINE FSLOPE



! no problem checked changed
!***********************************************************************
!
 SUBROUTINE ENHANC (I,EN,VLC,col,row)
!
!     THIS SUBROUTINE CALCULATES THE ENHANCEMENT COEFFICIENT FOR VAPOR
!    TRANSPORT DUE TO THERMAL GRADIENTS FOR NODE I
!
!***********************************************************************
    use soilproperty_mod,only:clay2d,SAT2d
    implicit none
!input
    integer(i4),intent(in)::I,COL,ROW
    real(r8),intent(inout)::EN
    real(r8),dimension(NSMAX),intent(in)::VLC
!temp
    real(r8)::EN1,EN2,EN4,EN5,EN3,EXPON

    real(r8),dimension(:),pointer::CLAY,SAT

    EN1=9.5
    EN2=3.0
    EN4=1.0
    EN5=4.0

    CLAY=>CLAY2d(col,row,:)
    SAT=>SAT2d(col,row,:)

!
      IF (CLAY(I) .LE. 0.02) THEN
!        EN3 BECOMES INFINITELY LARGE WHEN CLAY IS ZERO
!        (SMALLEST CLAY CONTENT IN DATA FROM CASS ET. AL. WAS 2%)
         EN3=SAT(I)*(1+2.6/SQRT(0.02))
        ELSE
         if(clay(I).lt.0.0) then
           print*,"Error in Enhanc, the arg for sqrt is less than 0"
           print*,col,row
           stop
         end if
         EN3=SAT(I)*(1+2.6/SQRT(CLAY(I)))
      END IF
      EXPON=-(EN3*VLC(I)/SAT(I))**EN5
!
!     EXP(-50.) APPROXIMATELY = 0.0, BUT UNDERFLOW MAY RESULT IF YOU TRY
!     TO USE LARGE NEGATIVE NUMBERS  --->  CHECK IF  < -50.
      IF (EXPON .LE. -50.) THEN
         EXPON=0.0
       ELSE
         EXPON=EXP(EXPON)
      END IF
!
!**** CALCULATE ENHANCEMENT FACTOR
      EN = EN1 + EN2*VLC(I)/SAT(I) - (EN1-EN4)*EXPON
      CONTINUE
      RETURN
    END SUBROUTINE ENHANC


!no problem changed
!***********************************************************************
!
 SUBROUTINE QVSOIL (NS,QSV,ZS,TS,MAT,VLC,VIC,CONC,DCONVP,PRESUR,nsalt,col,row)
!
!     THIS SUBROUTINE CALCULATES THE VAPOR FLUX BETWEEN SOIL NODES DUE
!     TO GRADIENTS IN POTENTIAL AND TEMPERATURE GRADIENTS. (FLUX DUE TO
!     TEMPERATURE GRADIENTS ARE MULTIPLIED BY AN ENHANCEMENT FACTOR.
!
!***********************************************************************
    use constvar_mod,only:G,P0,RHOM,RHOOM,UGAS,VDIFF
    use soilproperty_mod,only:om2d,rhob2d,VAPEXP2d,VAPCOF2d,ENTRY2d
    implicit none
!input
    integer(i4),intent(in)::NS,col,row,nsalt
    real(r8),dimension(NSMAX),intent(in)::ZS,TS,MAT,VLC,VIC
    real(r8),dimension(NSMAX),intent(inout)::QSV,DCONVP
    real(r8),dimension(NSALTMAX,NSMAX),intent(in)::CONC
    real(r8),intent(in)::PRESUR
!temp
    real(r8),dimension(NSMAX)::HUMID,DVT,DVP,CONVT,CONVP
    integer(i4)::I,J
    real(r8)::VAC,DV,EN,S,SATVAP,TLCONC,TOTPOT

    real(r8),dimension(:),pointer::OM,rhob,VAPEXP,VAPCOF,ENTRY

    OM=>om2d(col,row,:)
    RHOB=>RHOB2d(col,row,:)
    VAPEXP=>VAPEXP2d(col,row,:)
    VAPCOF=>VAPCOF2d(col,row,:)
    ENTRY=>ENTRY2d(col,row,:)

      DO 20 I=1,NS
         VAC=1. - RHOB(I)*((1.-OM(I))/RHOM+OM(I)/RHOOM) - VLC(I)- VIC(I)
         IF (MAT(I).LT.ENTRY(I) .AND. VAC.GT.0.0) THEN
            DV=VDIFF*(((TS(I)+273.16)/273.16)**2)*(P0/PRESUR)
            DV=DV*VAPCOF(I)*(VAC**VAPEXP(I))
            CALL ENHANC (I,EN,VLC,col,row)
            if(abs(TS(I)) .gt. 120.0)then
              print*,"In qvsoil,TS(I)",col,row,TS(I)
              stop
            end if
            CALL VSLOPE (S,SATVAP,TS(I))
!
!****       DETERMINE THE HUMIDITY FROM THE TOTAL WATER POTENTIAL
            TLCONC=0.0
            DO 10 J=1,NSALT
               TLCONC= TLCONC+CONC(J,I)
   10       CONTINUE
            TOTPOT=MAT(I) - TLCONC*UGAS*(TS(I)+273.16)/G
            HUMID(I)=EXP(.018*G/UGAS/(TS(I)+273.16)*TOTPOT)
            DVT(I)=DV*EN*HUMID(I)*S
            DVP(I)=DV*SATVAP
          ELSE
!           NO AIR POROSITY OR SATURATED --> NO VAPOR DIFFUSION
            DVT(I)=0.0
            DVP(I)=0.0
            HUMID(I)=1.0
         END IF
!
   20 CONTINUE
!
      CALL CONDUC (NS,ZS,DVT,CONVT,"QVSOIL DVT")
      CALL CONDUC (NS,ZS,DVP,CONVP,"QVSOIL DVP")
!
      DO 30 I=1,NS-1
         QSV(I)=CONVT(I)*(TS(I)-TS(I+1))+ CONVP(I)*(HUMID(I)-HUMID(I+1))
         DCONVP(I)=.018*G/UGAS/(TS(I)+273.16)*CONVP(I)
   30 CONTINUE
      RETURN
 END SUBROUTINE QVSOIL



! no problem changed
!***********************************************************************
!
 SUBROUTINE SOILHK (NS,HK,MAT,VLC,VIC,col,row)
!
!    THIS SUBROUTINE CALCULATES THE HYDRAULIC CONDUCTIVITY OF EACH SOIL
!     NODE
!
!***********************************************************************
    use dims_mod, only:NSMAX
    use constvar_mod, only:RHOM,RHOOM
    use soilproperty_mod, only:RHOB2d,SATK2d,OM2d,ENTRY2d,B2d
    implicit none
! Input
    integer(i4),intent(in)::NS,col,row
    real(r8),dimension(NSMAX),intent(inout)::HK
    real(r8),dimension(NSMAX),intent(in)::MAT,VLC,VIC
! temp variable
    real(r8),dimension(:),pointer::RHOB,SATK,OM,ENTRY,B

    integer(i4)::I,J
    real(r8)::FRACT,AlfaHN,POROS,limitice

    RHOB=>RHOB2d(col,row,:)
    OM=>OM2d(col,row,:)
    B=>B2d(col,row,:)
    ENTRY=>ENTRY2d(col,row,:)
    SATK=>SATK2d(col,row,:)

    limitice=0.13
    DO 10 I=1,NS
         IF (ENTRY(I) .GT. MAT(I)) THEN
            HK(I)=SATK(I)*(ENTRY(I)/MAT(I))**(2.+3/B(I))
          ELSE
            HK(I)=SATK(I)
         END IF
   10 CONTINUE
      RETURN
 END SUBROUTINE SOILHK





! changed
!***********************************************************************
!
      SUBROUTINE QLSOIL (NS,ZS,QSL,CON,MAT)
!
!     THIS SUBROUTINE CALCULATES THE LIQUID MOISTURE FLUX BETWEEN
!     SOIL NODES
!
!***********************************************************************
  implicit none
  !input
  real(r8),dimension(NSMAX),intent(in)::ZS,CON,MAT
  real(r8),dimension(NSMAX),intent(inout)::QSL
  integer(i4),intent(in)::NS

  !temp
  integer(i4)::I
!
      DO 10 I=1,NS-1
         QSL(I)=CON(I)*(MAT(I)-MAT(I+1)+ZS(I+1)-ZS(I))
   10 CONTINUE
      RETURN
      END SUBROUTINE QLSOIL


!no problem...checked changed
!***********************************************************************
 SUBROUTINE EBRES (N,NR,NSP,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT, &
    GMCMAX,RHOR,RESCOF,SR,UR,QVR,RHOSP,ITER,PRESUR,col,row)
!
!     THIS SUBOUTINE CALCULATES THE NEWTON-RAPHSON COEFFICIENTS FOR THE
!     ENERGY BALANCE OF THE RESIDUE LAYERS.
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d,DT2d
    use constvar_mod,only:LV
    use matrix_mod,only:A12d,B12d,C12d,D12d
    use residu_mod,only:EVAP2d,EVAPK2d
    use savedata_mod,only:CREST2d
    implicit none
!input
    integer(i4),intent(in)::NR,NSP,ITER,col,row
    integer(i4),intent(inout)::N
    real(r8),dimension(NRMAX),intent(in)::ZR,TR,TRDT,VAPR,VAPRDT,RHOR,SR,UR
    real(r8),dimension(NRMAX),intent(inout)::GMC,GMCDT,QVR
    real(r8),dimension(NSPMAX),intent(in)::RHOSP
    real(r8),intent(in)::GMCMAX,RESCOF,PRESUR
!temp
    real(r8),dimension(NRMAX)::CRES,CRESDT,TKRES,CONVEC,VAPCON,CON,CONV
    integer(i4)::NR1,I,J

    real(r8),dimension(:),pointer::A1,B1,C1,D1
    real(r8),dimension(:),pointer::EVAP,EVAPK
    real(r8),dimension(:),pointer::CREST
    real(r8),pointer::DT,WT,WDT

    A1=>A12d(col,row,:)
    B1=>B12d(col,row,:)
    C1=>C12d(col,row,:)
    D1=>D12d(col,row,:)
    EVAP=>EVAP2d(col,row,:)
    EVAPK=>EVAPK2d(col,row,:)
    CREST=>CREST2d(col,row,:)
    DT=>DT2d(col,row)
    WT=>WT2d(col,row)
    WDT=>WDT2d(col,row)

!**** DETERMINE THE EVAPORATION FROM RESIDUE ELEMENTS AT EACH NODE
      CALL RESVAP (NR,EVAP,EVAPK,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT, &
            GMCMAX,RHOR,RESCOF,UR,col,row)
!
!**** DETERMINE THE HEAT TRANSFER COEFFICIENT FOR EACH NODE
      CALL RESTK (NR,NSP,TKRES,CONVEC,TR,TRDT,GMC,GMCDT,RHOR,RHOSP,col,row)
      NR1=NR+1
      TKRES(NR1)=TKRES(NR)
!
!**** DETERMINE THE CONDUCTANCE TERM FOR HEAT TRANSPORT BETWEEN NODES
      CALL CONDUC (NR1,ZR,TKRES,CON,"EBRES TKRES")
!
!**** DETERMINE THE VAPOR TRANSPORT FROM THE THERMAL CONVECTION
      CALL RESVK (NR,NSP,TR,TRDT,CONVEC,VAPCON,presur,col,row)
      VAPCON(NR1)=VAPCON(NR)
!
!**** DETERMINE THE CONDUCTANCE TERM FOR CONVECTIVE VAPOR TRANSPORT
      CALL CONDUC (NR1,ZR,VAPCON,CONV,"EBRES VAPCON")
!
!**** DETERMINE THE VAPOR FLUX BETWEEN RESIDUE NODES
      CALL QVRES (NR,QVR,CONV,VAPR,VAPRDT,col,row)
!
!**** DETERMINE THE VOLUMETRIC HEAT CAPACITY OF EACH NODE
      IF (ITER .EQ. 1) CALL RESHT (NR,CREST,GMC,RHOR)
      CALL RESHT (NR,CRESDT,GMCDT,RHOR)
!
!**** CALCULATE THE AVERAGE HEAT CAPACITY OVER THE TIME STEP
      CALL WEIGHT (NR,CRES,CREST,CRESDT,col,row)
!
!
!**** DETERMINE THE MATRIX COEFFICIENTS FOR THE TOP LAYER
!
      A1(N+1)=A1(N+1) + WDT*CON(1)
      C1(N)=C1(N) + WDT*CON(1)
      IF (NR.GT.1) THEN
         B1(N)=B1(N) - WDT*CON(1) - (ZR(2)-ZR(1))/(2*DT)*CRES(1)
         D1(N)=D1(N) - CON(1)*(WT*(TR(1)-TR(2))+WDT*(TRDT(1)-TRDT(2))) &
            + LV*EVAP(1) + SR(1) &
            - (ZR(2)-ZR(1))/(2*DT)*CRES(1)*(TRDT(1)-TR(1))
       ELSE
         B1(N)=B1(N) - WDT*CON(1) - (ZR(2)-ZR(1))/DT*CRES(1)
         D1(N)=D1(N) - CON(1)*(WT*(TR(1)-TR(2))+WDT*(TRDT(1)-TRDT(2))) &
            + LV*EVAP(1) + SR(1) &
            - (ZR(2)-ZR(1))/DT*CRES(1)*(TRDT(1)-TR(1))
      END IF
!
!**** DETERMINE THE COEFFICIENTS FOR THE REST OF THE LAYERS
      DO 10 I=N+1,N+NR-1
         J=I-N+1
         A1(I+1)=A1(I+1) + WDT*CON(J)
         C1(I)=C1(I) + WDT*CON(J)
         IF (J .NE. NR) THEN
            B1(I)=B1(I) - WDT*(CON(J-1)+CON(J)) &
             - (ZR(J+1)-ZR(J-1))/(2*DT)*CRES(J)
            D1(I)=CON(J-1)*(WDT*(TRDT(J-1)-TRDT(J))+WT*(TR(J-1)-TR(J))) &
        -CON(J)*(WDT*(TRDT(J)-TRDT(J+1)) + WT*(TR(J)-TR(J+1))) &
        +LV*EVAP(J) + SR(J) &
        -(ZR(J+1)-ZR(J-1))/(2*DT)*CRES(J)*(TRDT(J)-TR(J))
        ELSE
            B1(I)=B1(I) - WDT*(CON(J-1)+CON(J)) &
            - (ZR(J+1)-ZR(J)+(ZR(J)-ZR(J-1))/2)/DT*CRES(J)
            D1(I)=CON(J-1)*(WDT*(TRDT(J-1)-TRDT(J))+WT*(TR(J-1)-TR(J))) &
            -CON(J)*(WDT*(TRDT(J)-TRDT(J+1)) + WT*(TR(J)-TR(J+1))) &
            +LV*EVAP(J) + SR(J)-(ZR(J+1)-ZR(J)+(ZR(J)-ZR(J-1))/2) &
            /DT*CRES(J)*(TRDT(J)-TR(J))
         END IF
   10 CONTINUE
!
!**** DETERMINE THE COEFFICIENTS FOR THE SOIL SURFACE
!
      N=N+NR
      B1(N)=B1(N) - WDT*CON(NR)
      D1(N)=CON(NR)*(WDT*(TRDT(NR)-TRDT(NR+1))+WT*(TR(NR)-TR(NR+1))) &
        + LV*QVR(NR)
      RETURN
 END SUBROUTINE EBRES



! no problem checked changed
!***********************************************************************
!
      SUBROUTINE RESHT (NR,CRES,GMC,RHOR)
!
!     THIS SUBOUTINE CALCULATES THE VOLUMETRIC SPECIFIC HEAT FOR EACH
!     RESIDUE NODE.
!
!***********************************************************************
    use constvar_mod, only:CR,CL
    implicit none

! input
    integer(i4),intent(in)::NR
    real(r8),dimension(NRMAX),intent(inout)::CRES
    real(r8),dimension(NRMAX),intent(in)::GMC,RHOR

!temp variable
    integer(i4)::I

    DO 10 I=1,NR
         CRES(I)=RHOR(I)*(CR + GMC(I)*CL)
   10 CONTINUE
      RETURN
 END SUBROUTINE RESHT


! no problem checked changed
!***********************************************************************
!
      SUBROUTINE RESVK (NR,NSP,TR,TRDT,CONVEC,VAPCON,PRESUR,col,row)
!
!     THIS SUBROUTINE CALCULATES THE VAPOR CONDUCTANCE TERM (CONV)
!     BETWEEN RESIDUE NODES USING THE CALCULATED THERMAL CONDUCTANCE.
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d
    use constvar_mod,only:RHOA,CA,P0
    use spwatr_mod,only:VDIFSP,VAPSPX
    implicit none
!input
    integer(i4),intent(in)::NR,NSP
    real(r8),intent(in)::PRESUR
    real(r8),dimension(NRMAX),intent(in)::TR,TRDT,CONVEC
    real(r8),dimension(NRMAX),intent(inout)::VAPCON
    integer(i4),intent(in)::col,row
!temp
    integer(i4)::I
    real(r8)::AVGTMP

    real(r8),pointer::WT,WDT

    WT=>WT2d(col,row)
    WDT=>WDT2d(col,row)

!
      DO 10 I=1,NR
         VAPCON(I)=CONVEC(I)/RHOA/CA
!
         IF (NSP .GT. 0) THEN
!           SNOW FILTERS DOWN THRU RESIDUE - VAPOR DIFFUSIVITY IS THAT
!           THROUGH SNOW
            AVGTMP = WT*TR(I) + WDT*TRDT(I)
            VAPCON(I) = VDIFSP*(P0/PRESUR)*(1.+AVGTMP/273.16)**VAPSPX
         END IF
   10 CONTINUE
      RETURN
END SUBROUTINE RESVK

! chagned
!***********************************************************************
!
      SUBROUTINE QVRES (NR,QVR,CONV,VAPR,VAPRDT,col,row)
!
!     THIS SUBOUTINE CALCULATES VAPOR FLUX BETWEEN RESIDUE NODES.
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d
    implicit none

!input
    integer(i4)::NR
    integer(i4),intent(in)::col,row
    real(r8),dimension(NRMAX),intent(in)::CONV,VAPR,VAPRDT
    real(r8),dimension(NRMAX),intent(inout)::QVR
!temp
    integer(i4)::I

    real(r8),pointer::WT,WDT
    WT=>WT2d(col,row)
    WDT=>WDT2d(col,row)

      DO 10 I=1,NR
         QVR(I) = CONV(I)*(WT*(VAPR(I)-VAPR(I+1)) &
            + WDT*(VAPRDT(I)-VAPRDT(I+1)))
   10 CONTINUE
      RETURN
 END SUBROUTINE QVRES

! changed
!***********************************************************************
!
 SUBROUTINE RESTK(NR,NSP,TKRES,CONVEC,TR,TRDT,GMC,GMCDT,RHOR,RHOSP,col,row)
!
!     THIS SUBROUTINE CALCULATES THE THERMAL CONDUCTANCE TERM BETWEEN
!     RESIDUE NODES USING THE CALCULATED WINDSPEED AT A NODE FOR THE
!     THE CONVECTIVE TRANSPORT TERM, AND THE THERMAL CONDUCTIVITIES OF
!     OF THE RESIDUE AND WATER THE THERMAL CONDUCTIVITY TERM.
!     CONVECTIVE AND CONDUCTIVE TERMS ARE THEN WEIGHTED ACCORDING TO
!     THE VOLUMETRIC FRACTIONS OF AIR, WATER AND RESIDUE.
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d!DT
    use constvar_mod,only:TKA,RHOL,TKL,TKR
    use rsparm_mod,only:RESTKA,RESTKB
    use windv_mod,only:WINDR2d
    implicit none


!input
    integer(i4),intent(in)::NR,NSP,col,row
    real(r8),dimension(NRMAX),intent(in)::TR,TRDT,GMC,GMCDT,RHOR
    real(r8),dimension(NSPMAX),intent(in)::RHOSP
    real(r8),dimension(NRMAX),intent(inout)::TKRES,CONVEC
!temp
    real(r8),dimension(NSPMAX)::TKSP
    real(r8)::RESDEN,TK,AIR,AVGGMC,AVGTMP,RESVOL,WATER
    integer(i4)::I

    real(r8),dimension(:),pointer::WINDR
    real(r8),pointer::WT,WDT


    RESDEN=170.0
    WINDR=>WINDR2d(col,row,:)
    WT=>WT2d(col,row)
    WDT=>WDT2d(col,row)

!
      IF (NSP .GT. 0) THEN
!        ASSUME THE SNOW FILTERS DOWN INTO THE RESIDUE - THEREFORE
!        REDEFINE THE THERMAL CONVECTION AS A CONDUCTION THRU SNOW
         CALL SNOWTK (NSP,TKSP,RHOSP)
      END IF
!
      DO 10 I=1,NR
!
!        AVERAGE THE TEMP AND WATER CONTENT OVER TIME
         AVGTMP = WDT*TRDT(I) + WT*TR(I)
         AVGGMC = WDT*GMCDT(I) + WT*GMC(I)
!
!        CALCULATE THE VOLUME FRACTION OF EACH MATERIAL
         RESVOL = RHOR(I)/RESDEN
         WATER = AVGGMC*RHOR(I)/RHOL
         AIR = 1. - RESVOL - WATER
!
!        IF SNOW IS PRESENT - DO NOT CALCULATE CONVECTION
         IF (NSP .GT. 0) THEN
            CONVEC(I) = TKSP(NSP)
           ELSE
!
!           CALCULATE THE CONVECTIVE HEAT TRANSFER COEFFICIENT
            CONVEC(I) = TKA*(1.+RESTKA*AVGTMP)*(1.+RESTKB*WINDR(I))
         END IF
!
!        ADJUST CONVECTIVE TRANSPORT FOR AIR POROSITY
!        (IN THE CASE OF SNOW, AIR POROSITY IS THE FRACTION OF SNOW)
         CONVEC(I) = CONVEC(I)*AIR
!
!        CALCULATE THERMAL CONDUCTIVITY THROUGH THE RESIDUE AND WATER
         TK =  WATER*TKL + RESVOL*TKR
!
!        CALCULATE THE EFFECTIVE HEAT TRANSFER COEFFICIENT
         TKRES(I) = CONVEC(I) + TK
   10 CONTINUE
      RETURN
END SUBROUTINE RESTK

! changed
!***********************************************************************
!
 SUBROUTINE RESVAP (NR,EVAP,EVAPK,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,&
    GMCMAX,RHOR,RESCOF,UR,col,row)
!     THIS SUBROUTINE IS TO DETERMINE THE EVAPORATION WITHIN THE
!     RESIDUE LAYERS ASSUMING THE AVERAGE MOISTURE CONTENT OF THE
!     RESIDUE IS THE MOISTURE CONTENT AT THE BEGINNING OF THE TIME STEP
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d,DT2d
    implicit none
!input
    integer(i4),intent(in)::NR,col,row
    real(r8),dimension(NRMAX),intent(in)::ZR,TR,TRDT,VAPR,VAPRDT,RHOR,UR
    real(r8),dimension(NRMAX),intent(inout)::GMC,GMCDT,EVAP,EVAPK
    real(r8),intent(in)::GMCMAX,RESCOF
!temp
    integer(i4)::I
    real(r8)::DHDW,DUMMY,DZ,HUMT,HUMDT,SATV,SATVDT
    real(r8),pointer::DT,WT,WDT
    DT =>DT2d(col,row)
    WT =>WT2d(col,row)
    WDT=>WDT2d(col,row)

      DO 20 I=1,NR
!
!****    DETERMINE THE RELATIVE HUMIDITY IN THE RESIDUE
         CALL RESHUM (1,HUMT,DUMMY,GMC(I),TR(I))
         CALL RESHUM (1,HUMDT,DHDW,GMCDT(I),TRDT(I))
!
!****    DETERMINE THE SATURATED VAPOR DENSITY OF THE RESIDUE
         if(TRDT(I).lt. -273.15 .or. TR(I).lt. -273.15)then
            Print*,"Error in RESVAP. The TRDT is less than -273.15"
            print*,TR(I),TRDT(I)
            stop
         end if
            if(abs(TR(I)) .gt. 120.0)then
              print*,"RESVAP TR",col,row,TR(I)
              stop
            end if
         CALL VSLOPE (DUMMY,SATV,TR(I))
            if(abs(TRDT(I)) .gt. 120.0)then
              print*,"RESVAP TRDT",col,row,TRDT(I)
              stop
            end if
         CALL VSLOPE (DUMMY,SATVDT,TRDT(I))
!
         IF (I .EQ. 1) THEN
            DZ=ZR(2)-ZR(1)
            IF (NR .NE. 1) DZ=DZ/2.
           ELSE
            IF (I .EQ. NR) THEN
               DZ=ZR(I+1)-ZR(I)+(ZR(I)-ZR(I-1))/2.
              ELSE
               DZ=(ZR(I+1)-ZR(I-1))/2.
            END IF
         END IF
!
!****    DEFINE EVAPK(I) -- FOR NOW, IT WILL SIMPLY BE SET EQUAL TO
!        THE INVERSE OF RESCOF, ASSUMING THAT RESCOF IS THE VAPOR
!        RESISTANCE.  LATER WORK MAY REQUIRE THAT THE VAPOR TRANSPORT
!        RESISTANCE BE CALCULATED OR OBTAINED FROM A SUBROUTINE.
         EVAPK(I)=1./RESCOF
!
         EVAP(I)=EVAPK(I)*(WDT*(VAPRDT(I)-HUMDT*SATV)&
            +WT*(VAPR(I)-HUMT*SATV))
         GMCDT(I)= GMC(I) + EVAP(I)*DT/DZ/RHOR(I) + UR(I)
         IF (GMCDT(I) .LT. 0.0) THEN
            GMCDT(I)=0.0
            EVAP(I)= DZ*RHOR(I)*(GMCDT(I)-GMC(I)-UR(I))/DT
         END IF
         IF (GMCDT(I) .GT. GMCMAX) THEN
            GMCDT(I)=GMCMAX
            EVAP(I)= DZ*RHOR(I)*(GMCDT(I)-GMC(I)-UR(I))/DT
         END IF
!
   20 CONTINUE
      RETURN
END SUBROUTINE RESVAP




!changed
!***********************************************************************
    SUBROUTINE EBSNOW (N,NSP,NR,ICESPT,TSP,TSPDT,DLW,DLWDT,RHOSP,ZSP, &
        DZSP,QVSP,VAPSP,VAPSPT,SSP,ITER,PRESUR,col,row)
!     THIS SUBROUTINE CALCULATES THE JACOBIAN MATRIX COEFFICIENTS FOR
!     THE SNOW PORTION OF THE NEWTON-RAPHSON SOLUTION OF THE ENERGY
!     BALANCE
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d,DT2d
    use constvar_mod,only:LF,LS,LV,RHOL
    use matrix_mod,only:A12d,B12d,C12d,D12d
    use savedata_mod,only:QVSPT2d,CON2d,CSPT2d
    implicit none
!input
    integer(i4),intent(in)::NSP,NR,ITER,col,row
    integer(i4),intent(inout)::N
    real(r8),dimension(NSPMAX),intent(in)::TSP,TSPDT,DLW,DLWDT,RHOSP,ZSP,DZSP,SSP
    real(r8),dimension(NSPMAX),intent(inout)::QVSP
    integer(i4),dimension(NSPMAX),intent(in)::ICESPT
    real(r8),intent(in)::VAPSP,VAPSPT,PRESUR
!temp
    real(r8),dimension(NSPMAX)::QVSPDT,TK,CONV,SLOPE,CSP,CSPDT
    integer(i4)::NSP1,I,J

    real(r8),dimension(:),pointer::A1,B1,C1,D1
    real(r8),dimension(:),pointer::QVSPT,CON,CSPT
    real(r8),pointer::DT,WT,WDT

    A1   =>A12d(col,row,:)
    B1   =>B12d(col,row,:)
    C1   =>C12d(col,row,:)
    D1   =>D12d(col,row,:)
    QVSPT=>QVSPT2d(col,row,:)
    CON  =>CON2d(col,row,:)
    CSPT =>CSPT2d(col,row,:)
    DT   =>DT2d(col,row)
    WT   =>WT2d(col,row)
    WDT  =>WDT2d(col,row)


!**** DETERMINE THE VAPOR FLUX BETWEEN NODES
      IF (ITER .EQ. 1) CALL QVSNOW (NSP,QVSPT,CONV,SLOPE,TSP,ZSP,VAPSP,PRESUR,col,row)
      CALL QVSNOW (NSP,QVSPDT,CONV,SLOPE,TSPDT,ZSP,VAPSPT,PRESUR,col,row)
!     IF UNDERLYING MATERIAL IS RESIDUE, VAPOR DENSITY IS NOT A FUNCTION
!     OF TEMPERATURE, I.E.  SLOPE(NSP+1) = 0.0
      IF (NR .GT. 0) SLOPE(NSP+1) = 0.0
!
!**** OBTAIN THE AVERAGE VAPOR FLUX OVER THE TIME STEP
      CALL WEIGHT (NSP,QVSP,QVSPT,QVSPDT,col,row)
!
      IF (ITER .EQ. 1) THEN
!        (DENSITY AND THERMAL CONDUCTIVITY ARE CONSTANT OVER TIME STEP)
!****    DETERMINE THE THERMAL CONDUCTIVITY OF EACH NODE
         CALL SNOWTK (NSP,TK,RHOSP)
         TK(NSP+1) = TK(NSP)
!
!****    DETERMINE THE CONDUCTANCE TERM BETWEEN NODES
         NSP1 = NSP+1
         CALL CONDUC (NSP1,ZSP,TK,CON,"EBSNOW TK")
      END IF
!
!**** CALCULATE THE SPECIFIC HEAT OF EACH NODE
      IF (ITER .EQ. 1) CALL SNOWHT (NSP,CSPT,TSP,RHOSP,col,row)
      CALL SNOWHT (NSP,CSPDT,TSPDT,RHOSP,col,row)
!
!**** OBTAIN THE AVERAGE SPECIFIC HEAT OVER THE TIME STEP
      CALL WEIGHT (NSP,CSP,CSPT,CSPDT,col,row)
!
!
!**** DETERMINE THE COEFFICIENTS FOR THE SURFACE LAYER
      D1(N)= D1(N) - CON(1)*(WT*(TSP(1)-TSP(2)) &
        + WDT*(TSPDT(1)-TSPDT(2))) + SSP(1) &
        - CSP(1)*(TSPDT(1)-TSP(1))*DZSP(1)/DT &
        - RHOL*LF*(DLWDT(1)-DLW(1))/DT - LS*QVSP(1)
!
      IF (ICESPT(1) .EQ. 0) THEN
!        LAYER IS NOT MELTING - ENERGY BUDGET BASED ON TEMPERATURE
         A1(N+1) = WDT*(CON(1) + CONV(1)*LS*SLOPE(1))
         B1(N) = B1(N) - WDT*(CON(1) + CONV(1)*LS*SLOPE(1)) &
            - DZSP(1)*CSP(1)/DT
        ELSE
!
!        LAYER IS MELTING - ENERGY BUDGET BASED ON WATER CONTENT
         A1(N+1)= 0.0
         B1(N)= B1(N) - RHOL*LF/DT
!        IF SNOW IS NOT FIRST MATERIAL, SET DERIV. FOR LAST NODE TO 0.0
         IF (N.GT.1) C1(N-1)=0.0
      END IF
!
!**** DETERMINE THE COEFFICIENTS FOR THE REMAINDER OF THE SNOWPACK
      DO 20 I=N+1,N+NSP-1
         J=I-N+1
         D1(I)= CON(J-1)*(WT*(TSP(J-1)-TSP(J)) &
            + WDT*(TSPDT(J-1)-TSPDT(J))) &
            - CON(J)*(WT*(TSP(J)-TSP(J+1)) &
            + WDT*(TSPDT(J)-TSPDT(J+1))) + SSP(J) &
            - CSP(J)*(TSPDT(J)-TSP(J))*DZSP(J)/DT &
            - RHOL*LF*(DLWDT(J)-DLW(J))/DT - LS*(QVSP(J)-QVSP(J-1))
!
         IF (ICESPT(J) .EQ. 0) THEN
!           LAYER IS NOT MELTING - ENERGY BUDGET BASED ON TEMPERATURE
            A1(I+1)= WDT*(CON(J) + CONV(J)*LS*SLOPE(J))
            B1(I)= -WDT*(CON(J-1)+CON(J)) &
                - WDT*LS*SLOPE(J)*(CONV(J-1) + CONV(J)) &
                - DZSP(J)*CSP(J)/DT
            C1(I-1)= WDT*(CON(J-1) + CONV(J-1)*LS*SLOPE(J))
           ELSE
!
!           LAYER IS MELTING - ENERGY BUDGET BASED ON WATER CONTENT
            A1(I+1)= 0.0
            B1(I)= -RHOL*LF/DT
            C1(I-1)= 0.0
         END IF
20    CONTINUE
!
!
!**** DETERMINE THE BOUNDARY CONDITIONS FOR TOP LAYER OF NEXT MATERIAL
      N=N+NSP
!
      IF (NR .GT. 0) THEN
!        SNOW OVERLYING RESIDUE
!        CHECK IF LAST SNOW NODE IS MELTING - IF SO, ENERGY BALANCE
!        IS BASED ON WATER CONTENT, NOT TEMPERATURE AND A1(N)=0.0
         IF (ICESPT(NSP) .EQ. 0) A1(N) =  WDT*CON(NSP)
         B1(N) = -WDT*CON(NSP)
         C1(N-1) = C1(N-1) + WDT*CON(NSP)
         D1(N) = CON(NSP)* (WT*(TSP(NSP) - TSP(NSP+1)) &
            + WDT*(TSPDT(NSP) - TSPDT(NSP+1)))
!
        ELSE
!        SNOW IS LYING ON BARE SOIL - INCLUDE LATENT HEAT TRANSFER
!        AND VAPOR FLUX DEPENDENCE ON TEMPERATURE OF SOIL SURFACE
!        CHECK IF LAST SNOW NODE IS MELTING - IF SO, ENERGY BALANCE
!        IS BASED ON WATER CONTENT, NOT TEMPERATURE AND A1(N)=0.0
         IF (ICESPT(NSP) .GT. 0) &
         A1(N)=WDT*(CON(NSP)+CONV(NSP)*LV*SLOPE(NSP))
         B1(N) =-WDT*(CON(NSP) + CONV(NSP)*LV*SLOPE(NSP+1))
         C1(N-1) = WDT*(CON(NSP) + CONV(NSP)*LS*SLOPE(NSP+1))
         D1(N) = CON(NSP)* (WT*(TSP(NSP) - TSP(NSP+1)) &
            + WDT*(TSPDT(NSP) - TSPDT(NSP+1))) + LV*QVSP(NSP)
      END IF
      RETURN
 END SUBROUTINE EBSNOW


! changed
!***********************************************************************
 SUBROUTINE SNOWHT (NSP,CSP,TSP,RHOSP,col,row)
!     THIS SUBOUTINE CALCULATES THE VOLUMETRIC SPECIFIC HEAT FOR EACH
!     RESIDUE NODE, INCLUDING LATENT HEAT EFFECTS.
!***********************************************************************
    use constvar_mod,only:LS,RHOI
    implicit none
!input
    integer(i4),intent(in)::NSP,col,row
    real(r8),dimension(NSPMAX),intent(inout)::CSP
    real(r8),dimension(NSPMAX),intent(in)::TSP,RHOSP
!temp
    integer(i4)::I
    real(r8)::SPHEAT,DUMMY,S

      DO 10 I=1,NSP
!        CALCULATE THE SPECIFIC HEAT OF THE SNOW
         SPHEAT = 92.96 + 7.37*(TSP(I)+273.16)
!        INCLUDE THE LATENT HEAT TERM
            if(abs(TSP(I)) .gt. 120.0)then
              print*,"In snowht,TSP(I)",col,row,TSP(I)
              stop
            end if
         CALL VSLOPE (S,DUMMY,TSP(I))
!        CALCULATE THE VOLUMETRIC SPECIFIC HEAT INCLUDING LATENT HEAT
         CSP(I) = RHOSP(I)*SPHEAT + (1.-RHOSP(I)/RHOI)*LS*S
   10 CONTINUE
      RETURN
 end subroutine SNOWHT



! changed
!***********************************************************************
!
 SUBROUTINE SNOWTK (NSP,TK,RHOSP)
!
!     THIS SUBROUTINE CALCULATES THE THERMAL CONDUCTIVITY OF THE
!     SNOWPACK LAYERS.  THE EQUATION USED IS OF THE FORM:
!
!           K = A + B*(RHOSP/RHOL)**C
!
!     WHERE:    A = TKSPA      B = TKSPB     C = TKSPEX
!
!***********************************************************************
    use spwatr_mod,only:TKSPA,TKSPB,TKSPEX
    use constvar_mod,only:RHOL
    implicit none
!input
    integer(i4),intent(in)::NSP
    real(r8),dimension(NSPMAX),intent(in)::RHOSP
    real(r8),dimension(NSPMAX),intent(inout)::TK
!temp
    integer(i4)::I


      DO 10 I=1,NSP
         TK(I) = TKSPA + TKSPB*(RHOSP(I)/RHOL)**TKSPEX
   10 CONTINUE
      RETURN
  END SUBROUTINE


!no problem changed
!***********************************************************************
!
      SUBROUTINE CONDUC (N,Z,RK,CON,mainname)
!
!     THIS SUBROUTINE CALCULATES THE CONDUCTANCE TERM BETWEEN TWO NODES
!     USING THE GEOMETRIC MEAN AND THE SPACE INCREMENT
!
!                                                     K(I)
!       K(I) = ( K(I)*K(I+1) )**1/2 =>  CON(I) = -------------
!                                                 Z(I+1) - Z(I)
!
!***********************************************************************
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4
    implicit none
    integer(i4),intent(in)::N
    real(r8),dimension(N),intent(in)::Z
    real(r8),dimension(N),intent(inout)::CON,RK
    integer(i4)::I
    character(*)::mainname
      DO 10 I=1,N-1
            if(RK(I)*RK(I+1).lt.0.0) then
              print*,"The call of CONDUC in "//trim(mainname)//", the arg in sqrt is le 0.0"
              print*,I,N-1,RK(I),RK(I+1)
              if(RK(I).lt.0.0) then
                RK(I)=amax1(RK(I),RK(I+1))
              end if
              if(RK(I+1).lt.0.0) then
                RK(I+1)=amax1(RK(I),RK(I+1))
              end if
            end if
            CON(I)= SQRT(RK(I)*RK(I+1))/(Z(I+1) - Z(I))
   10 CONTINUE
      RETURN
 END SUBROUTINE CONDUC


!!no problem changed
!***********************************************************************
!
      SUBROUTINE QVSNOW (NSP,QVSP,CONV,SLOPE,TSP,ZSP,VAPSP,PRESUR,col,row)
!
!    THIS SUBROUTINE CALCULATES THE VAPOR DIFFUSION IN THE SNOWPACK
!
!***********************************************************************
    use constvar_mod,only:LF,P0,UGAS
    use spwatr_mod
    implicit none
!input
    integer(i4),intent(in)::NSP,col,row
    real(r8),dimension(NSPMAX),intent(in)::TSP,ZSP
    real(r8),dimension(NSPMAX),intent(inout)::SLOPE,CONV,QVSP
    real(r8),intent(in)::VAPSP,PRESUR

    real(r8),dimension(NSPMAX)::VAPSNO,VAPICE
    real(r8)::HUMID,SATV,TMP
    integer(i4)::I,NSP1

!**** DETERMINE THE VAPOR DIFFUSIVITY, DENSITY AND SLOPE OF EACH NODE
      DO 10 I=1,NSP
         VAPSNO(I) = VDIFSP*(P0/PRESUR)*(1.+TSP(I)/273.16)**VAPSPX
            if(abs(TSP(I)) .gt. 120.0)then
              print*,"QVSNOW",col,row,TSP(I)
              stop
            end if
         CALL VSLOPE (SLOPE(I),SATV,TSP(I))
!        CALCULATE THE SATURATED VAPOR DENSITY OVER ICE
         TMP = TSP(I) + 273.16
         HUMID = EXP(0.018/(UGAS*TMP)*LF*TSP(I)/TMP)
         VAPICE(I) = HUMID*SATV
         SLOPE(I) = HUMID*SLOPE(I)
   10 CONTINUE
      VAPICE(NSP+1) = VAPSP
      VAPSNO(NSP+1) = VAPSNO(NSP)
            if(abs(TSP(NSP+1)) .gt. 120.0)then
              print*,"QVSNOW",col,row,TSP(NSP+1)
              stop
            end if
      CALL VSLOPE (SLOPE(NSP+1),SATV,TSP(NSP+1))
      SLOPE(NSP+1) = SLOPE(NSP+1)*VAPSP/SATV
!
!**** DETERMINE THE CONDUCTANCE TERM BETWEEN NODES
      NSP1 = NSP+1
      CALL CONDUC (NSP1,ZSP,VAPSNO,CONV,"qvsnow VAPSNO")
!
!**** CALCULATE THE VAPOR FLUX BETWEEN EACH NODE
      DO 20 I=1,NSP
         QVSP(I) = CONV(I)*(VAPICE(I)-VAPICE(I+1))
   20 CONTINUE
      RETURN
 END SUBROUTINE QVSNOW



!no problem--checked changed
!***********************************************************************
 SUBROUTINE EBCAN (N,NPLANT,NC,NSP,NR,NS,ZC,TC,TCDT,TLC,TLCDT, &
    VAPC,VAPCDT,WCAN,WCANDT,PCAN,PCANDT,MAT,MATDT,SWCAN,LWCAN, &
    CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,ITYPE,ITER,col,row,julian,&
    errorflag)
!     THIS SUBROUTINE CALCULATES THE JACOBIAN MATRIX COEFFICIENTS FOR
!     THE CANOPY PORTION OF THE NEWTON-RAPHSON SOLUTION OF THE ENERGY
!     BALANCE
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d,DT2d
    use constvar_mod,only:CA,LS,LV,RHOA
    !use windv_mod
    use savedata_mod,only:CONT2d
    use writeit_mod,only:hnc2d
    use matrix_mod,only:A12d,B12d,C12d,D12d
    implicit none
!input
    integer(i4),intent(inout)::N
    integer(i4),intent(in)::NPLANT,NC,NSP,NR,NS,ITER,col,row,julian
    real(r8),dimension(NCMAX),intent(in)::ZC,TC,TCDT,VAPC,VAPCDT
    real(r8),dimension(NPMAX,NCMAX-1),intent(in)::TLC
    real(r8),dimension(NPMAX,NCMAX-1),intent(inout)::TLCDT
    real(r8),dimension(NCMAX-1),intent(inout)::WCAN
    real(r8),dimension(NCMAX-1),intent(inout)::WCANDT
    real(r8),dimension(NPMAX),intent(in)::PCAN
    real(r8),dimension(NPMAX),intent(inout)::PCANDT
    real(r8),dimension(NSMAX),intent(in)::MAT,MATDT
    real(r8),dimension(NPMAX+1,NCMAX-1),intent(in)::SWCAN,LWCAN
    real(r8),dimension(NPMAX),intent(in)::DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0
    real(r8),intent(in)::CANMA,CANMB
    integer(i4),dimension(NPMAX),intent(in):: ITYPE
!temp
    real(r8),dimension(NCMAX-1)::CON,CONDT,HEATC,ETLYR,DETLYR,DHEATC,DTLDTC
    real(r8),dimension(NPMAX+1)::TRNSP
    real(r8),dimension(NSMAX)::XTRACT

    real(r8),dimension(:),pointer::A1,B1,C1,D1
    real(r8),dimension(:),pointer::CONT
    real(r8),pointer::hnc,dt,WT,WDT

    real(r8)::contk,CONV,QVCAN,HUMNC,HUMNC1,SATV,SLPNC,SLPNC1
    integer(i4)::I,J
    integer(i4),intent(inout)::errorflag

    A1  =>A12d(col,row,:)
    B1  =>B12d(col,row,:)
    C1  =>C12d(col,row,:)
    D1  =>D12d(col,row,:)
    hnc =>hnc2d(col,row)
    cont=>cont2d(col,row,:)
    dt  =>dt2d(col,row)
    WT  =>WT2d(col,row)
    WDT =>WDT2d(col,row)



!**** CALCULATE LEAF TEMPERATURE AND HEAT TRANSFER FROM CANOPY
      CALL LEAFT (NPLANT,NC,NS,ITER,ITYPE,TC,TCDT,TLC,TLCDT,VAPC,VAPCDT, &
        WCAN,WCANDT,PCAN,PCANDT,MAT,MATDT,TRNSP,XTRACT,SWCAN,LWCAN, &
        HEATC,ETLYR,DETLYR,CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0, &
        DHEATC,DTLDTC,col,row,julian,errorflag)
!
!**** DETERMINE THE EDDY CONDUCTANCE TERM BETWEEN NODES
      IF (ITER .EQ. 1) CALL CANTK (NC,CONT,TC,ZC,col,row)
      CALL CANTK (NC,CONDT,TCDT,ZC,col,row)
!
!**** CALCULATE THE AVERAGE CONDUCTANCE TERM OVER THE TIME STEP
      CALL WEIGHT (NC,CON,CONT,CONDT,col,row)
!xxxx
!XXX  if (nc .lt. 3) then
        contk = con(nc)
!XXX   else
!XXX    contk = con(3)
!XXX  end if
!
!**** DETERMINE THE MATRIX COEFFICIENTS FOR THE TOP LAYER
!
      A1(N+1)=A1(N+1) + WDT*CON(1)
      B1(N)=B1(N) - WDT*(CON(1)+DHEATC(1)*(1.-DTLDTC(1))) &
                    - (ZC(2)-ZC(1))/(2.*DT)*RHOA*CA
      C1(N)=C1(N) + WDT*CON(1)
      D1(N)=D1(N) - CON(1)*(WT*(TC(1)-TC(2))+WDT*(TCDT(1)-TCDT(2))) &
                    + HEATC(1) -(ZC(2)-ZC(1))/(2.*DT)*RHOA*CA*(TCDT(1)-TC(1))
!
!**** DETERMINE THE COEFFICIENTS FOR THE REST OF THE LAYERS
      DO 10 I=N+1,N+NC-1
         J=I-N+1
         A1(I+1)=A1(I+1) + WDT*CON(J)
         B1(I)=B1(I) - WDT*(CON(J-1) +CON(J) +DHEATC(J)*(1.-DTLDTC(J))) &
                - (ZC(J+1)-ZC(J-1))/(2.*DT)*RHOA*CA
         C1(I)=C1(I) + WDT*CON(J)
         D1(I)=CON(J-1)*(WDT*(TCDT(J-1)-TCDT(J)) + WT*(TC(J-1)-TC(J))) &
            -CON(J)*(WDT*(TCDT(J)-TCDT(J+1)) + WT*(TC(J)-TC(J+1))) &
            -(ZC(J+1)-ZC(J-1))/(2.*DT)*RHOA*CA*(TCDT(J)-TC(J)) &
            +HEATC(J)
   10 CONTINUE
      N=N+NC
!
!**** DETERMINE THE COEFFICIENTS FOR THE TOP LAYER OF THE NEXT MATERIAL
      IF (NSP.GT.0 .OR. NR.EQ.0) THEN
!        NEXT MATERIAL IS SNOW OR SOIL -- NEED TO CONSIDER LATENT HEAT
!        TRANSFER TO THE NEXT NODE.
!        DETERMINE THE CONVECTIVE VAPOR TRANSPORT FROM EDDY CONDUCTANCE
         CONV=CON(NC)/RHOA/CA
!        CALCULATE THE VAPOR TRANSFER BETWEEN NODES
         QVCAN=CONV*(WDT*(VAPCDT(NC)-VAPCDT(NC+1)) &
                + WT*(VAPC(NC)-VAPC(NC+1)))
!        CALCULATE HUMIDITY AND SLOPE OF SAT. VAPOR CURVE
            if(abs(TC(NC)) .gt. 120.0)then
              print*,"In EBCAN,TC(NC)",col,row,TC(NC)
              stop
            end if
            if(abs(TC(NC+1)) .gt. 120.0)then
              print*,"In EBCAN,TC(NC+1)",col,row,TC(NC+1)
              stop
            end if
         CALL VSLOPE (SLPNC,SATV,TC(NC))
         HUMNC=VAPCDT(NC)
         CALL VSLOPE (SLPNC1,SATV,TC(NC+1))
         HUMNC1=VAPCDT(NC+1)/SATV
      END IF
!
      IF (NSP.GT.0) THEN
!****    NEXT MATERIAL IS SNOW
         A1(N) = A1(N) + WDT*CONV*LS*SLPNC*HUMNC
         B1(N) = B1(N) - WDT*(CON(NC) + CONV*LS*SLPNC1*HUMNC1)
         D1(N) =CON(NC)*(WT*(TC(NC)-TC(NC+1))+WDT*(TCDT(NC)-TCDT(NC+1))) &
             + LS*QVCAN
!xxxx
!xxxx    xlenc = LS*QVCAN
        ELSE
!
      IF (NR.GT.0) THEN
!****    NEXT MATERIAL IS RESIDUE
         B1(N) = B1(N) - WDT*CON(NC)
         D1(N) =CON(NC)*(WT*(TC(NC)-TC(NC+1))+WDT*(TCDT(NC)-TCDT(NC+1)))
!xxxx
!xxxx    xlenc = 0.0
        ELSE
!
!**** NEXT MATERIAL IS SOIL
      A1(N) = A1(N) + WDT*CONV*LV*SLPNC*HUMNC
      B1(N) = B1(N) - WDT*(CON(NC) + CONV*LV*SLPNC1*HUMNC1)
      D1(N) = CON(NC)*(WT*(TC(NC)-TC(NC+1)) + WDT*(TCDT(NC)-TCDT(NC+1))) &
        + LV*QVCAN
!xxxx
!xxxx xlenc = LV*QVCAN
!
      END IF
      END IF
!xxxxx
      hnc = CON(NC)*(WT*(TC(NC)-TC(NC+1)) + WDT*(TCDT(NC)-TCDT(NC+1)))
!
      RETURN
 END SUBROUTINE EBCAN



!no problem checked changed
!***********************************************************************
      subroutine leaft (nplant,nc,ns,iter,itype,tc,tcdt,tlc,tlcdt, &
        vapc,vapcdt,wcan,wcandt,pcan,pcandt,mat,matdt,trnsp,xtract, &
        swcan,lwcan,heatc,etlyr,detlyr,canma,canmb,dchar,rstom0,rstexp, &
        pleaf0,rleaf0,dheatc,dtldtc,col,row,julian,errorflag)
!     this subroutine computes leaf temperature of each canopy type and
!     the total heat and water transferred from canopy to surrounding
!     air space in each canopy layer
!***********************************************************************
        use controlpara_mod,only:wt2d,wdt2d,dt2d
        use constvar_mod,only:ca,cl,cr,emitc,lv,rhoa,rhol,stefan
        use clayrs_mod,only:totlai2d,totrot2d,ievap2d,drycan2d,canlai2d,rleaf2d,rroot2d,rootdn2d,pxylem2d,rhcan2d,etcan2d,init2d
        use windv_mod,only:windc2d
        use radcan_mod,only:tlclwr2d,dirkl2d,difkl2d,tdircc2d,tdiffc2d
        use writeit_mod,only:tleaf2d
        implicit none
!input
        integer(i4),intent(in)::nplant,nc,ns,iter,col,row,julian
        real(r8),dimension(ncmax),intent(in)::tc,tcdt,vapc,vapcdt
        real(r8),dimension(npmax,ncmax-1),intent(in)::tlc
        real(r8),dimension(npmax,ncmax-1),intent(inout)::tlcdt
        real(r8),dimension(ncmax-1),intent(inout)::wcan
        real(r8),dimension(ncmax-1),intent(inout)::wcandt
        real(r8),dimension(npmax),intent(inout)::pcandt
        real(r8),dimension(nsmax),intent(inout)::xtract
        real(r8),dimension(npmax+1),intent(inout)::trnsp
        real(r8),dimension(ncmax-1),intent(inout)::heatc,etlyr,detlyr,dheatc,dtldtc
        integer(i4),intent(inout)::errorflag
        real(r8),dimension(npmax),intent(in):: pcan
        real(r8),dimension(nsmax),intent(in)::mat,matdt
        real(r8),dimension(npmax+1,ncmax-1),intent(in)::swcan,lwcan
        real(r8),dimension(npmax),intent(in)::dchar,rstom0,rstexp,pleaf0,rleaf0
        real(r8),intent(in)::canma,canmb
        integer(i4),dimension(npmax),intent(in)::itype

 !temp
        real(r8),dimension(ncmax-1)::avgtmp,avgvap,rstom,pleaf,pevap
        real(r8),dimension(ncmax-1)::humid,vtslop
        real(r8),dimension(nsmax)::avgmat
        real(r8),dimension(ncmax-1)::f1,df1dp,df1dt,df1dx,dp,aa1,cc1,f2,df2dp,df2dt,df2dx,dtlc,df3dp
        real(r8),dimension(npmax,ncmax-1)::a1lwr

        real(r8)::a1,b1,d1,a2,b2,d2,aneg,b1lv,cc3,compar,delmax,delta,deriv,det,df3dx,dummy,dx, &
           error,f3,etmax,ff1,ff2,humcan,hum,hwslop,pleaf1,resist,rhavg,rslog,rssoil,rstom1,satv, &
           soimat,srroot,vapdef,wchum,rsoil
        integer(i4)::i,j,iter0,iroot,iflag,iter1,max,min

        real(r8)::sumet,sumev,fractn

        real(r8),dimension(:),pointer::totlai,totrot
        real(r8),dimension(:,:),pointer::tleaf
        real(r8),dimension(:,:),pointer::tlclwr
        real(r8),dimension(:,:),pointer::dirkl,difkl
        real(r8),dimension(:),pointer::tdircc,tdiffc
        real(r8),dimension(:),pointer::windc
        integer(i4),dimension(:),pointer::ievap
        real(r8),dimension(:,:),pointer::drycan,canlai,rleaf
        real(r8),dimension(:,:),pointer::rroot,rootdn
        real(r8),dimension(:),pointer::pxylem
        real(r8),dimension(:,:),pointer::rhcan,etcan
        integer(i4),dimension(:),pointer::init
        real(r8),pointer::dt,wdt,wt
        
        real(r8)::tlcDelta


        totlai  =>totlai2d(col,row,:)
        totrot  =>totrot2d(col,row,:)
        tleaf   =>tleaf2d(col,row,:,:)
        tlclwr  =>tlclwr2d(col,row,:,:)
        difkl   =>difkl2d(col,row,:,:)
        dirkl   =>dirkl2d(col,row,:,:)
        tdircc  =>tdircc2d(col,row,:)
        tdiffc  =>tdiffc2d(col,row,:)
        windc   =>windc2d(col,row,:)
        ievap   =>ievap2d(col,row,:)
        canlai  =>canlai2d(col,row,:,:)
        drycan  =>drycan2d(col,row,:,:)
        rleaf   =>rleaf2d(col,row,:,:)
        rroot   =>rroot2d(col,row,:,:)
        rootdn  =>rootdn2d(col,row,:,:)
        pxylem  =>pxylem2d(col,row,:)
        rhcan   =>rhcan2d(col,row,:,:)
        etcan   =>etcan2d(col,row,:,:)
        init    =>init2d(col,row,:)
        dt      =>dt2d(col,row)
        wt      =>wt2d(col,row)
        wdt     =>wdt2d(col,row)



!     initialize root extraction and compute average matric potential
        do i=1,ns
          xtract(i)=0.0
          avgmat(i)=wt*mat(i)+wdt*matdt(i)
!         limit the water potentials that the plant sees
          if (avgmat(i) .gt. 0.0) avgmat(i) = -0.0001
        end do
!
!     initialize the total transp. from the entire transpiring canopy
        trnsp(nplant+1)=0.0
!
!     initialize canopy temp, vapor density, and heat and water fluxes.
!     heatc(i) and etlyr are heat and water fluxes; dheatc and detlyr
!     are derivatives of flux terms.
        do i=1,nc
          avgtmp(i)=wt*tc(i) + wdt*tcdt(i)
          avgvap(i)=wt*vapc(i) + wdt*vapcdt(i)
          heatc(i)=0.0
          dheatc(i)=0.0
          etlyr(i)=0.0
          detlyr(i)=0.0
          dtldtc(i)=0.0
!        initialize resistance to transport from canopy leaves
!        and long-wave emittance coefficient (for both sides of leaves)
          do j=1,nplant
            if(difkl(nplant+1,i).ne.0) then
              a1lwr(j,i) = 8.*(1.-tdiffc(i))*difkl(j,i)/difkl(nplant+1,i)*emitc*stefan*((tlclwr(j,i)+273.16)**3)
            else
              print*,col,row,difkl(nplant+1,i)
              stop
            end if
            if(iter .eq. 1)then
!             rh not necessary if already calculated this time step
              if (windc(i).ne.0.0 .and. dchar(j)*windc(i).gt. 0.0) then
                rhcan(j,i)=307.*sqrt(dchar(j)/windc(i))
              else
                print*,col,row,windc(i),dchar(j)
                stop
              end if
            init(j)=1              
            endif
          end do
        end do
!!!! level1
    do 60 j=1,nplant
      sumet=0.0
      sumev = 0.0
      fractn = 1.0
!
      if (totlai(j) .eq. 0.0) then
!           no leaf area for this plant(perhaps dormant or snow-covered)
        trnsp(j)=0.0
        go to 60
      end if
!!!! level2
      if (itype(j) .ne. 0) then
!
!********   transpiring plant - check if conditions are such that plant
!           will transpire
        pcandt(j)=0.0

        ! level3
        if (ievap(j) .eq. 0 .or. pcan(j) .gt. 0.0) then
!*****         plant isn't transpiring - perhaps no sunlight, too cold,
!              or intercepted precip available on plant leaves
!              "fractn" is fraction of time step plants will transpire
               fractn=0.0
!              calculate leaf temperature
               do 12 i=1,nc
!                 check if plant has any leaf area in layer
                  if (canlai(j,i) .le. 0.0) go to 12
                  if (canlai(j,i) .lt. 0.05) canlai(j,i)=0.05
                  rstom(i)=0.0
                  etcan(j,i)=0.0
                  if (pcan(j) .gt. 0.0) then
!***                 intercepted precip available for evaporation
                     humid(i)=1.0
!xxxx                call vslope (vtslop(i),satv,avgtmp(i))
                     if(abs(tlclwr(j,i)) .gt. 120.0)then
                          print*,"In leaft, 4169",col,row,tlclwr(j,i)
                          stop
                     end if
                     call vslope (vtslop(i),satv,tlclwr(j,i))
                     vapdef = canlai(j,i)*(satv-avgvap(i))
                     tlcDelta=(swcan(j,i)+lwcan(j,i)-lv*vapdef/rhcan(j,i)-rhoa*ca*canlai(j,i)*(tlclwr(j,i)- &
                     avgtmp(i))/rhcan(j,i)-drycan(j,i)*cr*(tlclwr(j,i)-tlc(j,i))/dt)/(rhoa*ca*canlai(j,i)/rhcan(j,i) &
                     + a1lwr(j,i) + drycan(j,i)*cr/dt+lv*canlai(j,i)*vtslop(i)/rhcan(j,i))
                     if(abs(tlcDelta) .gt. 50.0) tlcDelta=tlcDelta*0.1
                     tlcdt(j,i)=tlclwr(j,i)+tlcDelta
                     if(abs(tlcdt(j,i)) .gt. 70.0) tlcdt(j,i)=tlclwr(j,i)
!                    determine amount of intercepted precip that
!                    evaporated from plant surfaces
                     pevap(i)=(canlai(j,i)*vtslop(i)*(tlcdt(j,i)-tlclwr(j,i))+vapdef)/rhcan(j,i)
                     sumev=sumev+pevap(i)
                  else
!***                 no water available for evaporation
                     humid(i)=0.0
                     vtslop(i)=0.0
                     pevap(i)=0.0
                     sumev=0.0
                     tlcDelta=(swcan(j,i)+lwcan(j,i)-rhoa*ca*canlai(j,i)*(tlclwr(j,i)-avgtmp(i))/rhcan(j,i) &
                     -drycan(j,i)*cr*(tlclwr(j,i)-tlc(j,i))/dt)/(rhoa*ca*canlai(j,i)/rhcan(j,i)+a1lwr(j,i) + drycan(j,i)*cr/dt)
                     if(abs(tlcDelta) .gt. 50.0) tlcDelta=tlcDelta*0.1
                     tlcdt(j,i)=tlclwr(j,i)+tlcDelta
                     if(abs(tlcdt(j,i)) .gt. 70.0) tlcdt(j,i)=tlclwr(j,i)
                  end if
               12 continue
               if (pcan(j) .gt. 0.0) then
!*****            calculate water on plants at end of time step
                  pcandt(j)=pcan(j)-dt*sumev/rhol
                  if (pcandt(j) .lt. 0.0) then
!                    no water remaining on plants - compute fraction of
!                    time step plants will be transpiring and adjust
!                    amount evaporated
                     pcandt(j)=0.0
                     fractn= 1 - (rhol*pcan(j)/dt)/sumev
                     sumev = (1.-fractn)*sumev
                     do 13 i=1,nc
                        pevap(i)=pevap(i)*(1.-fractn)
   13                continue
                  end if
               end if
        ! level 3
        end if
!
        ! level 3
        if (pcandt(j) .le. 0.0 .and. ievap(j) .ne. 0) then
!********      plant is transpiring
!
!              find the extremes of matric potential seen by roots
!              taking into account root resistance for each layer
               max=0
               do 14 i=1,ns
                  if (rootdn(j,i).gt.0.0) then
                    if (max .eq. 0) then
                       max=i
                       min=i
                     else
                       if (avgmat(i).lt.avgmat(min)) min=i
                       if (avgmat(i)/rroot(j,i) .gt. &
                        avgmat(max)/rroot(j,max)) max=i
                    end if
                  end if
               14 continue
!
!              determine leaf potential and temp
               vapdef=0.0
               rhavg=0.0
!              initize variables for canopy
               do 15 i=1,nc
!                 check if plant has any leaf area in this layer
                  if (canlai(j,i) .gt. 0.0) then
                     if (fractn .gt. 0.999) pevap(i)=0.0
!xxxx                tlc(j,i) = avgtmp(i)
                     humid(i)=1.0
!                    calculate average conditions in canopy for an
!                    intiial approximation to pleaf and pxylem
                  !if(tlcdt(j,i).le.-273.0) then
                  !   errorflag=1
                  !   return
                  !end if
                     if( abs(tlcdt(j,i)) .gt. 120.0) then
                       print*,"In leaft, 4246",col,row,tlcdt(j,i)
                       stop
                     end if
                     call vslope (vtslop(i),satv,tlcdt(j,i))
                     vapdef = vapdef + canlai(j,i)*(satv-avgvap(i))
                     rhavg = rhavg + canlai(j,i)/rhcan(j,i)
                  end if
               15 continue
               vapdef=vapdef/totlai(j)
               if (vapdef .lt. 0.0) vapdef=0.0
               rhavg=totlai(j)/rhavg
!
!*****         begin iteration to find initial estimate for pxylem
               iter0=0
!              calculate initial guess for pxylem if this is first time
!              for this time step
   18          if (init(j) .eq. 1) pxylem(j)=2.*avgmat(max)
!***           compute sum of root conductance times matric potential
   19          rsoil=0.0
               srroot=0.0
               iroot=0
               do 20 i=1,ns
!                 do not consider if no plant j roots in soil layer
                  if (rootdn(j,i).gt. 0.0) then
!                    do not include soil layers dryer than plant xylem
                     if (avgmat(i) .ge. pxylem(j)) then
                        rsoil=rsoil + rroot(j,i)*avgmat(i)
                        srroot= srroot + rroot(j,i)
                       else
                        iroot=1
                     end if
                  end if
               20 continue
               if (init(j) .gt. 1) then
!                 use values for pxylem and etcan from previous
!                 calculations for this time step if available
!                 (calculate sumet and go directly to iterative scheme)
                  sumet=rsoil - pxylem(j)*srroot
                  if (sumet*srroot .le. 0.0) then
!                    previous value of pxylem will not work for updated
!                    end-of-time-step conditions
                     init(j) = 1
                     go to 18
                  end if
                  go to 24
               end if
!***           calc. effective matric pot. and total resistance of plant
               if (srroot .eq. 0.0) then
                  rsoil=rroot(j,max)*avgmat(max)
                  srroot=rroot(j,max)
               end if
               soimat=rsoil/srroot
               resist= 1/srroot + 1/rleaf0(j)
               if (iter0 .eq. 0) then
!                 estimate starting point for iteration
                  pleaf1=pxylem(j)
                  if (pleaf1/pleaf0(j) .gt. 40.) then
!                    likelihood of arithmetic overflow -- proceed with
!                    caution by taking logarithm
                     rslog=log10(rstom0(j)) &
                     + rstexp(j)*log10(pleaf1/pleaf0(j))
                     if (rslog .gt. 20) then
!                       log of stomatal resistance extremely large --
!                       transp is essentially zero - calc leaf temp
                        sumet=0.0
                        do 21 i=1,nc
                           if (canlai(j,i) .le. 0.0) go to 21
                           humid(i)=0.0
                           rstom(i)=1.0e20
                           etcan(j,i)=0.0
                           tlcDelta=(swcan(j,i)+lwcan(j,i)-rhoa*ca*canlai(j,i)*(tlclwr(j,i)-avgtmp(i))/rhcan(j,i) &
                           -drycan(j,i)*cr*(tlclwr(j,i)-tlc(j,i))/dt)/(rhoa*ca*canlai(j,i)/rhcan(j,i) &
                           + a1lwr(j,i) + drycan(j,i)*cr/dt)
                           if(abs(tlcDelta) .gt. 50.0) tlcDelta=tlcDelta*0.1
                           tlcdt(j,i)=tlclwr(j,i)+tlcDelta
                           if(abs(tlcdt(j,i)) .gt. 70.0) tlcdt(j,i)=tlclwr(j,i)
                        21 continue
                        go to 48
                     end if
                  end if
                  rstom1=rstom0(j)*(1. + (pleaf1/pleaf0(j))**rstexp(j))
                  sumet=totlai(j)*vapdef/(rstom1+rhavg)
                  pleaf1=soimat-sumet*resist/2.
               end if
!***           update stomatal resistance and transpiration
   22          rstom1=rstom0(j)*(1. + (pleaf1/pleaf0(j))**rstexp(j))
               sumet=totlai(j)*vapdef/(rstom1+rhavg)
!              calculate error in et estimate, derivative with respect
!              to leaf potential, and new approx. to leaf potential
               error = (soimat-pleaf1)/resist - sumet
               deriv = -1/resist + sumet*rstom0(j)*rstexp(j) &
                    *(pleaf1/pleaf0(j))**(rstexp(j)-1) &
                    /pleaf0(j)/(rstom1+rhavg)
               delta=error/deriv
!***           depending on magnitude of rstexp, a drastic point of
!              inflection occurs in the error function at pleaf0. if
!              updated pleaf1 crosses this point, cut delta in half
               aneg=(pleaf1-pleaf0(j))*(pleaf1-delta-pleaf0(j))
               if (aneg .lt. 0.0) then
                  pleaf1=pleaf1-delta/2.
                 else
                  pleaf1=pleaf1-delta
               end if
!              calculate updated et and xylem potential
               sumet = (soimat-pleaf1)/resist
               pxylem(j)=(rsoil-sumet)/srroot
               if (abs(pleaf1) .lt. 1.0) then
!                 avoid division by zero
                  compar = 1.0
                 else
                  compar = pleaf1
               end if
               if (abs(delta/compar).gt.0.01 .and. iter0.le.20) then
!***              pleaf and pxylem not close enough
                  iter0=iter0+1
!                 if pxylem > minimum soil potential, recalculate
!                 rsoil, srroot and pxylem to exclude dry layers
                  if (pxylem(j).gt.avgmat(min) .or. iroot.gt.0) go to 19
                  go to 22
               end if
!*****         estimate transp (etcan) within each layer to begin iter
               do 23 i=1,nc
                  etcan(j,i)=canlai(j,i)*sumet/totlai(j)
   23          continue
               init(j) = 2
!
!*****         begin iteration to find leaf temperature, leaf potential
!              and transpiration from each canopy layer for plant
   24          iter1 = 0
               do 25 i=1,nc
!                 intial estimate of leaf potential in each layer
!                 and define newton-raphson coeff. that are constant
                  if (canlai(j,i) .gt. 0.0) then
                     pleaf(i)=pxylem(j) - etcan(j,i)/rleaf(j,i)
                     df1dx(i) = -rleaf(j,i)
                     df2dt(i) = -canlai(j,i)*rhoa*ca/rhcan(j,i) &
                        -a1lwr(j,i) - drycan(j,i)*cr/dt
                  end if
   25          continue

!
   26          iflag= 0
               sumet=0.0
               df3dx=0.0
!              set up coefficients for newton raphson solution
               do 30 i=1,nc
                  if (canlai(j,i) .gt. 0.0) then
                  if(tlcdt(j,i).le.-273.15) then
                        print*,col,row,tlcdt(j,i),j,i
                  !   errorflag=1
                  !   return
                  end if
                     if( abs(tlcdt(j,i)) .gt. 120.0) then
                       print*,"In leaft, 4398",col,row,tlcdt(j,i)
                       stop
                     end if
                  call vslope (vtslop(i),satv,tlcdt(j,i))
                  rstom(i)=rstom0(j)*(1+(pleaf(i)/pleaf0(j))**rstexp(j))
                  f1(i) =canlai(j,i)*(satv-avgvap(i)) &
                    /(rstom(i)+rhcan(j,i))
                  if (f1(i) .lt. 0.0) then
!***                 no transpiration
                     vtslop(i)=0.0
                     etcan(j,i)=0.0
                     df1dp(i) = rleaf(j,i)
                     df1dt(i) = 0.0
                     df2dp(i) = 0.0
                     df2dx(i) = 0.0
                     df3dp(i) = 0.0
!                    force leaf potential equal to pxylem potential,
!                    i.e. force transpiration in layer to zero
                     f1(i) = 0.0
                     pleaf(i) = pxylem(j)
                    else
!***                 calculate tranpiration in each layer and set up
!***                 matrix for newton-raphson approx. of updated xylem
!                    potential, leaf temperature and leaf potential
                     etcan(j,i)=rleaf(j,i)*(pxylem(j)-pleaf(i))
                     sumet = sumet +etcan(j,i)
                     df1dp(i) = rleaf(j,i) - f1(i)*rstom0(j)*rstexp(j) &
                        *(pleaf(i)/pleaf0(j))**(rstexp(j)-1) &
                        /pleaf0(j)/(rstom(i)+rhcan(j,i))
                     df1dt(i) = canlai(j,i)*vtslop(i) &
                        /(rstom(i)+rhcan(j,i))
                     df2dp(i) = lv*rleaf(j,i)
                     df2dx(i) = -df2dp(i)
                     df3dp(i) = rleaf(j,i)
                     df3dx = df3dx - rleaf(j,i)
                     f1(i) = f1(i) - etcan(j,i)
                  end if
                  f2(i)= swcan(j,i) + lwcan(j,i) - lv*etcan(j,i) &
                        -a1lwr(j,i)*(tlcdt(j,i)-tlclwr(j,i)) &
                        -canlai(j,i)*rhoa*ca*(tlcdt(j,i)-avgtmp(i))/rhcan(j,i) &
                        -drycan(j,i)*cr*(tlcdt(j,i)-tlc(j,i))/dt
                  end if
   30          continue
               f3 = rsoil - pxylem(j)*srroot - sumet
               df3dx = df3dx - srroot
               cc3=df3dx
!
               if (sumet .le. 0.0) then
!                 set tranpiration to zero and calculate leaf temp
                  sumet=0.0
                  do 31 i=1,nc
                     if (canlai(j,i) .gt. 0.0) then
                     humid(i)=0.0
                     rstom(i)=1.0e20
                     etcan(j,i)=0.0
                     tlcDelta=(swcan(j,i)+lwcan(j,i)-rhoa*ca*canlai(j,i)*(tlclwr(j,i)-avgtmp(i))/rhcan(j,i) &
                -drycan(j,i)*cr*(tlclwr(j,i)-tlc(j,i))/dt)/(-df2dt(i))
                     if(abs(tlcDelta) .gt. 50.0) tlcDelta=tlcDelta*0.1
                     tlcdt(j,i)=tlclwr(j,i)+tlcDelta
                     if(abs(tlcdt(j,i)) .gt. 70.0) tlcdt(j,i)=tlclwr(j,i)
                     init(j)=1
                     end if
   31             continue
                  go to 48
               end if
!
!              solve matrix for change in pxylem(j) and tlcdt(j,i)
               do 32 i=1,nc
                  if (canlai(j,i) .gt. 0.0) then
                  aa1(i)=df1dp(i)-(df1dt(i)/df2dt(i))*df2dp(i)
                  cc1(i)=df1dx(i)-(df1dt(i)/df2dt(i))*df2dx(i)
                  f1(i) =  f1(i) -(df1dt(i)/df2dt(i))*f1(i)
                  cc3=cc3-(df3dp(i)/aa1(i))*cc1(i)
                  f3 = f3-(df3dp(i)/aa1(i))*f1(i)
                  end if
   32          continue
               dx=f3/cc3
               pxylem(j)=pxylem(j)-dx
!
!              solve matrix for change in pleaf(i) and tlcdt(j,i)
               do 33 i=1,nc
                  if (canlai(j,i) .gt. 0.0) then
                  dp(i)=(f1(i)-cc1(i)*dx)/aa1(i)
                  dtlc(i)=(f2(i)-df2dx(i)*dx-df2dp(i)*dp(i))/df2dt(i)
!                 adjust leaf temp & potential
                  pleaf(i)=pleaf(i)-dp(i)
                  if (pleaf(i) .gt. pxylem(j)) &
                    pleaf(i)=(pleaf(i)+dp(i)+pxylem(j))/2.
                  tlcdt(j,i)=tlcdt(j,i)-dtlc(i)
                     if(abs(tlcdt(j,i)) .gt. 70.0) tlcdt(j,i)=tlclwr(j,i)
                     !print*,"2222",f2(i),-df2dx(i),dx,-df2dp(i),dp(i)
!                 check if temperature change is within 0.01 c
                  if (abs(dtlc(i)) .gt. 0.01) iflag=iflag+1
                  end if
   33          continue
!
!*****         check if tolerances have been met
               if (iflag .gt. 0) then
                  iter1 = iter1 + 1
                  if (iter1 .lt. 20) go to 26
               end if
!
!*****         solution has been found for leaf temp and tranpiration
!              find final xylem potential for use in root extraction
               if (srroot*sumet .ne. 0.0) then
                  pxylem(j)=(rsoil-sumet)/srroot
                  if (pxylem(j).gt.avgmat(min) .or. iroot.gt.0) then
                     rsoil=0.0
                     srroot=0.0
                     do 34 i=1,ns
                     if (rootdn(j,i).gt. 0.0) then
!                       don't include soil layers dryer than plant xylem
                        if (avgmat(i) .ge. pxylem(j)) then
                           rsoil=rsoil + rroot(j,i)*avgmat(i)
                           srroot= srroot + rroot(j,i)
                        end if
                     end if
   34                continue
!                    recalculate water potential in xylem
                     pxylem(j)=(rsoil-sumet)/srroot
                  end if
!
!                 calculate root extraction from each soil layer
                  do 35 i=1,ns
                     if (avgmat(i) .gt. pxylem(j)) xtract(i) = xtract(i) &
                + fractn*totrot(j)*rroot(j,i)*(avgmat(i)-pxylem(j))
   35             continue
               else
!                 soil too dry for plants to extact water
                  sumet=0.0
               end if
!
!              sum transpiration from entire transpiring canopy
               trnsp(nplant+1)=trnsp(nplant+1)+fractn*sumet
        ! level 3
        end if
!!!   level2
      else
!********   dead plant material -- compute water content and temperature
            do 45 i=1,nc
               if (canlai(j,i) .le. 0.0) go to 45
!
               pevap(i)=0.0
               rstom(i)=0.0
               b1lv = drycan(j,i)/dt
               a1 = -a1lwr(j,i) - canlai(j,i)*rhoa*ca/rhcan(j,i) &
                    -drycan(j,i)*(cr+wcan(i)*cl)/dt
               b1=lv*b1lv
!              calculate water content based on vapor density of air
                  !if(tlcdt(j,i).le.-273.15) then
                  !   errorflag=1
                  !   return
                  !end if
                     if( abs(tlcdt(j,i)) .gt. 70.0) then
                       print*,"In leaft, 4550",col,row,tlcdt(j,i)
                       stop
                     end if
               call vslope (dummy,satv,tlcdt(j,i))
               humcan=avgvap(i)/satv
               call canhum (2,humcan,dummy,wchum,tcdt(i),canma,canmb)
!
               if (init(j) .eq. 1) then
!                 initialize variables for this time step
                  etcan(j,i)= -b1lv*(wcandt(i)-wcan(i))
                  if (etcan(j,i) .eq. 0.0) then
!xxx                 tlc(j,i) = avgtmp(i)
!                    compute humidity in plant material
                     call canhum (1,humid(i),dummy,wcandt(i),tcdt(i), &
                canma,canmb)
                     if( abs(tlcdt(j,i)) .gt. 70.0) then
                       print*,col,row,tlcdt(j,i)
                       stop
                     end if
                     call vslope (dummy,satv,tlcdt(j,i))
                     etcan(j,i)=canlai(j,i)* &
                (humid(i)*satv-avgvap(i))/rhcan(j,i)
                     wcandt(i) = wcan(i) - etcan(j,i)/b1lv
!                    check if water content is reasonable --
                     if ((wcandt(i)-wchum)*(wcan(i)-wchum).lt.0.0) then
!                       water content went beyond equilibruim with air
                        etcan(j,i)= -b1lv*(wchum-wcan(i))
                        humid(i)=(avgvap(i) + &
            etcan(j,i)*rhcan(j,i)/canlai(j,i))/satv
                        call canhum (2,humid(i),dummy,wcandt(i),tcdt(i), &
                canma,canmb)
                        etcan(j,i)= -b1lv*(wcandt(i)-wcan(i))
                     end if
                  end if
                  tlcdt(j,i)=-(swcan(j,i)+lwcan(j,i)-lv*etcan(j,i)+rhoa*ca*canlai(j,i)*avgtmp(i)/rhcan(j,i)&
                     +drycan(j,i)*(cr+wcan(i)*cl)*tlc(j,i)/dt) /a1
                     if(abs(tlcdt(j,i)) .gt. 70.0) tlcdt(j,i)=tlclwr(j,i)
                           if(abs(tlcdt(j,i)) .gt. 70.0) then 
                              print*,"The tlcdt is wrong 6",tlcdt(j,i)
                              print*,swcan(j,i)+lwcan(j,i)-lv*etcan(j,i)+rhoa*ca*canlai(j,i)*avgtmp(i)/rhcan(j,i)&
                     +drycan(j,i)*(cr+wcan(i)*cl)*tlc(j,i)/dt
                              print*,a1
                              print*,tlclwr(j,i),swcan(j,i),lwcan(j,i),canlai(j,i),avgtmp(i),rhcan(j,i)
                              stop
                           end if
               end if
!
                  !if(tlcdt(j,i).le.-273.15) then
                  !   errorflag=1
                  !   return
                  !end if
                     if( abs(tlcdt(j,i)) .gt. 120.0) then
                       print*,"In leaft, 4594",col,row,tlcdt(j,i)
                       stop
                     end if
               call vslope (vtslop(i),satv,tlcdt(j,i))
!
!*****         begin iterations to find leaf temp and water content
               iter1 = 0
   40          etcan(j,i)= -b1lv*(wcandt(i)-wcan(i))
!              compute humidity in plant material at end of time step
               call canhum(1,humid(i),hwslop,wcandt(i),tcdt(i), &
            canma,canmb)
!***           set up and solve 2x2 matric for newton-raphson approx.
!              for leaf temp and water content
               a2 = canlai(j,i)*humid(i)*vtslop(i)/rhcan(j,i)
               b2 = canlai(j,i)*satv*hwslop/rhcan(j,i) + b1lv
               ff1 = swcan(j,i) + lwcan(j,i) - lv*etcan(j,i) &
        -a1lwr(j,i)*(tlcdt(j,i)-tlclwr(j,i)) &
        -canlai(j,i)*rhoa*ca*(tlcdt(j,i)-avgtmp(i))/rhcan(j,i) &
        -drycan(j,i)*(cr+wcan(i)*cl)*(tlcdt(j,i)-tlc(j,i))/dt
               ff2=canlai(j,i)*(humid(i)*satv-avgvap(i))/rhcan(j,i) &
        - etcan(j,i)
               det = a1*b2 - a2*b1
               d1 = (ff1*b2 - ff2*b1)/det
               d2 = (a1*ff2 - a2*ff1)/det
!
!***           update values
               tlcdt(j,i)=tlcdt(j,i) - d1
                     if(abs(tlcdt(j,i)) .gt. 70.0) tlcdt(j,i)=tlclwr(j,i)
                           if(abs(tlcdt(j,i)) .gt. 70.0) then 
                              print*,"The tlcdt is wrong 7",tlcdt(j,i),d1
                              print*,tlclwr(j,i),swcan(j,i),lwcan(j,i),canlai(j,i),avgtmp(i),rhcan(j,i)
                              stop
                           end if
               wcandt(i) = wcandt(i) - d2
!              check if water content is reasonable --
                  !if(tlcdt(j,i).le.-273.15) then
                  !   errorflag=1
                  !   return
                  !end if
                     if( abs(tlcdt(j,i)) .gt. 120.0) then
                       print*,"In leaft, 4628",col,row,tlcdt(j,i)
                       stop
                     end if
               call vslope (vtslop(i),satv,tlcdt(j,i))
               humcan=avgvap(i)/satv
               call canhum (2,humcan,dummy,wchum,tcdt(i),canma,canmb)
               delmax=wcan(i)-wchum
               delta=wcandt(i)-wchum
!cccc          if(delta*delmax.lt.0.0 .or.abs(delmax).lt.abs(delta))then
               if(delta*delmax.lt.0.0)then
!                 water content went beyond equilibruim with humidity
                  etcan(j,i)= -b1lv*(wchum-wcan(i))
                  call canhum (1,hum,dummy,wcan(i),tcdt(i),canma,canmb)
                  etmax=canlai(j,i)*(hum*satv-avgvap(i))/rhcan(j,i)
                  if (abs(etcan(j,i)) .gt. abs(etmax)) etcan(j,i)=etmax
                  humid(i)=(avgvap(i) + &
            etcan(j,i)*rhcan(j,i)/canlai(j,i))/satv
                  call canhum (2,humid(i),dummy,wcandt(i),tcdt(i), &
            canma,canmb)
               end if
               if (abs(d1) .gt. 0.01) then
                  iter1 = iter1 +1
                  if (iter1 .lt. 10) go to 40
               end if
!              calculate evaporation from this layer
               etcan(j,i)= -b1lv*(wcandt(i)-wcan(i))
               sumet=sumet+etcan(j,i)
   45       continue
            init(j) = 2
      end if
!
!********store evaporation/transpiration from plant species
   48    trnsp(j)=fractn*sumet+sumev
!
!        sum heat and water transfer in each layer from all canopy types
      do 50 i=1,nc
            if (canlai(j,i) .le. 0.0) go to 50
            heatc(i)=heatc(i) &
        + canlai(j,i)*rhoa*ca*(tlcdt(j,i)-avgtmp(i))/rhcan(j,i)
            etlyr(i)=etlyr(i) + fractn*etcan(j,i) + pevap(i)
            dheatc(i)=dheatc(i) &
        + canlai(j,i)*rhoa*ca/rhcan(j,i)
            detlyr(i)=detlyr(i) + canlai(j,i)/(rstom(i)+rhcan(j,i))
            dtldtc(i)=dtldtc(i) + lv*humid(i)*vtslop(i)* &
        canlai(j,i)/(rstom(i)+rhcan(j,i))
      50 continue
    60 continue
!
!xxxx
      do i=1,nc
         tleaf(nplant+1,i)=0.0
         do j=1,nplant
            tleaf(j,i)=tlcdt(j,i)
            if (canlai(j,i) .le. 0.0) tleaf(j,i)=0.0
            tleaf(nplant+1,i)=tleaf(nplant+1,i) &
        + tlcdt(j,i)*difkl(j,i)/difkl(nplant+1,i)
         end do
      end do
!
!     calculate change of leaf temp with respect to canopy temp
      do 70 i=1,nc
         dtldtc(i)=dheatc(i)/(dheatc(i)+dtldtc(i))
   70 continue
      return
end subroutine leaft


!no problem checked changed
!***********************************************************************
      subroutine leaft2 (nplant,nc,ns,iter,itype,tc,tcdt,tlc,tlcdt, &
        vapc,vapcdt,wcan,wcandt,pcan,pcandt,mat,matdt,trnsp,xtract, &
        swcan,lwcan,heatc,etlyr,detlyr,canma,canmb,dchar,rstom0,rstexp, &
        pleaf0,rleaf0,dheatc,dtldtc,col,row,julian,errorflag)
!     this subroutine computes leaf temperature of each canopy type and
!     the total heat and water transferred from canopy to surrounding
!     air space in each canopy layer
!***********************************************************************
        use controlpara_mod,only:wt2d,wdt2d,dt2d
        use constvar_mod,only:ca,cl,cr,emitc,lv,rhoa,rhol,stefan
        use clayrs_mod,only:totlai2d,totrot2d,ievap2d,drycan2d,canlai2d,rleaf2d,rroot2d,rootdn2d,pxylem2d,rhcan2d,etcan2d,init2d
        use windv_mod,only:windc2d
        use radcan_mod,only:tlclwr2d,dirkl2d,difkl2d,tdircc2d,tdiffc2d
        use writeit_mod,only:tleaf2d
        implicit none
!input
        integer(i4),intent(in)::nplant,nc,ns,iter,col,row,julian
        real(r8),dimension(ncmax),intent(in)::tc,tcdt,vapc,vapcdt
        real(r8),dimension(npmax,ncmax-1),intent(in)::tlc
        real(r8),dimension(npmax,ncmax-1),intent(inout)::tlcdt
        real(r8),dimension(ncmax-1),intent(inout)::wcan
        real(r8),dimension(ncmax-1),intent(inout)::wcandt
        real(r8),dimension(npmax),intent(inout)::pcandt
        real(r8),dimension(nsmax),intent(inout)::xtract
        real(r8),dimension(npmax+1),intent(inout)::trnsp
        real(r8),dimension(ncmax-1),intent(inout)::heatc,etlyr,detlyr,dheatc,dtldtc
        integer(i4),intent(inout)::errorflag
        real(r8),dimension(npmax),intent(in):: pcan
        real(r8),dimension(nsmax),intent(in)::mat,matdt
        real(r8),dimension(npmax+1,ncmax-1),intent(in)::swcan,lwcan
        real(r8),dimension(npmax),intent(in)::dchar,rstom0,rstexp,pleaf0,rleaf0
        real(r8),intent(in)::canma,canmb
        integer(i4),dimension(npmax),intent(in)::itype

 !temp
        real(r8),dimension(ncmax-1)::avgtmp,avgvap,rstom,pleaf,pevap
        real(r8),dimension(ncmax-1)::humid,vtslop
        real(r8),dimension(nsmax)::avgmat
        real(r8),dimension(ncmax-1)::f1,df1dp,df1dt,df1dx,dp,aa1,cc1,f2,df2dp,df2dt,df2dx,dtlc,df3dp
        real(r8),dimension(npmax,ncmax-1)::a1lwr

        real(r8)::a1,b1,d1,a2,b2,d2,aneg,b1lv,cc3,compar,delmax,delta,deriv,det,df3dx,dummy,dx, &
           error,f3,etmax,ff1,ff2,humcan,hum,hwslop,pleaf1,resist,rhavg,rslog,rssoil,rstom1,satv, &
           soimat,srroot,vapdef,wchum,rsoil
        integer(i4)::i,j,iter0,iroot,iflag,iter1,max,min

        real(r8)::sumet,sumev,fractn

        real(r8),dimension(:),pointer::totlai,totrot
        real(r8),dimension(:,:),pointer::tleaf
        real(r8),dimension(:,:),pointer::tlclwr
        real(r8),dimension(:,:),pointer::dirkl,difkl
        real(r8),dimension(:),pointer::tdircc,tdiffc
        real(r8),dimension(:),pointer::windc
        integer(i4),dimension(:),pointer::ievap
        real(r8),dimension(:,:),pointer::drycan,canlai,rleaf
        real(r8),dimension(:,:),pointer::rroot,rootdn
        real(r8),dimension(:),pointer::pxylem
        real(r8),dimension(:,:),pointer::rhcan,etcan
        integer(i4),dimension(:),pointer::init
        real(r8),pointer::dt,wdt,wt
        
        real(r8)::tlcDelta


        totlai  =>totlai2d(col,row,:)
        totrot  =>totrot2d(col,row,:)
        tleaf   =>tleaf2d(col,row,:,:)
        tlclwr  =>tlclwr2d(col,row,:,:)
        difkl   =>difkl2d(col,row,:,:)
        dirkl   =>dirkl2d(col,row,:,:)
        tdircc  =>tdircc2d(col,row,:)
        tdiffc  =>tdiffc2d(col,row,:)
        windc   =>windc2d(col,row,:)
        ievap   =>ievap2d(col,row,:)
        canlai  =>canlai2d(col,row,:,:)
        drycan  =>drycan2d(col,row,:,:)
        rleaf   =>rleaf2d(col,row,:,:)
        rroot   =>rroot2d(col,row,:,:)
        rootdn  =>rootdn2d(col,row,:,:)
        pxylem  =>pxylem2d(col,row,:)
        rhcan   =>rhcan2d(col,row,:,:)
        etcan   =>etcan2d(col,row,:,:)
        init    =>init2d(col,row,:)
        dt      =>dt2d(col,row)
        wt      =>wt2d(col,row)
        wdt     =>wdt2d(col,row)



!     initialize root extraction and compute average matric potential
        do i=1,ns
          xtract(i)=0.0
          avgmat(i)=wt*mat(i)+wdt*matdt(i)
!         limit the water potentials that the plant sees
          if (avgmat(i) .gt. 0.0) avgmat(i) = -0.0001
        end do
!
!     initialize the total transp. from the entire transpiring canopy
        trnsp(nplant+1)=0.0
!
!     initialize canopy temp, vapor density, and heat and water fluxes.
!     heatc(i) and etlyr are heat and water fluxes; dheatc and detlyr
!     are derivatives of flux terms.
        do i=1,nc
          avgtmp(i)=wt*tc(i) + wdt*tcdt(i)
          avgvap(i)=wt*vapc(i) + wdt*vapcdt(i)
          heatc(i)=0.0
          dheatc(i)=0.0
          etlyr(i)=0.0
          detlyr(i)=0.0
          dtldtc(i)=0.0
!        initialize resistance to transport from canopy leaves
!        and long-wave emittance coefficient (for both sides of leaves)
          do j=1,nplant
            if(difkl(nplant+1,i).ne.0) then
              a1lwr(j,i) = 8.*(1.-tdiffc(i))*difkl(j,i)/difkl(nplant+1,i)*emitc*stefan*((tlclwr(j,i)+273.16)**3)
            else
              print*,col,row,difkl(nplant+1,i)
              stop
            end if
            if(iter .eq. 1)then
!             rh not necessary if already calculated this time step
              if (windc(i).ne.0.0 .and. dchar(j)*windc(i).gt. 0.0) then
                rhcan(j,i)=307.*sqrt(dchar(j)/windc(i))
              else
                print*,col,row,windc(i),dchar(j)
                stop
              end if
            init(j)=1              
            endif
          end do
        end do
!!!! level1
    do 60 j=1,nplant
      sumet=0.0
      sumev = 0.0
      fractn = 1.0
!
      if (totlai(j) .eq. 0.0) then
!           no leaf area for this plant(perhaps dormant or snow-covered)
        trnsp(j)=0.0
        go to 60
      end if
!!!! level2
      if (itype(j) .ne. 0) then
!
!********   transpiring plant - check if conditions are such that plant
!           will transpire
        pcandt(j)=0.0

        ! level3
        if (ievap(j) .eq. 0 .or. pcan(j) .gt. 0.0) then
!*****         plant isn't transpiring - perhaps no sunlight, too cold,
!              or intercepted precip available on plant leaves
!              "fractn" is fraction of time step plants will transpire
               fractn=0.0
!              calculate leaf temperature
               do 12 i=1,nc
!                 check if plant has any leaf area in layer
                  if (canlai(j,i) .le. 0.0) go to 12
                  if (canlai(j,i) .lt. 0.05) canlai(j,i)=0.05
                  rstom(i)=0.0
                  etcan(j,i)=0.0
                  if (pcan(j) .gt. 0.0) then
!***                 intercepted precip available for evaporation
                     humid(i)=1.0
!xxxx                call vslope (vtslop(i),satv,avgtmp(i))
                     if(abs(tlclwr(j,i)) .gt. 120.0)then
                          print*,"In leaft, 4169",col,row,tlclwr(j,i)
                          stop
                     end if
                     call vslope (vtslop(i),satv,tlclwr(j,i))
                     vapdef = canlai(j,i)*(satv-avgvap(i))
                     tlcdt(j,i)=avgtmp(i)
!                    determine amount of intercepted precip that
!                    evaporated from plant surfaces
                     pevap(i)=(canlai(j,i)*vtslop(i)*(tlcdt(j,i)-tlclwr(j,i))+vapdef)/rhcan(j,i)
                     sumev=sumev+pevap(i)
                  else
!***                 no water available for evaporation
                     humid(i)=0.0
                     vtslop(i)=0.0
                     pevap(i)=0.0
                     sumev=0.0
                     tlcdt(j,i)=avgtmp(i)
                  end if
               12 continue
               if (pcan(j) .gt. 0.0) then
!*****            calculate water on plants at end of time step
                  pcandt(j)=pcan(j)-dt*sumev/rhol
                  if (pcandt(j) .lt. 0.0) then
!                    no water remaining on plants - compute fraction of
!                    time step plants will be transpiring and adjust
!                    amount evaporated
                     pcandt(j)=0.0
                     fractn= 1 - (rhol*pcan(j)/dt)/sumev
                     sumev = (1.-fractn)*sumev
                     do 13 i=1,nc
                        pevap(i)=pevap(i)*(1.-fractn)
   13                continue
                  end if
               end if
        ! level 3
        end if
!
        ! level 3
        if (pcandt(j) .le. 0.0 .and. ievap(j) .ne. 0) then
!********      plant is transpiring
!
!              find the extremes of matric potential seen by roots
!              taking into account root resistance for each layer
               max=0
               do 14 i=1,ns
                  if (rootdn(j,i).gt.0.0) then
                    if (max .eq. 0) then
                       max=i
                       min=i
                     else
                       if (avgmat(i).lt.avgmat(min)) min=i
                       if (avgmat(i)/rroot(j,i) .gt. &
                        avgmat(max)/rroot(j,max)) max=i
                    end if
                  end if
               14 continue
!
!              determine leaf potential and temp
               vapdef=0.0
               rhavg=0.0
!              initize variables for canopy
               do 15 i=1,nc
!                 check if plant has any leaf area in this layer
                  if (canlai(j,i) .gt. 0.0) then
                     if (fractn .gt. 0.999) pevap(i)=0.0
!xxxx                tlc(j,i) = avgtmp(i)
                     humid(i)=1.0
!                    calculate average conditions in canopy for an
!                    intiial approximation to pleaf and pxylem
                  !if(tlcdt(j,i).le.-273.0) then
                  !   errorflag=1
                  !   return
                  !end if
                     if( abs(tlcdt(j,i)) .gt. 120.0) then
                       print*,"In leaft, 4246",col,row,tlcdt(j,i)
                       stop
                     end if
                     call vslope (vtslop(i),satv,tlcdt(j,i))
                     vapdef = vapdef + canlai(j,i)*(satv-avgvap(i))
                     rhavg = rhavg + canlai(j,i)/rhcan(j,i)
                  end if
               15 continue
               vapdef=vapdef/totlai(j)
               if (vapdef .lt. 0.0) vapdef=0.0
               rhavg=totlai(j)/rhavg
!
!*****         begin iteration to find initial estimate for pxylem
               iter0=0
!              calculate initial guess for pxylem if this is first time
!              for this time step
   18          if (init(j) .eq. 1) pxylem(j)=2.*avgmat(max)
!***           compute sum of root conductance times matric potential
   19          rsoil=0.0
               srroot=0.0
               iroot=0
               do 20 i=1,ns
!                 do not consider if no plant j roots in soil layer
                  if (rootdn(j,i).gt. 0.0) then
!                    do not include soil layers dryer than plant xylem
                     if (avgmat(i) .ge. pxylem(j)) then
                        rsoil=rsoil + rroot(j,i)*avgmat(i)
                        srroot= srroot + rroot(j,i)
                       else
                        iroot=1
                     end if
                  end if
               20 continue
               if (init(j) .gt. 1) then
!                 use values for pxylem and etcan from previous
!                 calculations for this time step if available
!                 (calculate sumet and go directly to iterative scheme)
                  sumet=rsoil - pxylem(j)*srroot
                  if (sumet*srroot .le. 0.0) then
!                    previous value of pxylem will not work for updated
!                    end-of-time-step conditions
                     init(j) = 1
                     go to 18
                  end if
                  go to 24
               end if
!***           calc. effective matric pot. and total resistance of plant
               if (srroot .eq. 0.0) then
                  rsoil=rroot(j,max)*avgmat(max)
                  srroot=rroot(j,max)
               end if
               soimat=rsoil/srroot
               resist= 1/srroot + 1/rleaf0(j)
               if (iter0 .eq. 0) then
!                 estimate starting point for iteration
                  pleaf1=pxylem(j)
                  if (pleaf1/pleaf0(j) .gt. 40.) then
!                    likelihood of arithmetic overflow -- proceed with
!                    caution by taking logarithm
                     rslog=log10(rstom0(j)) &
                     + rstexp(j)*log10(pleaf1/pleaf0(j))
                     if (rslog .gt. 20) then
!                       log of stomatal resistance extremely large --
!                       transp is essentially zero - calc leaf temp
                        sumet=0.0
                        do 21 i=1,nc
                           if (canlai(j,i) .le. 0.0) go to 21
                           humid(i)=0.0
                           rstom(i)=1.0e20
                           etcan(j,i)=0.0
                           tlcdt(j,i)=avgtmp(i)
                        21 continue
                        go to 48
                     end if
                  end if
                  rstom1=rstom0(j)*(1. + (pleaf1/pleaf0(j))**rstexp(j))
                  sumet=totlai(j)*vapdef/(rstom1+rhavg)
                  pleaf1=soimat-sumet*resist/2.
               end if
!***           update stomatal resistance and transpiration
   22          rstom1=rstom0(j)*(1. + (pleaf1/pleaf0(j))**rstexp(j))
               sumet=totlai(j)*vapdef/(rstom1+rhavg)
!              calculate error in et estimate, derivative with respect
!              to leaf potential, and new approx. to leaf potential
               error = (soimat-pleaf1)/resist - sumet
               deriv = -1/resist + sumet*rstom0(j)*rstexp(j) &
                    *(pleaf1/pleaf0(j))**(rstexp(j)-1) &
                    /pleaf0(j)/(rstom1+rhavg)
               delta=error/deriv
!***           depending on magnitude of rstexp, a drastic point of
!              inflection occurs in the error function at pleaf0. if
!              updated pleaf1 crosses this point, cut delta in half
               aneg=(pleaf1-pleaf0(j))*(pleaf1-delta-pleaf0(j))
               if (aneg .lt. 0.0) then
                  pleaf1=pleaf1-delta/2.
                 else
                  pleaf1=pleaf1-delta
               end if
!              calculate updated et and xylem potential
               sumet = (soimat-pleaf1)/resist
               pxylem(j)=(rsoil-sumet)/srroot
               if (abs(pleaf1) .lt. 1.0) then
!                 avoid division by zero
                  compar = 1.0
                 else
                  compar = pleaf1
               end if
               if (abs(delta/compar).gt.0.01 .and. iter0.le.20) then
!***              pleaf and pxylem not close enough
                  iter0=iter0+1
!                 if pxylem > minimum soil potential, recalculate
!                 rsoil, srroot and pxylem to exclude dry layers
                  if (pxylem(j).gt.avgmat(min) .or. iroot.gt.0) go to 19
                  go to 22
               end if
!*****         estimate transp (etcan) within each layer to begin iter
               do 23 i=1,nc
                  etcan(j,i)=canlai(j,i)*sumet/totlai(j)
   23          continue
               init(j) = 2
!
!*****         begin iteration to find leaf temperature, leaf potential
!              and transpiration from each canopy layer for plant
   24          iter1 = 0
               do 25 i=1,nc
!                 intial estimate of leaf potential in each layer
!                 and define newton-raphson coeff. that are constant
                  if (canlai(j,i) .gt. 0.0) then
                     pleaf(i)=pxylem(j) - etcan(j,i)/rleaf(j,i)
                     df1dx(i) = -rleaf(j,i)
                     df2dt(i) = -canlai(j,i)*rhoa*ca/rhcan(j,i) &
                        -a1lwr(j,i) - drycan(j,i)*cr/dt
                  end if
   25          continue

!
   26          iflag= 0
               sumet=0.0
               df3dx=0.0
!              set up coefficients for newton raphson solution
               do 30 i=1,nc
                  if (canlai(j,i) .gt. 0.0) then
                  if(tlcdt(j,i).le.-273.15) then
                        print*,col,row,tlcdt(j,i),j,i
                  !   errorflag=1
                  !   return
                  end if
                     if( abs(tlcdt(j,i)) .gt. 120.0) then
                       print*,"In leaft, 4398",col,row,tlcdt(j,i)
                       stop
                     end if
                  call vslope (vtslop(i),satv,tlcdt(j,i))
                  rstom(i)=rstom0(j)*(1+(pleaf(i)/pleaf0(j))**rstexp(j))
                  f1(i) =canlai(j,i)*(satv-avgvap(i)) &
                    /(rstom(i)+rhcan(j,i))
                  if (f1(i) .lt. 0.0) then
!***                 no transpiration
                     vtslop(i)=0.0
                     etcan(j,i)=0.0
                     df1dp(i) = rleaf(j,i)
                     df1dt(i) = 0.0
                     df2dp(i) = 0.0
                     df2dx(i) = 0.0
                     df3dp(i) = 0.0
!                    force leaf potential equal to pxylem potential,
!                    i.e. force transpiration in layer to zero
                     f1(i) = 0.0
                     pleaf(i) = pxylem(j)
                    else
!***                 calculate tranpiration in each layer and set up
!***                 matrix for newton-raphson approx. of updated xylem
!                    potential, leaf temperature and leaf potential
                     etcan(j,i)=rleaf(j,i)*(pxylem(j)-pleaf(i))
                     sumet = sumet +etcan(j,i)
                     df1dp(i) = rleaf(j,i) - f1(i)*rstom0(j)*rstexp(j) &
                        *(pleaf(i)/pleaf0(j))**(rstexp(j)-1) &
                        /pleaf0(j)/(rstom(i)+rhcan(j,i))
                     df1dt(i) = canlai(j,i)*vtslop(i) &
                        /(rstom(i)+rhcan(j,i))
                     df2dp(i) = lv*rleaf(j,i)
                     df2dx(i) = -df2dp(i)
                     df3dp(i) = rleaf(j,i)
                     df3dx = df3dx - rleaf(j,i)
                     f1(i) = f1(i) - etcan(j,i)
                  end if
                  f2(i)= swcan(j,i) + lwcan(j,i) - lv*etcan(j,i) &
                        -a1lwr(j,i)*(tlcdt(j,i)-tlclwr(j,i)) &
                        -canlai(j,i)*rhoa*ca*(tlcdt(j,i)-avgtmp(i))/rhcan(j,i) &
                        -drycan(j,i)*cr*(tlcdt(j,i)-tlc(j,i))/dt
                  end if
   30          continue
               f3 = rsoil - pxylem(j)*srroot - sumet
               df3dx = df3dx - srroot
               cc3=df3dx
!
               if (sumet .le. 0.0) then
!                 set tranpiration to zero and calculate leaf temp
                  sumet=0.0
                  do 31 i=1,nc
                     if (canlai(j,i) .gt. 0.0) then
                     humid(i)=0.0
                     rstom(i)=1.0e20
                     etcan(j,i)=0.0
                     tlcdt(j,i)=avgtmp(i)
                     init(j)=1
                     end if
   31             continue
                  go to 48
               end if
!
!              solve matrix for change in pxylem(j) and tlcdt(j,i)
               do 32 i=1,nc
                  if (canlai(j,i) .gt. 0.0) then
                  aa1(i)=df1dp(i)-(df1dt(i)/df2dt(i))*df2dp(i)
                  cc1(i)=df1dx(i)-(df1dt(i)/df2dt(i))*df2dx(i)
                  f1(i) =  f1(i) -(df1dt(i)/df2dt(i))*f1(i)
                  cc3=cc3-(df3dp(i)/aa1(i))*cc1(i)
                  f3 = f3-(df3dp(i)/aa1(i))*f1(i)
                  end if
   32          continue
               dx=f3/cc3
               pxylem(j)=pxylem(j)-dx
!
!              solve matrix for change in pleaf(i) and tlcdt(j,i)
               do 33 i=1,nc
                  if (canlai(j,i) .gt. 0.0) then
                  dp(i)=(f1(i)-cc1(i)*dx)/aa1(i)
                  dtlc(i)=(f2(i)-df2dx(i)*dx-df2dp(i)*dp(i))/df2dt(i)
!                 adjust leaf temp & potential
                  pleaf(i)=pleaf(i)-dp(i)
                  if (pleaf(i) .gt. pxylem(j)) &
                    pleaf(i)=(pleaf(i)+dp(i)+pxylem(j))/2.
                     tlcdt(j,i)=avgtmp(i)
                     !print*,"2222",f2(i),-df2dx(i),dx,-df2dp(i),dp(i)
!                 check if temperature change is within 0.01 c
                  if (abs(dtlc(i)) .gt. 0.01) iflag=iflag+1
                  end if
   33          continue
!
!*****         check if tolerances have been met
               if (iflag .gt. 0) then
                  iter1 = iter1 + 1
                  if (iter1 .lt. 20) go to 26
               end if
!
!*****         solution has been found for leaf temp and tranpiration
!              find final xylem potential for use in root extraction
               if (srroot*sumet .ne. 0.0) then
                  pxylem(j)=(rsoil-sumet)/srroot
                  if (pxylem(j).gt.avgmat(min) .or. iroot.gt.0) then
                     rsoil=0.0
                     srroot=0.0
                     do 34 i=1,ns
                     if (rootdn(j,i).gt. 0.0) then
!                       don't include soil layers dryer than plant xylem
                        if (avgmat(i) .ge. pxylem(j)) then
                           rsoil=rsoil + rroot(j,i)*avgmat(i)
                           srroot= srroot + rroot(j,i)
                        end if
                     end if
   34                continue
!                    recalculate water potential in xylem
                     pxylem(j)=(rsoil-sumet)/srroot
                  end if
!
!                 calculate root extraction from each soil layer
                  do 35 i=1,ns
                     if (avgmat(i) .gt. pxylem(j)) xtract(i) = xtract(i) &
                + fractn*totrot(j)*rroot(j,i)*(avgmat(i)-pxylem(j))
   35             continue
               else
!                 soil too dry for plants to extact water
                  sumet=0.0
               end if
!
!              sum transpiration from entire transpiring canopy
               trnsp(nplant+1)=trnsp(nplant+1)+fractn*sumet
        ! level 3
        end if
!!!   level2
      else
!********   dead plant material -- compute water content and temperature
            do 45 i=1,nc
               if (canlai(j,i) .le. 0.0) go to 45
!
               pevap(i)=0.0
               rstom(i)=0.0
               b1lv = drycan(j,i)/dt
               a1 = -a1lwr(j,i) - canlai(j,i)*rhoa*ca/rhcan(j,i) &
                    -drycan(j,i)*(cr+wcan(i)*cl)/dt
               b1=lv*b1lv
!              calculate water content based on vapor density of air
                  !if(tlcdt(j,i).le.-273.15) then
                  !   errorflag=1
                  !   return
                  !end if
                     if( abs(tlcdt(j,i)) .gt. 70.0) then
                       print*,"In leaft, 4550",col,row,tlcdt(j,i)
                       stop
                     end if
               call vslope (dummy,satv,tlcdt(j,i))
               humcan=avgvap(i)/satv
               call canhum (2,humcan,dummy,wchum,tcdt(i),canma,canmb)
!
               if (init(j) .eq. 1) then
!                 initialize variables for this time step
                  etcan(j,i)= -b1lv*(wcandt(i)-wcan(i))
                  if (etcan(j,i) .eq. 0.0) then
!xxx                 tlc(j,i) = avgtmp(i)
!                    compute humidity in plant material
                     call canhum (1,humid(i),dummy,wcandt(i),tcdt(i), &
                canma,canmb)
                     if( abs(tlcdt(j,i)) .gt. 70.0) then
                       print*,col,row,tlcdt(j,i)
                       stop
                     end if
                     call vslope (dummy,satv,tlcdt(j,i))
                     etcan(j,i)=canlai(j,i)* &
                (humid(i)*satv-avgvap(i))/rhcan(j,i)
                     wcandt(i) = wcan(i) - etcan(j,i)/b1lv
!                    check if water content is reasonable --
                     if ((wcandt(i)-wchum)*(wcan(i)-wchum).lt.0.0) then
!                       water content went beyond equilibruim with air
                        etcan(j,i)= -b1lv*(wchum-wcan(i))
                        humid(i)=(avgvap(i) + &
            etcan(j,i)*rhcan(j,i)/canlai(j,i))/satv
                        call canhum (2,humid(i),dummy,wcandt(i),tcdt(i), &
                canma,canmb)
                        etcan(j,i)= -b1lv*(wcandt(i)-wcan(i))
                     end if
                  end if
                     tlcdt(j,i)=avgtmp(i)
               end if
!
                  !if(tlcdt(j,i).le.-273.15) then
                  !   errorflag=1
                  !   return
                  !end if
                     if( abs(tlcdt(j,i)) .gt. 120.0) then
                       print*,"In leaft, 4594",col,row,tlcdt(j,i)
                       stop
                     end if
               call vslope (vtslop(i),satv,tlcdt(j,i))
!
!*****         begin iterations to find leaf temp and water content
               iter1 = 0
   40          etcan(j,i)= -b1lv*(wcandt(i)-wcan(i))
!              compute humidity in plant material at end of time step
               call canhum(1,humid(i),hwslop,wcandt(i),tcdt(i), &
            canma,canmb)
!***           set up and solve 2x2 matric for newton-raphson approx.
!              for leaf temp and water content
               a2 = canlai(j,i)*humid(i)*vtslop(i)/rhcan(j,i)
               b2 = canlai(j,i)*satv*hwslop/rhcan(j,i) + b1lv
               ff1 = swcan(j,i) + lwcan(j,i) - lv*etcan(j,i) &
        -a1lwr(j,i)*(tlcdt(j,i)-tlclwr(j,i)) &
        -canlai(j,i)*rhoa*ca*(tlcdt(j,i)-avgtmp(i))/rhcan(j,i) &
        -drycan(j,i)*(cr+wcan(i)*cl)*(tlcdt(j,i)-tlc(j,i))/dt
               ff2=canlai(j,i)*(humid(i)*satv-avgvap(i))/rhcan(j,i) &
        - etcan(j,i)
               det = a1*b2 - a2*b1
               d1 = (ff1*b2 - ff2*b1)/det
               d2 = (a1*ff2 - a2*ff1)/det
!
!***           update values
                     tlcdt(j,i)=avgtmp(i)
               wcandt(i) = wcandt(i) - d2
!              check if water content is reasonable --
                  !if(tlcdt(j,i).le.-273.15) then
                  !   errorflag=1
                  !   return
                  !end if
                     if( abs(tlcdt(j,i)) .gt. 120.0) then
                       print*,"In leaft, 4628",col,row,tlcdt(j,i)
                       stop
                     end if
               call vslope (vtslop(i),satv,tlcdt(j,i))
               humcan=avgvap(i)/satv
               call canhum (2,humcan,dummy,wchum,tcdt(i),canma,canmb)
               delmax=wcan(i)-wchum
               delta=wcandt(i)-wchum
!cccc          if(delta*delmax.lt.0.0 .or.abs(delmax).lt.abs(delta))then
               if(delta*delmax.lt.0.0)then
!                 water content went beyond equilibruim with humidity
                  etcan(j,i)= -b1lv*(wchum-wcan(i))
                  call canhum (1,hum,dummy,wcan(i),tcdt(i),canma,canmb)
                  etmax=canlai(j,i)*(hum*satv-avgvap(i))/rhcan(j,i)
                  if (abs(etcan(j,i)) .gt. abs(etmax)) etcan(j,i)=etmax
                  humid(i)=(avgvap(i) + &
            etcan(j,i)*rhcan(j,i)/canlai(j,i))/satv
                  call canhum (2,humid(i),dummy,wcandt(i),tcdt(i), &
            canma,canmb)
               end if
               if (abs(d1) .gt. 0.01) then
                  iter1 = iter1 +1
                  if (iter1 .lt. 10) go to 40
               end if
!              calculate evaporation from this layer
               etcan(j,i)= -b1lv*(wcandt(i)-wcan(i))
               sumet=sumet+etcan(j,i)
   45       continue
            init(j) = 2
      end if
!
!********store evaporation/transpiration from plant species
   48    trnsp(j)=fractn*sumet+sumev
!
!        sum heat and water transfer in each layer from all canopy types
      do 50 i=1,nc
            if (canlai(j,i) .le. 0.0) go to 50
            heatc(i)=heatc(i) &
        + canlai(j,i)*rhoa*ca*(tlcdt(j,i)-avgtmp(i))/rhcan(j,i)
            etlyr(i)=etlyr(i) + fractn*etcan(j,i) + pevap(i)
            dheatc(i)=dheatc(i) &
        + canlai(j,i)*rhoa*ca/rhcan(j,i)
            detlyr(i)=detlyr(i) + canlai(j,i)/(rstom(i)+rhcan(j,i))
            dtldtc(i)=dtldtc(i) + lv*humid(i)*vtslop(i)* &
        canlai(j,i)/(rstom(i)+rhcan(j,i))
      50 continue
    60 continue
!
!xxxx
      do i=1,nc
         tleaf(nplant+1,i)=0.0
         do j=1,nplant
            tleaf(j,i)=tlcdt(j,i)
            if (canlai(j,i) .le. 0.0) tleaf(j,i)=0.0
            tleaf(nplant+1,i)=tleaf(nplant+1,i) &
        + tlcdt(j,i)*difkl(j,i)/difkl(nplant+1,i)
         end do
      end do
!
!     calculate change of leaf temp with respect to canopy temp
      do 70 i=1,nc
         dtldtc(i)=dheatc(i)/(dheatc(i)+dtldtc(i))
   70 continue
      return
end subroutine leaft2


! no problem checked changed
!***********************************************************************
!
 SUBROUTINE WEIGHT (N,AVG,BEGIN,END,col,row)
!
!     THIS SUBROUTINE CALCULATES THE WEIGHTED AVERAGE OF VARIABLES AT
!     THE BEGINNING AND END OF THE TIME STEP
!
!       WT = WEIGHTING FOR VALUES AT THE BEGINNING OF THE TIME STEP
!        WDT = WEIGHTING FOR VALUES AT THE END OF THE TIME STEP
!
!***********************************************************************
    use controlpara_mod,only: wt2d,wdt2d
    implicit none
    integer(i4),intent(in)::N
    integer(i4),intent(in)::col,row
    integer(i4)::I
    real(r8),dimension(N),intent(in):: begin,end
    real(r8),dimension(N),intent(inout)::avg
    real(r8),pointer::WT,WDT
    WT=>WT2d(col,row)
    WDT=>WDT2d(col,row)


      DO 10 I=1,N
         AVG(I) = WT*BEGIN(I) + WDT*END(I)
   10 CONTINUE
      RETURN
 END SUBROUTINE WEIGHT




!no problem changed
!***********************************************************************
!
 SUBROUTINE CANTK (NC,CON,TC,ZC,col,row)
!
!     THIS SUBROUTINE COMPUTES THE EDDY CONDUCTANCE COEFFICIENT THROUGH
!     THE CANOPY
!
!***********************************************************************
    use constvar_mod,only:CA,RHOA,TKA,VONKRM
    use windv_mod,only:STABLE2d,WINDC2d,ZERSUB2d,ZERO2d,ZHSUB2d,ZMSUB2d
    implicit none
!input
    integer(i4),intent(in)::NC,col,row
    real(r8),dimension(NCMAX-1),intent(inout)::CON
    real(r8),dimension(NCMAX),intent(in)::TC,ZC
!temp
    real(r8)::RICHRD,dudz,TKCAN,TMPAIR,TMPSFC,HEIGHT,HFLUX,PHIH,PHIW,RH,stabl
    integer(i4)::I

    real(r8),pointer::STABLE,ZERSUB,ZERO,ZHSUB,ZMSUB
    real(r8),dimension(:),pointer::WINDC

    STABLE=>STABLE2d(col,row)
    WINDC =>WINDC2d(col,row,:)
    ZERSUB=>ZERSUB2d(col,row)
    ZERO  =>ZERO2d(col,row)
    ZHSUB =>ZHSUB2d(col,row)
    ZMSUB =>ZMSUB2d(col,row)


      RICHRD=STABLE
      IF (RICHRD.GT. 1.0) RICHRD= 1.0
      IF (RICHRD.LT. -2.) RICHRD= -2.

      stabl = richrd

      DO 10 I=1,NC-1
!
!XX      AVGZ=ZC(NC+1)-(ZC(I+1)+ZC(I))/2.
!XX      IF (AVGZ .GT. ZERO) THEN
!           USE LOGARITHRMIC WIND PROFILE CONDUCTANCE FROM ATSTAB
!XX         CON(I)=CONRH/(ZC(I+1)-ZC(I))
!XX        ELSE
            dudz = (windc(i)-windc(i+1))/(ZC(I+1)-ZC(I))
            IF (STABLE .GE. 0) THEN
!              STABLE CONDITIONS
               PHIW = 1. + 5.0*STABL
               PHIH = 1. + 5.0*STABL
              ELSE
!              UNSTABLE CONDITIONS
               PHIW = 1./(1 - 16.*STABL)**0.25
               if(1 - 16.*STABL.lt.0.0) then
                 print*,"In CANTK, the arg for sqrt is le 0.0"
                 print*,col,row
                 stop
               end if
               PHIH = 1./SQRT(1 - 16.*STABL)
            END IF
            TKCAN=RHOA*CA*dudz*(VONKRM*(ZC(NC+1)-ZERO))**2/(PHIW*PHIH)
            IF (TKCAN .LT. TKA) TKCAN=TKA  
            CON(I) = TKCAN/(ZC(I+1)-ZC(I))
!XX      END IF
   10 CONTINUE
!
!     Get conductance for surface
      TMPAIR = TC(NC)
      TMPSFC = TC(NC+1)
      HEIGHT = MAX(ZC(NC+1)-ZC(NC),ZERSUB+0.01)
      
!
      CALL STAB (TMPAIR,TMPSFC,WINDC(NC),HEIGHT,HFLUX,RH, &
        ZMSUB,ZHSUB,ZERSUB)
      CON(NC) = RHOA*CA/RH

      RETURN

 END SUBROUTINE CANTK



! no problem changed
!***********************************************************************
      SUBROUTINE STAB(TMPAIR,TMPSFC,WIND,HEIGHT,HFLUX,RH,ZM,ZH,ZERO)
!
!***********************************************************************
    use constvar_mod
    implicit none
!input
    real(r8),intent(in)::TMPAIR,TMPSFC,WIND,HEIGHT,ZM,ZH,ZERO
    real(r8),intent(inout)::HFLUX,RH
!temp
    real(r8)::HTOLER,HFLUX1,PSIM,PSIH,ZMLOG,ZHLOG,USTAR,STABLE,HFLUX2,ERROR
    integer(i4)::ITER1

!****    DEFINE INITIAL ASSUMPTIONS AND CONSTANTS FOR CURRENT TIME STEP
      HTOLER = 0.001
      HFLUX=0.0
      HFLUX1=0.0
      PSIM=0.0
      PSIH=0.0
      ZMLOG=LOG((HEIGHT+ZM-ZERO)/ZM)
      ZHLOG=LOG((HEIGHT+ZH-ZERO)/ZH)
      USTAR=WIND*VONKRM/(ZMLOG + PSIM)
!
      ITER1=0
!
!**** START ITERATIVE PROCESS TO OBTAIN HEAT FLUX
   10 ITER1=ITER1+1
      STABLE=VONKRM*(HEIGHT-ZERO)*G*HFLUX1/(RHOA*CA*(TMPAIR+273.16)*USTAR**3)
!
      IF (STABLE .GE. 0.0) THEN
!****    ATMOSPHERE IS STABLE
!        IF STABILITY IS GREATER THAN ZMLOG/9.4, COMPUTED HFLUX WILL
!        ACTUALLY DECREASE WITH INCREASING TEMPERATURE DIFFERENTIAL,
!        CAUSING NUMERICAL INSTABILITY AND CONVERGENCE PROBLEMS -- AVOID
!        THIS SITUATION BY LIMITING STABLE TO ZMLOG/9.4
         IF (STABLE .GT. ZMLOG/9.4) STABLE=ZMLOG/9.4
         PSIH=4.7*STABLE
         PSIM=PSIH
       ELSE
!****    ATMOSPHERE IS UNSTABLE
         PSIH=-2.*LOG((1 + (1-16.*STABLE)**0.5)/2.0)
         PSIM=0.6*PSIH
!        IF ATMOSPHERE IS SUFFICIENTLY UNSTABLE, EQUATIONS MAY RESULT
!        IN FRICTION VELOCITY LESS THEN ZERO !!??!!?? (LIMIT PSIM)
         IF (PSIM/ZMLOG .LT. -0.50) THEN
            PSIM=-0.50*ZMLOG
            PSIH=PSIM/0.6
            STABLE = -((2.*EXP(-PSIH/2.) - 1)**2 -1)/16.
         END IF
      END IF
      USTAR=WIND*VONKRM/(ZMLOG + PSIM)
      RH=(ZHLOG + PSIH)/(USTAR*VONKRM)
      HFLUX2=RHOA*CA*(TMPAIR-TMPSFC)/RH
      ERROR=ABS(HFLUX2-HFLUX1)
      HFLUX1=HFLUX2
      IF(ABS(HFLUX1).LT.ABS(TKA*(TMPAIR-TMPSFC)/(HEIGHT-ZERO))) GO TO 15
!        WIND IS SUFFICIENTLY SMALL, OR ATMOSPHERE IS SUFFICIENTLY
!        STABLE THAT TURBULENCE DOES NOT INCREASE FLUX -- NO NEED TO
!        ITERATE FURTHER
!
!**** CHECK IF TOLERANCE HAS BEEN MET
      IF (ITER1 .GT. 40) THEN
!        CONVERGENCE NOT MET, BUT PROBABLY CLOSE ENOUGH
         GO TO 15
      END IF
      IF (ERROR .GT. HTOLER) GO TO 10
!
   15 HFLUX = HFLUX1
      RETURN
 END SUBROUTINE STAB



!abs no problem changed
!***********************************************************************
 SUBROUTINE ATSTAB (NPLANT,NC,NSP,NR,TA,TADT,T,TDT,VAPA,VAPADT, &
        VAP,VAPDT,WIND,HEIGHT,HTOLER,HFLUX,VFLUX,ZC,ZR,RHOR, &
        ICE,ITER,PRESUR,col,row,soilbeta)
!     THIS SUBROUTINE CALCULATES THE ATMOSPHERIC STABILITY, THE TRANSFER
!     COEFFICIENTS FOR THE BOTH HEAT AND VAPOR FROM THE SURFACE,
!     DEFINES THE BOUNDARY CONDITIONS FOR THE MATRIX SOLUTIONS, AND
!     CALCULATES THE WINDSPEED PROFILES IN THE RESIDUE AND CANOPY
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d
    use constvar_mod
    use windv_mod,only:WINDC2d,WINDR2d,ZERO2d,ZERSUB2d,ZH2d,ZM2d,ZMSUB2d,CONRH2d,CONRV2d,STABLE2d,USTAR2d
    use writeit_mod,only:hflux12d,rh2d
    use matrix_mod,only:B12d,D12d,B22d,D22d
    use clayrs_mod,only:TOTLAI2d,CANLAI2d
    use savedata_mod,only:TMPAIR2d,VAPAIR2d,ZMLOG2d,ZHLOG2d,PSIM2d,PSIH2d
    use statevar_mod,only:LANDUSE2d
    implicit none
!input
    integer(i4),intent(in)::NPLANT,NC,NSP,NR,ITER,col,row
    real(r8),dimension(NCMAX),intent(in)::ZC
    real(r8),dimension(NRMAX),intent(in)::ZR,RHOR
    real(r8),intent(in)::TA,TADT,T,TDT,VAPA,VAPADT,VAP,VAPDT,WIND,HEIGHT,HTOLER,PRESUR
    real(r8),intent(inout)::HFLUX,VFLUX
    integer(i4),intent(in)::ICE
    real(r8),intent(in)::soilbeta
!temp
    real(r8)::AWIND,CON,CONV,DUMMY,DV,ERROR,hflux2,HMIN,HNODE1,HNODE2,RV,S,SUML,SUMLAI,TMPSFC,VAPSFC, &
            VMIN,WINDEXP,WINDLOG,WINDRH,WSTAR,ZMCAN1,ZMCAN2
    integer(i4)::I,ITER1,J

    real(r8),pointer::TMPAIR,VAPAIR,ZMLOG,ZHLOG,PSIM,PSIH
    real(r8),dimension(:),pointer::B1,D1,B2,D2
    real(r8),dimension(:),pointer::WINDC
    real(r8),dimension(:),pointer::WINDR
    real(r8),dimension(:),pointer::TOTLAI
    real(r8),dimension(:,:),pointer::CANLAI
    real(r8),pointer::CONRH,CONRV,hflux1,rh,STABLE,USTAR,ZERO,ZERSUB,ZH,ZM,ZMSUB
    real(r8),pointer::WT,WDT



    TMPAIR =>TMPAIR2d(col,row)
    VAPAIR =>VAPAIR2d(col,row)
    ZMLOG  =>ZMLOG2d(col,row)
    ZHLOG  =>ZHLOG2d(col,row)
    PSIM   =>PSIM2d(col,row)
    PSIH   =>PSIH2d(col,row)

    B1     =>B12d(col,row,:)
    D1     =>D12d(col,row,:)
    B2     =>B22d(col,row,:)
    D2     =>D22d(col,row,:)
    WINDC  =>WINDC2d(col,row,:)
    WINDR  =>WINDR2d(col,row,:)

    CONRH  =>CONRH2d(col,row)
    CONRV  =>CONRV2d(col,row)
    hflux1 =>hflux12d(col,row)
    rh     =>rh2d(col,row)
    STABLE =>STABLE2d(col,row)
    USTAR  =>USTAR2d(col,row)
    ZERO   =>ZERO2d(col,row)
    ZERSUB =>ZERSUB2d(col,row)
    ZH     =>ZH2d(col,row)
    ZM     =>ZM2d(col,row)
    ZMSUB  =>ZMSUB2d(col,row)
    TOTLAI =>TOTLAI2d(col,row,:)
    CANLAI =>CANLAI2d(col,row,:,:)

    WT     =>WT2d(col,row)
    WDT    =>WDT2d(col,row)

    ZERO   =MIN(ZERO,ZC(NC+1)-0.01)
!
      IF (ITER .LE. 2) THEN
         IF (ITER .LE. 1) THEN
            TMPAIR=WT*TA + WDT*TADT
            VAPAIR=WT*VAPA + WDT*VAPADT
            ZMLOG=LOG((HEIGHT+ZM-ZERO)/ZM)
            ZHLOG=LOG((HEIGHT+ZH-ZERO)/ZH)
         END IF
!
!****    DEFINE INITIAL ASSUMPTIONS FOR CURRENT TIME STEP
         HFLUX1=0.0
         PSIM=0.0
         PSIH=0.0
         USTAR=WIND*VONKRM/(ZMLOG + PSIM)
      END IF
!
      TMPSFC=WT*T + WDT*TDT
      VAPSFC=WT*VAP + WDT*VAPDT
      HFLUX=0.0
      VFLUX=0.0
      IF (WIND .LE. 0.0) GO TO 20
!
      ITER1=0
!
!**** START ITERATIVE PROCESS TO OBTAIN HEAT FLUX
   10 ITER1=ITER1+1
      STABLE=VONKRM*(HEIGHT-ZERO)*G*HFLUX1 &
        /(RHOA*CA*(TMPAIR+273.16)*USTAR**3)
!
      IF (STABLE .GE. 0.0) THEN
!****    ATMOSPHERE IS STABLE
!        IF STABILITY IS GREATER THAN ZMLOG/9.4, COMPUTED HFLUX WILL
!        ACTUALLY DECREASE WITH INCREASING TEMPERATURE DIFFERENTIAL,
!        CAUSING NUMERICAL INSTABILITY AND CONVERGENCE PROBLEMS -- AVOID
!        THIS SITUATION BY LIMITING STABLE TO ZMLOG/9.4
         IF (STABLE .GT. ZMLOG/9.4) STABLE=ZMLOG/9.4
         PSIH=4.7*STABLE
         PSIM=PSIH
       ELSE
!****    ATMOSPHERE IS UNSTABLE
         PSIH=-2.*LOG((1 + (1-16.*STABLE)**0.5)/2.0)
         PSIM=0.6*PSIH
!        IF ATMOSPHERE IS SUFFICIENTLY UNSTABLE, EQUATIONS MAY RESULT
!        IN FRICTION VELOCITY LESS THEN ZERO !!??!!?? (LIMIT PSIM)
         IF (PSIM/ZMLOG .LT. -0.50) THEN
            PSIM=-0.50*ZMLOG
            PSIH=PSIM/0.6
            STABLE = -((2.*EXP(-PSIH/2.) - 1)**2 -1)/16.
         END IF
      END IF
      USTAR=WIND*VONKRM/(ZMLOG + PSIM)
      if(USTAR.lt.0.0) then
         print*,USTAR,wind,VONKRM,ZMLOG,PSIM,HEIGHT,ZM,ZERO
         stop
      end if
      RH=(ZHLOG + PSIH)/(USTAR*VONKRM)
      if (nc .gt. 0) then
!        RH IS TO ZERO PLANE DISPLACEMENT - ADJUST TO TOP OF CANOPY
         if(ZC(NC+1)<ZERO) then
           print*,col,row,ZC(NC+1),ZERO
           stop
         end if
         RH=RH*(1.-((LOG((ZC(NC+1)-ZERO+ZH)/ZH)+PSIH)/(ZHLOG + PSIH)))
      end if
      HFLUX2=RHOA*CA*(TMPAIR-TMPSFC)/RH
      ERROR=ABS(HFLUX2-HFLUX1)
      HFLUX1=HFLUX2
      IF(ABS(HFLUX1).LT.ABS(TKA*(TMPAIR-TMPSFC)/(HEIGHT-ZC(NC+1)))) &
        GO TO 15
!        WIND IS SUFFICIENTLY SMALL, OR ATMOSPHERE IS SUFFICIENTLY
!        STABLE THAT TURBULENCE DOES NOT INCREASE FLUX -- NO NEED TO
!        ITERATE FURTHER
!
!**** CHECK IF TOLERANCE HAS BEEN MET
      IF (ITER1 .GT. 40) THEN
!        CONVERGENCE NOT MET, BUT PROBABLY CLOSE ENOUGH
         GO TO 15
      END IF
      IF (ERROR .GT. HTOLER) GO TO 10
!
!**** HEAT FLUX IS WITHIN TOLERABLE ERROR -- NOW CALCULATE VAPOR FLUX
   15 HFLUX=HFLUX1
      RV=RH
      VFLUX=soilbeta*(VAPAIR-VAPSFC)/RV
!
!**** COMPARE FLUX WITH MINIMUM FLUX AS CALCULATED FROM THERMAL
!     CONDUCTIVITY AND VAPOR DIFFUSIVITY OF STILL AIR -- HEAT FLUX FIRST
   20 if (nc .gt. 0) then
         HMIN=TKA*(TMPAIR-TMPSFC)/(HEIGHT-ZC(NC+1))
        else
         HMIN=TKA*(TMPAIR-TMPSFC)/(HEIGHT-ZERO)
      end if
      IF (ABS(HMIN) .LE. ABS(HFLUX)) THEN
!        WIND IS SUFFICIENT THAT HEAT FLUX IS ENHANCED
         CON=RHOA*CA/RH
         if (nc .gt. 0) CONRH=RHOA*CA/RH*(HEIGHT-ZC(NC+1))
        ELSE
!        WIND IS SUFFICIENTLY LOW TO BE CONSIDERED AS "STILL AIR"
         HFLUX=HMIN
         if (nc .gt. 0) then
            CON=TKA/(HEIGHT-ZC(NC+1))
          else
            CON=TKA/(HEIGHT-ZERO)
         end if
         CONRH=TKA
      END IF
!
!**** NOW COMPARE VAPOR FLUXES
      DV=VDIFF*(((TMPAIR+273.16)/273.16)**2)*(P0/PRESUR)
      if (nc .eq. 0) then
         VMIN=DV*(VAPAIR-VAPSFC)/(HEIGHT-ZERO)
        else
         VMIN=DV*(VAPAIR-VAPSFC)/(HEIGHT-ZC(NC+1))
      end if
      IF (ABS(VMIN) .LE. ABS(VFLUX)) THEN
!        WIND IS SUFFICIENT THAT VAPOR FLUX IS ENHANCED
         CONV=1./RV
         CONRV=1./RV
       ELSE
!        WIND IS SUFFICIENTLY LOW TO BE CONSIDERED AS "STILL AIR"
         VFLUX=VMIN
         CONV=DV
       END IF
!
!**** DEFINE MATRIX COEFFICIENTS FOR SURFACE MATERIAL PRESENT
!
      IF (NC .GT. 0) THEN
!****    SURFACE MATERIAL IS CANOPY
         B1(1)=B1(1) - WDT*CON
         D1(1)=HFLUX
         B2(1)=-WDT*CONV
         D2(1)=VFLUX
         WINDC(1)=max(USTAR*LOG((ZC(NC+1)+ZM-ZERO)/ZM)/VONKRM,0.01)
         if(windc(1).le.0.0) then
           print*,"the wind speed in the canopy is wrong"
           print*,"USTAR,ZM,ZC(NC+1)+ZM-ZERO)/ZM,VONKRM"
           print*,USTAR,ZM,ZC(NC+1)+ZM-ZERO,VONKRM
           print*,WIND
           stop
         end if
!        CALCULATE WINDSPEED AT THE CANOPY NODES ASSUMING AN
!        EXPONENTIAL DECREASE IN WINDSPEED WITH DEPTH
!        AWIND = 0  FOR SPARSE CANOPY: AWIND >= 4 FOR DENSE CANOPY
         SUMLAI=0.0
         DO 40 J=1,NPLANT
            SUMLAI = SUMLAI + TOTLAI(J)
   40    CONTINUE
!        USE EXTINCTION COEFF FROM NIKOLOV & ZELLER (2003) ENVIRON.POLL.
         AWIND = 2.879*(1-EXP(-SUMLAI))
         DO 45 I=2,NC+1
!           COMPUTE WIND BASED ON EXPONENTIAL DECREASE WITH LAI
            SUML = 0.0
            DO 42 J=1,NPLANT
               SUML = SUML+CANLAI(J,I-1)
   42       CONTINUE
!           u(L) = uh*EXP(-AWIND*L/Ltot)  L = Total LAI above layer
!           ==>  u(L+dL) = u(L)*exp(AWIND*dL/Ltot)
            WINDEXP = WINDC(I-1)*EXP(-AWIND*SUML/SUMLAI)
            IF (I .LT. NC+1) THEN
!             COMPUTE WIND BASED ON LOGARITHMIC WIND PROFILE
              HNODE1 = ZC(NC+1) - ZC(I-1)
              HNODE2 = ZC(NC+1) - ZC(I)
              if(ZERSUB .ge. HNODE1) HNODE1=ZERSUB+0.01
              if(ZERSUB .ge. HNODE2) HNODE2=ZERSUB+0.01
              ZMCAN1 = LOG((HNODE1+ZMSUB-ZERSUB)/ZMSUB)
              ZMCAN2 = LOG((HNODE2+ZMSUB-ZERSUB)/ZMSUB)
              WSTAR = WINDC(I-1)*VONKRM/ZMCAN1
              WINDLOG = WSTAR*ZMCAN2/VONKRM
!             WIND AT NODE IS MINIMUM OF THE TWO APPROACHES
!             (THIS RESOLVES PROBLEM WITH LITTLE OR NO LAI IN A LAYER)
              WINDC(I) = max(0.01,MIN(WINDLOG,WINDEXP))
              if(windc(I).le.0.0) then
                print*,"the wind speed in the canopy is wrong"
                print*,"WINDLOG,WINDEXP in line 5022",ZMCAN2,WSTAR,WINDLOG,WINDEXP,ZERSUB,HNODE1,HNODE2
                print*,WIND
                stop
              end if
             ELSE
!               WIND AT TOP OF RESIDUE CANNOT BE ZERO
                WINDC(I) = max(0.01,WINDEXP)
              if(windc(I).le.0.0) then
                print*,"the wind speed in the canopy is wrong"
                print*,"WINDLOG,WINDEXPin line 5031",WINDLOG,WINDEXP
                print*,WIND
                stop
              end if
            END IF
   45    CONTINUE
         IF (NR.GT.0 .AND. NSP.EQ.0) THEN
!           RESIDUE LIES BENEATH THE CANOPY -- CALCULATE WINDRH
            WINDRH = WINDC(NC+1)
           ELSE
!           NO RESIDUE LIES BENEATH CANOPY -- PERHAPS IT IS SNOW-COVERED
            WINDRH=0.0
         END IF
      ELSE
!
      IF (NSP .GT. 0) THEN
!****    SURFACE MATERIAL IS SNOW
         if(abs(TDT) .gt. 120.0)then
            print*,"In ATSTAB,TDT",col,row,TDT
            stop
         endif
         CALL VSLOPE (S,DUMMY,TDT)
         IF (ICE .EQ. 0) THEN
!           SNOW IS NOT MELTING - ENERGY BALANCE BASED ON TEMPERATURE
            B1(1) = B1(1) - WDT*CON - WDT*LS*S*CONV
          ELSE
!           SNOW IS MELTING - ENERGY BALANCE BASED ON LIQUID WATER
            B1(1) = B1(1)
         END IF
         D1(1)=HFLUX + LS*VFLUX
         B2(1)=0.0
         D2(1)=0.0
!        IF RESIDUE LIES BENEATH THE SNOW, WINDRH = 0.0
         WINDRH=0.0
      ELSE
!
      IF (NR .GT. 0) THEN
!****    SURFACE MATERIAL IS RESIDUE
         B1(1)=B1(1) - WDT*CON
         D1(1)=HFLUX
         B2(1)=-WDT*CONV
         D2(1)=VFLUX
!        CALCULATE WINDSPEED AT THE TOP OF THE RESIDUE LAYER
         WINDRH=USTAR*LOG((ZR(NR+1)+ZM-ZERO)/ZM)/VONKRM
      ELSE
!
!**** SURFACE MATERIAL IS BARE SOIL
         if(abs(TDT) .gt. 120.0)then
            print*,"In ATSTAB,TDT 2",col,row,TDT
            stop
         endif
         CALL VSLOPE (S,DUMMY,TDT)
         D1(1)=HFLUX+LV*VFLUX
         D2(1)=VFLUX/RHOL
!        CHECK IF ICE IS PRESENT -- IF SO, THE WATER BALANCE IS WRITTEN
!        IN TERMS OF ICE CONTENT.  THE ENERGY BALANCE IS WRITTEN WITH
!        LIQUID WATER AS THE REFERENCE STATE, THEREFORE LATENT HEAT OF
!        SUBLIMATION IS NEVER USED -- ONLY HEAT OF VAPORIZATION.
         IF (ICE .EQ. 0) THEN
!           NO ICE IS PRESENT AT SOIL SURFACE
            B1(1) = B1(1) - WDT*CON - WDT*LV*S*CONV
            B2(1) = -WDT*VAPDT*CONV*0.018*G/(UGAS*(TDT+273.16))/RHOL
           ELSE
!           ICE IS PRESENT AT SOIL SURFACE
            B1(1) = B1(1)
            B2(1) = 0.0
         END IF
!
      END IF
      END IF
      END IF
!
      IF (NR .GT. 0) THEN
!        CALCULATE THE WINDSPEED AT THE MID-POINT OF EACH RESIDUE
!        LAYER ASSUMING AN EXPONENTIAL DECREASE WITH DEPTH
!        AWIND = 0  FOR SPARSE RESIDUE: AWIND >= 4 FOR DENSE RESIDUE
         AWIND = 0.1*(RHOR(1) - 20.)
         IF (AWIND .LT. 0.0) AWIND=0.0
         WINDR(1)=WINDRH*EXP(-AWIND*(ZR(1)+ZR(2))/4/ZR(NR+1))
         DO 65 I=2,NR+1
            WINDR(I) = WINDRH*EXP(-AWIND*ZR(I)/ZR(NR+1))
   65    CONTINUE
      END IF
!
      RETURN
!
 END SUBROUTINE ATSTAB


!
!abs no problem changed
!***********************************************************************
!
 SUBROUTINE SOURCE (NC,NSP,NR,NS,NSALT,NPLANT,UC,SC,USP,SSP,UR,SR, &
    US,SS,SINK,SWCAN,SWSNOW,SWRES,SWSOIL,LWCAN,LWSNOW,LWRES,LWSOIL, &
    SOILXT)
!     THIS SUBROUTINE SUMS UP THE SOURCE-SINK TERMS FOR EACH NODE.
!***********************************************************************
!input
    integer(i4),intent(in)::NC,NSP,NR,NS,NSALT,NPLANT
    real(r8),dimension(NCMAX-1),intent(inout)::UC,SC
    real(r8),dimension(NSPMAX),intent(inout)::USP,SSP
    real(r8),dimension(NRMAX),intent(inout)::UR,SR,SWRES,LWRES
    real(r8),dimension(NSMAX),intent(inout)::US,SS,SOILXT
    real(r8),dimension(NSALTMAX,NSMAX),intent(inout)::SINK
    real(r8),dimension(NPMAX+1,NCMAX-1),intent(inout)::SWCAN,LWCAN
    real(r8),dimension(NSPMAX),intent(inout)::SWSNOW
    real(r8),intent(inout)::SWSOIL,LWSNOW,LWSOIL
!temp
    integer(i4)::I,J

!**** CANOPY NODES
      DO 10 I=1,NC
         UC(I)=0.0
         SC(I)=SWCAN(NPLANT+1,I)+LWCAN(NPLANT+1,I)
   10 CONTINUE
!
!**** SNOWPACK NODES
      IF (NSP .GT. 0) THEN
         DO 20 I=1,NSP
            USP(I)=0.0
            SSP(I)=SWSNOW(I)
   20    CONTINUE
         SSP(1)=SSP(1)+LWSNOW
      END IF
!
!**** RESIDUE NODES
      DO 30 I=1,NR
         UR(I)=0.0
         SR(I)=SWRES(I)+LWRES(I)
   30 CONTINUE
!**** SOIL NODES
      DO 50 I=1,NS
         US(I)=-SOILXT(I)
         SS(I)=0.0
         DO 40 J=1,NSALT
            SINK(J,I)=0.0
   40    CONTINUE
   50 CONTINUE
      SS(1)=SS(1)+SWSOIL+LWSOIL
!
      RETURN
 END SUBROUTINE SOURCE


! no problem changed
!***********************************************************************
!
 SUBROUTINE LWRMAT (NC,NSP,NR,TC,TSP,TR,TS,ICESPT,col,row)
!     THIS SUBROUTINE CALCULATES THE CONTRIBUTION TO THE ENERGY BALANCE
!     MATRIX BY LONG-WAVE RADIATION.
!***********************************************************************
    use controlpara_mod,only:WDT2d
    use constvar_mod,only:STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
    use matrix_mod,only:A12d,B12d,C12d
    use radcan_mod,only:TDIFFC2d
    use radres_mod,only:TDIFFU2d
    implicit none

!input
    integer(i4),intent(in)::NC,NSP,NR,COL,ROW
    real(r8),dimension(NCMAX),intent(in)::TC
    real(r8),dimension(NSPMAX),intent(in)::TSP
    real(r8),dimension(NRMAX),intent(in)::TR
    real(r8),dimension(NSMAX),intent(in)::TS
    integer(i4),dimension(NSPMAX),intent(in)::ICESPT
!temp
    integer(i4)::N,MATERL,I
    real(r8)::EMITNC,TSP3,TR3,TS3

    real(r8),dimension(:),pointer::B1,A1,C1
    real(r8),dimension(:),pointer::TDIFFC
    real(r8),dimension(:),pointer::TDIFFU
    real(r8),pointer::WDT

    A1    =>A12d(col,row,:)
    B1    =>B12d(col,row,:)
    C1    =>C12d(col,row,:)
    TDIFFC=>TDIFFC2d(col,row,:)
    TDIFFU=>TDIFFU2d(col,row,:)
    WDT   =>WDT2d(col,row)


      N = 1
      MATERL=2
      B1(N) = 0.0
!
!     RADIATION BALANCE MATRIX COEFFICIENTS FOR CANOPY NODES
      IF (NC .GT. 0) THEN
         EMITNC=(1.-TDIFFC(1))*EMITC
!       TC3=STEFAN*4.*WDT*(TC(1)+273.16)**3
!        B1(N) =-2.*EMITNC*TC3
!        C1(N) = EMITNC
!        A1(N+1) = EMITNC*TC3
         B1(N) =0.0
         C1(N) =0.0
         A1(N+1) =0.0

         N = N + 1
!
         DO 10 I=2,NC
            EMITNC=(1.-TDIFFC(I))*EMITC
!           TC3=STEFAN*4.*WDT*(TC(I)+273.16)**3
!           C1(N-1)=C1(N-1)*EMITNC*TC3
!           A1(N) = A1(N)*EMITNC
!           B1(N) = -2.*EMITNC*TC3
!           C1(N)= EMITNC
!           A1(N+1) = EMITNC*TC3
            C1(N-1)=0.0
            A1(N) = 0.0
            B1(N) = 0.0
            C1(N)= 0.0
            A1(N+1) = 0.0
            N = N + 1
   10    CONTINUE
         MATERL=MATERL + 1
      END IF
!
!     RADIATION BALANCE MATRIX COEFFICIENTS FOR SNOW PACK
      IF (NSP .GT. 0) THEN
         IF (ICESPT(1) .GT. 0) THEN
!           TEMP IS KNOWN - ENERGY BALANCE BASED ON LIQUID CONTENT
            IF (MATERL .GT. 2) THEN
               A1(N) = A1(N)*EMITSP
               C1(N-1)=0.0
            END IF
            B1(N) = 0.0
          ELSE
!           ENERGY BALANCE BASED ON TEMPERATURE
            TSP3=STEFAN*4.*WDT*(TSP(1)+273.16)**3
            IF (MATERL .GT. 2) THEN
!              ADJUST COEFFICIENTS FOR MATERIAL ABOVE SNOWPACK
               A1(N) = A1(N)*EMITSP
               C1(N-1)=C1(N-1)*EMITSP*TSP3
            END IF
            B1(N) = -EMITSP*TSP3
         END IF
         N = N + NSP
         MATERL=MATERL + 1
!        INITIALIZE COEFFICIENTS FOR BOTTOM NODE AND UNDERLYING MATERIAL
         IF (NSP .GT. 1) B1(N-1) = 0.0
         C1(N-1) = 0.0
         DO 5 I=1,NR+1
!           SNOW OVERLYING RESIDUE AND SOIL - INITIALIZE MATRIX COEFF.
            A1(N) = 0.0
            B1(N) = 0.0
            C1(N) = 0.0
            N = N + 1
    5    CONTINUE
         GO TO 50
      END IF
!
!     RADIATION BALANCE MATRIX COEFFICIENTS FOR RESIDUE NODES
!
      IF (NR .GT. 0) THEN
         EMITNC=(1-TDIFFU(1))*EMITR
         TR3=STEFAN*4.*WDT*(TR(1)+273.16)**3
         IF (MATERL .GT. 2) THEN
!           ADJUST COEFFICIENTS FOR MATERIAL ABOVE RESIDUE
            A1(N) = A1(N)*EMITNC
            C1(N-1)=C1(N-1)*EMITNC*TR3
         END IF
         B1(N) =-2.*EMITNC*TR3
         C1(N) = EMITNC
         A1(N+1) = EMITNC*TR3
         N = N + 1
!
         DO 20 I=2,NR
            EMITNC=(1-TDIFFU(I))*EMITR
            TR3=STEFAN*4.*WDT*(TR(I)+273.16)**3
            C1(N-1)=C1(N-1)*EMITNC*TR3
            A1(N) = A1(N)*EMITNC
            B1(N) = -2.*EMITNC*TR3
            C1(N) = EMITNC
            A1(N+1) =EMITNC*TR3
            N = N + 1
   20    CONTINUE
         MATERL = MATERL + 1
      END IF
!
!     RADIATION BALANCE MATRIX COEFFICIENTS FOR SOIL SURFACE
      TS3=STEFAN*4.*WDT*(TS(1)+273.16)**3
      IF (MATERL .GT. 2) THEN
!        ADJUST COEFFICIENTS FOR MATERIAL ABOVE SOIL
         A1(N) = A1(N)*EMITS
         C1(N-1)=C1(N-1)*EMITS*TS3
      END IF
      B1(N)=-EMITS*TS3
   50 RETURN
 END SUBROUTINE LWRMAT


!no problem changed
!***********************************************************************
 SUBROUTINE LWRBAL (NC,NSP,NR,NPLANT,TA,TADT,TLC,TLCDT,TSP,TSPDT, &
    TR,TRDT,TS,TSDT,VAPA,VAPADT,CLOUDS,LWCAN,LWSNOW,LWRES,LWSOIL,col,row)
!     THIS SUBROUTINE CALLS SUBROUTINES TO SET UP THE LONG-WAVE
!     RADIATION BALANCE FOR SYSTEM.
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d
    use constvar_mod,only:STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
    use radcan_mod,only:tlclwr2d
    use radout_mod,only:dwnlwr2d,uplwr2d
    implicit none
!input
    integer(i4),intent(in)::NC,NSP,NR,NPLANT,col,row
    real(r8),dimension(NPMAX,NCMAX-1),intent(in)::TLC,TLCDT
    real(r8),dimension(NSPMAX),intent(in)::TSP,TSPDT
    real(r8),dimension(NRMAX),intent(in)::TR,TRDT
    real(r8),dimension(NSMAX),intent(in)::TS,TSDT
    real(r8),dimension(NPMAX+1,NCMAX-1),intent(inout)::LWCAN
    real(r8),dimension(NRMAX),intent(inout)::LWRES
    real(r8),intent(in)::TA,TADT,VAPA,VAPADT,CLOUDS
    real(r8),intent(inout)::LWSNOW,LWSOIL
!temp
    real(r8)::ABOVE,TK,EMITSFC,dwnlwr1,dwnlwr2
    integer(i4)::I,J
    real(r8),dimension(:),pointer::dwnlwr,uplwr
    real(r8),dimension(:,:),pointer::tlclwr
    real(r8),pointer::WT,WDT

    dwnlwr=>dwnlwr2d(col,row,:)
    tlclwr=>tlclwr2d(col,row,:,:)
    uplwr =>uplwr2d(col,row,:)
    WT    =>WT2d(col,row)
    WDT   =>WDT2d(col,row)



!**** DETERMINE THE LONGWAVE RADIATION FLUX ABOVE EACH NODE.  START WITH
!     ATMOSPHERE AND WORK DOWNWARDS.
!
      CALL LWRATM (ABOVE,CLOUDS,TA,TADT,VAPA,VAPADT,col,row)
!xxxxx
      dwnlwr(1) = above
!
      IF (NSP .GT. 0) THEN
!        SURFACE IS SNOW
         TK = TSP(1)+273.16
         EMITSFC=EMITSP*STEFAN *(TK**4 + 4.*WDT*(TK**3) &
                *(TSPDT(1)-TSP(1)))
        ELSE
!        SURFACE IS SOIL
         TK=TS(1) + 273.16
         EMITSFC=EMITS*STEFAN*((TK**4)+4.*WDT*(TK**3)*(TSDT(1)-TS(1)))
      END IF
!
      IF (NC .GT. 0) THEN
         DO 25 I=1,NC
           DO 20 J=1,NPLANT
!             save temperatures used to compute canopy lwr balance to be
!             used in LEAFT
              tlclwr(j,i) = tlcdt(j,i)
   20      CONTINUE
   25    CONTINUE
         IF (NSP.EQ.0 .AND. NR.GT.0) THEN
!           RESIDUE LIES BENEATH CANOPY - DETERMINE LONGWAVE RADIATION
!           BALANCE THROUGH CANOPY AND RESIDUE TOGETHER
            CALL LWRCAN (NPLANT,NC,NSP,NR,TLC,TLCDT,TR,TRDT, &
                ABOVE,EMITSFC,LWCAN,LWRES,col,row)
           ELSE
!           CALCULATE RADIATION THROUGH CANOPY DOWN TO SNOW OR SOIL
            CALL LWRCAN (NPLANT,NC,NSP,0,TLC,TLCDT,TR,TRDT, &
                ABOVE,EMITSFC,LWCAN,LWRES,col,row)
         END IF
      END IF
!
      IF (NSP .GT. 0) THEN
!xxxxx
         if (nc .eq. 0) then
            uplwr(1) = emitsfc
            uplwr(2) = emitsfc
            dwnlwr(1) = above
            dwnlwr(2) = above
         end if

!        RADIATION BALANCE FOR SNOW PACK
         LWSNOW = EMITSP*ABOVE - EMITSFC
!        NO LWR BENEATH THE SNOW
         LWSOIL = 0.0
         DO 30 I=1,NR
            LWRES(I) = 0.0
   30    CONTINUE
        ELSE
!
!        RADIATION BALANCE OF RESIDUE AND SOIL
         IF (NR .GT. 0) THEN
!           CHECK IF ALREADY COMPUTED WITH CANOPY RADIATION BALANCE
            IF (NC.EQ.0) THEN
!xxxxx
               dwnlwr1 = above
               dwnlwr2 = above
!              RADIATION BALANCE FOR THE RESIDUE LAYERS WITH NO CANOPY
               CALL LWRCAN (NPLANT,0,NSP,NR,TLC,TLCDT,TR,TRDT, &
                    ABOVE,EMITSFC,LWCAN,LWRES,col,row)
!xxxx
               uplwr(2) = uplwr(1)
               dwnlwr(1) = dwnlwr1
               dwnlwr(2) = dwnlwr2
            END IF
         END IF
!        RADIATION BALANCE OF SOIL
         LWSOIL = EMITS*ABOVE - EMITSFC
      END IF
!
      RETURN
 END SUBROUTINE LWRBAL



!no problem changed
!***********************************************************************
 SUBROUTINE LWRCAN (NPLANT,NC,NSP,NR,TLC,TLCDT,TR,TRDT, &
    ABOVE,EMITSFC,LWCAN,LWRES,col,row)
!     THIS SUBROUTINE CALCULATES THE AMOUNT OF LONGWAVE RADIATION
!     ABSORBED BY EACH CANOPY AND RESIDUE NODE
!
!
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d
    use constvar_mod,only:STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
    use radcan_mod,only:TDIFFC2d,DIFKL2d,FDDU2d
    use radres_mod,only:TDIFFU2d
    use radout_mod,only:uplwr2d,dwnlwr2d
    implicit none
!input
    integer(i4),intent(in)::NPLANT,NC,NSP,NR,col,row
    real(r8),dimension(NPMAX,NCMAX-1),intent(in)::TLC,TLCDT
    real(r8),dimension(NRMAX),intent(in)::TR,TRDT
    real(r8),dimension(NRMAX),intent(inout)::LWRES
    real(r8),dimension(NPMAX+1,NCMAX-1),intent(inout)::LWCAN
    real(r8),intent(inout):: ABOVE
    real(r8),intent(in)::EMITSFC
!temp
    real(r8),dimension(2*(NCMAX+NRMAX))::A,B,C,D
    real(r8),dimension(NCMAX+NRMAX)::RD,RU
    real(r8),dimension(NCMAX-1)::DDU

    real(r8),dimension(:),pointer::uplwr,dwnlwr
    real(r8),dimension(:),pointer::TDIFFC
    real(r8),dimension(:,:),pointer::DIFKL
    real(r8),dimension(:),pointer::TDIFFU
    real(r8),dimension(:),pointer::FDDU
    real(r8),pointer::WT,WDT


    integer(i4)::N,I,J,NCLST,NLAYR
    real(r8)::TK,EMIT,FRACTN,ABSORB

    dwnlwr=>dwnlwr2d(col,row,:)
    uplwr =>uplwr2d(col,row,:)
    TDIFFC=>TDIFFC2d(col,row,:)
    DIFKL =>DIFKL2d(col,row,:,:)
    TDIFFU=>TDIFFU2d(col,row,:)
    FDDU  =>FDDU2d(col,row,:)
    WT    =>WT2d(col,row)
    WDT   =>WDT2d(col,row)

!     INITIALIZE MATRIX ELEMENTS
      N = 1
      A(1) = 0.0
      B(1) = 1.0
      C(1) = 0.0
      D(1) = ABOVE
!
      DO 20 I=1,NC
         N = N + 1
         D(N) = 0.0
         DDU(I) = 0.0
!        COMPUTE RADIATION EMITTED BY EACH PLANT WITH LAYER AND TOTAL
!        EMITTED BY ENTIRE LAYER
         DO 10 J=1,NPLANT
            TK=TLC(J,I)+273.16
            EMIT=EMITC*STEFAN*(TK**4 + 4.*(TK**3)*(TLCDT(J,I)-TLC(J,I)))
!           COMPUTE TOTAL RADIATION EMITTED FROM BOTH SIDES OF LAYER
!           BASED ON FRACTION OF TRANSMISSIVITY OF LAYER FOR PLANT
            FRACTN = (1.-TDIFFC(I))*DIFKL(J,I)/DIFKL(NPLANT+1,I)
            LWCAN(J,I)= -2.*FRACTN*EMIT
            D(N) = D(N) + FRACTN*EMIT
!           WEIGHT THE FRACTION OF UPWARD SCATTERED RADIATION
            DDU(I) = DDU(I) + FDDU(J)*DIFKL(J,I)
   10    CONTINUE
         DDU(I) = DDU(I)/DIFKL(NPLANT+1,I)
!
!        COMPUTE COEFFICIENTS FOR UPWELLING RADIATION ABOVE LAYER
         A(N) = -DDU(I)*(1.-EMITC)*(1.-TDIFFC(I))
         B(N) = 1.0
         C(N) = -(TDIFFC(I) + (1.-DDU(I))*(1.-EMITC)*(1.-TDIFFC(I)))
!        COMPUTE COEFFICIENTS FOR DOWNWELLING RADIATION BELOW LAYER
         N = N + 1
         A(N) = C(N-1)
         B(N) = 1.0
         C(N) = A(N-1)
         D(N) = D(N-1)
   20 CONTINUE
!
      DO 30 I=1,NR
          N = N + 1
         TK=TR(I)+273.16
         EMIT=(1.-TDIFFU(I))*EMITR*STEFAN &
                *(TK**4 + 4.*WDT*(TK**3)*(TRDT(I)-TR(I)))
         LWRES(I)=-2.*EMIT
!        COMPUTE COEFFICIENTS FOR UPWELLING RADIATION ABOVE LAYER
!        HORIZONTAL LEAF ORIENTATION IS ASSUMED FOR SCATTERING.
!        THEREFORE ONLY BACKSCATTERING; NO FORWARD SCATTERING
         A(N) = -(1.-EMITR)*(1.-TDIFFU(I))
         B(N) = 1.0
         C(N) = -TDIFFU(I)
         D(N) = EMIT
!        COMPUTE COEFFICIENTS FOR DOWNWELLING RADIATION BELOW LAYER
         N = N + 1
         A(N) = C(N-1)
         B(N) = 1.0
         C(N) = A(N-1)
         D(N) = D(N-1)
   30 CONTINUE
!
!     COMPUTE MATRIX ELEMENTS FOR LOWER SURFACE BOUNDARY
      N = N + 1
      B(N) = 1.0
      C(N) = 0.0
      IF (NSP .GT. 0) THEN
         A(N) = -(1.-EMITSP)
        ELSE
         A(N) = -(1.-EMITS)
      END IF
      D(N) = EMITSFC
!
!     CALL MATRIX SOLVER
      NLAYR = NC + NR
      CALL SOLVRAD (NLAYR,A,B,C,D,RD,RU)
!
      N = 0
      DO 40 I=1,NC
          N = N + 1
         LWCAN(NPLANT+1,I) = 0.0
!        COMPUTE THE TOTAL RADIATION ABSORBED BY CANOPY LAYER
         ABSORB = (1.-TDIFFC(I))*EMITC*(RD(N)+RU(N+1))
         DO 35 J=1,NPLANT
!           ADD THE FRACTION OF TOTAL RADIATION ABSORBED BY EACH PLANT
!          TO THAT EMITTED BY EACH PLANT
            LWCAN(J,I)=LWCAN(J,I) + ABSORB*DIFKL(J,I)/DIFKL(NPLANT+1,I)
            LWCAN(NPLANT+1,I) = LWCAN(NPLANT+1,I) + LWCAN(J,I)
   35    CONTINUE
   40 CONTINUE
!
      DO 50 I=1,NR
         N = N + 1
         LWRES(I)=LWRES(I) + (1.-TDIFFU(I))*EMITR*(RD(N)+RU(N+1))
   50 CONTINUE
!
      ABOVE = RD(N+1)
!xxxx
      uplwr(1) = ru(1)
!xxx  change to ABOVE to last canopy node to be above grass
      NClst = NC
      if (NClst .eq. 0) NClst = 1
      uplwr(2) = ru(NClst)
      dwnlwr(2) = RD(NClst)
!xxx  uplwr(2) = ru(NC+1)
!xxx  dwnlwr(2) = RD(NC+1)
!
      RETURN
 END SUBROUTINE LWRCAN



!no problem changed
!***********************************************************************
 SUBROUTINE LWRATM (ABOVE,CLOUDS,TA,TADT,VAPA,VAPADT,col,row)
!     THIS SUBROUTINE DETERMINES THE LONG-WAVE RADIATIVE FLUX FROM THE
!     ATMOSPHERE.
!***********************************************************************
    use controlpara_mod,only:WT2d,WDT2d
    use constvar_mod,only:STEFAN,EMATM1,EMATM2,EMITC,EMITR,EMITSP,EMITS
    implicit none
!input
    real(r8),intent(in)::CLOUDS,TA,TADT,VAPA,VAPADT
    real(r8),intent(inout)::ABOVE
    integer(i4),intent(in)::col,row

!temp
    real(r8)::AVGTMP,AVGVAP,CLEAR,BLACK,EMITAT

    real(r8),pointer::WT,WDT
    WT=>WT2d(col,row)
    WDT=>WDT2d(col,row)

!**** DETERMINE THE EMISSIVITY OF CLEAR-SKY ATMOSPHERE BASED ON
!     DILLEY & O'BRIEN (1998), QJR MET SOC.
      AVGTMP = WT*TA + WDT*TADT
      AVGVAP = WT*VAPA + WDT*VAPADT
!     CONVERT TO VAPOR PRESSURE (kPa)
      AVGVAP = 0.4619*AVGVAP*(AVGTMP+273.16)
      CLEAR =59.38+113.7*((AVGTMP+273.16)/273.16)**6+96.96*SQRT((465.*AVGVAP/(AVGTMP+273.16))/2.5)
      BLACK = STEFAN*(AVGTMP+273.16)**4
!     BACK-CALCULATE EMISSIVITY FROM CLEAR-SKY LONG-WAVE RADIATION
      EMITAT=CLEAR/BLACK
!
!**** ADJUST CLEAR-SKY EMISSIVITY FOR FRACION OF CLOUD COVER
!     USE EQUATION BY UNSWORTH & MONTEITH (1975) QJR MET. SOC.
!     (GIVEN IN CAMPBELL, 1985, SOIL PHYSICS WITH BASIC)
      EMITAT=(1.-0.84*CLOUDS)*EMITAT + 0.84*CLOUDS
!
!**** NOW CALCULATE THE LONG-WAVE RADIATIVE FLUX
      ABOVE = EMITAT*BLACK
!
      RETURN
 END SUBROUTINE LWRATM


! no problem changed
!***********************************************************************
 SUBROUTINE UPDATE (NS,NR,NSP,NC,NPLANT,NSALT,ICES,ICESDT,TS,TSDT,&
    MAT,MATDT,CONC,CONCDT,VLC,VLCDT,VIC,VICDT,SALT,SALTDT,&
    TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,TC,TCDT,TLC,TLCDT,VAPC,VAPCDT,&
    WCAN,WCANDT,PCAN,PCANDT,TSP,TSPDT,DLW,DLWDT,ICESP,ICESPT)
!     THIS SUBROUTINE UPDATES THE BEGINNING OF TIME STEP VALUES FOR THE
!     NEW TIME STEP.
!***********************************************************************
    implicit none
!input
    integer(i4),intent(in)::NS,NR,NSP,NC,NPLANT,NSALT
    integer(i4),dimension(NSMAX),intent(inout)::ICES
    integer(i4),dimension(NSMAX),intent(in)::ICESDT
    real(r8),dimension(NSALTMAX,NSMAX),intent(inout)::CONC,SALT
    real(r8),dimension(NSALTMAX,NSMAX),intent(in)::CONCDT,SALTDT
    real(r8),dimension(NSMAX),intent(inout)::TS,MAT,VLC,VIC
    real(r8),dimension(NSMAX),intent(in)::TSDT,MATDT,VLCDT,VICDT
    real(r8),dimension(NRMAX),intent(inout)::TR,VAPR,GMC
    real(r8),dimension(NRMAX),intent(in)::TRDT,VAPRDT,GMCDT
    real(r8),dimension(NCMAX),intent(inout)::TC,VAPC
    real(r8),dimension(NCMAX),intent(in)::TCDT,VAPCDT
    real(r8),dimension(NPMAX,NCMAX-1),intent(inout)::TLC
    real(r8),dimension(NPMAX,NCMAX-1),intent(in)::TLCDT
    real(r8),dimension(NCMAX-1),intent(inout)::WCAN
    real(r8),dimension(NCMAX-1),intent(in)::WCANDT
    real(r8),dimension(NPMAX),intent(inout)::PCAN
    real(r8),dimension(NPMAX),intent(in)::PCANDT
    real(r8),dimension(NSPMAX),intent(inout)::TSP,DLW
    real(r8),dimension(NSPMAX),intent(in)::TSPDT,DLWDT
    integer(i4),dimension(NSPMAX),intent(inout)::ICESP
    integer(i4),dimension(NSPMAX),intent(in)::ICESPT
!temp
    integer(i4)::I,J

!     SOIL PROPERTIES
      DO 15 I=1,NS
         TS(I)=TSDT(I)
         MAT(I)=MATDT(I)
         VLC(I)=VLCDT(I)
         VIC(I)=VICDT(I)
         ICES(I)=ICESDT(I)
         DO 10 J=1,NSALT
            SALT(J,I)=SALTDT(J,I)
            CONC(J,I)=CONCDT(J,I)
   10    CONTINUE
   15 CONTINUE
!
!     RESIDUE PROPERTIES
      DO 20 I=1,NR
         TR(I)=TRDT(I)
         VAPR(I)=VAPRDT(I)
         GMC(I)=GMCDT(I)
   20 CONTINUE
!
!     SNOWPACK PROPERTIES
      DO 30 I=1,NSP
         ICESP(I) = ICESPT(I)
         TSP(I) = TSPDT(I)
         DLW(I) = DLWDT(I)
   30 CONTINUE
!
!     CANOPY PROPERTIES
      DO 40 I=1,NC
         TC(I)=TCDT(I)
         VAPC(I)=VAPCDT(I)
         WCAN(I)=WCANDT(I)
         DO 35 J=1,NPLANT
            TLC(J,I)=TLCDT(J,I)
   35    CONTINUE
   40 CONTINUE
      DO 45 J=1,NPLANT
         PCAN(J)=PCANDT(J)
   45 CONTINUE
!
      RETURN
 END SUBROUTINE UPDATE


! no problem changed
!***********************************************************************
 SUBROUTINE SWRBAL (NPLANT,NC,NSP,NR,XANGLE,CLUMPNG,SWCAN,SWSNOW,&
        SWRES,SWSOIL,SUNHOR,CANALB,ZSP,DZSP,RHOSP,ZR,RHOR,ALBRES,DIRRES,&
        ALBDRY,ALBEXP,VLC,ALATUD,longtitude,SLOPE,ASPECT,HRNOON,HAFDAY,DECLIN,&
        HOUR,NHRPDT,shadow,skyview,shadeef,col,row,julian,year,InDirect,InDiffuse)
!       THIS SUBROUTINE DETERMINES THE SHORT-WAVE RADIATION BALANCE FOR
!       FOR EACH NODE ABOVE THE GROUND SURFACE.
!***********************************************************************
    use radcan_mod,only:TDIRCC2d,TDIFFC2d,DIRKL2d,DIFKL2d,FBDU2d,FDDU2d
    use radres_mod,only:TDIREC2d,TDIFFU2d
    use radout_mod,only:dwnswr2d,upswr2d,dir2d,down2d,up2d,albsrf2d,sw_on_soil2d,sw_on_snow2d
    implicit none
!   input
    integer(i4),intent(in)::NPLANT,NC,NSP,NR,col,row,julian,year
    real(r8),dimension(NPMAX+1,NCMAX-1),intent(inout)::SWCAN
    real(r8),dimension(NSPMAX),intent(inout)::SWSNOW
    real(r8),dimension(NRMAX),intent(inout)::SWRES
    real(r8),dimension(NPMAX),intent(in)::XANGLE,CLUMPNG,CANALB
    real(r8),intent(inout)::SWSOIL,InDirect,InDiffuse
    real(r8),intent(in)::SUNHOR
    real(r8),dimension(NSPMAX),intent(in)::ZSP,DZSP,RHOSP
    real(r8),dimension(NRMAX),intent(in)::ZR,RHOR
    real(r8),intent(in)::ALBRES,DIRRES,ALBDRY,ALBEXP,VLC,ALATUD,SLOPE,ASPECT,HRNOON,HAFDAY,DECLIN
    real(r8),intent(in)::longtitude
    integer(i4),intent(in)::HOUR,NHRPDT,shadeef
    real(r8),intent(in)::shadow,skyview
!temp
    real(r8),dimension(NCMAX)::dirsave,upsave,downsave
    integer(i4)::I,J
    real(r8)::ALBSOI,ALBNXT,ALBSNO,ALTITU,DIRECT,DIFFUS,SUNSLP


    real(r8),dimension(:),pointer::dwnswr,upswr
    real(r8),dimension(:),pointer::dir,down,up
    real(r8),pointer::albsrf,sw_on_snow,sw_on_soil
    real(r8),dimension(:),pointer::TDIRCC,TDIFFC
    real(r8),dimension(:,:),pointer::DIRKL,DIFKL
    real(r8),dimension(:),pointer::FBDU,FDDU
    real(r8),dimension(:),pointer::TDIREC,TDIFFU


    dwnswr=>dwnswr2d(col,row,:)
    upswr=>upswr2d(col,row,:)
    dir=>dir2d(col,row,:)
    down=>down2d(col,row,:)
    up=>up2d(col,row,:)
    albsrf=>albsrf2d(col,row)
    TDIRCC=>TDIRCC2d(col,row,:)
    TDIFFC=>TDIFFC2d(col,row,:)
    DIRKL=>DIRKL2d(col,row,:,:)
    DIFKL=>DIFKL2d(col,row,:,:)
    FBDU=>FBDU2d(col,row,:)
    FDDU=>FDDU2d(col,row,:)
    TDIREC=>TDIREC2d(col,row,:)
    TDIFFU=>TDIFFU2d(col,row,:)
    sw_on_soil=>SW_on_soil2d(col,row)
    sw_on_snow=>SW_on_snow2d(col,row)

      sw_on_snow=0.0
      sw_on_soil=0.0
!**** SEPARATE THE TOTAL SOLAR RADIATION INTO DIRECT AND DIFFUSE
      call SOLAR2(DIRECT,DIFFUS,SUNSLP,ALTITU,SUNHOR,ALATUD,longtitude,SLOPE,&
            ASPECT,HOUR,NHRPDT,shadow,skyview,shadeef,julian,year)
      !CALL SOLAR (DIRECT,DIFFUS,SUNSLP,ALTITU,SUNHOR,ALATUD,SLOPE, &
      !     ASPECT,HRNOON,HAFDAY,DECLIN,HOUR,NHRPDT,shadow,skyview,shadeef,julian)
      InDirect=DIRECT
      InDiffuse=DIFFUS

      IF ((DIRECT+DIFFUS) .LE. 0.0) GO TO 60
!xxxx
      dwnswr(1) = direct+diffus
      dir(1) = direct
      down(1) = diffus
      dirsave(1) = direct
      downsave(1) = diffus
!
!**** DETERMINE THE SOIL AND SNOW ALBEDO
      ALBSOI=ALBDRY*EXP(-ALBEXP*VLC)
      IF (NSP.GT.0) THEN
         IF (NR .GT. 0) THEN
            ALBNXT=ALBRES
           ELSE
            ALBNXT=ALBSOI
         END IF
         CALL SNOALB (NSP,ALBSNO,ZSP,RHOSP,SUNSLP,ALBNXT)
!xxxx
         albsrf=albsno
      END IF
!
!
!**** DETERMINE THE SOLAR RADIATION BALANCE FOR EACH MATERIAL, STARTING
!     FROM THE SURFACE MATERIAL AND WORKING DOWNWARD.
!
!**** SOLAR RADIATION BALANCE OF THE CANOPY
      IF (NC .GT. 0) THEN
         IF (NSP.EQ.0 .AND. NR.GT.0) THEN
!           RESIDUE LIES BENEATH CANOPY - DETERMINE SOLAR RADIATION
!           BALANCE THROUGH CANOPY AND RESIDUE TOGETHER
            CALL SWRCAN (NPLANT,NC,NSP,NR,XANGLE,CLUMPNG,SWCAN,SWRES, &
                DIRECT,DIFFUS,SUNSLP,ALTITU,CANALB,ZR,RHOR,ALBRES,DIRRES, &
                ALBSOI,col,row)
!xxxxx
              dwnswr(2) = dir(1+nc)+ down(1+nc)
              upswr(1)  = up(1)
              upswr(2) = up(1+nc)
              albsrf = albres
           ELSE
!           CALCULATE RADIATION THROUGH CANOPY DOWN TO SNOW OR SOIL
            IF (NSP .GT. 0) THEN
               ALBNXT=ALBSNO
              ELSE
               ALBNXT=ALBSOI
            END IF
            CALL SWRCAN (NPLANT,NC,NSP,0,XANGLE,CLUMPNG,SWCAN,SWRES, &
                DIRECT,DIFFUS,SUNSLP,ALTITU,CANALB,ZR,RHOR,ALBRES,DIRRES, &
                ALBNXT,col,row)
!xxxxx
              albsrf = albnxt
              dwnswr(2) = dir(1+nc)+ down(1+nc)
              upswr(1)  = up(1)
              upswr(2) = up(1+nc)
              do i=1,nc+1
                 dirsave(i)=dir(i)
                 upsave(i)=up(i)
                 downsave(i)=down(i)
              end do
!
         END IF
      END IF
!
!**** SOLAR RADIATION BALANCE OF THE SNOWPACK
      IF (NSP.GT.0) CALL SWRSNO (NSP,SWSNOW,DZSP,RHOSP,DIRECT,DIFFUS, &
            SUNSLP,ALBSNO)
!xxxx if (nsp.gt.0 .and. nr.eq.0 .and. nc.eq.0)
      if (nsp.gt.0 .and.  nc.eq.0) &
            upswr(1)  = dwnswr(1)*albsno

!
!
!**** SOLAR RADIATION BALANCE OF THE RESIDUE IF NO CANOPY PRESENT OR IF
!     RESIDUE IS COVERED BY SNOW
      IF (NR.GT.0) THEN
         IF (NSP.GT.0 .OR. NC.EQ.0) THEN
            CALL SWRCAN (NPLANT,0,NSP,NR,XANGLE,CLUMPNG,SWCAN,SWRES, &
                DIRECT,DIFFUS,SUNSLP,ALTITU,CANALB,ZR,RHOR,ALBRES,DIRRES, &
                ALBSOI,col,row)
!xxxx
            if (nsp .eq. 0) then
               albsrf=albres
               upswr(1) = up(1)
              else
            end if
            do i=1,nc+1
             dir(i)=dirsave(i)
             up(i)=upsave(i)
             down(i)=downsave(i)
            end do
         END IF
      END IF
!
!**** SOLAR RADIATION BALANCE FOR SOIL SURFACE
      SWSOIL=(DIRECT+DIFFUS)*(1.-ALBSOI)
!xxxxx
      if (nc.eq.0 .and.nsp.eq.0 .and.nr.eq.0) then
        albsrf=albsoi
        upswr(1)=dwnswr(1)*albsoi
      end if
!
      RETURN
!
!**** NO SOLAR RADIATION -- SET ABSORBED SOLAR RADIATION TO ZERO
!
!     CANOPY NODES
   60 DO 70 I=1,NC
         DO 65 J=1,NPLANT+1
            SWCAN(J,I)=0.0
   65    CONTINUE
!xxxx
         dir(i)=0.0
         down(i)=0.0
         up(i)=0.0
   70 CONTINUE
!xxxx
      dir(nc+1)=0.0
      down(nc+1)=0.0
      up(nc+1)=0.0
!     DEFINE THE DIFFUSE RADIATION TRANSMISSION FOR LONG-WAVE BALANCE
      IF (NC.GT.0) CALL TRANSC (NPLANT,NC,XANGLE,CLUMPNG,TDIRCC,TDIFFC, &
            DIFKL,DIRKL,FBDU,FDDU,0.0_r8,0.0_r8,col,row)
!
!     SNOWPACK NODES
      DO 80 I=1,NSP
         SWSNOW(I)=0.0
   80 CONTINUE
!
!     RESIDUE NODES
      DO 90 I=1,NR
         SWRES(I)=0.0
   90 CONTINUE
!     DEFINE THE DIFFUSE RADIATION TRANSMISSION FOR LONG-WAVE BALANCE
      IF (NR .GT. 0) CALL TRANSR (NR,TDIREC,TDIFFU,ZR,RHOR,0.785_r8,DIRRES)
!
!     SOIL SURFACE
      SWSOIL=0.0
!xxxx
      dwnswr(1) = 0.0
      dwnswr(2) = 0.0
      upswr(1)  = 0.0
      upswr(2)  = 0.0
!
      RETURN
 END SUBROUTINE SWRBAL




!no problem changed
!***********************************************************************
 SUBROUTINE SWRSNO (NSP,SWSNOW,DZSP,RHOSP,DIRECT,DIFFUS, &
        SUNSLP,ALBSNO)
!     THIS SUBROUTINE COMPUTES THE AMOUNT OF OF SOLAR RADIATION
!     ABSORBED BY EACH LAYER WITHIN THE SNOWPACK
!***********************************************************************

    use constvar_mod,only:RHOL
    use spwatr_mod,only:EXTSP,G1,G2,G3
    implicit none

!input
    integer(i4),intent(in)::NSP
    real(r8),dimension(NSPMAX),intent(inout)::SWSNOW
    real(r8),dimension(NSPMAX),intent(in)::DZSP,RHOSP
    real(r8),intent(inout)::DIRECT,DIFFUS
    real(r8),intent(in)::SUNSLP,ALBSNO
!temp
    real(r8)::SPGRAV,GRAIN,dzv,dzir,albdir,TOTAL,EXC,PERCAB
    integer(i4)::I
    real(r8)::vis

!   ASSUME 58% OF SOLAR RADIATION IS IN VISIBLE SPECTRAL
    vis=0.58
!     COMPUTE THE TOTAL ABSORBED FOR THE ENTIRE SNOWPACK
!xxx  TOTAL = (DIRECT + DIFFUS)*(1. - ALBSNO)
      SPGRAV = RHOSP(1)/RHOL
      GRAIN = G1 + G2*SPGRAV*SPGRAV + G3*SPGRAV**4
!     include influence of sun elevation angle
      if (sunslp .gt. 0.0) then
         dzv = (1.375e-3*sqrt(grain*1000.))*(1.-sin(sunslp))
         dzir= (2.0e-3*sqrt(grain*1000.) + 0.1)*(1.-sin(sunslp))
        else
!        SUN IS NOT ABOVE THE HORIZON
         DZV = 0.0
         DZIR = 0.0
      END IF
      albdir = albsno + dzv*vis + dzir*(1-vis)
      TOTAL = DIRECT*(1.-albdir) + DIFFUS*(1. - ALBSNO)
!
!     CALCULATE EXTINCTION COEFF. AND RADIATION ABSORBED BY EACH LAYER
      DO 10 I= 1,NSP
         SPGRAV = RHOSP(I)/RHOL
         GRAIN = G1 + G2*SPGRAV*SPGRAV + G3*SPGRAV**4
!        EXC = EXTINCTION COEFFICIENT --> CONVERT FROM 1/CM TO 1/M
         EXC = EXTSP*SPGRAV*SQRT(1.0/GRAIN)*100.
         PERCAB = 1.0 - EXP(-EXC*DZSP(I))
         SWSNOW(I) = TOTAL*PERCAB
         TOTAL = TOTAL - SWSNOW(I)
   10 CONTINUE
!
!     DEFINE THE DIFFUSE RADIATION OUT THE BOTTOM OF THE SNOWPACK
!     (THIS WILL BE ABSORBED BY EITHER THE RESIDUE OR THE SOIL.)
      DIFFUS = TOTAL
      DIRECT = 0.0
      RETURN
 END SUBROUTINE SWRSNO


! no problem changed
!***********************************************************************
 SUBROUTINE SWRCAN (NPLANT,NC,NSP,NR,XANGLE,CLUMPNG,SWCAN,SWRES,&
    DIRECT,DIFFUS,SUNSLP,ALTITU,CANALB,ZR,RHOR,ALBRES,DIRRES,ALBNXT,col,row)
!     THIS SUBROUTINE CALCULATES THE AMOUNT OF SOLAR RADIATION ABSORBED
!     BY EACH CANOPY AND RESIDUE NODE
!***********************************************************************
    use radcan_mod,only:FBDU2d,FDDU2d,DIRKL2d,DIFKL2d,TDIRCC2d,TDIFFC2d
    use radres_mod,only:TDIREC2d,TDIFFU2d
    use radout_mod,only:dir2d,down2d,up2d
    implicit none
!input
    integer(i4),intent(in)::NPLANT,NC,NSP,NR,col,row
    real(r8),dimension(NPMAX),intent(in)::XANGLE,CLUMPNG,CANALB
    real(r8),dimension(NPMAX+1,NCMAX-1),intent(inout)::SWCAN
    real(r8),dimension(NRMAX),intent(inout)::SWRES
    real(r8),intent(inout)::DIRECT,DIFFUS
    real(r8),intent(in)::SUNSLP,ALTITU,ALBRES,DIRRES,ALBNXT
    real(r8),dimension(NRMAX),intent(in)::ZR,RHOR
!temp
    real(r8),dimension(NCMAX-1)::DIRALB,DIFALB,DIRTRN,DIFTRN,FBREFU,FDREFU,FBTRNU,FDTRNU,FBREFD,FDREFD,FBTRND,FDTRND
    real(r8),dimension(2*(NRMAX+NCMAX))::A,B,C,D
    real(r8),dimension(NPMAX)::CANTRN
    integer(i4)::I,J,NR1,K,N,NLAYR
    real(r8)::absorb,fract1,fract2,totdif,totdir

    real(r8),dimension(:),pointer::FBDU,FDDU
    real(r8),dimension(:,:),pointer::DIRKL,DIFKL
    real(r8),dimension(:),pointer::TDIREC,TDIFFU
    real(r8),dimension(:),pointer::dir
    real(r8),dimension(:),pointer::TDIRCC,TDIFFC
    real(r8),dimension(:),pointer::down,up


    CANTRN=0.20
    FBDU=>FBDU2d(col,row,:)
    FDDU=>FDDU2d(col,row,:)
    DIRKL=>DIRKL2d(col,row,:,:)
    DIFKL=>DIFKL2d(col,row,:,:)
    TDIREC=>TDIREC2d(col,row,:)
    TDIFFU=>TDIFFU2d(col,row,:)
    dir=>dir2d(col,row,:)
    TDIRCC=>TDIRCC2d(col,row,:)
    TDIFFC=>TDIFFC2d(col,row,:)
    down=>down2d(col,row,:)
    up=>up2d(col,row,:)

!**** OBTAIN TRANSMISSIVITY TO DIRECT AND DIFFUSE RADIATION THROUGH
!     CANOPY AND/OR RESIDUE DEPENDING ON LAYERING OF MATERIALS
      IF (NC .GT. 0) THEN
!        GET TRANSMISSION COEFFICIENTS THROUGH LAYERS
         CALL TRANSC (NPLANT,NC,XANGLE,CLUMPNG,TDIRCC,TDIFFC, &
            DIFKL,DIRKL,FBDU,FDDU,SUNSLP,ALTITU,col,row)
!        CALCULATE WEIGHTED ALBEDO, TRANSMISSION AND SCATTERING FOR
!        DIRECT AND DIFFUSE RADIATION OR EACH OF THE CANOPY LAYERS.
!        THE FRACTION OF DOWNWARD RADIATION TRANSMITTED THROUGH THE
!        LEAVES AND SCATTERED DOWNWARD EQUALS DOWNWARD REFLECTED UPWARD.
!        UPWARD SCATTERED REFLECTED AND TRANSMITTED (FBREFU AND FBTRNU)
!        MUST BE WEIGHTED BASED ON ALBEDO AND TRANSMISSION RESPECTIVELY.
         DO 15 I=1,NC
            DIRALB(I)=0.0
            DIFALB(I)=0.0
            DIRTRN(I)=0.0
            DIFTRN(I)=0.0
            FBREFU(I)=0.0
            FDREFU(I)=0.0
            FBTRNU(I)=0.0
            FDTRNU(I)=0.0
            DO 10 J=1,NPLANT
               DIRALB(I) = DIRALB(I) + CANALB(J)*DIRKL(J,I)
               DIFALB(I) = DIFALB(I) + CANALB(J)*DIFKL(J,I)
               DIRTRN(I) = DIRTRN(I) + CANTRN(J)*DIRKL(J,I)
               DIFTRN(I) = DIFTRN(I) + CANTRN(J)*DIFKL(J,I)
               FBREFU(I) = FBREFU(I) + FBDU(J)*CANALB(J)*DIRKL(J,I)
               FDREFU(I) = FDREFU(I) + FDDU(J)*CANALB(J)*DIFKL(J,I)
!              FOR A SINGLE PLANT SPECIES, FRACTION TRANSMITTED UPWARD
!              IS EQUAL TO 1 - REFLECTED UPWARD (BUT NOT NECESSARILY FOR
!              CANOPY LAYER BECAUSE THEY ARE WEIGHTED BY LAI, K, ETC.)
               FBTRNU(I) = FBTRNU(I) + (1.-FBDU(J))*CANTRN(J)*DIRKL(J,I)
               FDTRNU(I) = FDTRNU(I) + (1.-FDDU(J))*CANTRN(J)*DIFKL(J,I)
   10       CONTINUE
            IF (DIRALB(I) .GT. 0.0) THEN
               FBREFU(I)=FBREFU(I)/DIRALB(I)
               FDREFU(I)=FDREFU(I)/DIFALB(I)
            END IF
            IF (DIRTRN(I) .GT. 0.0) THEN
               FBTRNU(I)=FBTRNU(I)/DIRTRN(I)
               FDTRNU(I)=FDTRNU(I)/DIFTRN(I)
            END IF
!           FRACTION SCATTERED DOWNWARD IS 1 - SCATTERED UPWARD
            FBREFD(I)=1. - FBREFU(I)
            FDREFD(I)=1. - FDREFU(I)
            FBTRND(I)=1. - FBTRNU(I)
            FDTRND(I)=1. - FDTRNU(I)
            DIRALB(I)=DIRALB(I)/DIRKL(NPLANT+1,I)
            DIFALB(I)=DIFALB(I)/DIFKL(NPLANT+1,I)
            DIRTRN(I)=DIRTRN(I)/DIRKL(NPLANT+1,I)
            DIFTRN(I)=DIFTRN(I)/DIFKL(NPLANT+1,I)
   15    CONTINUE
      END IF
!
      IF (NR .GT. 0) &
            CALL TRANSR (NR,TDIREC,TDIFFU,ZR,RHOR,SUNSLP,DIRRES)
      NR1=NR
      IF (NR.EQ.1) THEN
!        RADIATION EXCHANGE IS NOT ACCURATE WITH ONLY ONE RESIDUE LAYER.
!        SPLIT RESIDUE LAYER IN TWO AND REDEFINE TRANSMISSION COEFF.
!        FOR TWO RESIDUE LAYERS.
         NR1=2
         if(TDIREC(1) .lt. 0.0 .or. TDIFFU(1).lt. 0.0)then
           print*,"In SWRCAN, the arg for sqrt is less than 0.0"
           print*,col,row,TDIFFU(1),TDIREC(1)
           stop
         end if
         TDIREC(1)=SQRT(TDIREC(1))
         TDIFFU(1)=SQRT(TDIFFU(1))
         TDIREC(2)=TDIREC(1)
         TDIFFU(2)=TDIFFU(1)
      END IF
!
!**** RADIATION EXCHANGE WITHIN THE CANOPY AND RESIDUE
!     INITIALIZE MATRIX ELEMENTS
      N = 1
      A(1) = 0.0
      B(1) = 1.0
      C(1) = 0.0
      D(1) = DIFFUS
      DIR(1)=DIRECT
!
      DO 20 I=1,NC
         N = N + 1
!
!        COMPUTE COEFFICIENTS FOR UPWELLING RADIATION ABOVE LAYER
         A(N) =-(1.-TDIFFC(I))*(FDREFU(I)*DIFALB(I)+FDTRNU(I)*DIFTRN(I))
         B(N) = 1.0
         C(N) = -TDIFFC(I) - &
            (1.-TDIFFC(I))*(FDREFD(I)*DIFALB(I)+FDTRND(I)*DIFTRN(I))
         D(N) = (1.-TDIRCC(I))*(FBREFU(I)*DIRALB(I)+FBTRNU(I)*DIRTRN(I)) &
            *DIR(I)
!        COMPUTE COEFFICIENTS FOR DOWNWELLING RADIATION BELOW LAYER
         N = N + 1
         A(N) = C(N-1)
         B(N) = 1.0
         C(N) = A(N-1)
         D(N) = (1.-TDIRCC(I))*(FBREFD(I)*DIRALB(I)+FBTRND(I)*DIRTRN(I)) &
            *DIR(I)
!        COMPUTE DIRECT RADIATION PASSING THROUGH THIS LAYER
         DIR(I+1)=DIR(I)*TDIRCC(I)
   20 CONTINUE
!
      DO 30 I=1,NR1
          N = N + 1
!        COMPUTE COEFFICIENTS FOR UPWELLING RADIATION ABOVE LAYER
!        HORIZONTAL LEAF ORIENTATION IS ASSUMED FOR SCATTERING.
!        THEREFORE ONLY BACKSCATTERING; NO FORWARD SCATTERING
         A(N) = -ALBRES*(1.-TDIFFU(I))
         B(N) = 1.0
         C(N) = -TDIFFU(I)
         D(N) = ALBRES*(1.-TDIREC(I))*DIR(NC+I)
!        COMPUTE COEFFICIENTS FOR DOWNWELLING RADIATION BELOW LAYER
         N = N + 1
         A(N) = C(N-1)
         B(N) = 1.0
         C(N) = A(N-1)
         D(N) = 0.0
!        COMPUTE DIRECT RADIATION PASSING THROUGH THIS LAYER
         DIR(NC+I+1)=DIR(NC+I)*TDIREC(I)
   30 CONTINUE
!
!     COMPUTE MATRIX ELEMENTS FOR LOWER SURFACE BOUNDARY
      N = N + 1
      B(N) = 1.0
      C(N) = 0.0
      A(N) = -ALBNXT
      D(N) = ALBNXT*DIR(NC+NR1+1)
!
!     CALL MATRIX SOLVER
      NLAYR = NC + NR1
      CALL SOLVRAD (NLAYR,A,B,C,D,DOWN,UP)
!
!**** DETERMINE THE AMOUNT OF RADIATION ABSORBED AT EACH CANOPY NODE BY
!     WEIGHTING ACCORDING TO ABSORBTION, LAI, AND ATTENUATION FACTOR
      DO 65 I=1,NC
         TOTDIR=0.0
         TOTDIF=0.0
         DO 63 J=1,NPLANT
            TOTDIR = TOTDIR + (1.-CANALB(J)-CANTRN(J))*DIRKL(J,I)
            TOTDIF = TOTDIF + (1.-CANALB(J)-CANTRN(J))*DIFKL(J,I)
   63    CONTINUE
         SWCAN(NPLANT+1,I)=0.0
         DO 64 J=1,NPLANT
            FRACT1=(1.-CANALB(J)-CANTRN(J))*DIRKL(J,I)/TOTDIR
            FRACT2=(1.-CANALB(J)-CANTRN(J))*DIFKL(J,I)/TOTDIF
            SWCAN(J,I)= FRACT1*(1.-DIRALB(I)-DIRTRN(I)) &
                *(1.-TDIRCC(I))*DIR(I) &
                + FRACT2*(1.-DIFALB(I)-DIFTRN(I)) &
                *(1.-TDIFFC(I))*(DOWN(I)+UP(I+1))
            SWCAN(NPLANT+1,I)=SWCAN(NPLANT+1,I)+SWCAN(J,I)
   64    CONTINUE
   65 CONTINUE
!**** DETERMINE THE AMOUNT OF RADIATION ABSORBED AT EACH RESIDUE NODE
      ABSORB=0.0
      DO 70 I=1,NR1
         K=I+NC
         SWRES(I)=DIR(K)+DOWN(K)+UP(K+1)-DIR(K+1)-DOWN(K+1)-UP(K)
         ABSORB=ABSORB + SWRES(I)
   70 CONTINUE
      IF (NR .EQ. 1) THEN
         SWRES(1)=SWRES(1)+SWRES(2)
         SWRES(2)=0.0
      END IF
!
!**** IF SNOW LAYER LIES ABOVE THE RESIDUE, ALL RADIATION LEAVING THE
!     THE RESIDUE IS ASSUMED TO BE REFLECTED BACK AND ABSORBED.
!     SPLIT UP THIS RADIATION PROPORTIONATELY BETWEEN THE LAYERS
      IF (NR.GT.0 .AND. NSP.GT.0) THEN
         IF (ABSORB .GT. 0.0) THEN
         DO 80 I=1,NR
            SWRES(I)=SWRES(I) + UP(1)*SWRES(I)/ABSORB
   80    CONTINUE
         END IF
      END IF
!
!**** DEFINE THE AMOUNT OF RADIATION AVAILABLE AT THE SOIL SURFACE
      DIRECT=DIR(NC+NR1+1)
      DIFFUS=DOWN(NC+NR1+1)
!
      RETURN
 END SUBROUTINE SWRCAN






! no problem checked changed
!***********************************************************************
!
      SUBROUTINE SOLVRAD (NLAYR,A,B,C,D,RD,RU)
!
!     THIS SUBROUTINE SOLVES THE MATRIX OF EQUATIONS FOR RADIATION
!     TRANSFER WITHIN THE CANOPY.  THE FORM OF THE MATRIX IS:
!
!     THIS SUBROUTINE SOLVES A TRI-DIAGONAL MATRIX OF THE FORM :
!
!              | B1  C1   0   0  .  .  .   0 | |RD1|   | D1|
!              | A2  B2   0  C2  .  .  .   0 | |RU1|   | D2|
!              | A3   0  B3  C3  .  .  .   0 | |RD2| = | D3|
!              |  .   .   .   .     .  .   . | | . |   |  .|
!              |  0   . AI1 BI1    0 CI1   0 | |RU |   |DI1|
!              |  0   . AI+1  0 BI+1 CI+1  0 | |RD | = |DI2|
!              |  .   .   .   .     .  .   . | | . |   |  .|
!              |  0   0   0   0  .  . AN  BN | |RUN|   | DN|
!
!
!     RD(I) IS THE DOWN FLUX ABOVE LAYER I
!     RU(I) IS THE UPWARD FLUX ABOVE LAYER I
!     (C1 IS NORMALLY ZERO, BUT IT DOES NOT HAVE TO BE)
!
!***********************************************************************
    implicit none
!input
    integer(i4),intent(in)::NLAYR
    real(r8),dimension(2*(NRMAX+NCMAX)),intent(inout)::A,B,C,D
    real(r8),dimension(NRMAX+NCMAX),intent(inout)::RD,RU

!temp
    integer(i4)::N,I,J,K
    real(r8)::QUOT

      DO 40 N=1,NLAYR
!        ELIMINATE EACH COEFFICIENT BELOW DIAGONAL IN COLUMN K
         K = 2*N - 1
         I = K + 1
!        DETERMINE MULTIPLICATION FACTOR TO ELIMATE A(I)
         QUOT = A(I)/B(K)
!        ADJUST ELEMENTS IN ROW (EQUATION) I
         B(I)=B(I)-QUOT*C(K)
         D(I)=D(I)-QUOT*D(K)
         A(I) = 0.0
!
!        DETERMINE MULTIPLICATION FACTOR TO ELIMATE A(I), BUT ZERO
!        ELEMENT IN ARRAY NOW TEMPORARILY ASSIGNED TO A(I)
         I = K + 2
         QUOT = A(I)/B(K)
!        ADJUST ELEMENTS IN ROW (EQUATION) I
         A(I)=-QUOT*C(K)
         D(I)=D(I)-QUOT*D(K)
!
         K = K + 1
!        DETERMINE MULTIPLICATION FACTOR TO ELIMATE A(I)
         QUOT = A(I)/B(K)
!        ADJUST ELEMENTS IN ROW (EQUATION) I
         C(I)=C(I)-QUOT*C(K)
         D(I)=D(I)-QUOT*D(K)
         A(I) = 0.0
!
   40 CONTINUE
!
      K = 2*(NLAYR+1) - 1
      I = K + 1
      QUOT = A(I)/B(K)
      B(I)=B(I)-QUOT*C(K)
      D(I)=D(I)-QUOT*D(K)
      A(I)=0.0
!
!     COMPUTE UPWELLING RADIATION FROM SNOW OR SOIL
      RU(NLAYR+1)=D(2*(NLAYR+1))/B(2*(NLAYR+1))
!
!     Remainder of back-substitution
      DO 60 N=NLAYR,1,-1
         I = 2*N + 1
!        COMPUTE DOWNWELLING RADIATION BELOW LAYER N
         RD(N+1) = (D(I) - C(I)*RU(N+1))/B(I)
         I = 2*N
!        COMPUTE UPWELLING RADIATION ABOVE LAYER N
         RU(N) = (D(I) - C(I)*RU(N+1))/B(I)
   60 CONTINUE
      RD(1) = D(1)/B(1)
!
      RETURN
 END SUBROUTINE SOLVRAD


! no problem checked changed
!***********************************************************************
!
 SUBROUTINE TRANSR (NR,TDIREC,TDIFFU,ZR,RHOR,ANGLE,DIRRES)
!
!THIS SUBROUTINE CALCULATES THE DIRECT AND DIFFUSE TRANSMISSION
!COEFFICIENTS FOR THE RESIDUE LAYERS, GIVEN THE ANGLE BETWEEN THE
!SURFACE AND THE LINE OF DIRECT APPROACH (SUN-SLOPE ANGLE)
!
!TRANSMISSION COEFFICIENTS ARE BASED ON THE TONNES/HA OF RESIDUE
!***********************************************************************
    use swrcoe_mod,only:DIFRES
    implicit none

!input
    real(r8),dimension(NRMAX),intent(inout)::TDIREC,TDIFFU
    real(r8),dimension(NRMAX),intent(in)::ZR,RHOR
    real(r8),intent(in)::ANGLE,DIRRES
    integer(i4),intent(in)::NR
!temp
    integer(i4)::I
    real(r8)::W


!**** DETERMINE THE TRANSMISSIVITY TO DIRECT AND DIFFUSE RADIATION
!     BETWEEN RESIDUE NODES BASED ON FRACTION OF SURFACE COVERED
!     CALCULATED FROM THE WEIGHT OF RESIDUE (TON/HA)
      DO 10 I=1,NR
!        DETERMINE THE WEIGHT OF THE RESIDUE AT EACH NODE
         IF (I .EQ. 1) THEN
            W=RHOR(I)*(ZR(2)-ZR(1))*10.
            IF (NR .NE. 1) W=W/2.
           ELSE
            IF (I .EQ. NR) THEN
               W=RHOR(I)*(ZR(NR+1)-ZR(NR)+(ZR(NR)-ZR(NR-1))/2)*5.
              ELSE
               W=RHOR(I)*(ZR(I+1)-ZR(I-1))*5.
            END IF
         END IF
!        NOW DETERMINE THE DIRECT AND DIFFUSE TRANMISSIVITIES
         IF (ANGLE .LE. 0.0) THEN
!           SUN IS NOT ABOVE HORIZON OF LOCAL SLOPE -- SET DIRECT
!           TRANSMISSIVITY TO ZERO
            TDIREC(I)=0.0
           ELSE
            TDIREC(I)=EXP(-DIRRES*W)*SIN(ANGLE)
         END IF
         TDIFFU(I)=DIFRES*EXP(-DIRRES*W)
   10 CONTINUE
!
      RETURN
 END SUBROUTINE TRANSR

! no problem checked changed
!***********************************************************************
 SUBROUTINE TRANSC (NPLANT,NC,XANGLE,CLUMPNG,TDIRCC,TDIFFC,&
    DIFKL,DIRKL,FBDU,FDDU,SUNSLP,ALTITU,col,row)
!THIS SUBROUTINE CALCULATES THE DIRECT AND DIFFUSE RADIATION
!TRANSMISSION COEFFICIENTS AND THE SUM OF THE ATTENUATION FACTOR
!TIMES LEAF AREA INDEX FOR EACH CANOPY LAYER. (THIS SUM IS USED IN
!CALCULATING AN EFFECTIVE ALBEDO FOR THE CANOPY LAYER AND IN
!PROPORTIONING ABSORBED RADIATION INTO EACH OF THE PLANT SPECIES.
!***********************************************************************
    use clayrs_mod,only:TOTLAI2d,CANLAI2d
    implicit none
!input
    integer(i4),intent(in)::NPLANT,NC,col,row
    real(r8),dimension(NPMAX),intent(in)::XANGLE,CLUMPNG
    real(r8),intent(in)::SUNSLP,ALTITU
    real(r8),dimension(NCMAX-1),intent(inout)::TDIRCC,TDIFFC
    real(r8),dimension(NPMAX+1,NCMAX-1),intent(inout)::DIRKL,DIFKL
    real(r8),dimension(NPMAX),intent(inout)::FBDU,FDDU
!temp
    real(r8),dimension(NPMAX)::DIRCAN,DIFCAN
    real(r8)::AA,BB,CC,DD,EE,FF,GG
    integer(i4)::I,J
    real(r8)::DIFINF,zen
    real(r8),dimension(:),pointer::TOTLAI
    real(r8),dimension(:,:),pointer::CANLAI
    TOTLAI=>TOTLAI2d(col,row,:)
    CANLAI=>CANLAI2d(col,row,:,:)
    AA=0.65
    BB=1.9
    CC=1.46
    DD=0.585
    EE=0.569
    FF=1.09
    GG=1.585
      DO 5 J=1,NPLANT
!        COMUPUTE ATTENUATION FACTOR BASED ON SOLAR AND LEAF ANGLE
         IF (SUNSLP .GT. 0.0) THEN
!           SUN IS ABOVE HORIZON
            IF (XANGLE(J) .LE. 0.0) THEN
!              LEAF ANGLE INCLINATION IS VERTICLE -- USE SOLAR ALTITUDE
!              INSTEAD OF ANGLE ON LOCAL SLOPE AND USE SIMPLIFIED
!              EXPRESSION FOR K
               IF (ALTITU .GE. 1.57) THEN
                  DIRCAN(J) = 0.0
                 ELSE
                  DIRCAN(J)=1./(1.5708*TAN(ALTITU))
               END IF
             ELSE
!                LEAF ORIENTATION OTHER THAN VERTICLE -- USE ANGLE ON
!                LOCAL SLOPE -- COMPUTE ZENITH ANGLE
               ZEN = 1.571 - SUNSLP
               DIRCAN(J)=(XANGLE(J)*XANGLE(J) + (TAN(ZEN))**2)**0.5 &
                /(XANGLE(J) + 1.774*(XANGLE(J)+1.182)**(-0.733))
            END IF
          ELSE
!          SUN IS BELOW HORIZON - JUST SET DIRCAN(J) TO 1.0
            DIRCAN(J) = 1.0
         END IF
!        COMPUTE Kd AT INFINITE LEAF AREA
         IF (XANGLE(J) .LE. 1.0) THEN
            DIFINF = ATAN(XANGLE(J))/1.5708
          ELSE
            DIFINF = XANGLE(J)**CC/(XANGLE(J)**CC + 1.)
         END IF
         DIFCAN(J)= (DIFINF*CLUMPNG(J)*TOTLAI(J)**AA + BB)/ &
            (CLUMPNG(J)*TOTLAI(J)**AA + BB)
!        DEFINE SCATTERING COEFF. FOR INDIVIDUAL PLANT TYPES COMPUTED
!        FROM FLERCHINGER & YU (2007, AG & FOREST MET).  FBDU IS
!        FRACTION OF DOWNWARD REFLECTED BEAM RADIATION SCATTERED
!        UPWARD. FDDU IS THE SAME FOR DIFFUSE RADIATION.
         FBDU(J) = 0.5 + 0.5*(ATAN(XANGLE(J))/1.5708)**DD &
            *XANGLE(J)**((COS(ALTITU))**(EE+FF*XANGLE(J)))*SIN(ALTITU)
         FDDU(J) = 0.5 + 0.5*(ATAN(XANGLE(J))/1.5708)**GG
    5 CONTINUE
!
!**** DETERMINE THE TRANSMISSIVITY TO DIRECT AND DIFFUSE RADIATION
!     BETWEEN CANOPY NODES
      DO 20 I=1,NC
         DIRKL(NPLANT+1,I)=0.0
         DIFKL(NPLANT+1,I)=0.0
!        CALCULATE PRODUCT OF ATTENUATION COEFFICIENT AND LAI
         DO 10 J=1,NPLANT
            DIRKL(J,I)=CLUMPNG(J)*CANLAI(J,I)*DIRCAN(J)
            DIFKL(J,I)=CLUMPNG(J)*CANLAI(J,I)*DIFCAN(J)
            DIRKL(NPLANT+1,I)=DIRKL(NPLANT+1,I)+DIRKL(J,I)
            DIFKL(NPLANT+1,I)=DIFKL(NPLANT+1,I)+DIFKL(J,I)
   10    CONTINUE
!        TRANSMISSIVITY OF LAYER IS EQUAL TO EXPONENT OF THE SUM OF
!        ATTENUATION FACTOR TIMES LEAF AREA INDEX
         IF (SUNSLP .LE. 0.0) THEN
!           SUN IS ON OR BELOW HORIZON OF LOCAL SLOPE -
!           SET DIRECT TRANSMISSIVITY TO ZERO.
            TDIRCC(I)=0.0
           ELSE
            TDIRCC(I)=EXP(-DIRKL(NPLANT+1,I))
         END IF
         TDIFFC(I)=EXP(-DIFKL(NPLANT+1,I))
   20 CONTINUE
      RETURN
 END SUBROUTINE TRANSC


! no problem checked changed
 !***********************************************************************
 SUBROUTINE SOLAR2 (DIRECT,DIFFUS,SUNSLP,ALTITU,SUNHOR,ALATUD,longtitude,SLOPE,&
    ASPECT,HOUR,NHRPDT,shadow,skyview,shadeef,julian,year)
!     THIS SUBROUTINE SEPARATES THE TOTAL RADIATION MEASURED ON THE
!     HORIZONTAL INTO THE DIRECT AND DIFFUSE ON THE LOCAL SLOPE.
!***********************************************************************
      !use controlpara_mod,only:WT,WDT,DT
      use controlpara_mod,only:timezone
      use calsun_mod
      use swrcoe_mod,only:SOLCON,DIFATM
      implicit none

!input
      integer(i4),intent(in)::HOUR,NHRPDT,julian,year
      real(r8),intent(in)::SUNHOR,SLOPE,ASPECT
      real(r8),intent(in)::ALATUD,longtitude
      real(r8),intent(inout)::DIRECT,DIFFUS,ALTITU,SUNSLP
      real(r8),intent(in)::shadow,skyview
      integer(i4),intent(in)::shadeef
!temp
      real(r8)::SUNRIS,SUNSET,HRWEST,SINAZM,COSAZM,SUMALT,COSALT,SUNMAX,HRANGL,SUN,TTOTAL,AZM,AZMUTH,DIRHOR,SINALT,TDIFFU
      integer(i4)::IHR,THOUR
      real(r8)::temp
      REAL(r8),parameter::Pi=3.1415926
      real(r8),parameter::DegToRad=Pi/180.0
      integer(i4)::mm,dd
      real(r8)::azimuth,solarZen,eqTime,solarDec
      real(r8)::zone,lat,lon

!**** CHECK IF SUN HAS RISEN YET (OR IF IT HAS ALREADY SET)
      IF (SUNHOR .LE. 0.0) THEN
         DIRECT=0.0
         DIFFUS=0.0
         RETURN
      END IF
      call GetDateFromJulianDay(year,julian,mm,dd)
!**** SUM UP VALUES AND FIND AVERAGE SUN POSITION FOR TIME STEP
      SINAZM=0.0
      COSAZM=0.0
      SUMALT=0.0
      COSALT=0.0
      SUNMAX=0.0
      lat=ALATUD/DegToRad
      lon=-longtitude/DegToRad
      !zone=-floor((longtitude/DegToRad+7.5)/15.0)
      zone=timezone
      do IHR=HOUR-NHRPDT,HOUR
        THOUR=IHR
        call calcsun(lat,lon,year,mm,dd,real(IHR),0.0,0.0,zone,azimuth,solarZen,eqTime,solarDec)
        ALTITU=(90.0-solarZen)*DegToRad
        AZM=azimuth*DegToRad
        if(ALTITU .gt. 0.0)then
            SUN=SOLCON*sin(ALTITU)
            SUMALT=SUMALT+SUN*sin(ALTITU)
            COSALT=COSALT+SUN*COS(ALTITU)
            SINAZM=SINAZM+SUN*SIN(AZM)
            COSAZM=COSAZM+SUN*COS(AZM)
            SUNMAX=SUNMAX+SUN
        end if
      end do
!**** DETERMINE AVERAGE SOLAR RADIATION, AVERAGE ALTITUDE AND AZIMUTH OF
!     THE SUN AND ANGLE ON LOCAL SLOPE
      IF (SUNMAX .EQ. 0) THEN
         ALTITU=0.0
         SUNSLP=0.0
      else
         ALTITU=ATAN(SUMALT/COSALT)
         AZMUTH=ATAN2(SINAZM,COSAZM)
         SUNMAX=SUNMAX/(NHRPDT+1)
         !SUNSLP=acos(SIN(slope)*COS(ALTITU)*(cos(ASPECT)*cos(AZMUTH)+sin(ASPECT)*sin(AZMUTH)) &
         !   +COS(slope)*SIN(ALTITU))
         SUNSLP=ASIN( SIN(ALTITU)*COS(SLOPE) &
            + COS(ALTITU)*SIN(SLOPE)*COS(AZMUTH-ASPECT))
      END IF
!
!**** SEPARATE THE SOLAR RADIATION INTO DIRECT AND DIFFUSE COMPONENTS
      IF (ALTITU .LE. 0.0) THEN
!     SUN IS BELOW THE HORIZON - ALL RADIATION MUST BE DIFFUSE
         DIFFUS=SUNHOR
         SUNSLP=0.0   !!!!!!!!!!!!!!!!!!bug
         DIRECT=0.0
         RETURN
      END IF
      TTOTAL=SUNHOR/SUNMAX
!     LIMIT TOTAL TRANSMISSIVITY TO MAXIMUM (DIFATM) WHICH WILL
!     CAUSE TDIFFU TO BE 0.0
      IF (TTOTAL .GT. DIFATM) TTOTAL = DIFATM
      TDIFFU=TTOTAL*(1. - EXP(0.6*(1.-DIFATM/TTOTAL)/(DIFATM-0.4)))
      DIFFUS=TDIFFU*SUNMAX
      DIRHOR=SUNHOR-DIFFUS
      if(shadeef.eq.1)then
        DIRHOR=DIRHOR*(1.0-shadow)
        DIFFUS=DIFFUS*skyview
      end if
!
!**** NOW CALCULATE THE DIRECT SOLAR RADIATION ON THE LOCAL SLOPE
      IF (SUNSLP .LE. 0.0) THEN
!        SUN HAS NOT RISEN ON THE LOCAL SLOPE -- NO DIRECT RADIATION
         DIRECT=0.0
       ELSE
         DIRECT=DIRHOR*SIN(SUNSLP)/SIN(ALTITU)
!        IF THE SUN'S ALTITUDE IS NEAR ZERO, THE CALCULATED DIRECT
!        RADIATION ON THE SLOPING SURFACE MAY BE UNREALISTICALLY LARGE --
!        LIMIT DIRECT TO 5*DIRHOR (THIS IS APPROXIMATELY THE CASE WHEN
!        THE SUN IS 10 DEGREES ABOVE THE HORIZON AND THE SLOPING SURFACE
!        IS PERPENDICULAR TO THE SUN`S RAYS
         IF (DIRECT .GT. 5.*DIRHOR) DIRECT=5.*DIRHOR
      END IF
      RETURN
  END SUBROUTINE SOLAR2


! no problem checked changed
 !***********************************************************************
 SUBROUTINE SOLAR (DIRECT,DIFFUS,SUNSLP,ALTITU,SUNHOR,ALATUD,SLOPE,&
    ASPECT,HRNOON,HAFDAY,DECLIN,HOUR,NHRPDT,shadow,skyview,shadeef,julian)
!     THIS SUBROUTINE SEPARATES THE TOTAL RADIATION MEASURED ON THE
!     HORIZONTAL INTO THE DIRECT AND DIFFUSE ON THE LOCAL SLOPE.
!***********************************************************************
    !use controlpara_mod,only:WT,WDT,DT
    use swrcoe_mod,only:SOLCON,DIFATM
    implicit none

!input
    integer(i4),intent(in)::HOUR,NHRPDT,julian
    real(r8),intent(in)::SUNHOR,ALATUD,SLOPE,ASPECT,HRNOON,HAFDAY,DECLIN
    real(r8),intent(inout)::DIRECT,DIFFUS,ALTITU,SUNSLP
    real(r8),intent(in)::shadow,skyview
    integer(i4),intent(in)::shadeef
!temp
    real(r8)::SUNRIS,SUNSET,HRWEST,SINAZM,COSAZM,SUMALT,COSALT,SUNMAX,HRANGL,SUN,TTOTAL,AZM,AZMUTH,DIRHOR,SINALT,TDIFFU
    integer(i4)::IHR,THOUR
    real(r8)::temp



!**** CHECK IF SUN HAS RISEN YET (OR IF IT HAS ALREADY SET)
      IF (SUNHOR .LE. 0.0) THEN
         DIRECT=0.0
         DIFFUS=0.0
         RETURN
      END IF
      SUNRIS=HRNOON - HAFDAY/0.261799
      SUNSET=HRNOON + HAFDAY/0.261799
!

!**** CALCULATE HOUR ANGLE AT WHICH THE SUN WILL BE DUE EAST/WEST IN
!     ORDER TO ADJUST AZIMUTH ANGLE FOR SOUTHERN AZIMUTHS
!     -- SIN(AZIMUTH) TELLS YOU ONLY THE EAST/WEST DIRECTION - NOT
!     WHETHER THE SUN IS NORTH/SOUTH.
      IF (ABS(DECLIN).GE.ABS(ALATUD)) THEN
!        LATITUDE IS WITHIN THE TROPICS (EQUATION WON'T WORK)
         HRWEST=3.14159
        ELSE
         HRWEST=ACOS(TAN(DECLIN)/TAN(ALATUD))
      END IF
!
!**** SUM UP VALUES AND FIND AVERAGE SUN POSITION FOR TIME STEP
      SINAZM=0.0
      COSAZM=0.0
      SUMALT=0.0
      COSALT=0.0
      SUNMAX=0.0
      DO 10 IHR=HOUR-NHRPDT,HOUR
         THOUR=IHR
!
!****    DETERMINE THE GEOMETRY OF THE SUN'S RAYS AT CURRENT TIME
         HRANGL=0.261799*(IHR-HRNOON)
         IF (THOUR .GT. SUNRIS  .AND.  THOUR .LT. SUNSET) THEN
!           SUN IS ABOVE HORIZON -- CALCULATE ITS ALTITUDE ABOVE THE
!           HORIZON (ALTITU) AND ANGLE FROM DUE NORTH (AZMUTH)
            SINALT=SIN(ALATUD)*SIN(DECLIN) &
                + COS(ALATUD)*COS(DECLIN)*COS(HRANGL)
            ALTITU=ASIN(SINALT)
            temp=-COS(DECLIN)*SIN(HRANGL)/COS(ALTITU)
            if(temp.lt.-1.0) temp=-1.0
            if(temp.gt.1.0)  temp=1.0
            AZM = ASIN(temp)

!           CORRECT AZIMUTH FOR SOUTHERN ANGLES
            IF (ALATUD-DECLIN .GT. 0.0) THEN
!              NORTHERN LATITUDES   (HRANGL=0.0 AT NOON)
               IF (ABS(HRANGL).LT.HRWEST) AZM=3.14159-AZM
              ELSE
!              SOUTHERN LATITUDES
               IF (ABS(HRANGL).GE.HRWEST) AZM=3.14159-AZM
            END IF
!           SUM CONDITIONS TO GET AVERAGE ALTITUDE AND AZMUTH
!           (OBTAIN AVERAGE BY SUMMING VECTOR COMPONENTS)
            SUN=SOLCON*SINALT
            SUMALT=SUMALT+SUN*SINALT
            COSALT=COSALT+SUN*COS(ALTITU)
            SINAZM=SINAZM+SUN*SIN(AZM)
            COSAZM=COSAZM+SUN*COS(AZM)
            SUNMAX=SUNMAX+SUN
         END IF
!
   10 CONTINUE
!
!**** DETERMINE AVERAGE SOLAR RADIATION, AVERAGE ALTITUDE AND AZIMUTH OF
!     THE SUN AND ANGLE ON LOCAL SLOPE
      IF (SUNMAX .EQ. 0) THEN
         ALTITU=0.0
         SUNSLP=0.0
        ELSE
         ALTITU=ATAN(SUMALT/COSALT)
         AZMUTH=ATAN2(SINAZM,COSAZM)
         SUNMAX=SUNMAX/(NHRPDT+1)
         SUNSLP=ASIN( SIN(ALTITU)*COS(SLOPE) &
            + COS(ALTITU)*SIN(SLOPE)*COS(AZMUTH-ASPECT))
      END IF
      print*,"solar",ALTITU,AZMUTH
!
!**** SEPARATE THE SOLAR RADIATION INTO DIRECT AND DIFFUSE COMPONENTS
      IF (ALTITU .LE. 0.0) THEN
!     SUN IS BELOW THE HORIZON - ALL RADIATION MUST BE DIFFUSE
         DIFFUS=SUNHOR
         SUNSLP=0.0   !!!!!!!!!!!!!!!!!!bug
         DIRECT=0.0
         RETURN
      END IF
      TTOTAL=SUNHOR/SUNMAX
!     LIMIT TOTAL TRANSMISSIVITY TO MAXIMUM (DIFATM) WHICH WILL
!     CAUSE TDIFFU TO BE 0.0
      IF (TTOTAL .GT. DIFATM) TTOTAL = DIFATM
      TDIFFU=TTOTAL*(1. - EXP(0.6*(1.-DIFATM/TTOTAL)/(DIFATM-0.4)))
      DIFFUS=TDIFFU*SUNMAX
      DIRHOR=SUNHOR-DIFFUS
      if(shadeef.eq.1)then
        DIRHOR=DIRHOR*(1.0-shadow)
        DIFFUS=DIFFUS*skyview
      end if
!**** NOW CALCULATE THE DIRECT SOLAR RADIATION ON THE LOCAL SLOPE
      IF (SUNSLP .LE. 0.0) THEN
!        SUN HAS NOT RISEN ON THE LOCAL SLOPE -- NO DIRECT RADIATION
         DIRECT=0.0
       ELSE
         DIRECT=DIRHOR*SIN(SUNSLP)/SIN(ALTITU)
!        IF THE SUN'S ALTITUDE IS NEAR ZERO, THE CALCULATED DIRECT
!        RADIATION ON THE SLOPING SURFACE MAY BE UNREALISTICALLY LARGE --
!        LIMIT DIRECT TO 5*DIRHOR (THIS IS APPROXIMATELY THE CASE WHEN
!        THE SUN IS 10 DEGREES ABOVE THE HORIZON AND THE SLOPING SURFACE
!        IS PERPENDICULAR TO THE SUN`S RAYS
         IF (DIRECT .GT. 5.*DIRHOR) DIRECT=5.*DIRHOR
      END IF
!
      RETURN
  END SUBROUTINE SOLAR

! changed
!***********************************************************************
 SUBROUTINE SNOALB (NSP,ALBSNO,ZSP,RHOSP,SUNSLP,ALBNXT)
!     THIS SUBROUTINE COMPUTES THE ALBEDO, THE EXTINCTION COEFFICIENT
!     AND THE AMOUNT OF SOLAR RADIATION ABSORBED BY EACH LAYER
!***********************************************************************
    use constvar_mod,only:RHOL
    use swrcoe_mod,only:SNOCOF,SNOEXP
    use spwatr_mod,only:EXTSP,G1,G2,G3
    implicit none
!input
    integer(i4),intent(in)::NSP
    real(r8),dimension(NSPMAX),intent(in)::ZSP,RHOSP
    real(r8),intent(inout)::ALBSNO
    real(r8),intent(in)::SUNSLP,ALBNXT
!temp
    real(r8)::SPGRAV,GRAIN,ABS,W,EXC,Y


!     DETERMINE THE ALBEDO OF THE SNOW
      SPGRAV = RHOSP(1)/RHOL
      GRAIN = G1 + G2*SPGRAV*SPGRAV + G3*SPGRAV**4
      ALBSNO = 1. - 0.206*EXTSP*SQRT(GRAIN)
      IF (ALBSNO .LT. 0.35) ALBSNO = 0.35
!
      IF (ZSP(NSP+1) .LE. 0.04) THEN
!        SNOWPACK IS LESS THAN 4.0 CM -- ALBEDO IS AFFECTED BY
!        UNDERLYING MATERIAL
         ABS = 1. - ALBNXT
         W = 2.*(1. - ALBSNO)/(1. + ALBSNO)
!        EXC = EXTINCTION COEFFICIENT --> CONVERT FROM 1/CM TO 1/M
         EXC = EXTSP*SPGRAV*SQRT(1.0/GRAIN)*100.
         Y = EXP(-EXC*ZSP(NSP+1))*(ABS + W*(ABS/2.-1.)) &
            / (W*(ABS/2.-1.)*COSH(EXC*ZSP(NSP+1))-ABS*SINH(EXC*ZSP(NSP+1)))
         ALBSNO = (1. - W*(1.-Y)/2)/(1. + W*(1.-Y)/2)
      END IF
!
      IF (ZSP(NSP+1) .LE. SNOCOF) THEN
!        SNOCOF IS MINIMUM DEPTH OF SNOW FOR COMPLETE GROUND COVER --
!        SNOW IS NOT COMPLETELY COVERED WITH SNOW
         ALBSNO = ALBNXT + (ALBSNO-ALBNXT)*(ZSP(NSP+1)/SNOCOF)**SNOEXP
      END IF
      RETURN
 END SUBROUTINE SNOALB


! no problem checked changed
!***********************************************************************
 subroutine systm (nc,nplant,nsp,zc,zsp,zmsrf,zhsrf,zersrf,zmsp,&
    zhsp,height,solr,ta,tccrit,col,row)
!   this subroutine defines the roughness parameters for the system,
!   and determines whether the canopy will transpire
!***********************************************************************
    !xxxx common /windv/ zh,zm,zero,ustar,stable,windc(11),windr(10)
    use windv_mod,only:zm2d,zh2d,zero2d,zersub2d,zmsub2d,zhsub2d
    use clayrs_mod,only:ievap2d,canlai2d,totlai2d
    implicit none
!input
    integer(i4),intent(in)::nc,nplant,nsp,col,row
    real(r8),dimension(ncmax),intent(in)::zc
    real(r8),dimension(nspmax),intent(in)::zsp
    real(r8),dimension(npmax),intent(in)::tccrit
    real(r8),intent(in)::zmsrf,zhsrf,zersrf,zmsp,zhsp,height,solr,ta

!temp variable
    real(r8),pointer::zm,zh,zero,zmsub,zhsub,zersub
    integer(i4),dimension(:),pointer::ievap
    real(r8),dimension(:,:),pointer::canlai
    real(r8),dimension(:),pointer::totlai

    real(r8)::zerfull,zmfull,sumlai,v
    integer(i4)::i,j,np


    zm=>zm2d(col,row)
    zh=>zh2d(col,row)
    zero=>zero2d(col,row)
    zersub=>zersub2d(col,row)
    zmsub=>zmsub2d(col,row)
    zhsub=>zhsub2d(col,row)
    ievap=>ievap2d(col,row,:)
    canlai=>canlai2d(col,row,:,:)
    totlai=>totlai2d(col,row,:)


      zm=zmsrf
      zh=zhsrf
      zero=zersrf
      if (nc .gt. 0) then
         zero = 0.77*zc(nc+1)
         zm = 0.13*zc(nc+1)
         if (nsp .gt. 0) then
            zmsub=zmsp
            zhsub=zhsp
            zersub=0.0
           else
            zmsub=zmsrf
            zhsub=zhsrf
            zersub=zersrf
         end if
!xxxx
!        compute lai and reduce roughness parameters for sparse canopies
!        based on work of zeng and wang (2007; j hydromet 8:730-737)
!        (assumes lai = 2.0 is a full canopy)
!
!        start at top of canopy and compute roughness parameters based
!        plants present - stop when all plants are considered or when
!        it gets below computed zero displacement.
!        this will compensate for a sparse tree canopy over a dense
!        understory if zero is above the understory
         zerfull = zero
         zmfull = zm
         do i=1,nc
            np = 0
            sumlai=0.0
            do j=1,nplant
               if (canlai(j,i) .gt. 0.0) sumlai = sumlai + totlai(j)
               np = np + 1
            end do
            if (sumlai .lt. 2.0) then
               v = (1.-exp(-sumlai))/(1.-exp(-2.))
               zero = v*zerfull + (1.-v)*zersub
!xxxx          zm = exp(v*alog(zm) + (1.-v)*alog(zmsub))
               zm = v*zmfull + (1.-v)*zmsub
              else
               zero = zerfull
               zm = zmfull
            end if
            if (np .eq. nplant .or. zc(nc+1)-zc(i+1).le.zero) go to 15
        end do
   15    zh = 0.2*zm

!
!        no tranpiration at night
         if (solr .lt. 10.0) then
            do j=1,nplant
               ievap(j)=0
            end do
          else
            do j=1,nplant
               ievap(j)=1
!              check if temperature is too cold for transpiration
               if (ta.le.tccrit(j)) ievap(j)=0
            end do
         end if
        else
!
         if (nsp .gt. 0) then
!           snow is surface material -- define roughness elements and
!           zero plane of displacement
            zm=zmsp
            zh=zhsp
            zero=zsp(nsp+1)
         end if
         if(zero>zc(nc+1))zero=zc(nc+1)-0.01
      end if
!     do not let displacement plane be above height of instruments
      if (zero+zm .gt. height) zero=height/2.
      return
 end subroutine systm


! no problem checked changed
!***********************************************************************
 SUBROUTINE SOILHT (NS,CS,VLC,VIC,TS,MAT,CONC,nsalt,col,row)
    !THIS SUBROUTINE CALCULATES THE SPECIFIC HEAT THE SOIL LAYERS
    !THE LATENT HEAT OF VAPORIZATION IS INCLUDED.
!***********************************************************************
    use constvar_mod,only:CA,CI,CL,CM,COM,G,LV,RHOA,RHOI,RHOL,RHOM,RHOOM,UGAS
    use soilproperty_mod,only:OM2D,RHOB2d,SAND2d,CLAY2d,SILT2d,ENTRY2d
    implicit none
! input
    integer(i4),intent(in)::NS,col,row,nsalt
    real(r8),dimension(NSMAX),intent(in)::VLC,VIC,TS,MAT
    real(r8),dimension(NSMAX),intent(inout)::CS
    real(r8),dimension(NSALTMAX,NSMAX),intent(in)::CONC
! temp
    integer(i4)::I,J
    real(r8)::XSM,VAC,TLCONC,TOTPOT,HUMID,S,RHOV
    real(r8),dimension(:),pointer::OM,RHOB,SAND,CLAY,SILT,ENTRY



    OM=>OM2D(col,row,:)
    RHOB=>RHOB2d(col,row,:)
    SAND=>SAND2d(col,row,:)
    CLAY=>CLAY2d(col,row,:)
    SILT=>SILT2d(col,row,:)
    ENTRY=>ENTRY2d(col,row,:)

      DO 20 I=1,NS
!        CORRECTION OF MINERAL FRACTION FOR ORGANIC MATTER
         XSM = 1.0 - OM(I)
!
         VAC = 1. -XSM*RHOB(I)/RHOM -OM(I)*RHOB(I)/RHOOM -VLC(I) -VIC(I)
         IF (VAC .LT. 0.0) VAC=0.0
         CS(I)=XSM*RHOB(I)*CM + OM(I)*RHOB(I)*COM + VLC(I)*RHOL*CL &
            + VIC(I)*RHOI*CI + VAC*RHOA*CA
!
!****    INCLUDE THE LATENT HEAT OF VAPORIZATION IN THE HEAT CAPACITY
!        IF LAYER IS NOT SATURATED
         IF (MAT(I).LT.ENTRY(I) .AND. VAC.GT.0.0) THEN
!           DETERMINE HUMIDITY OF LAYER FROM TOTAL WATER POTENTIAL
            TLCONC=0.0
            DO 10 J=1,NSALT
               TLCONC= TLCONC+CONC(J,I)
   10       CONTINUE
            TOTPOT=MAT(I) - TLCONC*UGAS*TS(I)/G
            HUMID=EXP(.018*G/UGAS/(TS(I)+273.16)*TOTPOT)
!           OBTAIN SLOPE OF VAPOR DENSITY CURVE
         if(abs(TS(I)) .gt. 120.0)then
            print*,"In SOILHT TS",col,row,TS(I)
            stop
         endif
            CALL VSLOPE (S,RHOV,TS(I))
            CS(I)=CS(I)+VAC*LV*HUMID*S
         END IF
   20 CONTINUE
      RETURN
 END SUBROUTINE SOILHT

!no problem checked changed
!***********************************************************************
 SUBROUTINE SOILTK (NS,TK,VLC,VIC,col,row)
    !THIS SUBROUTINE CALCULATES THE THERMAL CONDUCTIVITY OF THE SOIL
    !LAYERS USING DE VRIES'S METHOD.
!***********************************************************************
    use constvar_mod,only:GAC,GAOM,GASA,GASI,RHOM,RHOOM,TKA,TKCL,TKI,TKL,TKOM,TKSA,TKSI
    use soilproperty_mod,only:OM2D,RHOB2d,SAND2d,CLAY2d,SILT2d
    use savedata_mod,only:WFAIRD2d,WFSAD2d,WFSID2d,WFCLD2d,WFOMD2d,WFICED2d,TKMA2d,&
        WFL2d,WFSA2d,WFSI2d,WFCL2d,WFOM2d,WFICE2d,IFIRST2d
    implicit none
! input
    integer(i4),intent(in)::NS,col,row
    real(r8),dimension(NSMAX),intent(in)::VLC,VIC
    real(r8),dimension(NSMAX),intent(inout)::TK

! temp
    integer(i4)::I
    real(r8)::XSM,XSAND,XSILT,XCLAY,VSAND,VSILT,VCLAY,VOM,PORO,VAC,VLMT, &
        GAAIR,WFAIR,TK0,TKLMT

    real(r8),pointer::WFAIRD,WFSAD,WFSID,WFCLD,WFOMD,WFICED,TKMA,&
        WFL,WFSA,WFSI,WFCL,WFOM,WFICE
    integer(i4),pointer ::IFIRST

    real(r8),dimension(:),pointer::OM,RHOB,SAND,CLAY,SILT

!   WEIGHTING FACTOR FOR ICE IS NOW SET TO 1.0 TO CORRECT PROBLEMS
!   POINTED OUT BY KENNEDY & SHARRATT, SOIL SCI. 163(8):636-645.
!   DATA GAI, GAOM, GASA, GASI, GAC /0.33, 0.5, 0.144, 0.144, 0.125/
!   See in constvar_mod

    WFAIRD=>WFAIRD2d(col,row)
    WFSAD=>WFSAD2d(col,row)
    WFSID=>WFSID2d(col,row)
    WFCLD=>WFCLD2d(col,row)
    WFOMD=>WFOMD2d(col,row)
    WFICED=>WFICED2d(col,row)
    TKMA=>   TKMA2d(col,row)
    WFL=>WFL2d(col,row)
    WFSA=>WFSA2d(col,row)
    WFSI=>WFSI2d(col,row)
    WFCL=>WFCL2d(col,row)
    WFOM=>WFOM2d(col,row)
    WFICE=>WFICE2d(col,row)
    IFIRST=>IFIRST2d(col,row)

    OM=>OM2D(col,row,:)
    RHOB=>RHOB2d(col,row,:)
    SAND=>SAND2d(col,row,:)
    CLAY=>CLAY2d(col,row,:)
    SILT=>SILT2d(col,row,:)


      IF (IFIRST .EQ. 0) THEN
!        FIRST TIME INTO SUBROUTINE -- CALCULATE WEIGHTING FACTORS FOR
!        SOIL COMPONENTS WHEN SOIL IS COMPLETELY DRY AND WHEN WET
         WFAIRD = 1.0
         WFSAD = DeVWF(TKA,TKSA,GASA)
         WFSID = DeVWF(TKA,TKSI,GASI)
         WFCLD = DeVWF(TKA,TKCL,GAC)
         WFOMD = DeVWF(TKA,TKOM,GAOM)
         WFICED = 1.0
!
!        SET THERMAL CONDUCTIVITY OF MOIST AIR EQUAL TO THAT OF DRY AIR
!        (VAPOR TRANSFER IS CALCULATED SEPARATELY AND INCLUDED IN HEAT
!        TRANSFER THROUGH SOIL; THEREFORE IT SHOULD NOT BE COMPENSATED
!        FOR IN COMPUTATION FOR THERMAL CONDUCTIVITY AS IS USUALLY DONE
!        WHEN USING DEVRIES METHOD.)
         TKMA = TKA
         WFL = 1.0
         WFSA = DeVWF(TKL,TKSA,GASA)
         WFSI = DeVWF(TKL,TKSI,GASI)
         WFCL = DeVWF(TKL,TKCL,GAC)
         WFOM = DeVWF(TKL,TKOM,GAOM)
         WFICE = 1.0
         IFIRST=1
      END IF
!
      DO 10 I=1,NS
!       CORRECTION OF TEXTURE FOR ORGANIC MATTER
        XSM = 1.0 - OM(I)
        XSAND=XSM*SAND(I)
        XSILT=XSM*SILT(I)
        XCLAY=XSM*CLAY(I)
!
!       CONVERT WEIGHT FRACTION TO VOLUMETRIC FRACTION
        VSAND = XSAND*RHOB(I)/RHOM
        VSILT = XSILT*RHOB(I)/RHOM
        VCLAY = XCLAY*RHOB(I)/RHOM
        VOM = OM(I)*RHOB(I)/RHOOM
        PORO = 1.-VSAND-VSILT-VCLAY-VOM
        VAC = PORO-VIC(I)-VLC(I)
        IF (VAC .LT. 0.0) VAC=0.0
!
!       DETERMINE LIMIT OF DEVRIES METHOD FROM TEXTURE
!       SAND ==> 0.05;  LOAM ==> 0.10;  CLAY ==> 0.15
        VLMT = 0.10 + 0.2*CLAY(I) - 0.1*SAND(I)
        IF (VLMT .LT. 0.05) VLMT = 0.05
        IF (VLMT .GT. 0.15) VLMT = 0.15
!
        IF (VLC(I) .GT. VLMT) THEN
!         MOIST SOIL -- USE DEVRIES METHOD DIRECTLY
          GAAIR = 0.035+0.298*(VLC(I)-VLMT)/(PORO-VLMT)
          WFAIR= DeVWF(TKL,TKMA,GAAIR)
          TK(I)=(WFSA*VSAND*TKSA + WFL*VLC(I)*TKL + WFICE*VIC(I)*TKI &
            + WFAIR*VAC*TKMA + WFSI*VSILT*TKSI + WFCL*VCLAY*TKCL &
            + WFOM*VOM*TKOM)  /  (WFSA*VSAND + WFL*VLC(I) + WFICE*VIC(I) &
            + WFAIR*VAC + WFSI*VSILT + WFCL*VCLAY + WFOM*VOM)
        ELSE
!         INTERPOLATE THERMAL CONDUCTIVITY BETWEEN WATER CONTENT AT 0
!         AND LIMIT OF DEVRIES METHOD
          TK0=1.25*(WFSAD*VSAND*TKSA + WFICED*VIC(I)*TKI &
            + WFAIRD*(VAC+VLC(I))*TKA &
            + WFSID*VSILT*TKSI + WFCLD*VCLAY*TKCL + WFOMD*VOM*TKOM) &
            /(WFSAD*VSAND + WFICED*VIC(I) + WFAIRD*(VAC+VLC(I)) &
            + WFSID*VSILT + WFCLD*VCLAY + WFOMD*VOM)
!
          WFAIR= DeVWF(TKL,TKMA,0.035_r8)
          TKLMT =(WFSA*VSAND*TKSA + WFL*VLMT*TKL + WFICE*VIC(I)*TKI &
            + WFAIR*VAC*TKMA + WFSI*VSILT*TKSI + WFCL*VCLAY*TKCL &
            + WFOM*VOM*TKOM) / (WFSA*VSAND + WFL*VLMT + WFICE*VIC(I) &
            + WFAIR*VAC + WFSI*VSILT + WFCL*VCLAY + WFOM*VOM)
!
          TK(I) = TK0 + (TKLMT-TK0)*VLC(I)/VLMT
        ENDIF
   10 CONTINUE
      RETURN
 END SUBROUTINE SOILTK

 real(r8) FUNCTION DeVWF(TK0,TK1,GA)
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    real(r8),intent(in)::TK0,TK1,GA
    DeVWF = 2.0/3.0/(1.0+(TK1/TK0-1.0)*GA)&
            + 1.0/3.0/(1.0+(TK1/TK0-1.0)*(1.0-2.0*GA))
    return
 END FUNCTION



!***********************************************************************
!
 SUBROUTINE WBALNC (NPLANT,NC,NSP,NR,NS,JULIAN,HOUR,YEAR,&
    ITYPE,INITAL,ZC,WCAN,WCANDT,PCAN,PCANDT,VAPC,VAPCDT,RHOSP,DZSP,&
    DLWDT,WLAG,STORE,ZR,GMC,GMCDT,VAPR,VAPRDT,RHOR,ZS,VLC,VLCDT,&
    VIC,VICDT,TOTFLO,PRECIP,RUNOFF,POND,EVAP1,ETSUM,col,row)
!
!     THIS SUBROUTINE SUMS THE EVAPORATION AND DEEP PERCOLATION AT THE
!     END OF EACH HOUR, THEN PRINTS THE SUMMARY AT THE DESIRED OUTPUT
!     INTERVAL
!
!***********************************************************************
    use constvar_mod
    use clayrs_mod
    use controlpara_mod,only:LVLOUT
    use waterbal_mod
    implicit none
!input
    integer(i4),intent(in)::NPLANT,NC,NSP,NR,NS,JULIAN,HOUR,YEAR,INITAL,col,row
    integer(i4),dimension(NPMAX),intent(in)::ITYPE
    real(r8),dimension(NCMAX),intent(in)::ZC
    real(r8),dimension(NCMAX),intent(in)::VAPC,VAPCDT
    real(r8),dimension(NCMAX-1),intent(in)::WCAN,WCANDT
    real(r8),dimension(NPMAX),intent(in)::PCAN,PCANDT
    real(r8),dimension(NSPMAX),intent(in)::RHOSP,DZSP,DLWDT
    real(r8),dimension(11),intent(in)::WLAG
    real(r8),intent(in)::STORE
    real(r8),dimension(NRMAX),intent(in)::ZR,GMC,GMCDT,VAPR,VAPRDT,RHOR
    real(r8),dimension(NSMAX),intent(in)::ZS,VLC,VLCDT,VIC,VICDT,TOTFLO
    real(r8),intent(in)::PRECIP,RUNOFF,POND,EVAP1,ETSUM

!temp variable
    real(r8),pointer::RAIN,DPCAN,DCAN,DSNOW,DRES,DSOIL,TRUNOF,POND2,TPERC,&
    TETSUM,TEVAP,CUMVAP,SWE
    real(r8),dimension(:,:),pointer::DRYCAN
    integer(i4)::I,J
    real(r8)::SWEDT,DPOND,DZ,ERROR,TPCAN



    RAIN=>RAIN2d(col,row)
    DPCAN=>DPCAN2d(col,row)
    DCAN=>DCAN2d(col,row)
    DSNOW=>DSNOW2d(col,row)
    DRES=>DRES2d(col,row)
    DSOIL=>DSOIL2d(col,row)
    TRUNOF=>TRUNOF2d(col,row)
    POND2=>POND22d(col,row)
    TPERC=>TPERC2d(col,row)
    TETSUM=>TETSUM2d(col,row)
    TEVAP=>TEVAP2d(col,row)
    CUMVAP=>CUMVAP2d(col,row)
    SWE=>SWE2d(col,row)

    DRYCAN=>DRYCAN2d(col,row,:,:)


    IF (INITAL .EQ. 0) THEN
    !CALCULATE THE SNOW WATER EQUIVALENT FOR THE INITIAL CONDITIONS
        SWE=0.0
        DO I=1,NSP
            SWE=SWE + RHOSP(I)*DZSP(I)/RHOL + DLWDT(I)
        ENDDO
        RETURN
    END IF

    !END OF THE HOUR -- DETERMINE THE CHANGE IN STORAGE
    !CHANGE IN STORAGE OF CANOPY
    DO I=1,NC
        DO J=1,NPLANT
        !CALCULATE CHANGE IN WATER CONTENT OF DEAD PLANT MATERIAL
            IF (ITYPE(J) .EQ. 0) DCAN = DCAN +(WCANDT(I)-WCAN(I))*DRYCAN(J,I)*1000./RHOL
        ENDDO
        !INCLUDE CHANGE IN VAPOR DENSITY OF AIR SPACE
        IF (I .EQ. 1)  THEN
            DZ = (ZC(2) - ZC(1))/2.
        ELSE
            DZ = (ZC(I+1) - ZC(I-1))/2.
        END IF
        DCAN = DCAN + DZ*(VAPCDT(I)-VAPC(I))*1000./RHOL
    ENDDO

    !CHANGE IN SNOWPACK
    SWEDT = 0.0
    IF (NSP .GT. 0) THEN
        DO I=1,NSP
            SWEDT = SWEDT + RHOSP(I)*DZSP(I)/RHOL + DLWDT(I)
        ENDDO
        !WATER IN PROCESS OF BEING LAGGED THROUGH THE SNOWPACK
        SWEDT = SWEDT + STORE
        DO I=1,11
            SWEDT = SWEDT + WLAG(I)
        ENDDO
    END IF
    !CALCULATE CHANGE IN WATER CONTENT OF SNOWPACK IN MILLIMETERS
    DSNOW = DSNOW + (SWEDT - SWE)*1000.
    SWE = SWEDT
    !CHANGE IN STORAGE OF RESIDUE
    DO I=1,NR
        IF (I .EQ. 1 .OR.  I .EQ. NR) THEN
            IF (NR .EQ.1) THEN
                DZ = ZR(NR+1)
            ELSE
                IF (I .EQ. 1)  DZ = (ZR(2) - ZR(1))/2.
                IF (I .EQ. NR) DZ = ZR(NR+1)-ZR(NR)+(ZR(NR)-ZR(NR-1))/2.
            END IF
        ELSE
            DZ = (ZR(I+1) - ZR(I-1))/2.
        END IF
        DRES = DRES + DZ*((GMCDT(I)-GMC(I))*RHOR(I)&
                +(VAPRDT(I)-VAPR(I)))*1000./RHOL
    ENDDO

    !CHANGE IN STORAGE FOR THE SOIL
    DO I=1,NS-1
        IF (I .EQ. 1) THEN
            DZ = (ZS(2) - ZS(1))/2.
        ELSE
            DZ = (ZS(I+1) - ZS(I-1))/2.
        END IF
        DSOIL= DSOIL + DZ*1000*(VLCDT(I) - VLC(I)&
            + (VICDT(I) - VIC(I))*RHOI/RHOL)
    ENDDO

    !COMPUTE THE CHANGE IN PRECIP INTERCEPTED ON PLANT LEAVES
    TPCAN=0.0
    DO J=1,NPLANT
        IF (ITYPE(J).NE.0) TPCAN=TPCAN + PCANDT(J)
        IF (ITYPE(J).NE.0) DPCAN=DPCAN + PCANDT(J)-PCAN(J)
    ENDDO

    RAIN = RAIN + PRECIP
    TRUNOF = TRUNOF + RUNOFF
    TETSUM = TETSUM + ETSUM
    TEVAP = TEVAP + EVAP1
    TPERC = TPERC + TOTFLO(NS-1)

    !RETURN IF HOURLY OUTPUT IS NOT REQUIRED AND IT IS NOT END OF DAY
    IF (MOD(HOUR,LVLOUT(7)) .NE. 0) RETURN

    !CONVERT TO MILLIMETERS
    TPCAN=TPCAN*1000.
    DPCAN=DPCAN*1000.
    DPOND = POND*1000. - POND2
    POND2 = POND*1000.
    RAIN = RAIN*1000.
    TRUNOF = TRUNOF*1000.
    TETSUM = TETSUM*1000.
    TEVAP = TEVAP*1000.
    TPERC = TPERC*1000.
    CUMVAP = CUMVAP + TEVAP
    ERROR = RAIN+TEVAP-TPERC-TRUNOF-DPOND-DPCAN-DCAN-DSNOW-DRES-DSOIL
!          WRITE (26,110)JULIAN,HOUR,YEAR,RAIN,TPCAN,-TEVAP,TETSUM,DCAN,
!         >              DSNOW,DRES,DSOIL,TPERC,TRUNOF,POND2,-CUMVAP,ERROR
!    call output_waterbal()
    RAIN=0.0
    DPCAN=0.0
    TRUNOF=0.0
    TETSUM=0.0
    TEVAP=0.0
    TPERC=0.0
    DCAN=0.0
    DSNOW=0.0
    DRES=0.0
    DSOIL=0.0
    RETURN
 END SUBROUTINE WBALNC



! no problem checked changed
!***********************************************************************
!
 SUBROUTINE RTDIST (NPLANT,ITYPE,NS,ZS,ROOTDP,RROOT0,col,row)
!THIS PROGRAM CALCULATES THE ROOT RESISTANCES AND FRACTION OF ROOTS
!WITHIN EACH SOIL LAYER GIVEN THE MAXIMUM ROOTING DEPTH.  A
!TRIANGULAR ROOTING DENSITY IS ASSUMED HAVING A MAXIMUM DENSITY AT
!A DEPTH OF ZMXDEN AND ZERO DENSITY AT THE SURFACE AND THE MAXIMUM
!ROOTING DEPTH.  ZMXDEN IS ASSUMED TO BE A CONSTANT
!FRACTION (RMXDEN) OF THE MAXIMUM ROOTING DEPTH.
!(CURRENTLY RMXDEN IS SET IS DATA STATEMENT)
!***********************************************************************
    use clayrs_mod,only:TOTROT2d,TOTLAI2d,ROOTDN2d,RROOT2d
    implicit none
!input
    integer(i4),intent(in)::NPLANT,NS,col,row
    real(r8),dimension(NSMAX),intent(in)::ZS
    real(r8),dimension(NPMAX),intent(in)::ROOTDP,RROOT0
    integer(i4),dimension(NPMAX),intent(in)::ITYPE
!temp
    real(r8)::RMXDEN,ZMXDEN,ZMID1,ZMID2,AREA1,AREA2
    integer(i4)::I,J

    real(r8),dimension(:),pointer::TOTROT,TOTLAI
    real(r8),dimension(:,:),pointer::ROOTDN,RROOT

    TOTROT=>TOTROT2d(col,row,:)
    TOTLAI=>TOTLAI2d(col,row,:)
    ROOTDN=>ROOTDN2d(col,row,:,:)
    RROOT=>RROOT2d(col,row,:,:)


    RMXDEN=0.1
      DO 30 J=1,NPLANT
        TOTROT(J)=0.0
        IF (ITYPE(J).NE.0 .AND. TOTLAI(J).NE.0.0) THEN
!         TRANSPIRING PLANT -- CALCULATE FRACTION OF ROOTS IN EACH
!         SOIL LAYER; START BY COMPUTING DEPTH OF MAXIMUM ROOT DENSITY
          ZMXDEN=RMXDEN*ROOTDP(J)
!
          ZMID1=0.0
          DO 10 I=1,NS
!           CALCULATE MID-POINT BETWEEN THIS AND NEXT NODE, I.E. THE
!           LOWER BOUNDARY OF THIS LAYER
            IF (I.LT.NS) THEN
               ZMID2=(ZS(I+1)+ZS(I))/2.
              ELSE
               ZMID2=ZS(NS)
            END IF
!
            IF (ZMID2.LT.ZMXDEN) THEN
!             BOTTOM OF LAYER IS LESS THAN DEPTH OF MAXIMUM DENSITY
              ROOTDN(J,I)=(ZMID2-ZMID1)*(ZMID2+ZMID1)/ZMXDEN &
                /ROOTDP(J)
             ELSE
!             BOTTOM OF LAYER IS BEYOND DEPTH OF MAXIMUM DENSITY
              IF (ZMID2.LT.ROOTDP(J)) THEN
!               BOTTOM OF LAYER IS WITHIN ROOTING DEPTH
                IF (ZMID1.LT.ZMXDEN) THEN
!                 LAYER STRATTLES DEPTH OF MAXIMUM DENSITY
                  AREA1=(ZMXDEN-ZMID1)*(ZMXDEN+ZMID1)/ZMXDEN
                  AREA2=(ZMID2-ZMXDEN) &
                    *(2.-(ZMID2-ZMXDEN)/(ROOTDP(J)-ZMXDEN))
                  ROOTDN(J,I)=(AREA1+AREA2)/ROOTDP(J)
                 ELSE
!                 LAYER IS BEYOND DEPTH OF MAXIMUM ROOTING DENSITY BUT
!                 IS FULLY WITHIN THE ROOTING DEPTH
                  ROOTDN(J,I)=(ZMID2-ZMID1)*(2.*ROOTDP(J)-ZMID2-ZMID1) &
                    /(ROOTDP(J)-ZMXDEN)/ROOTDP(J)
                END IF
               ELSE

!               BOTTOM OF LAYER IS BEYOND ROOTING DEPTH
                IF (ZMID1.LT.ROOTDP(J)) THEN
!                 TOP OF LAYER IS STILL WITHIN ROOTING DEPTH; THE
!                 REMAINING FRACTION OF ROOTS ARE WITHIN THIS LAYER
                  ROOTDN(J,I)=1.0-TOTROT(J)
                 ELSE
!                 LAYER IS BEYOND THE ROOTING DEPTH -- NO ROOTS
                  ROOTDN(J,I)=0.0
                END IF
              END IF
            END IF
!           SUM THE TOTAL FRACTION OF ROOTS
            TOTROT(J)=TOTROT(J) + ROOTDN(J,I)
            ZMID1=ZMID2
   10     CONTINUE
!         CALCULATE EFFECTIVE ROOT CONDUCTANCE FOR EACH SOIL LAYER
          DO 20 I=1,NS
             RROOT(J,I)=RROOT0(J)*ROOTDN(J,I)/TOTROT(J)
   20     CONTINUE
        END IF
   30 CONTINUE
!
      RETURN
 END SUBROUTINE RTDIST



! no problem checked changed
!***********************************************************************
!
 SUBROUTINE CANLAY (NC,NPLANT,ITYPE,ZC,WCANDT,TCDT,TLCDT,VAPCDT,PLANTZ,PLANTW,PLTLAI,RLEAF0,col,row)
!
!      This subroutine splits the canopy into layers.  Each layer
!      will have the same total leaf area index, but the actual layer
!      dimension will vary.  The different plant variables are
!      dimensioned by plant type and layer, and their values are
!      apportioned according to their representation in the layer.
!
!***********************************************************************
    use clayrs_mod,only:TOTLAI2d,DRYCAN2d,CANLAI2d,RLEAF2d
    implicit none
!input
    integer(i4),intent(in)::NPLANT,col,row
    integer(i4),intent(inout)::NC
    real(r8),dimension(NCMAX),intent(inout)::ZC
    real(r8),dimension(NCMAX-1),intent(inout)::WCANDT
    real(r8),dimension(NCMAX),intent(inout)::TCDT
    real(r8),dimension(NCMAX),intent(inout)::VAPCDT
    real(r8),dimension(NPMAX,NCMAX-1),intent(inout)::TLCDT
    real(r8),dimension(NPMAX),intent(in)::PLANTZ,PLANTW,PLTLAI,RLEAF0
    integer(i4),dimension(NPMAX),intent(in)::ITYPE
!temp
    real(r8)::LAISUM,LAILYR,LYRTOP,LYRTOT,DZMAX,A, FACTOR,TMPLAI
    real(r8),dimension(NCMAX-1)::TMP,DZ
    integer(i4)::HTINDX,J,LAYER,OLDNC,NZEROS
    integer(i4),dimension(NPMAX)::HTORDR
!-- Local variables:
!A         Leaf area index per layer.
!DZMAX     Maximum theoretical layer thickness for current values.
!HTINDX    Index for ordered plant height - in increasing order.
!HTORDR    Array of plant numbers ordered by height (1 = smallest)
!LAILYR    Current accumulated LAI for current layer.
!LAISUM    The delta LAI/unit distance for current calculations.
!LAYER     Current layer number, used for loop control.
!LYRTOP    Height of top of current layer.
!LYRTOT    Current accumulated layer thickness for current layer.
!TMPLAI    Sum of total leaf area index for all plants.

    real(r8),dimension(:),pointer::TOTLAI
    real(r8),dimension(:,:),pointer::DRYCAN,CANLAI,RLEAF

    TOTLAI=>TOTLAI2d(col,row,:)
    DRYCAN=>DRYCAN2d(col,row,:,:)
    CANLAI=>CANLAI2d(col,row,:,:)
    RLEAF=>RLEAF2d(col,row,:,:)

!**   --Initialization
      LAILYR = 0.0
      LYRTOT = 0.0
      LAISUM = 0.0
      TMPLAI = 0.0
      NZEROS = 0
      LYRTOP = 0.0
      OLDNC = NC

!     -- Initialize temporary variable used for sorting,
!        sum of total leaf area index for all plants,
!        and sum of delta LAI per unit distance.
      DO 10 J=1, NPLANT
         IF (PLANTZ(J)*PLTLAI(J) .NE. 0.0) THEN
            TMP(J) = PLANTZ(J)
            TOTLAI(J) = PLTLAI(J)
            if(TOTLAI(J) .lt. 0.1) TOTLAI(J)=0.1
            TMPLAI = TMPLAI + TOTLAI(J)
            LAISUM = LAISUM + TOTLAI(J) / PLANTZ(J)
           ELSE
            TMP(J) = 9999.
            NZEROS=NZEROS+1
            HTORDR(NZEROS)=J
            TOTLAI(J) = 0.0
         END IF
   10 CONTINUE

!**   Arrange the plants by height - in increasing order
      DO 40 HTINDX = NZEROS+1, NPLANT
         HTORDR(HTINDX) = 1
         DO 20 J = 1, NPLANT
            IF (TMP(HTORDR(HTINDX)).GT.TMP(J)) HTORDR(HTINDX) = J
   20    CONTINUE
         TMP(HTORDR(HTINDX)) = 9999.0
   40 CONTINUE
!
!xxxx  NC = NINT (TMPLAI / 0.75)
       NC = NINT (TMPLAI / 0.50)
!     -- Set upper limit for number of canopy layers.
      IF (NC.EQ.0) THEN
         NC=1
        ELSE
         IF (NC.GT.10) NC = 10
      END IF
!
      HTINDX = NZEROS+1
!
      A  = TMPLAI / NC
      LAYER  = NC
      ZC(NC+1) = PLANTZ(HTORDR(NPLANT))


!**   --Main loop   ***g

!     --Calculate theoretical maximum layer based on current LAISUM.
  100 IF (LAYER .NE. 1) THEN
         DZMAX = (A - LAILYR) / LAISUM
        ELSE
         DZMAX = PLANTZ(HTORDR(NPLANT)) - LYRTOP
      END IF
!
      IF(LYRTOP+DZMAX.LE.PLANTZ(HTORDR(HTINDX)).OR.LAYER.EQ.1)THEN
!**   -- Top of layer is below or equal to next tallest plant.

         LYRTOP = LYRTOP + DZMAX
         DZ(LAYER) = LYRTOT + DZMAX

!**      --Apportioning routine start.
         DO 110 J = 1, NPLANT
           IF (PLANTZ(J).GE.LYRTOP) THEN
!               -- Plant fully contained in layer, full apportioning.
                FACTOR = DZ(LAYER) / PLANTZ(J)
           ELSE IF (PLANTZ(J).GT.LYRTOP-DZ(LAYER)) THEN
!          -- Plant partially contained in layer, partial apportioning.
                FACTOR = (PLANTZ(J)-(LYRTOP-DZ(LAYER))) / PLANTZ(J)
           ELSE
!          -- Plant is not contained in layer, set values to zero.
                FACTOR = 0.0
           END IF

           CANLAI(J,LAYER) = TOTLAI(J) * FACTOR
           DRYCAN(J,LAYER) = PLANTW(J) * FACTOR
           IF (ITYPE(J).NE.0) THEN
!          -- If plant is not dead.
              IF (FACTOR*TOTLAI(J).NE.0) THEN
!                -- Plant is within layer.
                 RLEAF(J,LAYER) = RLEAF0(J)*(CANLAI(J,LAYER)/TOTLAI(J))
               ELSE
!                -- Plant is outside of layer.
                 RLEAF(J,LAYER) = 0.0
              END IF
            ELSE
!             -- Plant is dead.
              RLEAF (J,LAYER) = 0.0
           END IF

!**        -- End apportioning routine.
  110    CONTINUE

!          --Update variables if more layers were added.
           IF (OLDNC.LT.LAYER) THEN
                TCDT(LAYER) = TCDT(OLDNC)
                VAPCDT(LAYER) = VAPCDT(OLDNC)
                WCANDT(LAYER) = WCANDT(OLDNC)
                DO 120 J=1,NPLANT
                   TLCDT(J,LAYER)=TLCDT(J,OLDNC)
  120           CONTINUE
           END IF

!          --Calc midpoint of the layer
           IF (LAYER.EQ.1) THEN
               ZC(1) = 0.0
           ELSE
               ZC(LAYER) = PLANTZ(HTORDR(NPLANT))-(LYRTOP-DZ(LAYER)/2.)
           END IF

!          Check if done.  If so, drop out of loop
           IF (LAYER.EQ.1) THEN
                IF (ABS(PLANTZ(HTORDR(NPLANT))-LYRTOP).GT.0.0001) THEN
                     WRITE (6,*)'ERROR MATCHING TOP OF CANOPY'
                     WRITE (6,*)'  LYRTOP =  ', LYRTOP
                END IF
                GO TO 190
           END IF

!**        -- If top of layer equaled top of the next tallest plant
           IF (LYRTOP .GE. PLANTZ(HTORDR(HTINDX))-.0001) THEN
                LAISUM = LAISUM - &
                    TOTLAI(HTORDR(HTINDX))/PLANTZ(HTORDR(HTINDX))
                HTINDX = HTINDX + 1
           END IF

           LYRTOT = 0.0
           LAILYR = 0.0
           LAYER = LAYER - 1

      ELSE
!**   -- Top of layer extends beyond top of plant being considered.
!          -- Update calcs based on limitation of top of next plant.
           LYRTOT = LYRTOT + PLANTZ(HTORDR(HTINDX))-LYRTOP
           LAILYR = LAILYR + (A-LAILYR) * &
                ((PLANTZ(HTORDR(HTINDX))-LYRTOP)/DZMAX)
           LAISUM = LAISUM - &
                    TOTLAI(HTORDR(HTINDX))/PLANTZ(HTORDR(HTINDX))
           LYRTOP = PLANTZ(HTORDR(HTINDX))
           HTINDX = HTINDX + 1
      END IF
      GO TO 100
!     -- End of loop

  190 CONTINUE
      RETURN
 END SUBROUTINE CANLAY


! no problem checked changed
!***********************************************************************
 SUBROUTINE CANOPY (NPLANT,NC,NSP,NSPLST,MZCINP,ITYPE,INITAL,PLTHGT,PLTWGT,PLTLAI,RLEAF0,ZSP,ZC,TC,TCDT,TLC,&
                    TLCDT,VAPC,VAPCDT,WCAN,WCANDT,WCMAX,PCANDT,CANMA,CANMB,TA,HUMID,col,row)
! THIS SUBROUTINE DETERMINES WHICH LAYERS OF THE CANOPY ARE COVERED
! WITH SNOW
!***********************************************************************
    use clayrs_mod,only:TOTLAI2d,TOTROT2d,CANLAI2d,ZZC2d
    use writeit_mod,only:nnc2d
    implicit none
!      COMMON /CLAYRS/ DRYCAN(8,10),CANLAI(8,10),TOTLAI(8),IEVAP(8),
!     >                RLEAF(8,10),RROOT(8,50),ROOTDN(8,50),TOTROT(8)
!      common /writeit/ hflux1,rh,hnc,xlenc,contk,tleaf(8,10),nnc

!input
    integer(i4),intent(in)::NPLANT,NSP,NSPLST,col,row,MZCINP,INITAL
    integer(i4),intent(inout)::NC
    integer(i4),dimension(NPMAX),intent(in)::ITYPE
    real(r8),dimension(NPMAX),intent(in)::PLTHGT,PLTWGT,PLTLAI,RLEAF0
    real(r8),dimension(NPMAX),intent(inout)::PCANDT
    real(r8),dimension(NSPMAX),intent(in)::ZSP
    real(r8),dimension(NCMAX),intent(inout)::ZC
    real(r8),dimension(NCMAX),intent(inout)::TCDT,TC
    real(r8),dimension(NPMAX,NCMAX-1),intent(inout)::TLC,TLCDT
    real(r8),dimension(NCMAX),intent(inout)::VAPC,VAPCDT
    real(r8),dimension(NCMAX-1),intent(inout)::WCAN,WCANDT
    real(r8),intent(in)::WCMAX,CANMA,CANMB,TA
    real(r8),intent(inout)::HUMID
!temp
    integer(i4)::I,J,K,NCCHK,NCLAST,NC1
    real(r8)::DUMMY,ZMID1,ZMID2,SATV
    real(r8),dimension(NPMAX)::PLANTZ,PLANTW

    real(r8),dimension(:,:),pointer::CANLAI
    integer(i4),pointer::nnc
    real(r8),dimension(:),pointer::TOTLAI,TOTROT
    real(r8),dimension(:),pointer::ZZC


    nnc=>nnc2d(col,row)
    TOTLAI=>TOTLAI2d(col,row,:)
    TOTROT=>TOTROT2d(col,row,:)
    CANLAI=>CANLAI2d(col,row,:,:)
    ZZC=>ZZC2d(col,row,:)

      IF (INITAL .EQ. 0) THEN
         IF (MZCINP .EQ. 1) THEN
!           LAYERING OF CANOPY NODES SPECIFIED BY USER -- SAVE SPACING
            NNC=NC
            ZZC(1)=0.0
            DO 5 I=2,NC
               ZZC(I)=ZC(I)
               TCDT(I)=TCDT(I-1)
               VAPCDT(I)=VAPCDT(I-1)
               DO 4 J=1,NPLANT
                  TLCDT(J,I)=TCDT(I)
    4          CONTINUE
    5       CONTINUE
            ZZC(NC+1)=ZC(NC+1)
            TCDT(NC+1)=TCDT(NC)
            VAPCDT(NC+1)=VAPCDT(NC)
!           IF NO NEED TO ADJUST FOR SNOWPACK, RETURN TO MAIN PROGRAM
            IF (NSP .EQ. 0) GO TO 80
         END IF
      END IF
!
!     CHECK IF CANOPY HAS EMERGED
      NCCHK=0
      DO 8 J=1,NPLANT
         IF (PLTLAI(J).GT.0) NCCHK=1
    8 CONTINUE
      IF (NCCHK.EQ.0) THEN
         NC=0
         GO TO 80
      END IF
!
!     CANOPY IS PRESENT - CHECK IF COVERED BY SNOW
      NCLAST=NC
!
      IF (NCLAST .EQ. 0) THEN
!        NO CANOPY LAST TIME STEP
!        DEFINE STATE VARIABLES FOR TOP OF CANOPY
         TC(1)=TA
         TCDT(1)=TA
         DO 9 J=1,NPLANT
            TLCDT(J,1)=TCDT(1)
    9    CONTINUE
         if(abs(TA) .gt. 120.0 )then
            print*,"In CANOPY,TA",col,row,TA
            stop
         endif
         CALL VSLOPE (DUMMY,SATV,TA)
         VAPC(1)=HUMID*SATV
         VAPCDT(1)=HUMID*SATV
         IF (NSPLST.GT.0) THEN
!           SNOW MELTED AND EXPOSED TOP OF CANOPY--SET MAX WATER CONTENT
            WCAN(1)=WCMAX
            WCANDT(1)=WCMAX
           ELSE
!           SET WATER CONTENT TO EQUILIBRIUM WITH HUMIDITY
            CALL CANHUM (2,HUMID,DUMMY,WCANDT(1),TCDT(1),CANMA,CANMB)
            WCAN(1)=WCANDT(1)
         END IF
         NC=1
      END IF
!
      IF (NSP .GT. 0) THEN
!        DETERMINE IF CANOPY IS COVERED WITH SNOW
         DO 10 J=1,NPLANT
            IF (PLTHGT(J).GT.ZSP(NSP+1) .AND. PLTLAI(J).GT.0.0) GO TO 15
   10    CONTINUE
!        SNOW COVERS CANOPY
         NC=0
         GO TO 80
!
   15    IF (MZCINP .EQ. 1) GO TO 25
!        ALLOW MODEL TO REDEFINE LAYERING AND NODES WITHIN CANOPY
         DO 20 J=1,NPLANT
            PLANTZ(J)=PLTHGT(J)-ZSP(NSP+1)
            IF (PLANTZ(J).LT.0.0) THEN
               PLANTZ(J)=0.0
               PLANTW(J)=0.0
               TOTLAI(J)=0.0
              ELSE
               PLANTW(J)=PLTWGT(J)*PLANTZ(J)/PLTHGT(J)
               TOTLAI(J)=PLTLAI(J)*PLANTZ(J)/PLTHGT(J)
            END IF
   20    CONTINUE
!
!        DETERMINE LAYERING OF CANOPY ACCOUNTING FOR SNOW
         CALL CANLAY (NC,NPLANT,ITYPE,ZC,WCANDT,TCDT,TLCDT,VAPCDT, &
            PLANTZ,PLANTW,TOTLAI,RLEAF0,col,row)
         GO TO 55
!
!        CANOPY LAYERING WAS INPUT - DO NOT ALLOW MOVING OF NODES
   25    DO 30 I=NNC,1,-1
            IF ((ZZC(NNC+1)-ZZC(I)) .GT. ZSP(NSP+1)) THEN
               NC=I
               ZC(NC+1) = ZZC(NNC+1)-ZSP(NSP+1)
               IF (NC .EQ. 1) THEN
                  ZC(NC)=0.0
                  GO TO 35
               END IF
               ZMID1=(ZZC(NC)+ZZC(NC-1))/2
               ZMID2=(ZZC(NC+1)+ZZC(NC))/2
               IF (ZC(NC+1) .GT. ZMID2) THEN
                  ZC(NC)=ZZC(NC)
                 ELSE
                  ZC(NC)=ZZC(NC) - (ZZC(NC)-ZMID1)*(ZMID2-ZC(NC+1)) &
                         /(ZMID2-ZZC(NC))
               END IF
               GO TO 35
            END IF
   30    CONTINUE
        ELSE
!
!        NO SNOW ON GROUND
         IF (MZCINP .EQ. 1) THEN
!           RESET NUMBER OF NODES AND POSITION OF BOTTOM TWO NODES
!           (POSITION OF OTHER NODES RESET BELOW)
            NC=NNC
            ZC(NC)=ZZC(NC)
            ZC(NC+1)=ZZC(NC+1)
           ELSE
!           DETERMINE LAYERING OF CANOPY WITH NO SNOWPACK (ADJUSTING FOR
!           ANY CHANGES IN CANOPY CHARACTERISTICS)
            CALL CANLAY (NC,NPLANT,ITYPE,ZC,WCANDT,TCDT,TLCDT,VAPCDT, &
                    PLTHGT,PLTWGT,PLTLAI,RLEAF0,col,row)
         END IF
      END IF
!
   35 IF (MZCINP.EQ.1) THEN
!        CALCULATE TOTAL LEAF AREA INDEX ABOVE SNOW FOR EACH PLANT TYPE
         IF (NC .NE. NCLAST) THEN
            DO 45 J=1,NPLANT
               TOTLAI(J)=0.0
               DO 40 I=1,NC
                  TOTLAI(J)=TOTLAI(J)+CANLAI(J,I)
   40          CONTINUE
   45       CONTINUE
!
            IF (NC .GT. NCLAST) THEN
!              SNOW HAS MELTED AND EXPOSED ADDITIONAL CANOPY LAYERS --
!              DEFINE STATE VARIABLES OF NEWLY EXPOSED LAYERS
               NC1=NCLAST
               IF (NCLAST .EQ. 0) NC1=1
               DO 50 I=NC1+1,NC
                  TCDT(I)=TCDT(I-1)
                  VAPCDT(I)=VAPCDT(I-1)
                  ZC(I-1)=ZZC(I-1)
                  DO 49 J=1,NPLANT
                     TLCDT(J,I)=TLCDT(J,I-1)
   49             CONTINUE
   50          CONTINUE
            END IF
!
         END IF
      END IF
   55 IF (NC .GT. NCLAST) THEN
!        SET WATER CONTENT, TEMPERATURE AND VAPOR OF NEW LAYERS
!        FOR BEGINNING OF TIME STEP
         NC1=NCLAST
         IF (NCLAST .EQ. 0) NC1=1
         DO 60 I=NC1+1,NC
            IF (NSPLST .GT. 0) THEN
!              SNOW HAS MELTED EXPOSING LAYERS - WATER CONTENT = MAXIMUM
               WCANDT(I)=WCMAX
              ELSE
               WCANDT(I)=WCANDT(I-1)
            END IF
            WCAN(I)=WCANDT(I)
            TC(I)=TCDT(I)
            VAPC(I)=VAPCDT(I)
            DO 59 J=1,NPLANT
               TLC(J,I)=TLCDT(J,I)
   59       CONTINUE
   60    CONTINUE
      END IF
!
   80 DO 85 J=1,NPLANT
         IF (PLTLAI(J).LE. 0.0 .OR. TOTLAI(J) .LE. 0.0) PCANDT(J)=0.0
   85 CONTINUE
!xxxx
      if (mzcinp .eq. 0) nnc=nc
      RETURN
 END SUBROUTINE CANOPY


! no problem changed
!***********************************************************************
!
 SUBROUTINE CANHUM (NCALC,HUM,DHDW,WCAN,TC,CANMA,CANMB)
!
!     THIS SUBROUTINE DEFINES THE RELATION BETWEEN THE HUMIDITY AND THE
!     MOISTURE CONTENT OF THE DEAD CANOPY MATERIAL USING THE CONVERSION
!     BETWEEN WATER POTENTIAL AND HUMIDITY.  DEPENDING ON THE VALUE OF
!     NCALC, THE SUBROUTINE WILL EITHER CALCULATE HUMIDITY AND THE
!     DERIVATIVE OF HUMIDTY WITH RESPECT TO WATER CONTENT FOR A GIVEN
!     WATER CONTENT, OR AN EQUILIBRIUM WATER CONTENT GIVEN THE HUMIDITY.
!     (NCALC = 1: CALCULATE HUMIDITY; OTHERWISE CALCULATE WATER CONTENT)
!***********************************************************************
    use constvar_mod,only:G,UGAS
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4
    implicit none

    integer(i4),intent(in)::NCALC
    real(r8),intent(in)::TC,CANMA,CANMB
    real(r8),intent(inout)::HUM,DHDW,WCAN


      IF (NCALC .EQ. 1) THEN
!        CALCULATE HUMIDITY AND SLOPE BASED ON WATER CONTENT
         IF (WCAN .LE. 0.0) THEN
!           EQUATIONS CANNOT HANDLE CASE WHERE WCAN=0.0 -- SET HUM=0.0
            HUM=0.0
            DHDW=0.0
           ELSE
!
            HUM=EXP(0.018*G/(UGAS*(TC+273.16))*CANMA*WCAN**(-CANMB))
            DHDW=-0.018*G*CANMA*CANMB*HUM*WCAN**(-CANMB-1.)&
                /(UGAS*(TC+273.16))
         END IF
        ELSE
!        CALCULATE WATER CONTENT BASED ON HUMIDITY
         IF (HUM .LT. 0.999) THEN
          IF (HUM .LE. 0.0) THEN
            WCAN = 0.0
           ELSE
            WCAN=(log(HUM)*UGAS*(TC+273.16)/0.018/G/CANMA)**(-1./CANMB)
          END IF
         ELSE
          WCAN=(log(0.999)*UGAS*(TC+273.16)/0.018/G/CANMA)**(-1./CANMB)
         END IF
      END IF
      RETURN
 END SUBROUTINE


! no problem  changed
!***********************************************************************
!
 SUBROUTINE RESHUM (NCALC,HUM,DHDW,GMC,TR)
!
!     THIS SUBROUTINE DEFINES THE RELATION BETWEEN THE HUMIDITY AND THE
!     MOISTURE CONTENT OF THE RESIDUE USING THE CONVERSION
!     BETWEEN WATER POTENTIAL AND HUMIDITY.  DEPENDING ON THE VALUE OF
!     NCALC, THE SUBROUTINE WILL EITHER CALCULATE HUMIDITY AND THE
!     DERIVATIVE OF HUMIDTY WITH RESPECT TO WATER CONTENT FOR A GIVEN
!     WATER CONTENT, OR AN EQUILIBRIUM WATER CONTENT GIVEN THE HUMIDITY.
!     (NCALC = 1: CALCULATE HUMIDITY; OTHERWISE CALCULATE WATER CONTENT)
!***********************************************************************
    use constvar_mod,only:G,RESMA,RESMB,UGAS
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4
!    use rsparm_mod
    implicit none
    integer(i4),intent(in)::NCALC
    real(r8),intent(inout)::HUM,DHDW,GMC
    real(r8),intent(in)::TR


!**** THIS SUBROUTINE USES FOLLOWING RELATION:
!                TENSION = A*W**B
!                H=EXP(M*TENSION/(UGAS*TEMP))
!     THEREFORE  H=EXP(M*A*W**B/(UGAS*TEMP))
!
      IF (NCALC .EQ. 1) THEN
!        CALCULATE HUMIDITY AND SLOPE BASED ON WATER CONTENT
         IF (GMC .LE. 0.0) THEN
!           EQUATIONS CANNOT HANDLE CASE WHERE GMC=0.0 -- SET HUM=0.0
            HUM=0.0
            DHDW=0.0
           ELSE
            HUM=EXP(0.018*G/(UGAS*(TR+273.16))*RESMA*GMC**(-RESMB))
            DHDW=-0.018*G*RESMA*RESMB*HUM*GMC**(-RESMB-1.) &
                /(UGAS*(TR+273.16))
         END IF
        ELSE
!        CALCULATE WATER CONTENT BASED ON HUMIDITY
         IF (HUM .GT. 0.999) THEN
           GMC=(log(0.999)*UGAS*(TR+273.16)/0.018/G/RESMA)**(-1./RESMB)
          ELSE
           GMC=(log(HUM)*UGAS*(TR+273.16)/0.018/G/RESMA)**(-1./RESMB)
         END IF
      END IF
      RETURN
 END SUBROUTINE


!no problem changed
!***********************************************************************
!
 SUBROUTINE VSLOPE (S,SATV,T)
!
!     THIS SUBROUTINE CALCULATES THE SATURATED VAPOR DENSITY AND THE
!     SLOPE OF THE VAPOR DENSITY CURVE  (SATV IS IN KG/M**3)
!
!***********************************************************************
    use constvar_mod,only:UGAS
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4
    implicit none
    real(r8),intent(in)::T
    real(r8),intent(inout)::S,SATV
    real(r8)::TMP

      TMP=T+273.16
      SATV=EXP(52.57633-6790.4985/TMP-5.02808*log(TMP))*1000./(UGAS*TMP/.018)
      S=.0000165 + 4944.43*SATV/(TMP**2)
      RETURN
 END SUBROUTINE VSLOPE





! no problem changed
!***********************************************************************
!
 SUBROUTINE FROZEN (I,VLC,VIC,MAT,CONC,TS,SALT,ICES,nsalt,col,row)
!
!     THIS SUBROUTINE DETERMINES THE LIQUID AND ICE WATER CONTENT OF THE
!     SOIL FOR TMP < 0 C  AND DETERMINES IF ANY ICE MAY BE PRESENT.
!     IT THEN UPDATES THE MATRIC POTENTIAL AND SOLUTE CONCENTRATION.
!
!***********************************************************************
    use dims_mod
    use constvar_mod,only:G,LF,RHOI,RHOL,UGAS
    use soilproperty_mod,only:SALTKQ2d,RHOB2d,SAT2d
    implicit none
    integer(i4),intent(in)::I,col,row,nsalt
    real(r8),dimension(NSMAX),intent(inout)::VLC,VIC,MAT
    real(r8),dimension(NSALTMAX,NSMAX),intent(inout)::CONC,SALT
    real(r8),dimension(NSMAX)::TS
    integer(i4),dimension(NSMAX),intent(inout)::ICES
    real(r8)::MAT1

!   soil property pointer
    real(r8),dimension(:,:),pointer::SALTKQ
    real(r8),dimension(:),pointer::RHOB,SAT

    integer(i4)::ICONV,ITER,J
    real(r8)::tmp,TOTPOT,VLC1,OSMPOT,DERIV,TERM,ERROR,DELTA,DLDM,VLC2


    SALTKQ=>SALTKQ2d(col,row,:,:)
    RHOB=>RHOB2d(col,row,:)
    SAT=>SAT2d(col,row,:)
!     ITERATE TO FIND THE MAXIMUM WATER CONTENT FROM THE TOTAL, MATRIC
!     AND OSMOTIC POTENTIALS:
!               TOTAL - MATRIC - OSMOTIC = 0 => F
!     MATRIC AND OSMOTIC POTENTIALS ARE FUNCTIONS OF WATER CONTENT.
!     SOLVE BY NEWTON-RAPHSON ITERATION  =>  D(VLC) = -F/(DF/D(VLC)
!
      ICONV=0
    5 TMP=TS(I)+273.16
      TOTPOT = LF*TS(I)/TMP/G
      CALL MATVL2 (I,TOTPOT,VLC1,col,row)
!
!**** BEGIN ITERATIONS
      ITER=0
!     CALCULATE OSMOTIC POTENTIAL AND DERIVATIVE WITH RESPECT TO LIQUID
   10 OSMPOT=0.0
      DERIV=0.0
      DO 20 J=1,NSALT
         TERM = SALT(J,I)*UGAS*TMP/(SALTKQ(J,I)+VLC1*RHOL/RHOB(I))
         OSMPOT= OSMPOT-TERM
         DERIV= DERIV+TERM*RHOL/RHOB(I)/(SALTKQ(J,I)+VLC1*RHOL/RHOB(I))
   20 CONTINUE
      OSMPOT= OSMPOT/G
      DERIV= DERIV/G
!
!     CALCULATE MATRIC POTENTIAL AND DERIVATIVE WITH RESPECT TO LIQUID
      CALL MATVL1 (I,MAT1,VLC1,col,row)
      CALL MATVL3 (I,MAT1,VLC1,DLDM,col,row)
!
!     DERTERMINE ERROR IN WATER POTENTIAL AND ADJUST WATER CONTENT
      ERROR = TOTPOT - OSMPOT - MAT1
      DELTA = ERROR/(DERIV + 1/DLDM)
      VLC2 = VLC1 + DELTA
      IF (VLC2 .GT. SAT(I)) VLC2=SAT(I)
      IF (VLC2 .LE. VLC1/2) VLC2=VLC1/2.
      DELTA = VLC2 - VLC1
      VLC1=VLC2
      IF (ABS(DELTA) .LT. .00001) GOTO 30
      ITER=ITER+1
      IF (ICONV .GT. 0) &
        WRITE (21,*) ITER,OSMPOT,MAT1,VLC1,DERIV,1/DLDM,ERROR
      IF (ITER .GT. 30) THEN
         ICONV = ICONV+1
         WRITE (21,100)
         WRITE (21,*) I,TS(I),VLC(I),VIC(I),TOTPOT
         IF (ICONV .EQ. 1) GO TO 5
           print*,"stop from frozen"
         STOP
      END IF
      GO TO 10
!
!**** IF ACTUAL LIQUID PLUS ICE CONTENT IS LESS THAN THAT CALCULATED
!     FROM ABOVE, THERE IS NO ICE PRESENT IN THE LAYER.  (THE TOTAL
!     POTENTIAL IS SUFFICIENTLY NEGATIVE THAT THE WATER WILL NOT
!     FREEZE.)
   30 IF ((VLC(I)+VIC(I)*RHOI/RHOL).LT.VLC2) THEN
!*       NO ICE IS PRESENT AND LAYER IS UNSATURATED
         ICES(I)=0
         VLC(I)=VLC(I) + VIC(I)*RHOI/RHOL
         VIC(I)=0.0
         CALL MATVL1 (I,MAT(I),VLC(I),col,row)
        ELSE
!*       ICE MAY BE PRESENT; CONVERT ANY WATER PRESENT ABOVE THE MAXIMUM
!        LIQUID WATER CONTENT TO ICE
         VIC(I)=VIC(I)+(VLC(I)-VLC2)*RHOL/RHOI
         VLC(I)=VLC2
         IF (VIC(I) .GT. 0.0) THEN
!           ICE IS PRESENT --> VLC AND MAT ARE FUNCTION OF TEMP
            ICES(I)=1
!           COMPUTE MATRIC POTENTIAL FROM LIQUID WATER CONTENT
            CALL MATVL1 (I,MAT(I),VLC(I),col,row)
          ELSE
!           NO ICE IS PRESENT AND WATER CONTENT IS AT SATURATION
            VIC(I)=0.0
            ICES(I)=0
         END IF
      END IF
!     REDEFINE SOLUTE CONCENTRATIONS
      DO 40 J=1,NSALT
         CONC(J,I)=SALT(J,I)/(SALTKQ(J,I) + VLC(I)*RHOL/RHOB(I))
   40 CONTINUE
      RETURN
  100 FORMAT (//10X,'*****  PROGRAM WAS STOPPED IN SUBROUTINE FROZEN DUE &
        TO CONVERGENCE PROBLEMS *****')
 END SUBROUTINE


!no problem changed
!***********************************************************************
!
  SUBROUTINE MATVL1(I,MAT,VLC,col,row)
!
!     THIS SUBROUTINE DEFINES THE RELATION BETWEEN THE VOLUMETRIC
!     LIQUID CONTENT AND THE MATRIC POTENTIAL. THE SUBROUTINE IS DIVIDED
!     INTO THREE PART, AND THE OUTPUT DEPENDS ON WHICH PART IS CALLED
!           MATVL1 : THE MATRIC POTENTIAL IS CALCULATED FROM MOISTURE
!           MATVL2 : THE MOISTURE CONTENT IS CALCULATED FROM MATRIC
!           MATVL3 : THE DERIVATIVE OF MOISTURE CONTENT WITH RESPECT
!                    TO MATRIC POTENTIAL IS CALCULTED.
!                I = NODE NUMBER
!
!***********************************************************************
    use soilproperty_mod,only:B2d,SAT2d,ENTRY2d
    implicit none
    integer(i4),intent(in)::I,col,row
    real(r8),intent(in)::VLC
    real(r8),intent(inout)::MAT
!   pointer to soil property variables
    real(r8),dimension(:),pointer::B,SAT,ENTRY
    B=>B2d(col,row,:)
    SAT=>SAT2d(col,row,:)
    ENTRY=>ENTRY2d(col,row,:)
!   DETERMINE THE MATRIC POTENTIAL FROM THE MOISTURE CONTENT
      IF (VLC .LT. SAT(I)) THEN
         MAT=ENTRY(I)*(VLC/SAT(I))**(-B(I))
        ELSE
         MAT=ENTRY(I)
      END IF
      RETURN
  END SUBROUTINE


!no problem changed
!***********************************************************************
!
 SUBROUTINE MATVL2 (I,MAT,VLC,col,row)
!
!     THIS SUBROUTINE DEFINES THE RELATION BETWEEN THE VOLUMETRIC
!     LIQUID CONTENT AND THE MATRIC POTENTIAL. THE SUBROUTINE IS DIVIDED
!     INTO THREE PART, AND THE OUTPUT DEPENDS ON WHICH PART IS CALLED
!           MATVL1 : THE MATRIC POTENTIAL IS CALCULATED FROM MOISTURE
!           MATVL2 : THE MOISTURE CONTENT IS CALCULATED FROM MATRIC
!           MATVL3 : THE DERIVATIVE OF MOISTURE CONTENT WITH RESPECT
!                    TO MATRIC POTENTIAL IS CALCULTED.
!                I = NODE NUMBER
!***********************************************************************
    use soilproperty_mod,only:ENTRY2d,SAT2d,B2d
    implicit none
    integer(i4),intent(in)::I,col,row
    real(r8),intent(in)::MAT
    real(r8),intent(inout)::VLC
    real(r8),dimension(:),pointer::ENTRY,SAT,B

    ENTRY=>ENTRY2d(col,row,:)
    SAT=>SAT2d(col,row,:)
    B=>B2d(col,row,:)
!     DETERMINE THE MOISTURE CONTENT FROM THE MATRIC POTENTIAL
      IF (ENTRY(I) .GT. MAT) THEN
         VLC=SAT(I)*(MAT/ENTRY(I))**(-1./B(I))
       ELSE
         VLC=SAT(I)
      END IF
      RETURN
 END SUBROUTINE


!no problem changed
!***********************************************************************
!
 SUBROUTINE MATVL3 (I,MAT,VLC,DLDM,col,row)
!
!     THIS SUBROUTINE DEFINES THE RELATION BETWEEN THE VOLUMETRIC
!     LIQUID CONTENT AND THE MATRIC POTENTIAL. THE SUBROUTINE IS DIVIDED
!     INTO THREE PART, AND THE OUTPUT DEPENDS ON WHICH PART IS CALLED
!           MATVL1 : THE MATRIC POTENTIAL IS CALCULATED FROM MOISTURE
!           MATVL2 : THE MOISTURE CONTENT IS CALCULATED FROM MATRIC
!           MATVL3 : THE DERIVATIVE OF MOISTURE CONTENT WITH RESPECT
!                    TO MATRIC POTENTIAL IS CALCULTED.
!                I = NODE NUMBER
!
!***********************************************************************
    use soilproperty_mod,only:ENTRY2d,B2d
    implicit none
    integer(i4),intent(in)::I,col,row
    real(r8),intent(in)::MAT,VLC
    real(r8),intent(inout)::DLDM
    real(r8),dimension(:),pointer::ENTRY,B

    ENTRY=>ENTRY2d(col,row,:)
    B=>B2d(col,row,:)
!
!   DETERMINE THE DERIVATIVE OF THE MOISTURE CONTENT WITH RESPECT
!   TO MATRIC POTENTIAL
      IF (MAT .GT. ENTRY(I)) THEN
!        LAYER IS SATURATED
         DLDM=0.0
       ELSE
!        UNSATURATED CONDITIONS
         DLDM=-VLC/B(I)/MAT
      END IF
      RETURN
 END SUBROUTINE
 
 
!***********************************************************************
!
      SUBROUTINE FROST (NSP,NS,JULIAN,HOUR,YEAR,INITAL,ZSP,RHOSP,&
       DZSP,DLWDT,WLAG,STORE,ZS,VLCDT,VICDT,TSDT,ICESDT,col,row)
!
!     THIS SUBROUTINE INTERPOLATES BETWEEN NODES TO DETERMINE THE FROST
!     DEPTH, THEN PRINT THE FROST DEPTH AND SNOW DEPTH AT THE DESIRED
!     OUTPUT INTERVAL
!
!***********************************************************************
      use constvar_mod,only:G,LF,RHOI,RHOL,UGAS
      use savedata_mod,only:LAST2d,SWE2d,FDEPTH2d,TDEPTH2d
      implicit none
      integer(i4),intent(in)::NSP,NS,JULIAN,HOUR,YEAR,INITAL,col,row
      real(r8),dimension(NSPMAX),intent(in)::ZSP,RHOSP,DZSP,DLWDT
      real(r8),dimension(11),intent(in)::WLAG
      real(r8),intent(in)::STORE
      real(r8),dimension(NSMAX),intent(in)::ZS,VLCDT,VICDT,TSDT
      integer(r8),dimension(NSMAX),intent(in)::ICESDT
	  
	  
      integer(i4),pointer::NPRINT,LAST,LDAY,LHOUR,LYR
      real(r8),pointer::SWE,FDEPTH,TDEPTH
	  
      integer(i4)::NTHAW,NFROST,I,NF,NT
      real(r8)::FRACTN

      LAST=>LAST2d(col,row) 
      SWE=>SWE2d(col,row)	  
      FDEPTH=>FDEPTH2d(col,row)
      TDEPTH=>TDEPTH2d(col,row)
	 
!
      IF (INITAL .EQ. 0) THEN
!        LIMIT PRINT OUT OF ICE CONTENT TO TOP 15 SOIL NODES -- ANY MORE
!        COULD NOT FIT ON ONE LINE.

!        LAST => FLAG STATING WHETHER THERE WAS FROST LAST TIME STEP
!        LDAY,LHOUR,LYR => LAST TIME THAT FROST CHECKED
         LAST = 0
         FDEPTH = 0.0
         TDEPTH = 0.0
         GO TO 5
      END IF
!
!      IF (MOD(HOUR,LVLOUT) .NE. 0) RETURN
!
!     FIND LAYERS OF MAXIMUM THAW AND FROST
    5 NTHAW = 0
      NFROST = 0
      DO 10 I=NS,1,-1
         IF (TSDT(I) <0.0 ) THEN
            IF (NFROST .EQ. 0) NFROST = I
         ELSE
            IF (NFROST .GT. 0 .AND. NTHAW .EQ. 0) THEN
               NTHAW = I+1
               GO TO 15
            END IF
         END IF
   10 CONTINUE
!
   15 IF (NFROST .EQ. 0  .AND.  LAST .EQ. 0  .AND.  NSP .LE. 0) THEN
!        SAVE TIME ON WHICH FROST WAS LAST CHECKED AND RETURN
      END IF
!
!     IF LAST = 0, SNOW OR FROST WAS NOT PRESENT LAST TIME STEP BUT IS
!     PRESENT THIS TIME STEP --> PRINT OUT ZEROES FOR LAST TIME STEP
!
      IF (NFROST .GT. 0) THEN
!        CALCULATE DEPTH OF FROST
         LAST = 1
         NF = NFROST
         FRACTN = VICDT(NF)/(VLCDT(NF) + VICDT(NF))
         IF (NF .EQ. 1  .OR.  NF .EQ. NS) THEN
            IF (NF .EQ. 1) THEN
              FDEPTH=ZS(1)+TSDT(1)/(TSDT(1)-TSDT(2))*(ZS(2)-ZS(1))/2.0
              IF(FDEPTH>(ZS(2)+ZS(1))/2.) FDEPTH=(ZS(2)+ZS(1))/2.
             ELSE
              FDEPTH= ZS(NS)
            END IF
          ELSE
		    FDEPTH=(ZS(NF)+ZS(NF-1))/2.0-TSDT(NF)/(TSDT(NF)-TSDT(NF+1))*(ZS(NF+1)-ZS(NF-1))/2.
            IF(FDEPTH>(ZS(NF)+ZS(NF+1))/2.0) FDEPTH=(ZS(NF)+ZS(NF+1))/2.0
            IF(FDEPTH<(ZS(NF)+ZS(NF-1))/2.0) FDEPTH=(ZS(NF)+ZS(NF-1))/2.0
          END IF
!
         TDEPTH = 0.0
         IF (NTHAW .GT. 0) THEN
!           CALCULATE DEPTH OF THAW
            NT = NTHAW
            FRACTN = VLCDT(NT)/(VLCDT(NT) + VICDT(NT))
            IF (NT .EQ. 1) THEN
              TDEPTH=ZS(1)+TSDT(1)/(TSDT(1)-TSDT(2))*(ZS(2)-ZS(1))/2.0
              IF(TDEPTH>(ZS(2)+ZS(1))/2.) TDEPTH=(ZS(2)+ZS(1))/2.
             ELSE
              IF (NT .EQ. NFROST) THEN
               TDEPTH=(ZS(NT)+ZS(NT-1))/2. + &
                     FRACTN*(FDEPTH-(ZS(NT)+ZS(NT-1))/2.)
              ELSE
               TDEPTH=(ZS(NT)+ZS(NT-1))/2.+(TSDT(NT-1)/(TSDT(NT-1)-TSDT(NT)))*(ZS(NT+1)-ZS(NT-1))/2.
			   IF(TDEPTH<(ZS(NT)+ZS(NT-1))/2.) TDEPTH=(ZS(NT)+ZS(NT-1))/2.
			   IF(TDEPTH>(ZS(NT)+ZS(NT+1))/2.) TDEPTH=(ZS(NT)+ZS(NT+1))/2.
              END IF
            END IF
         END IF
!
       ELSE
!        NO FROST EXISTS IN THE PROFILE - IF FROST EXISTED AT THE LAST
!        OUTPUT, INDICATE THAT FROST HAS LEFT
         IF (LAST .GT. 0) THEN
!           FROST EXISTED IN PROFILE AT LAST OUTPUT -- ENSURE THAT
!           FROST DEPTH AND THAW DEPTH LINES COME TOGETHER
            IF (TDEPTH .EQ. 0.0) THEN
!              GROUND WAS THAWING FROM THE BOTTOM UP
               FDEPTH = 0.0
              ELSE
!              THAW AND FROST DEPTH COME TOGETHER SOMEWHERE BETWEEN THE
!              DEPTHS AT LAST OUTPUT - USUALLY 2/3 BETWEEN
               FDEPTH = (2.*FDEPTH + TDEPTH)/3.
               TDEPTH = FDEPTH
            END IF
          ELSE
!           NO FROST IN PROFILE AT LAST OUTPUT
            FDEPTH = 0.0
            TDEPTH = 0.0
         END IF
         LAST = 0
      END IF
!
      SWE = 0.0
      IF (NSP .LE. 0) THEN
        ELSE
!        SAVE SNOW DEPTH AND CALCULATE SNOW WATER EQUIVALENT
         DO 20 I=1,NSP
            SWE = SWE + RHOSP(I)*DZSP(I)/RHOL + DLWDT(I)
   20    CONTINUE
!        WATER IN PROCESS OF BEING LAGGED THROUGH THE SNOWPACK
         SWE = SWE + STORE
         DO 25 I=1,11
            SWE = SWE + WLAG(I)
   25    CONTINUE
      END IF
      RETURN
      END SUBROUTINE
!***********************************************************************

end module shaw27_mod
