module goshaw_mod
  use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
  use dims_mod
  use shaw27_mod
!  use statevar_mod,only:NC2d,NSP2d,NR2d,NS2d,INBASIN2d,PLTHGT2d,PLTWGT2d,PLTLAI2d,ROOTDP2d,DCHAR2d,TCCRIT2d,&
!	RSTOM02d,RSTEXP2d,PLEAF02d,RLEAF02d,RROOT02d,CANALB2d,XANGLE2d,CLUMPNG2d,PCANDT2d,ITYPE2d,ZC2d,TCDT2d,&
!	VAPCDT2d,TLCDT2d,WCANDT2d,ZSP2d,DZSP2d,RHOSP2d,TSPDT2d,DLWDT2d,ICESPT2d,WLAG2d,STORE2d,SNOWEX2d,RUNOFF2d,&
!	HFLUX2d,VFLUX2d,ZR2d,RHOR2d,TRDT2d,VAPRDT2d,GMCDT2d,RLOAD2d,ZRTHIK2d,COVER2d,ALBRES2d,RESCOF2d,GMCMAX,&
!	DIRRES2d,ZS2d,TSDT2d,VLCDT2d,VICDT2d,MATDT2d,CONCDT2d,SALTDT2d,ICESDT2d,ALBDRY2d,ALBEXP2d,DGRADE2d,SLTDIF2d,&
!	ASALT2d,DISPER2d,ZMSRF2d,ZHSRF2d,ZERSRF2d,ZMSP2d,ZHSP2d,POND2d
  implicit none
contains

  subroutine goshaw(col,row,julian,hour,year,nsalt,nc,nsp,nr,ns,mzcinp,inph2o,mwatrxt,ivlcbc,itmpbc,nplant,plthgt,&
    pltwgt,pltlai,rootdp,dchar,tccrit,rstom0,rstexp,pleaf0,rleaf0,rroot0,pcandt,canalb,xangle,clumpng,itype,&
    zc,tcdt,tlcdt,vapcdt,wcandt,zsp,dzsp,rhosp,tspdt,dlwdt,icespt,wlag,store,snowex,zr,rhor,trdt,vaprdt,gmcdt,&
    gmcmax,rload,zrthik,cover,albres,rescof,dirres,zs,tsdt,vlcdt,vicdt,matdt,concdt,icesdt,saltdt,albdry,albexp,&
    dgrade,sltdif,asalt,disper,zmsrf,zhsrf,zersrf,zmsp,zhsp,pond,pondmx,alatud,lontitud,elevation,slope,aspect,&
    hrnoon,clouds,declin,hafday,sunhor,tmpday,winday,humday,presur,shadow,skyview,precip,snoden,soitmp,vlcday,soilxt,inbasin,&
    thflux,vflux,runoff,tsavg,evap1,etsum,totflo,lwcan,swcan,lwsnow,lwsoil,lwres,swres,swsnow,swsoil,wwdt,maskflag,&
    maskgrid,bot,dgl,infiltration,AbsorbedSW,AbsorbedLW,snowmelt,InDirect,InDiffuse)
    use controlpara_mod,only:DT2d,LVLOUT,CANMA,CANMB,HEIGHT,WT2d,WDT2d,DTIME,INITAL,NHRPDT,SNOTMP,WCMAX,TOLER,MAXSTEP,shadeef,&
        maxiter,GroundTempGradients
    use constvar_mod,only:G,UGAS,LF,RHOL,RHOI
    use soilproperty_mod,only:B2d,SAT2d,RHOB2d,ENTRY2d,SAND2d,SILT2d,CLAY2d,OM2d,SATK2d,SALTKQ2d,thfc2d
    use matrix_mod
    use savedata_mod,only:MINSTP2d,MAXSTP2d
    implicit none
    integer(i4),intent(in)::col,row,JULIAN,HOUR,YEAR,NSALT,maskflag
    integer(i4),intent(inout)::NC,NSP
    integer(i4),intent(in)::NR,NS
    integer(i4),intent(in)::MZCINP,INPH2O,MWATRXT,IVLCBC,ITMPBC,NPLANT
    integer(i4),intent(inout)::INBASIN
! cannopy
    real(r8),dimension(NPMAX),intent(in)::PLTHGT,PLTWGT,PLTLAI,ROOTDP,DCHAR,TCCRIT,RSTOM0,RSTEXP,PLEAF0,RLEAF0,RROOT0,&
                                          CANALB,XANGLE,CLUMPNG
    integer(i4),dimension(NPMAX),intent(in)::ITYPE
    real(r8),dimension(NPMAX),intent(inout)::PCANDT
    real(r8),dimension(NCMAX),intent(inout)::ZC,TCDT,VAPCDT
    real(r8),dimension(NPMAX,NCMAX-1),intent(inout)::TLCDT
    real(r8),dimension(NCMAX-1),intent(inout)::WCANDT
!snow
    real(r8),dimension(NSPMAX),intent(inout)::ZSP,DZSP,RHOSP,TSPDT,DLWDT
    integer(r4),dimension(NSPMAX),intent(inout)::ICESPT
    real(r8),dimension(11),intent(inout)::WLAG
!hydrology
    real(r8),intent(inout)::STORE,SNOWEX,RUNOFF,infiltration,AbsorbedSW,AbsorbedLW,snowmelt,InDirect,InDiffuse
    real(r8),intent(inout)::THFLUX,VFLUX
!residue
    real(r8),dimension(NRMAX),intent(inout)::ZR,RHOR
    real(r8),dimension(NRMAX),intent(inout)::TRDT,VAPRDT,GMCDT
    real(r8),intent(in)::GMCMAX,RLOAD,ZRTHIK,COVER,ALBRES,RESCOF
    real(r8),intent(inout)::DIRRES
!soil
    real(r8),dimension(NSMAX),intent(in)::ZS
    real(r8),dimension(NSMAX),intent(inout)::TSDT,VLCDT,VICDT,MATDT
    real(r8),dimension(NSALTMAX,NSMAX),intent(inout)::CONCDT,SALTDT
    integer(i4),dimension(NSMAX),intent(inout)::ICESDT
    real(r8),intent(in)::ALBDRY,ALBEXP
    real(r8),dimension(NSALTMAX),intent(inout)::DGRADE,SLTDIF
    real(r8),dimension(NSMAX),intent(inout)::ASALT,DISPER
!surface
    real(r8),intent(in)::ZMSRF,ZHSRF,ZERSRF,ZMSP,ZHSP
    real(r8),intent(inout)::POND
    real(r8),intent(in)::PONDMX
!location
    real(r8),intent(in)::ALATUD,LONTITUD,ELEVATION,SLOPE,ASPECT,HRNOON,CLOUDS,DECLIN,HAFDAY
!forcing
    real(r8),intent(inout)::SUNHOR,TMPDAY,WINDAY,PRECIP,SNODEN,PRESUR,shadow,skyview
    real(r8),intent(inout)::HUMDAY
    real(r8),intent(in)::SOITMP,VLCDAY
    real(r8),dimension(NSMAX),intent(inout)::SOILXT
    real(r8),intent(in)::TSAVG
!water balance
    real(r8),intent(inout)::ETSUM,EVAP1
    real(r8),dimension(NSMAX),intent(inout)::TOTFLO
!temp variable
    real(r8),dimension(NPMAX+1,NCMAX-1),intent(inout)::LWCAN,SWCAN
    real(r8),intent(inout)::LWSNOW,LWSOIL
    real(r8),dimension(NRMAX),intent(inout)::LWRES,SWRES
    real(r8),dimension(NSPMAX),intent(inout)::SWSNOW
    real(r8),intent(inout)::swsoil
    real(r8),intent(in)::WWDT
    real(r8),intent(inout)::maskgrid
    real(r8),intent(inout)::BOT
    real(r8),intent(inout)::DGL



!soil state variable
    real(r8),dimension(NSMAX)::TS,MAT,VLC,VIC,QSL,QSV,XTRACT,US,SS,TK,CS,ROOTXT
    real(r8),dimension(NSALTMAX,NSMAX)::CONC,SALT,SINK
    integer(i4),dimension(NSMAX)::IBICES,ICES
    real(r8),dimension(NSMAX)::BTS,BMAT,BVLC,BVIC
    real(r8),dimension(NSALTMAX,NSMAX)::BCONC,BSALT
!residue variable
    real(r8),dimension(NRMAX)::TR,VAPR,GMC,UR,SR,QVR
    real(r8),dimension(NRMAX)::BTR,BVAPR,BGMC
!cannopy variable
    real(r8),dimension(NCMAX)::TC,VAPC,BTC,BVAPC
    real(r8),dimension(NPMAX,NCMAX-1)::TLC,BTLC
    real(r8),dimension(NCMAX-1)::WCAN,UC,SC,QVC,BWCAN
    real(r8),dimension(NPMAX)::PCAN,BPCAN
    real(r8),dimension(NPMAX+1)::TRNSP
!snow variable
    real(r8),dimension(NSPMAX)::TSP,DLW,USP,SSP,QVSP,TQVSP
    real(r8),dimension(NSPMAX):: BTSP,BDLW
    integer(i4),dimension(NSPMAX)::ICESP
    real(r8),dimension(maxlayer)::DELTA,DELNRG,DELWTR

    real(r8),dimension(:),pointer::B,SAT,RHOB,ENTRY,SAND,CLAY,SILT,OM,SATK,thfc
    real(r8),dimension(:,:),pointer::SALTKQ
    real(r8),pointer::DT,WT,WDT
    real(r8),dimension(:),pointer::A1,B1,C1,D1,A2,B2,C2,D2

!
    integer(i4)::I,J,HRSTRT
    real(r8)::DUMMY,RESDEN,HUM,SATV,satvdt,TLCONC,TOTPOT,VAPDT,C1DDT,DAMPNG,DZS,TSNS,VLCNS, &
               VLC1,AVAILA
    real(r8)::HUMDT,TA,TADT,VAPA,VAPADT,WIND
    integer(i4)::MAXNDT,MAXDBL,NDT,NTIMES
    integer(i4)::ITER,ITRSLT,N,MATERL
    real(r8)::HTOLER,TLCNDT,TMP,TMPDT,TOTPDT,VAP
    real(r8):: VAPSPT,VAPSP
    real(r8)::GFLUX
    integer(i4)::IEFLAG,IWFLAG,ICE
    real(r8)::SEEP

    integer(i4),pointer::MAXSTP,MINSTP
    real(r8)::TGFLUX,TLWCAN,TLWRES,TLWSNO,TLWSOI,TOPSNO,TSEEP,TSWCAN,TSWRES,TSWSNO,TSWSOI,HFLUX

    integer(i4)::NSPLST,MAXTRY,extempflag
    real(r8)::RAIN
    real(r8)::CHKMAT
    real(r8)::VLC2

    real(r8)::CONG,DVLC,DZDGL,GDZ,GWS,hhh,HK1,HK2,PROSAT,QMAX,QSLG,RE,REACH
    real(r8),dimension(NSMAX)::HK,mat1,mat2
    integer(i4)::errorflag
    real(r8)::soilbeta
    integer(i4)::mark



    B=>B2d(col,row,:)
    SAT=>SAT2d(col,row,:)
    RHOB=>RHOB2d(col,row,:)
    ENTRY=>ENTRY2d(col,row,:)
    SAND=>SAND2d(col,row,:)
    SILT=>SILT2d(col,row,:)
    CLAY=>CLAY2d(col,row,:)
    OM=>OM2d(col,row,:)
    SATK=>SATK2d(col,row,:)
    SALTKQ=>SALTKQ2d(col,row,:,:)
    thfc=>thfc2d(col,row,:)

    DT=>DT2d(col,row)
    WT=>WT2d(col,row)
    WDT=>WDT2d(col,row)

    A1=>A12d(col,row,:)
    B1=>B12d(col,row,:)
    C1=>C12d(col,row,:)
    D1=>D12d(col,row,:)
    A2=>A22d(col,row,:)
    B2=>B22d(col,row,:)
    C2=>C22d(col,row,:)
    D2=>D22d(col,row,:)

    MAXSTP=>MAXSTP2d(col,row)
    MINSTP=>MINSTP2d(col,row)


    !input check
    if(HUMDAY<0)  then
      print*,col,row,"the input relative humidity is wrong",HUMDAY
      HUMDAY=0.10
    end if
    if(HUMDAY>1)  then
      print*,col,row,"the input relative humidity is wrong",HUMDAY
      HUMDAY=1.0
    end if

    !input check
    if(SUNHOR<0)  then
      print*,col,row,"the input relative SUNHOR is wrong",SUNHOR
      SUNHOR=0.0
    end if
    if(SUNHOR>1300.0)  then
      print*,col,row,"the input relative SUNHOR is wrong",SUNHOR
      SUNHOR=1300
    end if

    if(abs(TMPDAY)>80.0)  then
      print*,col,row,"the input relative TMPDAY is wrong",TMPDAY
      stop
    end if
    
    if(WINDAY .lt. 0.0)  then
      print*,col,row,"the input relative WINDAY is wrong",WINDAY
      stop
    end if
    
    if(PRECIP .lt. 0.0)  then
      print*,col,row,"the input relative PRECIP is wrong",PRECIP
      stop
    end if


    IF (INITAL .EQ. 0) THEN
!   FIRST TIME INTO SUBROUTINE FOR CURRENT PROFILE -- INITIALIZE
!   POORLY DEFINED STATE VARIABLES
!
!   DEFINE MATRIC POTENTIAL OR WATER CONTENT OF SOIL AND DETERMINE
!   WHETHER SOIL IS FROZEN
        DO 10 I=1,NS
            IF (INPH2O.NE.1) THEN
!           INPUT SOIL MOISTURE IS WATER CONTENT
                CALL MATVL1 (I,MATDT(I),VLCDT(I),col,row)
            ELSE
!           INPUT SOIL MOISTURE IS MATRIC POTENTIAL
                CALL MATVL2 (I,MATDT(I),VLCDT(I),col,row)
                VICDT(I)=0.0
                ICESDT(I)=0
            END IF
            IF (TSDT(I) .LE. 0.0) CALL FROZEN (I,VLCDT,VICDT,MATDT,&
                CONCDT,TSDT,SALTDT,ICESDT,nsalt,col,row)
        10 CONTINUE
!
!        INITIALIZE PONDING AND SNOWPACK LAG AND STORAGE
         POND=0.0
         STORE=0.0
         SNOWEX=0.0
         DO 15 I=1,11
            WLAG(I)=0.0
   15    CONTINUE
         IF (NR .GT. 0) THEN
!           SET INITIAL VAPOR DENSITY AND TEMPERATURE OF RESIDUE EQUAL
!           TO SOIL SURFACE NODE.
            TLCONC=0.0
            DO 20 J=1,NSALT
               TLCONC=TLCONC+CONCDT(J,1)
   20       CONTINUE
            TOTPOT=MATDT(1)-TLCONC*UGAS*(TSDT(1)+273.16)/G
            if(abs(TSDT(1)) .gt. 90.0) then
              print*,"GOSHAW 245",col,row,TSDT(1)
              stop
            end if
            CALL VSLOPE (DUMMY,SATVDT,TSDT(1))
            VAPDT=SATVDT*EXP(.018*G/UGAS/(TSDT(1)+273.16)*TOTPOT)
!           IF SOIL IT TOO TOO DRY, RESET VAPOR DENSITY TO 1% REL. HUM.
            IF (VAPDT .LT. 0.01*SATVDT) VAPDT = 0.01*SATVDT
            HUM=VAPDT/SATVDT
!
!           INITIALIZE CONDITIONS AND PARAMETERS FOR EACH LAYER
            RESDEN=RLOAD/ZRTHIK
            IF (COVER .LT. 0.999) THEN
               DIRRES=-LOG(1.-COVER)/RLOAD/10.
              ELSE
!              AVOID NUMERICAL OVERLOAD
               DIRRES = 6.9/RLOAD/10.
            END IF
!           INITIALIZE RESIDUE WATER CONTENT DEPENDING ON INPUT
            IF(GMCDT(1).LE.0.0) &
                CALL RESHUM(2,HUM,DUMMY,GMCDT(1),TSDT(1))
            DO 25 I=1,NR
               ZR(I)=ZRTHIK*I/(NR+1)
               GMCDT(I)=GMCDT(1)
               RHOR(I)=RESDEN
               VAPRDT(I)=VAPDT
               TRDT(I)= TSDT(1)
   25       CONTINUE
            ZR(1)=0.0
            ZR(NR+1)=ZRTHIK
         END IF
!
         IF (NPLANT.NE.0) THEN
!           IF NODE SPACING NOT USER SPECIFIED, INITIALIZE NC
            IF (MZCINP.EQ.0) NC=1
!           INITIALIZE TEMPERATURE AND VAPOR DENSITY OF TOP LAYER
!           (REMAINDER WILL BE SET IN SUBROUTINE CANLAY)
            if(abs(TMPDAY) .gt. 120)then
              print*,"GOSHAW 282,TMPDAY",col,row,tmpday
              stop
            end if
            CALL VSLOPE (DUMMY,SATV,TMPDAY)
            TCDT(1)=TMPDAY
            VAPCDT(1)=SATV*HUMDAY
            IF (WCANDT(1) .LE. 0.0) THEN
!             INITIALIZE WATER CONTENT BASED ON ATMOSPHERIC HUMIDITY
              CALL CANHUM (2,HUMDAY,DUMMY,WCANDT(1),TCDT(1),CANMA,CANMB)
            END IF
!           INITIALIZE PRECIP INTERCEPTION FOR CANOPY
            DO 30 J=1,NPLANT
               TLCDT(J,1)=TCDT(1)
               PCANDT(J)=0.0
   30       CONTINUE
!
!           DEFINE LAYERING OF CANOPY AND ROOT DISTRIBUTION
            CALL CANOPY (NPLANT,NC,NSP,0,MZCINP,ITYPE,INITAL, &
                PLTHGT,PLTWGT,PLTLAI,RLEAF0,ZSP,ZC,TC,TCDT,TLC,TLCDT,VAPC, &
                VAPCDT,WCAN,WCANDT,WCMAX,PCANDT,CANMA,CANMB,TMPDAY,HUMDAY,col,row)
!           DETERMINE ROOT DISTRIBUTION (IF NOT DEFINED BY USER)
            IF (MZCINP.NE.1 .AND. NC.GT.0) &
                CALL RTDIST (NPLANT,ITYPE,NS,ZS,ROOTDP,RROOT0,col,row)
           ELSE
            NC=0
         END IF
!
         HRSTRT=HOUR-NHRPDT
!        INITIALIZE THE WATER BALANCE SUMMARY
         IF (LVLOUT(7) .NE. 0) &
            CALL WBALNC(NPLANT,NC,NSP,NR,NS,JULIAN,HRSTRT,YEAR,&
            ITYPE,INITAL,ZC,WCAN,WCANDT,PCAN,PCANDT,VAPC,VAPCDT,RHOSP,DZSP,&
            DLWDT,WLAG,STORE,ZR,GMC,GMCDT,VAPR,VAPRDT,RHOR,ZS,VLC,VLCDT,&
            VIC,VICDT,TOTFLO,PRECIP,RUNOFF,POND,EVAP1,ETSUM,col,row)
!
!        PRINT OUT INITIAL CONDITIONS
!         CALL OUTPUT (NPLANT,NC,NSP,NR,NS,LVLOUT,0,INPH2O,
!     >    JULIAN,HRSTRT,YEAR,INITAL,ZC,TCDT,VAPCDT,WCANDT,ROOTXT,RHOSP,
!     >    ZSP,TSPDT,DLWDT,ZR,TRDT,VAPRDT,GMCDT,ZS,TSDT,VLCDT,VICDT,
!     >    MATDT,TOTFLO,CONCDT,SALTDT,ndt)
!xxxx
!    >    MATDT,TOTFLO,CONCDT,SALTDT)
!
!        PRINT OUT INITIAL FROST AND SNOW DEPTH
!	     CALL FROST (NSP,NS,JULIAN,&
!            HRSTRT,YEAR,INITAL,ZSP,RHOSP,DZSP,DLWDT,WLAG,STORE,&
!            ZS,VLCDT,VICDT,ICESDT,col,row)

      ELSE
        IF (NPLANT .GT. 0) THEN
!       DEFINE LAYERING OF CANOPY AND ROOT DISTRIBUTION
         CALL CANOPY (NPLANT,NC,NSP,0,MZCINP,ITYPE,1, &
            PLTHGT,PLTWGT,PLTLAI,RLEAF0,ZSP,ZC,TC,TCDT,TLC,TLCDT,VAPC,&
            VAPCDT,WCAN,WCANDT,WCMAX,PCANDT,CANMA,CANMB,TMPDAY,HUMDAY,col,row)
!        DETERMINE ROOT DISTRIBUTION (IF NOT DEFINED BY USER)
         IF (MZCINP.NE.1 .AND. NC.GT.0) &
            CALL RTDIST (NPLANT,ITYPE,NS,ZS,ROOTDP,RROOT0,col,row)
        ELSE
         NC=0
        END IF
!-----------------------------------------------------------------------
!        IF (LEVEL(1) .GE. 1) THEN
!        PRINT CONDITIONS AT BEGINNING OF TIME STEP FOR DEBUGGING
!         CALL OUTPUT (NPLANT,NC,NSP,NR,NS,LVLOUT,1,INPH2O,JULIAN,
!     >   HOUR-NHRPDT,YEAR,1,ZC,TCDT,VAPCDT,WCANDT,ROOTXT,RHOSP,
!     >   ZSP,TSPDT,DLWDT,ZR,TRDT,VAPRDT,GMCDT,ZS,TSDT,VLCDT,VICDT,
!     >   MATDT,TOTFLO,CONCDT,SALTDT,ndt)
!xxxx
!    >   MATDT,TOTFLO,CONCDT,SALTDT)
!        END IF
!-----------------------------------------------------------------------
      END IF
!
!
      WDT=WWDT
      WT=1.-WDT
      DT=DTIME
      NTIMES=1
      NDT=1
      MAXNDT=NHRPDT*MAXSTEP
      MAXDBL=MAXNDT
!
!
!
!**** UPDATE BOUNDARY CONDITIONS FOR NEXT TIME STEP
      TA = TMPDAY
      TADT=TMPDAY
      HUM = HUMDAY
      HUMDT=HUMDAY
      WIND=WINDAY
      if(abs(TMPDAY).gt.120)then
        print*,"GOSHAW 373 TMPDAY",col,row,tmpday
        stop
      end if
      CALL VSLOPE (DUMMY,SATV,TMPDAY)
      VAPA=HUM*SATV
      VAPADT=VAPA
      Bot=0.0
!
!
!**** SET LOWER BOUNDARY CONDITIONS
      IF (ITMPBC .eq. 0) THEN
!....    ESTIMATE TEMP AT LOWER BOUNDARY BY METHOD OF HIROTA ET AL 2002,
!        JGR 107, NO. d24, 4767, doi:10.1029/2001JD001280, 2002
!         IF (INITAL .EQ. 0) THEN
!            WRITE (6,*) ' Enter the average soil temperature at depth: '
!            READ (5,*) TSAVG
!         END IF
         CALL SOILTK (NS,TK,VLCDT,VICDT,col,row)
         CALL SOILHT (NS,CS,VLCDT,VICDT,TSDT,MATDT,CONCDT,nsalt,col,row)
         IF (TSDT(NS).LE.0.0) THEN
!           ADJUST HEAT CAPACITY FOR SOIL FREEZING - INCLUDE LATENT
!           HEAT OF FUSION OVER THE NEXT 1.0C INTERVAL
            TSNS = TSDT(NS)
            VLCNS = VLCDT(NS)
            IF (TSDT(NS-1) .GT. TSDT(NS)) THEN
               TSDT(NS) = TSDT(NS) + 1.0
              ELSE
               TSDT(NS) = TSDT(NS) - 1.0
            END IF
            CALL FROZEN (NS,VLCDT,VICDT,MATDT,CONCDT,TSDT,SALTDT,ICESDT,nsalt,col,row)
            CS(NS) = CS(NS) + ABS(RHOL*LF*(VLCNS-VLCDT(NS)))
!           RESET BOTTOM TEMPERATURE
            TSDT(NS) = TSNS
            CALL FROZEN (NS,VLCDT,VICDT,MATDT,CONCDT,TSDT,SALTDT,ICESDT,nsalt,col,row)
         END IF
!        CALCULATE DAMPING DEPTH, DELTA Z, AND C1(DELTA)/DT
         DAMPNG=SQRT(2.*TK(NS)/CS(NS)/1.99238E-07)
         DZS = ZS(NS)-ZS(NS-1)
         C1DDT = (1. + 2.*DZS/DAMPNG)/DT
         TSDT(NS) = TSDT(NS) + (-2.*TK(NS)/CS(NS)/DAMPNG/DZS/C1DDT) &
            *(TSDT(NS)-TSDT(NS-1)) - 1.99238E-07*(TSDT(NS) - TSAVG)/C1DDT
!        CALCULATE WEIGHTING COEFFICIENT BASED ON DAILY TIME STEP
!cc      AA=(-0.00082+0.00983957*DAMPNG/(ZS(NS)-ZS(NS-1)))
!cc  >       *(ZS(NS)/DAMPNG)**(-0.381266)
!cc      IF (AA .LT. 0.0) AA=0.0
!        ADJUST FOR ACTUAL TIME STEP
!cc      AA=AA*NHRPDT/24
!cc      TSDT(NS)=(1.-AA)*TSDT(NS) + AA*TSDT(NS-1)
        ELSEif(ITMPBC .eq. 1) then   !non heat flux boundary
          TSDT(NS)=TSDT(NS-1)+(ZS(NS)-ZS(NS-1))*GroundTempGradients
        ELSE
!        LOWER TEMPERATURE BOUNDARY SPECIFIED FROM INPUT FILE
         TSDT(NS)=SOITMP
      END IF

      IF (IVLCBC .GT. 0) THEN
!....    UNIT GRADIENT IS ASSUMED FOR WATER FLUX AT BOTTOM BOUNDARY;
!        SET MATRIC POTENTIAL FOR LAST NODE EQUAL TO SECOND TO LAST NODE
!        IF FREEZING FRONT IS IN NS-1, DO NOT CHANGE POTENTIAL AT BOTTOM
         IF(VICDT(NS-1).LE.0. .OR. VICDT(NS).GT.0.)MATDT(NS)=MATDT(NS-1)
         VLC1=VLCDT(NS)+VICDT(NS)*RHOI/RHOL
         VLC2=VLCDT(NS)
         VICDT(NS)=0.0
         CALL MATVL2 (NS,MATDT(NS),VLCDT(NS),col,row)
         IF (TSDT(NS) .LT. 0.0) VICDT(NS)=(VLC1-VLCDT(NS))*RHOL/RHOI
         IF (VICDT(NS) .LT. 0.0) VICDT(NS)=0.0
         Bot=VLC1-(VICDT(NS)+VLCDT(NS))
        ELSE
!        INPUT WATER CONTENT SPECIFIED FOR WATER FLUX AT BOTTOM BOUNDARY
         IF (INPH2O.NE.1) THEN
!           INPUT SOIL MOISTURE IS WATER CONTENT
            VLCDT(NS)=VLCDAY
            CALL MATVL1 (NS,MATDT(NS),VLCDT(NS),col,row)
          ELSE
!           INPUT SOIL MOISTURE IS MATRIC POTENTIAL
            MATDT(NS)=VLCDAY
            CALL MATVL2 (NS,MATDT(NS),VLCDT(NS),col,row)
         END IF
         VICDT(NS)=0.0
      END IF

      ICESDT(NS)=0
      DO 105 J=1,NSALT
        CONCDT(J,NS)=SALTDT(J,NS)/(SALTKQ(J,NS)+VLCDT(NS)*RHOL/RHOB(NS))
  105 CONTINUE
      IF (TSDT(NS) .LE. 0.0) CALL FROZEN (NS,VLCDT,VICDT,MATDT, &
                CONCDT,TSDT,SALTDT,ICESDT,nsalt,col,row)
!

!
      IF (MWATRXT.EQ.1) THEN
!****    SOIL SINK TERM INPUT FOR EACH LAYER
!        CHECK IF SOIL LAYERS CAN SATIFY SINK TERM
         AVAILA=VLCDT(1)*(ZS(2)-ZS(1))/2
         IF (SOILXT(1)*DT.GT.AVAILA) SOILXT(1)=0.0
         AVAILA=VLCDT(NS)*(ZS(NS)-ZS(NS-1))/2
         IF (SOILXT(NS)*DT.GT.AVAILA) SOILXT(NS)=0.0
         DO 1105 I=2,NS-1
            AVAILA=VLCDT(I)*(ZS(I+1)-ZS(I-1))/2
            IF (SOILXT(I)*DT.GT.AVAILA) SOILXT(I)=0.0
 1105    CONTINUE
        ELSE
!        SET SOIL SINK TERM TO ZERO
         DO 1106 I=1,NS
            SOILXT(I)=0.0
 1106    CONTINUE
      END IF
!
!
!**** DEFINE THE ORDER IN WHICH THE MATERIALS ARE LAYERED
      CALL SYSTM (NC,NPLANT,NSP,ZC,ZSP,ZMSRF,ZHSRF,ZERSRF,ZMSP, &
            ZHSP,HEIGHT,SUNHOR,TMPDAY,TCCRIT,col,row)
!
!
!**** CALCULATE THE SHORT-WAVE ENERGY BALANCE
      CALL SWRBAL (NPLANT,NC,NSP,NR,XANGLE,CLUMPNG,SWCAN,SWSNOW,SWRES,&
            SWSOIL,SUNHOR,CANALB,ZSP,DZSP,RHOSP,ZR,RHOR,ALBRES,DIRRES,&
            ALBDRY,ALBEXP,VLCDT(1),ALATUD,LONTITUD,SLOPE,ASPECT,HRNOON,HAFDAY,DECLIN,&
            HOUR,NHRPDT,shadow,skyview,shadeef,col,row,julian,year,InDirect,InDiffuse)
!
!
!**** SAVE VALUES AT THE BEGINNING OF THE TIME STEP
      CALL UPDATE (NS,NR,NSP,NC,NPLANT,NSALT,IBICES,ICESDT,BTS,TSDT, &
        BMAT,MATDT,BCONC,CONCDT,BVLC,VLCDT,BVIC,VICDT,BSALT,SALTDT, &
        BTR,TRDT,BVAPR,VAPRDT,BGMC,GMCDT,BTC,TCDT,BTLC,TLCDT, &
        BVAPC,VAPCDT,BWCAN,WCANDT,BPCAN,PCANDT,BTSP,TSPDT,BDLW,DLWDT,&
        ICESP,ICESPT)
!
!**** UPDATE PARAMETERS FOR NEXT TIME STEP
  110 CALL UPDATE (NS,NR,NSP,NC,NPLANT,NSALT,ICES,ICESDT,TS,TSDT,&
        MAT,MATDT,CONC,CONCDT,VLC,VLCDT,VIC,VICDT,SALT,SALTDT,&
        TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,TC,TCDT,TLC,TLCDT,VAPC,VAPCDT,&
        WCAN,WCANDT,PCAN,PCANDT,TSP,TSPDT,DLW,DLWDT,ICESP,ICESPT)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



!**** START ITERATIVE PROCEDURE TO SOLVE ENERGY AND MOISTURE BALANCE
!
  120 ITER = 0
      ITRSLT = 0
      extempflag=0
      errorflag=0
  200 ITER = ITER + 1

      IF (WDT.LE.1.0) THEN
!        IF SATURATED OR EXTREMELY DRY CONDITIONS EXIST,
!        SPECIFY FULLY IMPLICIT SOLUTION
         DO 100 I=NS-1,1,-1
            IF (MATDT(I).GT.ENTRY(I) .OR. VLCDT(I).LT.0.001) THEN
               WDT=1.0
               WT=0.0
               GO TO 101
            END IF
  100    CONTINUE
  101    CONTINUE
      END IF
!
!**** DETERMINE LONG-WAVE RADIATION BALANCE FOR EACH NODE
      CALL LWRBAL (NC,NSP,NR,NPLANT,TA,TADT,TLC,TLCDT,TSP,TSPDT, &
            TR,TRDT,TS,TSDT,VAPA,VAPADT,CLOUDS,LWCAN,LWSNOW,LWRES,LWSOIL,col,row)
!
!**** CALCULATE LONG-WAVE RAD. CONTRIBUTION TO ENERGY BALANCE MATRIX
      CALL LWRMAT (NC,NSP,NR,TC,TSP,TR,TS,ICESPT,col,row)
!**** SUM THE SOURCE-SINK TERMS FOR EACH NODE
      CALL SOURCE (NC,NSP,NR,NS,NSALT,NPLANT,UC,SC,USP,SSP,UR,SR,US,SS, &
            SINK,SWCAN,SWSNOW,SWRES,SWSOIL,LWCAN,LWSNOW,LWRES,LWSOIL, &
            SOILXT)
!
      N = 1
!
!**** DETERMINE THE BOUNDARY CONDITION FOR THE SURFACE MATERIAL
      HTOLER=0.001
      MATERL=2
!
      if(VLCDT(1).gt.thfc(1)) then
	soilbeta=1.0
      else
	soilbeta=0.25*(1.0-cos(3.14159*VLCDT(1)/thfc(1)))**2.0
      endif
      IF (NC .GT. 0) THEN
!        CANOPY IS THE SURFACE MATERIAL
         CALL ATSTAB(NPLANT,NC,NSP,NR,TA,TADT,TC(1),TCDT(1),VAPA,VAPADT, &
            VAPC(1),VAPCDT(1),WIND,HEIGHT,HTOLER,HFLUX,VFLUX,ZC,ZR,RHOR, &
            0,ITER,PRESUR,col,row,soilbeta)
!        GO THROUGH CALCULATION OF THE VAPOR DENSITY AT THE SOIL SURFACE
         GO TO 205
      END IF
!
      IF (NSP .GT. 0) THEN
!        SNOW IS THE SURFACE MATERIAL
         if(abs(tsp(1)) .gt. 70.0)then
            print*,"GOSHAW 593 TSP(1)",col,row,tsp(1)
            stop
         end if
         CALL VSLOPE (DUMMY,SATV,TSP(1))
         if(abs(TSPDT(1)) .gt. 70.0)then
            print*,"GOSHAW LINE 572",col,row,TSPDT(1)
            stop
         end if
         CALL VSLOPE (DUMMY,SATVDT,TSPDT(1))
!        CALCULATE THE SATURATED VAPOR DENSITY OVER ICE
         TMP = TSP(1) + 273.16
         TMPDT = TSPDT(1) + 273.16
         VAP = SATV*EXP(0.018/(UGAS*TMP)*LF*TSP(1)/TMP)
         VAPDT = SATVDT*EXP(0.018/(UGAS*TMPDT)*LF*TSPDT(1)/TMPDT)
         soilbeta=1.0
         CALL ATSTAB(NPLANT,NC,NSP,NR,TA,TADT,TSP(1),TSPDT(1),VAPA, &
            VAPADT,VAP,VAPDT,WIND,HEIGHT,HTOLER,HFLUX,VFLUX,ZC,ZR,RHOR, &
            ICESPT(1),ITER,PRESUR,col,row,soilbeta)
!        GO THROUGH CALCULATION OF THE VAPOR DENSITY AT THE SOIL SURFACE
         GO TO 205
      END IF
!
      IF (NR .GT. 0) THEN
!        RESIDUE IS THE SURFACE MATERIAL
         CALL ATSTAB(NPLANT,NC,NSP,NR,TA,TADT,TR(1),TRDT(1),VAPA,VAPADT, &
            VAPR(1),VAPRDT(1),WIND,HEIGHT,HTOLER,HFLUX,VFLUX,ZC,ZR,RHOR, &
            0,ITER,PRESUR,col,row,soilbeta)
      END IF
!
!     SOIL IS THE SURFACE MATERIAL
!     ---- FIRST CALCULATE THE TOTAL WATER POTENTIAL, THEN THE HUMIDITY
!          MAY BE DETERMINED.  CALL VSLOPE TO OBTAIN THE SATURATED
!          VAPOR DENSITY.  USING HUMIDITY AND SATURATED VAPOR DENSITY,
!          THE VAPOR PRESSURE AT THE SOIL SURFACE MAY BE CALCULATED.
  205 TLCONC=0.0
      TLCNDT=0.0
      DO 210 J=1,NSALT
         TLCONC=TLCONC+CONC(J,1)
         TLCNDT=TLCNDT+CONCDT(J,1)
  210 CONTINUE
      TOTPOT=MAT(1)-TLCONC*UGAS*(TS(1)+273.16)/G
      TOTPDT=MATDT(1)-TLCNDT*UGAS*(TSDT(1)+273.16)/G
      if(abs(TS(1)) .gt. 90.0)then
         print*,"GOSHAW,TS(1),636",col,row,TS(1)
         stop
      end if
      CALL VSLOPE (DUMMY,SATV,TS(1))
      if(abs(TSDT(1)) .gt. 90.0)then
         print*,"GOSHAW,641,TSDT(1)",col,row,TSDT(1)
         stop
      end if
      CALL VSLOPE (DUMMY,SATVDT,TSDT(1))
      VAP=SATV*EXP(.018*G/UGAS/(TS(1)+273.16)*TOTPOT)
      VAPDT=SATVDT*EXP(.018*G/UGAS/(TSDT(1)+273.16)*TOTPDT)
      IF (NC.EQ.0 .AND. NR.EQ.0 .AND. NSP.EQ.0) then
            CALL ATSTAB(NPLANT,NC,NSP,NR,TA,TADT,TS(1),TSDT(1),VAPA,VAPADT,&
            VAP,VAPDT,WIND,HEIGHT,HTOLER,HFLUX,VFLUX,ZC,ZR,RHOR,&
            ICESDT(1),ITER,PRESUR,col,row,soilbeta)
      endif
!
!**** DETERMINE ENERGY BALANCE FOR THE CANOPY LAYERS
      IF (NC .GT. 0) THEN
!        DEFINE VAPOR DENSITY FOR THE LOWER BOUNDARY OF CANOPY
         IF (NSP .GT. 0) THEN
!           CANOPY IS OVERLYING SNOWPACK
            if(abs(TSP(1)) .gt. 70.0)then
              print*,"GOSHAW,659,TSP(1)",col,row,TSP(1)
              stop
            end if
            CALL VSLOPE (DUMMY,SATV,TSP(1))
            if(abs(TSPDT(1)) .gt. 70.0)then
              print*,"GOSHAW,664,TSPDT(1)",col,row,TSPDT(1)
              stop
            end if
            CALL VSLOPE (DUMMY,SATVDT,TSPDT(1))
            TMP = TSP(1) + 273.16
            TMPDT = TSPDT(1) + 273.16
            VAPC(NC+1) = SATV*EXP(0.018/(UGAS*TMP)*LF*TSP(1)/TMP)
           VAPCDT(NC+1)=SATVDT*EXP(0.018/(UGAS*TMPDT)*LF*TSPDT(1)/TMPDT)
            TC(NC+1)=TSP(1)
            TCDT(NC+1)=TSPDT(1)
          ELSE
            IF (NR .GT. 0) THEN
!              CANOPY IS OVERLYING RESIDUE
               VAPC(NC+1) = VAPR(1)
               VAPCDT(NC+1) = VAPRDT(1)
               TC(NC+1) = TR(1)
               TCDT(NC+1) = TRDT(1)
             ELSE
!              CANOPY IS OVERLYING BARE SOIL
               VAPC(NC+1) = VAP
               VAPCDT(NC+1) = VAPDT
               TC(NC+1) = TS(1)
               TCDT(NC+1) = TSDT(1)
            END IF
         END IF

         CALL EBCAN (N,NPLANT,NC,NSP,NR,NS,ZC,TC,TCDT,TLC,TLCDT, &
            VAPC,VAPCDT,WCAN,WCANDT,PCAN,PCANDT,MAT,MATDT,SWCAN,LWCAN, &
            CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,ITYPE,ITER,col,row,julian,errorflag)
        !if(errorflag.gt.0)then
        ! maskgrid=1
        ! return	
        !end if
         MATERL = MATERL + 1
      END IF
!
!**** DETERMINE ENERGY BALANCE FOR THE SNOW LAYERS
      IF (NSP .GT. 0) THEN
!        DEFINE VAPOR DENSITY AND TEMPERATURE FOR LOWER BOUNDARY OF
!        SNOWPACK
         IF (NR .GT. 0) THEN
!           SNOW IS OVERLYING RESIDUE
            VAPSP = VAPR(1)
            VAPSPT = VAPRDT(1)
            TSP(NSP+1) = TR(1)
            TSPDT(NSP+1) = TRDT(1)
           ELSE
!           SNOW IS OVERLYING BARE SOIL
            VAPSP = VAP
            VAPSPT = VAPDT
            TSP(NSP+1) = TS(1)
            TSPDT(NSP+1) = TSDT(1)
         END IF
         CALL EBSNOW (N,NSP,NR,ICESPT,TSP,TSPDT,DLW,DLWDT,RHOSP,ZSP, &
             DZSP,QVSP,VAPSP,VAPSPT,SSP,ITER,presur,col,row)
         MATERL = MATERL + 1
      END IF
!
!**** DETERMINE ENERGY BALANCE FOR THE RESIDUE LAYERS
      IF (NR .GT. 0) THEN
         VAPR(NR+1) = VAP
         VAPRDT(NR+1)=VAPDT
         TR(NR+1)=TS(1)
         TRDT(NR+1)=TSDT(1)
         CALL EBRES (N,NR,NSP,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT, &
            GMCMAX,RHOR,RESCOF,SR,UR,QVR,RHOSP,ITER,presur,col,row)
         MATERL = MATERL + 1
      END IF
!
!**** SOLVE FOR ENERGY BALANCE OF SOIL
      CALL EBSOIL (N,NS,ZS,TS,TSDT,MAT,MATDT,CONC,CONCDT, &
        VLC,VLCDT,VIC,VICDT,ICES,ICESDT,QSL,QSV,SS,GFLUX,ITER,PRESUR,nsalt,col,row)
!
!**** SOLVE THE ENERGY BALANCE MATRIX
      mark=0
      CALL TDMA2 (N,A1,B1,C1,D1,DELTA,mark)
      if(mark.eq.0)then
        DELTA=TOLER*4.0
      end if
      DO I=1,N
         IF(abs(DELTA(I)).gt.60) extempflag=1
      ENDDO

      IF(extempflag.eq.1) then
        extempflag=0
        DELTA=DELTA/100.0
      endif



!**** SORT OUT THE SOLUTION OF THE MATRIX INTO THE PROPER MATERIALS
      MATERL=2
      N=1
      IEFLAG=0
      IWFLAG=0
!
      IF (NC .GT. 0) THEN
!**** CANOPY LAYERS
      DO 220 I=1,NC
         IF (ABS(DELTA(N)) .GT. 25.  .AND.  NDT .LT. MAXNDT) THEN
!           TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
!           BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
            IEFLAG=1
            ITER=11
            GO TO 350
         END IF
         IF (ABS(DELTA(N)) .GT. TOLER) IEFLAG = IEFLAG+1
         IF (ITER .GT. 3 .AND. DELNRG(N)*DELTA(N) .LT. 0) THEN
!           DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
!           TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
            DELTA(N)=DELTA(N)/2.
         END IF
         DELNRG(N)=DELTA(N)
         TCDT(I)=TCDT(I)-DELTA(N)
!         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),TCDT(I),VAPCDT(I),
!     >                                    WCANDT(I)
         N=N+1
  220 CONTINUE
      MATERL=MATERL+1
      END IF
!
      IF (NSP .GT. 0) THEN
!**** SNOW PACK LAYERS
      DO 230 I=1,NSP
         IF (ABS(DELTA(N)) .GT. 25.  .AND.  NDT .LT. MAXNDT) THEN
!           TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
!           BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
            IEFLAG=1
            ITER=11
            GO TO 350
         END IF
         IF (ICESPT(I) .EQ. 0) THEN
!           NO LIQUID WATER IN CURRENT LAYER AT END OF TIME STEP
            IF (ABS(DELTA(N)) .GT. TOLER) IEFLAG = IEFLAG + 1
            IF (ITER .GT. 3 .AND. DELNRG(N)*DELTA(N) .LT. 0) THEN
!              DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
!              TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
               DELTA(N)=DELTA(N)/2.
            END IF
            DELNRG(N)=DELTA(N)
            TSPDT(I) = TSPDT(I) - DELTA(N)
!           CHECK IF LAYER HAS BEGUN MELTING
            IF (TSPDT(I) .GT. 0) THEN
!              LAYER HAS GONE ABOVE 0 C - ADJUST WATER CONTENT IN LAYER
               ICESPT(I) = 1
!CCC           DLWDT(I) = RHOI*CI*DZSP(I)*TSPDT(I)/(RHOL*LF)
               TSPDT(I) = 0.0
            END IF
!
          ELSE
!           LAYER CONTAINS LIQUID WATER
!           CONVERT TOLERANCE FOR TEMPERATURE TO LIQUID EQUIVALENT
!           (DELTA LIQUID FRACTION) => 0.001*(DELTA TEMP.)
            IF (ABS(DELTA(N))/DZSP(I) .GT. 0.001*TOLER) THEN
!              IF DELTA < E-07*DLWDT, PRECISION OF COMPUTER IS EXCEEDED,
!              AND  ADDING DELTA TO DLWDT WILL NOT CHANGE DLWDT
               IF (ABS(DELTA(N)) .GT. DLWDT(I)*10E-07) IEFLAG=IEFLAG+1
            END IF
            IF (ITER .GT. 3 .AND. DELNRG(N)*DELTA(N) .LT. 0) THEN
!              DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
!              TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
               DELTA(N)=DELTA(N)/2.
            END IF
            DELNRG(N)=DELTA(N)
            DLWDT(I) = DLWDT(I) - DELTA(N)
!           CHECK IF ALL THE LIQUID HAS FROZEN
            IF (DLWDT(I) .LT. 0.0) THEN
!              LAYER HAS FROZEN COMPLETELY - ADJUST TEMPERATURE
               ICESPT(I) = 0
!CCC           TSPDT(I) = RHOL*LF*DLWDT(I)/(RHOI*CI)
               DLWDT(I) = 0.0
            END IF
         END IF
!         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),TSPDT(I),DLWDT(I),
!     >                                  DZSP(I)
         N=N+1
  230 CONTINUE
      MATERL=MATERL+1
      END IF
!
      IF (NR .GT. 0) THEN
!**** RESIDUE LAYERS
      DO 240 I=1,NR
         IF (ABS(DELTA(N)) .GT. 25.  .AND.  NDT .LT. MAXNDT) THEN
!           TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
!           BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
            IEFLAG=1
            ITER=11
            GO TO 350
         END IF
         IF (ABS(DELTA(N)) .GT. TOLER) IEFLAG = IEFLAG+1
         IF (ITER .GT. 3 .AND. DELNRG(N)*DELTA(N) .LT. 0) THEN
!           DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
!           TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
            DELTA(N)=DELTA(N)/2.
         END IF
         DELNRG(N)=DELTA(N)
         TRDT(I)=TRDT(I)-DELTA(N)
!         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),TRDT(I),VAPRDT(I),
!     >                                    GMCDT(I)
         N=N+1
  240 CONTINUE
      MATERL=MATERL+1
      END IF
!
!**** SOIL LAYERS
      DO 250 I=1,NS-1
         IF (ABS(DELTA(N)) .GT. 25.  .AND.  NDT .LT. MAXNDT) THEN
!           TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
!           BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
            IEFLAG=1
            ITER=11
            GO TO 350
         END IF
         IF (ABS(DELTA(N)) .GT. TOLER) IEFLAG = IEFLAG+1
         IF (ITER .GT. 3 .AND. DELNRG(N)*DELTA(N) .LT. 0) THEN
!           DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
!           TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
            DELTA(N)=DELTA(N)/2.
         END IF
         DELNRG(N)=DELTA(N)
         TSDT(I)=TSDT(I)-DELTA(N)
!
!        CHECK IF LAYER IS BELOW 0 C
         IF (TSDT(I) .LE. 0.0) THEN
!           ICE MAY BE PRESENT - CALL FROZEN TO DETERMINE IF ICE PRESENT
            ICE=ICESDT(I)
            CALL FROZEN (I,VLCDT,VICDT,MATDT,CONCDT,TSDT,SALTDT,ICESDT,nsalt,col,row)
!           CHECK IF LAYER HAS CROSSED THE FREEZING POINT AND ADJUST THE
!           TEMPERATURE FOR LATENT HEAT IF SO
            IF (ICE .EQ. 0  .AND.  ICESDT(I) .EQ. 1) &
            CALL ADJUST (I,VICDT,VLCDT,MATDT,CONCDT,TSDT,SALTDT,ICESDT,nsalt,col,row)
!
           ELSE
!           NO ICE IS PRESENT
            IF (ICESDT(I) .EQ. 1) THEN
!              CONVERT ANY REMAINING ICE TO WATER
               VLCDT(I) = VLCDT(I) + VICDT(I)*RHOI/RHOL
               IF (VLCDT(I) .GT. SAT(I)) VLCDT(I)=SAT(I)
               VICDT(I) = 0.0
               ICESDT(I) = 0
               CALL MATVL1 (I,MATDT(I),VLCDT(I),col,row)
            END IF
         END IF
!
!         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),TSDT(I),VLCDT(I),
!     >                           VICDT(I),MATDT(I),CONCDT(1,I),ICESDT(I)
         N=N+1
  250 CONTINUE
!
      IF (NC .GT. 0) THEN
!        DEFINE TEMPERATURE OF BOTTOM BOUNDARY OF CANOPY
         IF (NSP .GT. 0) THEN
!           SNOW IS UNDERLYING CANOPY
            TCDT(NC+1) = TSPDT(1)
           ELSE
            IF (NR .GT. 0) THEN
!              RESIDUE IS UNDERLYING CANOPY
               TCDT(NC+1) = TRDT(1)
              ELSE
!              SOIL IS UNDERLYING CANOPY
               TCDT(NC+1) = TSDT(1)
            END IF
         END IF
      END IF
!
      IF (NSP .GT. 0) THEN
!        DEFINE TEMPERATURE OF BOTTOM BOUNDARY OF SNOWPACK
         IF (NR .GT. 0) THEN
!           SNOW IS OVERLYING RESIDUE
            TSPDT(NSP+1) = TRDT(1)
           ELSE
!           SNOW IS OVERLYING BARE SOIL
            TSPDT(NSP+1) = TSDT(1)
         END IF
      END IF

!     DEFINE TEMPERATURE OF BOTTOM BOUNDARY OF RESIDUE
      IF (NR .GT. 0) TRDT(NR+1) = TSDT(1)


!{{{{{{{{{{{{{{{{{{
!Insert the code of water balance here after
!
!**** BEGIN CALCULATIONS FOR THE MOISTURE BALANCE OF THE SYSTEM
      N = 1
      do I=1,NS
        IF(abs(TSDT(I))>50.0) then
           TSDT(1:NS)=TS(1:NS)
           exit
        end if
      enddo
      do I=1,NC
        IF(abs(TCDT(I))>50.0) then
           TCDT(1:NC)=TC(1:NC)
           exit
        end if
      enddo
      
      do I=1,NR
        IF(abs(TRDT(I))>50.0) then
           TRDT(1:NR)=TR(1:NR)
           exit
        end if        
      enddo
      do I=1,NSP
        IF(abs(TSPDT(I))>50.0) then
           TSPDT(1:NSP)=TSP(1:NSP)
           exit
        end if
      enddo
      !IF(maskgrid.eq.1 .and. maskflag.eq.1) THEN
      !  DO I=1,NS
      !    TSDT(I)=-9999.0
      !  END DO
      !  INBASIN=0
      !  RETURN
      !ENDIF
!
!**** SOLVE WATER BALANCE FOR CANOPY LAYERS
    IF (NC .GT. 0) THEN
        CALL WBCAN(N,NPLANT,NC,NSP,NR,NS,ZC,TC,TCDT,TLC,TLCDT,VAPC,&
            VAPCDT,WCAN,WCANDT,PCAN,PCANDT,QVC,TRNSP,MAT,MATDT,XTRACT,&
            SWCAN,LWCAN,CANMA,CANMB,DCHAR,RSTOM0,RSTEXP,PLEAF0,RLEAF0,&
            ITYPE,ICESDT(1),ITER,col,row,julian,errorflag)
        !if(errorflag.gt.0)then
        ! maskgrid=1
        ! return	
        !end if
    ELSE
!        SET ROOT EXTRACTION AND PLANT TRANSPIRATION TO ZERO
        DO 300 I=1,NS
          XTRACT(I)=0.0
        300 CONTINUE
!        TRNSP(NPLANT+1) IS TOTAL TRANSPIRATION FOR ALL PLANTS
         DO 302 I=1,NPLANT+1
            TRNSP(I)=0.0
  302    CONTINUE
    END IF
!
!**** SET BOUNDARY CONDITIONS IF THERE IS A SNOWPACK PRESENT
!     (SNOWPACK IS NOT PART OF THE WATER BALANCE MATRIX - THE WATER
!     BALANCE FOR THE SNOWPACK IS DONE AT THE END OF THE HOUR)
      IF (NSP .GT. 0.0) THEN
         A2(N) = 0.0
         CALL SNOWBC (N,NSP,NR,ZSP,QVSP,VAPSPT,TSDT,TSPDT,ICESDT(1),PRESUR,col,row)
      END IF
!
!**** SOLVE FOR WATER BALANCE OF THE RESIDUE LAYERS
      IF (NR .GT. 0) &
         CALL WBRES (N,NR,NSP,ZR,TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,RHOR,&
                     QVR,RHOSP,ICESDT(1),ITER,presur,col,row)
!
!**** SOLVE FOR THE WATER BALANCE OF SOIL LAYERS
      CALL WBSOIL (N,NS,ZS,TS,TSDT,MAT,MATDT,VLC,VLCDT,VIC,VICDT,&
           CONC,CONCDT,ICESDT,QSL,QSV,XTRACT,SEEP,US,ITER,PRESUR,nsalt,col,row)
!
!-----------------------------------------------------------------------
!      IF (LEVEL(1) .GE. 2) THEN
!         WRITE (21,*) 'HFLUX (W/M2), VFLUX (KG/S): ', HFLUX,VFLUX
!         WRITE (21,*) ' VAPOR FLUXES (KG/S)'
!         WRITE (21,*) (QSV(K), K=1,NS-1)
!         WRITE (21,*) ' LIQUID FLUXES (M/S)'
!         WRITE (21,*) (QSL(K), K=1,NS-1)
!         WRITE (21,*) ' JACOBIAN MATRIX'
!         DO 305 I=1,N
!            WRITE (21,*) A2(I),B2(I),C2(I),D2(I)
!  305    CONTINUE
!         WRITE (21,*)
!         WRITE (21,*) ' VALUES AT WATER BALANCE ITERATION ', ITER
!      END IF
!-----------------------------------------------------------------------
!
!**** SOLVE THE WATER BALANCE MATRIX
      mark=0
      CALL TDMA2 (N,A2,B2,C2,D2,DELTA,mark)
      if(mark.eq.0)then
        DELTA=0.00001
      end if

!**** SORT OUT THE SOLUTION OF THE MATRIX INTO THE PROPER MATERIALS
      MATERL=2
      N=1
      IWFLAG=0
!
      IF (NC .GT. 0) THEN
!**** CANOPY LAYERS
      DO 310 I=1,NC
         IF (ABS(DELTA(N)/VAPCDT(I)) .GT. TOLER) IWFLAG = IWFLAG+1
         IF (ITER .GT. 3 .AND. DELWTR(N)*DELTA(N) .LT. 0) THEN
!           DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
!           TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
            DELTA(N)=DELTA(N)/2.
         END IF
         DELWTR(N)=DELTA(N)
         VAPCDT(I)=VAPCDT(I)-DELTA(N)
         IF (VAPCDT(I).le.0)VAPCDT(I)=0.000001
            if(abs(TC(I)) .gt. 70.0)then
              print*,"GOSHAW,1055,TC(I)",col,row,TC(I)
              stop
            end if
         CALL VSLOPE (DUMMY,SATV,TC(I))
         IF (VAPCDT(I).GT.SATV) VAPCDT(I)=SATV
!         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),VAPCDT(I),TCDT(I),&
!                                          WCANDT(I)
         N=N+1
  310 CONTINUE
      MATERL=MATERL+1
      END IF
!
      IF (NSP .GT. 0) THEN
!**** SNOW PACK LAYERS
!     SNOW IS NOT PART OF WATER BALANCE SOLUTION
      MATERL=MATERL+1
      END IF
!
      IF (NR .GT. 0) THEN
!**** RESIDUE LAYERS
      DO 320 I=1,NR
         IF (ABS(DELTA(N)/VAPRDT(I)) .GT. TOLER) IWFLAG = IWFLAG+1
         IF (ITER .GT. 3 .AND. DELWTR(N)*DELTA(N) .LT. 0) THEN
!           DELTA IS JUMPING BETWEEN NEG. AND POS. -- CUT DELTA IN HALF
!           TO SUPPRESS TENDENCY TO JUMP BACK AND FORTH AROUND SOLUTION
            DELTA(N)=DELTA(N)/2.
         END IF
         DELWTR(N)=DELTA(N)
         VAPRDT(I)=VAPRDT(I)-DELTA(N)
         IF (VAPRDT(I).le.0)VAPRDT(I)=0.000001
            if(abs(TR(I)) .gt. 70.0)then
              print*,"GOSHAW 1086 TR(I)",col,row,TR(I)
              stop
            end if
         CALL VSLOPE (DUMMY,SATV,TR(I))
         IF (VAPRDT(I).GT.SATV)VAPRDT(I)=SATV
!         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),VAPRDT(I),TRDT(I),
!     >                                    GMCDT(I)
         N=N+1
  320 CONTINUE
      MATERL=MATERL+1
      END IF
!
!**** SOIL LAYERS
      DO 330 I=1,NS-1
         IF (ICESDT(I) .EQ. 1) THEN
            IF (ABS(DELTA(N)) .GT. 1.0  .AND.  NDT .LT. MAXNDT) THEN
!              TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
!              BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
               IWFLAG=1
               ITER=11
               GO TO 350
            END IF
            IF (ABS(DELTA(N)) .GT. TOLER) IWFLAG = IWFLAG+1
            IF (ITER .GT. 3 .AND. DELWTR(N)*DELTA(N) .LT. 0) THEN
!              DELTA IS JUMPING BETWEEN NEG. AND POS.--CUT DELTA IN HALF
!              TO SUPPRESS TENDENCY TO JUMP AROUND SOLUTION
               DELTA(N)=DELTA(N)/2.
            END IF
            DELWTR(N)=DELTA(N)
            VICDT(I)=VICDT(I)-DELTA(N)
            IF (VICDT(I) .LT. 0.0) THEN
!              CHANGE IN ICE IS GREATER THAN ICE CONTENT -- ADJUST WATER
!              CONTENT FOR THE DIFFERENCE, IF NOT GREATER THAN VLCDT(I)
               IF ((VLCDT(I)+VICDT(I)*RHOI/RHOL) .GT. VLCDT(I)/2.) THEN
                  VLCDT(I) = VLCDT(I) + VICDT(I)*RHOI/RHOL
                  CALL MATVL1 (I,MATDT(I),VLCDT(I),col,row)
                 ELSE
                  VLCDT(I) = VLCDT(I)/2.
                  CALL MATVL1 (I,MATDT(I),VLCDT(I),col,row)
               END IF
               VICDT(I)=0.0
               ICESDT(I)=0
            END IF
          ELSE
!
!           SAVE VALUE TO COMPARE RELATIVE CHANGE IN MATRIC POTENTIAL;
!           IF MATRIC POTENTIAL IS NEAR ZERO, SET CHECK VALUE TO 1.0
            CHKMAT=MATDT(I)
            IF (ABS(CHKMAT).LT.1.0) CHKMAT=1.0
!
            IF (ABS(DELTA(N)/CHKMAT).GT.100. .AND. NDT.LT.MAXNDT) THEN
!              TOO LARGE OF A CHANGE -- CONVERGENCE WILL NOT LIKELY BE
!              BE MET WITH CURRENT TIME STEP.  CUT TIME STEP IN HALF
               IWFLAG=1
               ITER=11
               GO TO 350
            END IF
            IF (ABS(DELTA(N)/CHKMAT) .GT. TOLER) IWFLAG = IWFLAG+1
            IF (ITER .GT. 3 .AND. DELWTR(N)*DELTA(N) .LT. 0) THEN
!              DELTA IS JUMPING BETWEEN NEG. AND POS.--CUT DELTA IN HALF
!              TO SUPPRESS TENDENCY TO JUMP AROUND SOLUTION
               DELTA(N)=DELTA(N)/2.
            END IF
            DELWTR(N)=DELTA(N)
            !IF(MATDT(I)-DELTA(N).gt. 10.0 .and. DELTA(N).le.0.0 )DELTA(N)=DELTA(N)/2.0
            MATDT(I)=MATDT(I)-DELTA(N)

         END IF
         CALL MATVL2 (I,MATDT(I),VLCDT(I),col,row)
!         IF (LEVEL(1) .GE. 2) WRITE (21,*) DELTA(N),MATDT(I),VLCDT(I),
!     >                                    VICDT(I)
         N=N+1
  330 CONTINUE
!     LIMIT WATER POTENTIAL AT SURFACE TO FREE WATER CONDITION 
!     TO ALLOW SEEPAGE TO OCCUR

  331 continue
      IF (MATDT(1).GT.0.0) THEN
!       CHANGE ANY SUBSEQUENT SATURATED NODES BY SAME AMOUNT
        DO 335 I=2,NS
          IF (MATDT(I).GT.0.0) THEN
            MATDT(I) = MATDT(I) - MATDT(1)
            CALL MATVL2 (I,MATDT(I),VLCDT(I),col,row)
           ELSE
!           DONE WITH SATURATED NODES
            GO TO 340
          END IF
  335   CONTINUE
        !POND=POND+MATDT(1)
        !RUNOFF=RUNOFF+MATDT(1)
  340   MATDT(1)=0
        CALL MATVL2 (1,MATDT(1),VLCDT(1),col,row)
      END IF


!Insert the code of water balance here after
!}}}}}}}}}}}}}}}}}}





















!-----------------------------------------------------------------------
  350 IF (IWFLAG .GT. 0 .OR. IEFLAG .GT. 0) THEN
!        CONVERGENCE HAS NOT BEEN MET - IF ITERATIONS ARE UNDER 10, GO
!        BACK AND START NEXT ITERATION
         IF (ITER .LE. maxiter) GO TO 200
!
!        HAVING PROBLEMS REACHING CONVERGENCE WITH CURRENT TIME STEP
         IF (NDT .LT. MAXNDT) THEN
!           CUT TIME STEP IN HALF AND TRY AGAIN
            NDT=NDT*2
            NTIMES=NTIMES*2 - 1
            DT=DTIME/NDT
!           REDEFINE END-OF-TIME-STEP VALUES
            CALL BACKUP (NS,NR,NSP,NC,NPLANT,NSALT,ICES,ICESDT,TS,TSDT, &
                MAT,MATDT,CONC,CONCDT,VLC,VLCDT,VIC,VICDT,SALT,SALTDT, &
                TR,TRDT,VAPR,VAPRDT,GMC,GMCDT,TC,TCDT,TLC,TLCDT,VAPC, &
                VAPCDT,WCAN,WCANDT,PCAN,PCANDT,TSP,TSPDT,DLW,DLWDT, &
                ICESP,ICESPT)
!           GO BACK AND CALCULATE BOUNDARY CONDITIONS FOR HALF TIME-STEP
            GO TO 120
         END IF

!
!        NDT > MAXNDT  --> INDICATE CONVERGENCE PROBLEMS AND CONTINUE
!         IF (LVLOUT(13).NE.0) WRITE (6,535) JULIAN,HOUR,YEAR,NTIMES,NDT
!         WRITE (21,*) ' CONVERGENCE PROBLEMS AT : ',
!     >            JULIAN,HOUR,YEAR,NTIMES
!         IF (IEFLAG.GT.0) WRITE (21,*)
!     >      ' ENERGY BALANCE WILL NOT CONVERGE'
!         IF (IWFLAG.GT.0) WRITE (21,*)
!     >      ' WATER BALANCE WILL NOT CONVERGE'
         ITER = 0
!         IF (LEVEL(1) .EQ. 1 .OR. LEVEL(1) .EQ. 2) STOP
      END IF
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!{{{{{{{{{{{{{{{{{{
!Insert the code of salt solute here




!Insert the code of salt solute here
!}}}}}}}}}}}}}}}}}}







!**** END OF ITERATION FOR TIME STEP ***********************************
!
      IF (MINSTP.EQ.0) THEN
         MAXSTP=NDT
         MINSTP=NDT
      END IF
      IF (MAXSTP .LT. NDT) MAXSTP=NDT
      IF (MINSTP .GT. NDT) MINSTP=NDT
!
!     SUM THE NECESSARY FLUXES OCCURRING OVER THE TIME STEP
      CALL SUMDT (NC,NPLANT,NSP,NR,NS,HFLUX,VFLUX,GFLUX,LWCAN,LWSNOW,&
          LWRES,LWSOIL,SWCAN,SWSNOW,SWRES,SWSOIL,QVC,QVR,QVSP,QSL,QSV,&
          TRNSP,XTRACT,SEEP,TSWSNO,TLWSNO,TSWCAN,TLWCAN,TSWRES,TLWRES,&
          TSWSOI,TLWSOI,THFLUX,TGFLUX,EVAP1,ETSUM,TSEEP,ROOTXT,TOTFLO,&
          NTIMES,TOPSNO,TQVSP,col,row,AbsorbedSW,AbsorbedLW)

!
!
      IF (NTIMES .NE. NDT) THEN
!        NOT REACHED END OF TIME STEP -- GO BACK THROUGH SUB-TIME-STEP
         IF (MOD(NTIMES,2).EQ.0 .AND. (IEFLAG+IWFLAG).EQ.0) THEN
!           CHECK IF IT'S WORTH IT TO TRY TO DOUBLE THE TIME STEP
            IF (NDT .LE. MAXDBL) THEN
               IF (NDT .LT. MAXDBL) MAXTRY=0
               IF (NDT .LE. 8) MAXTRY=MAXTRY+2
               MAXTRY=MAXTRY+1
               MAXDBL=NDT
            END IF
            IF (NDT.GT.MAXDBL .OR. MAXTRY.LE.3 .OR. NDT.GE.MAXNDT) THEN
!              ATTEMPT TO INCREASE LENGTH OF SUB-TIME-STEP
               NDT=NDT/2
               NTIMES=NTIMES/2
               DT=DTIME/NDT
            END IF
         END IF
         NTIMES=NTIMES + 1
         GO TO 110
        ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!ground water balance
!------------------------------------------------------       
          PROSAT=SAT(NS)
          GWS=0.30
          DGL=DGL-TOTFLO(NS-1)/GWS

 

          hhh =0 
          if (hhh.eq.1)then
            IF (DGL .LT. ZS(NS)) THEN
              RE=(ZS(NS)-DGL)*GWS/(0.5*(ZS(NS)-ZS(NS-1)))
              IF(VLCDT(NS)+VICDT(NS)+RE .GT. PROSAT) THEN
                DGL=ZS(NS-1)
                VLCDT(NS)=PROSAT*VLCDT(NS)/(VLCDT(NS)+VICDT(NS))
                VICDT(NS)=PROSAT*VLCDT(NS)/(VLCDT(NS)+VICDT(NS))
              ELSE
                DGL = ZS(NS)
                VLCDT(NS)=(VLCDT(NS)+VICDT(NS)+RE)*VLCDT(NS)/(VLCDT(NS)+VICDT(NS))
                VICDT(NS)=(VLCDT(NS)+VICDT(NS)+RE)*VLCDT(NS)/(VLCDT(NS)+VICDT(NS)) 
              END IF
            END IF

            PROSAT=SAT(NS)
            GWS=0.30
            DVLC=TOTFLO(NS-1)/(ZS(NS)-ZS(NS-1))/2.0
            VLCDT(NS)=(VLC(NS)+VLCDT(NS))/2.0+DVLC         
            GDZ=DGL-(ZS(NS)+ZS(NS-1))/2.0
            CALL MATVL1(NS,MAT1(NS),VLCDT(NS),COL,ROW)         
            CALL MATVL1(NS,MAT2(NS),PROSAT,COL,ROW)
            IF(ENTRY(NS) .GT. MAT2(NS)) THEN
              HK2=SATK(NS)*(ENTRY(NS)/MAT2(NS))**(2.+3/B(NS))
            ELSE
              HK2=SATK(NS)
            END IF
            IF (VICDT(NS).GT. 0.0)THEN
              QMAX=(PROSAT-VLCDT(NS)-VICDT(NS))*(ZS(NS)-ZS(NS-1))/2./DTIME
              IF (QMAX.LT.0.0) QMAX=0.0
              CALL SOILHK(NS,HK,MAT1(NS),VLCDT,VICDT,COL,ROW)
              HK1=HK(NS)
              CONG=SQRT(HK1*HK2)/GDZ
              QSLG=CONG*(MAT1(NS)-MAT2(NS)+GDZ)
              IF (QSLG .LT.0 ) THEN
                IF (QSLG+QMAX .LT. 0.0) QSLG=-QMAX
              END IF
            ELSE
              IF(ENTRY(NS) .GT. MAT1(NS)) THEN
                HK1=SATK(NS)*(ENTRY(NS)/MAT1(NS))**(2.+3.0/B(NS))
              ELSE
                HK1=SATK(NS)
              END IF      
              QMAX=(PROSAT-VLCDT(NS))*(ZS(NS)-ZS(NS-1))/2./DTIME 
              IF (QMAX .LT. 0.0 ) QMAX=0.0
              CONG=SQRT(HK1*HK2)/GDZ
              QSLG=CONG*(MAT1(NS)-MAT2(NS)+GDZ)
              IF (QSLG .LT.0 ) THEN
                IF (QSLG+QMAX .LT. 0.0) QSLG=-QMAX
              END IF
            END IF
            VLCDT(NS)=VLCDT(NS)-QSLG*DTIME/(ZS(NS)-ZS(NS-1))/2.0
            REACH=QSLG*DTIME
!           IF(ISNAN(REACH)) THEN
!              PRINT*,COL,ROW,QSLG,VLCDT(NS),GDZ
!              STOP
!           END IF
            DGL=DGL-REACH/GWS
            DZDGL=REACH/GWS
            IF (DGL .LT. ZS(NS)) THEN
              RE=(ZS(NS)-DGL)*GWS/(0.5*(ZS(NS)-ZS(NS-1)))
              IF(VLCDT(NS)+VICDT(NS)+RE .GT. PROSAT) THEN
                DGL=ZS(NS)
                VLCDT(NS)=PROSAT*VLCDT(NS)/(VLCDT(NS)+VICDT(NS))
                VICDT(NS)=PROSAT*VLCDT(NS)/(VLCDT(NS)+VICDT(NS))
              ELSE
                DGL = ZS(NS)
                VLCDT(NS)=(VLCDT(NS)+VICDT(NS)+RE)*VLCDT(NS)/(VLCDT(NS)+VICDT(NS)) 
                VICDT(NS)=(VLCDT(NS)+VICDT(NS)+RE)*VLCDT(NS)/(VLCDT(NS)+VICDT(NS)) 
              END IF
           END IF

           
      end if
      !IF (DGL.GE.2.5) DGL=2.5
!         PRINT*,VLCDT(NS),ts(1:NS)
!         PRINT*,"11" 
!         WRITE(20,130),(VLCDT(I),I=2,NS-1)
!         WRITE(21,130),(VICDT(I),I=2,NS-1)

!         PRINT*,JTOP,DGL
!  130 FORMAT (50F7.4)
!------------------------------------------------------
!         WRITE(20,*),VLCDT(9),VLCDT(NS-1),TSDT(9),
!     >               TSDT(NS)
!
!         PRINT*,VLCDT(NS-1),VICDT(NS-1),VLCDT(NS),VICDT(NS),
!     >          DGL
!        CALCULATE THE RECHARGE TO GROUNDWATER   
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC







!        END OF TIME STEPS - INFILTRATE RAIN, DETERMINE OUTFLOW FROM
!        SNOW,AND CALCULATE WATER BALANCE.  SET TIME STEP TO ONE HOUR
!        FOR PRECIP CALULATIONS
         DT = DTIME
         RAIN=PRECIP
         TADT=TMPDAY
         HUMDT=HUMDAY
         TA=TMPDAY
         HUM=HUMDAY
         NSPLST=NSP
         CALL PRECP(NPLANT,NC,NSP,NR,NS,TA,TADT,HUM,HUMDT,&
         ZC,XANGLE,CLUMPNG,ITYPE,WCANDT,WCMAX,PCANDT,&
         ZSP,DZSP,TSPDT,BDLW,DLWDT,RHOSP,TQVSP,TOPSNO,WLAG,STORE,SNOWEX,&
         ZR,TRDT,GMCDT,GMCMAX,RHOR,&
         ZS,TSDT,VLCDT,VICDT,MATDT,TOTFLO,SALTDT,CONCDT,ICESPT,ICESDT,&
         RAIN,POND,RUNOFF,EVAP1,PONDMX,SNOTMP,SNODEN,DIRRES,&
         SLOPE,PRESUR,nsalt,col,row,infiltration,snowmelt)
!
!        ADD SEEPAGE TO RUNOFF
         RUNOFF=RUNOFF+TSEEP
!
!        ADJUST DEPTHS OF CANOPY LAYERS FOR ANY CHANGE IN SNOWPACK
         IF (NSP.GT.0 .OR. NSPLST.GT.0) THEN
            IF (NPLANT.NE.0) CALL CANOPY (NPLANT,NC,NSP,NSPLST,MZCINP,&
             ITYPE,1,PLTHGT,PLTWGT,PLTLAI,RLEAF0,ZSP,ZC,TC,TCDT,&
             TLC,TLCDT,BVAPC,VAPCDT,BWCAN,WCANDT,WCMAX,PCANDT,CANMA,&
             CANMB,TMPDAY,HUMDAY,col,row)
         END IF
!
!        PRINT OUT WATER BALANCE FOR THIS HOUR
!         IF (LVLOUT(7) .NE. 0)
!     >    CALL WBALNC (NPLANT,NC,NSP,NR,NS,LVLOUT(7),JULIAN,HOUR,YEAR,
!     >     ITYPE,1,ZC,BWCAN,WCANDT,BPCAN,PCANDT,BVAPC,VAPCDT,RHOSP,DZSP,
!     >     DLWDT,WLAG,STORE,ZR,BGMC,GMCDT,BVAPR,VAPRDT,RHOR,ZS,BVLC,
!     >     VLCDT,BVIC,VICDT,TOTFLO,PRECIP,RUNOFF,POND,EVAP1,ETSUM)
!
!        PRINT OUT ENERGY BALANCE FOR THIS HOUR
!         IF (LVLOUT(6) .NE. 0)
!     >     CALL ENERGY (NSPLST,LVLOUT(6),JULIAN,HOUR,YEAR,INITAL,DTIME,
!     >      TSWSNO,TLWSNO,TSWCAN,TLWCAN,TSWRES,TLWRES,TSWSOI,TLWSOI,
!     >      THFLUX,TGFLUX,EVAP1)
!
!        PRINT OUT FROST AND SNOW DEPTH
      CALL FROST (NSP,NS,JULIAN,&
         HOUR,YEAR,1,ZSP,RHOSP,DZSP,DLWDT,WLAG,STORE,&
         ZS,VLCDT,VICDT,TSDT,ICESDT,col,row)
         
         if(NSP .gt. 0 .and. NSP.lt.NSPMAX)then
           ZSP(NSPMAX)=ZSP(NSP+1)
         else
           ZSP(NSPMAX)=0.0
         end if
!        PRINT OUTPUT FOR THIS HOUR
!         CALL OUTPUT (NPLANT,NC,NSP,NR,NS,LVLOUT,0,INPH2O,
!     >    JULIAN,HOUR,YEAR,1,ZC,TCDT,VAPCDT,WCANDT,ROOTXT,RHOSP,
!     >    ZSP,TSPDT,DLWDT,ZR,TRDT,VAPRDT,GMCDT,ZS,TSDT,VLCDT,VICDT,
!     >    MATDT,TOTFLO,CONCDT,SALTDT,ndt)
!xxxx
!    >    MATDT,TOTFLO,CONCDT,SALTDT)
!
!        PRINT OUTPUT TO SCREEN
!         IF (LVLOUT(13) .NE. 0) THEN
!            NPRINT=NPRINT+1
!            IF (MOD(NPRINT,LVLOUT(13)).EQ.0) THEN
!               WRITE (6,530) JULIAN,HOUR,YEAR,MINSTP,MAXSTP
!               NPRINT=0
!               MAXSTP=0
!               MINSTP=0
!            END IF
!         END IF
!
         RETURN
!
      END IF












158 FORMAT (5X,F8.4,F8.4,F10.4,F10.4,F10.4,2X,4F7.3,F7.3,F8.3)
159 FORMAT ('JDAY,JH,HYR,TMP,WIND,HUM,SUNHOR:',I4,I3,I5,2X,F8.2,F8.2,F8.2,F8.2,F8.2,F8.2)
522 FORMAT (I5,F5.2,F8.3,F6.3,F6.3,F8.2,E14.5,F8.3,F8.3,E14.5)
530 FORMAT('+ Completed :      ',I4,3X,I4,3X,I4,5X,I4,6X,I4)
535 FORMAT('+ Converg. Prob. : ',I4,3X,I4,3X,I4,5X,I4,6X,I4)

END SUBROUTINE GOSHAW

end module
