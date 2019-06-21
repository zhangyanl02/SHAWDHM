module statevar_mod
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
    use dims_mod
    implicit none

    public
    save
    integer(i4),dimension(NX,NY)::NC2d,NSP2d,NR2d,NS2d
    integer(i4),dimension(NX,NY)::INBASIN2d,LANDUSE2d,RIVERNET
!cannopy
    real(r8),dimension(NX,NY,NPMAX)::PLTHGT2d,PLTWGT2d,PLTLAI2d,ROOTDP2d,DCHAR2d,TCCRIT2d,RSTOM02d,RSTEXP2d,PLEAF02d,&
        RLEAF02d,RROOT02d,CANALB2d,XANGLE2d,CLUMPNG2d,PCANDT2d
    integer(i4),dimension(NX,NY,NPMAX)::ITYPE2d
    real(r8),dimension(NX,NY,NCMAX)::ZC2d,TCDT2d,VAPCDT2d
    real(r8),dimension(NX,NY,NPMAX,NCMAX-1)::TLCDT2d
    real(r8),dimension(NX,NY,NCMAX-1)::WCANDT2d
    integer(i4),dimension(NX,NY)::NPLANT2d
!snow
    real(r8),dimension(NX,NY,NSPMAX)::ZSP2d,DZSP2d,RHOSP2d,TSPDT2d,DLWDT2d
    integer(r4),dimension(NX,NY,NSPMAX)::ICESPT2d
    real(r8),dimension(NX,NY,11)::WLAG2d
!hydrology
    real(r8),dimension(NX,NY)::STORE2d,SNOWEX2d,RUNOFF2d,RUNOFF12d
    real(r8),dimension(NX,NY)::THFLUX2d,VFLUX2d
!residue
    real(r8),dimension(NX,NY,NRMAX)::ZR2d,RHOR2d
    real(r8),dimension(NX,NY,NRMAX)::TRDT2d,VAPRDT2d,GMCDT2d
    real(r8),dimension(NX,NY)::RLOAD2d,ZRTHIK2d,COVER2d,ALBRES2d,RESCOF2d,DIRRES2d
    real(r8)::GMCMAX
!soil
    real(r8),dimension(NX,NY,NSMAX)::ZS2d,TSDT2d,VLCDT2d,VICDT2d,MATDT2d,HKDT2d
    real(r8),dimension(NX,NY,NSALTMAX,NSMAX)::CONCDT2d,SALTDT2d
    integer(i4),dimension(NX,NY,NSMAX)::ICESDT2d
    real(r8),dimension(NX,NY)::ALBDRY2d,ALBEXP2d
    real(r8),dimension(NX,NY,NSALTMAX)::DGRADE2d,SLTDIF2d
    real(r8),dimension(NX,NY,NSMAX)::ASALT2d,DISPER2d
!surface
    real(r8),dimension(NX,NY)::ZMSRF2d,ZHSRF2d,ZERSRF2d,ZMSP2d,ZHSP2d,POND2d
!location
    real(r8),dimension(NX,NY)::ALATUD2d,LONTITUD2d,ELEVATION2d,SLOPE2d,ASPECT2d,CLOUDS2d

    real(r8),dimension(NX,NY)::TSAVG2d
    real(r8),dimension(NX,NY)::HRNOON2d
    integer(i4)::maskflag

!forcing
    real(r8),dimension(NX,NY,24)::SUNHOR2d,TMPDAY2d,WINDAY2d,PRECIP2d,SNODEN2d,HUMDAY2d,shadow2d
    real(r8),dimension(NX,NY,24)::SOITMP2d,VLCDAY2d
    real(r8),dimension(NX,NY,24,NSMAX)::SOILXT2d
    real(r8),dimension(NX,NY,24)::PRESUR2d

!water balance
    real(r8),dimension(NX,NY)::ETSUM2d,EVAP12d,AbsorbedSW2d,AbsorbedLW2d,snowmelt2d,InDirect2d,InDiffuse2d
    real(r8),dimension(NX,NY,NSMAX)::TOTFLO2d

!radiation
    real(r8),dimension(NX,NY,NPMAX+1,NCMAX-1)::LWCAN2d,SWCAN2d
    real(r8),dimension(NX,NY)::LWSNOW2d,LWSOIL2d,SWSOIL2d
    real(r8),dimension(NX,NY,NRMAX)::LWRES2d,SWRES2d
    real(r8),dimension(NX,NY,NSPMAX)::SWSNOW2d


!time
    integer(i4)::JULIAN,HOUR,YEAR

    real(r8),dimension(NX,NY)::GRIDAREA2d,SLOPELEN2d
    real(r8),dimension(NX,NY)::SKYVIEW2d

    real(r8),dimension(NX,NY)::maskgrid
    real(r8),dimension(NX,NY)::runoff_inter,runoff_g,waterflx
    real(r8),dimension(NX,NY)::Dr2d,Drw2d,DGL2D,DS2D,DG2D,BOT2d
!daily output
    real(r8),allocatable::TSDT2d_day(:,:,:,:),VLCDT2d_day(:,:,:,:),VICDT2d_day(:,:,:,:)
    real(r8),allocatable::TRDT2d_day(:,:,:,:)
    real(r8),allocatable::TSPDT2d_day(:,:,:),zsp2d_day(:,:,:)
    real(r8),allocatable::EVAP12d_day(:,:,:),ETSUM2d_day(:,:,:),THFLUX2d_day(:,:,:)
    integer(i4),allocatable::NSP2d_day(:,:,:)
    real(r8),allocatable::RUNOFF2d_day(:,:,:),DGL2D_day(:,:,:),RUNOFFDIS_day(:,:,:)
    real(r8),allocatable::infill_day(:,:,:)
    real(r8),allocatable::snowmelt_day(:,:,:)

end module statevar_mod
