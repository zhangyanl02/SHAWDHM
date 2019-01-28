module savedata_mod
    use shr_kind_mod, only: r8 => shr_kind_r8
    use dims_mod
    implicit none
!   used in ebsnow
    real(r8),dimension(NX,NY,NSPMAX),target::QVSPT2d,CON2d,CSPT2d
    real(r8),dimension(NX,NY,NRMAX),target::CREST2d
    real(r8),dimension(NX,NY,NSMAX),target::QSLT2d,TKT2d,CST2d,QSVT2d
    !used in WBSOIL
    real(r8),dimension(NX,NY,NSMAX),target::QSLT_WB2d,QSVT_WB2d

!   used in subroutine ATSTAB
    real(r8),dimension(NX,NY),target::TMPAIR2d,VAPAIR2d,ZMLOG2d,ZHLOG2d,PSIM2d,PSIH2d


    integer(i4)::NSTEP,NSTART,LTMPDY,LTMPHR,LTMPYR,LVLCDY,LVLCHR,LVLCYR,LWTRDY,LWTRHR,LWTRYR,LWTRMX
    integer(i4),dimension(NPMAX)::LCANDY,LCANYR
    real(r8)::TMPLST,TMP,VLCLST,VLC
    real(r8),dimension(NSMAX)::TMP1,VLC1,FLUX
    real(r8),dimension(NPMAX)::ZCLST,DCHLST,WLST,TLALST,RDPLST
    integer(i4),dimension(NPMAX)::IFLAGC

!   used in SUBROUTINE EBCAN
    real(r8),dimension(NX,NY,NCMAX-1),target::CONT2d,CONT_WB2d

!   used in SUBROUTINE EBSOIL  in common block spheat
    real(r8),dimension(NX,NY,NSMAX),target::CS2d

!   used in soilTK function
    real(r8),dimension(NX,NY),target::WFAIRD2d,WFSAD2d,WFSID2d,WFCLD2d,WFOMD2d,WFICED2d,TKMA2d,&
        WFL2d,WFSA2d,WFSI2d,WFCL2d,WFOM2d,WFICE2d
    integer(i4),dimension(NX,NY),target ::IFIRST2d

!   used in snowmelt subroutine
    integer(i4),dimension(NX,NY),target::IFIRST_snowm2d,NLAG2d

!   used in frost subroutine
      integer(i4),dimension(NX,NY),target::LAST2d
      real(r8),dimension(NX,NY),target::SWE2d,FDEPTH2d,TDEPTH2d



!   used in sumdt
    real(r8),dimension(NX,NY,NPMAX),target::TRANSP2d

    integer(i4),dimension(NX,NY),target::MINSTP2d,MAXSTP2d


    DATA IFIRST2d /TOTGRID * 0/
    DATA MAXSTP2d/ TOTGRID *0/ MINSTP2d/ TOTGRID *0/
    data IFIRST_snowm2d/TOTGRID*0/
end module savedata_mod



