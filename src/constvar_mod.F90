module constvar_mod
 use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8
 implicit none

!BLOCK //CONSTN// except PRESUR
 real(r8),parameter::LF=335000.0
 real(r8),parameter::LV=2500000.0
 real(r8),parameter::LS=2835000.0
 real(r8),parameter::G=9.81
 real(r8),parameter::UGAS=8.3143
 real(r8),parameter::RHOL=1000.
 real(r8),parameter::RHOI=920.0
 real(r8),parameter::RHOM=2650.0
 real(r8),parameter::RHOOM=1300.0
 real(r8),parameter::RHOA=1.25
 real(r8),parameter::CL=4200.0
 real(r8),parameter::CI=2100.0
 real(r8),parameter::CM=900.0
 real(r8),parameter::COM=1920.0
 real(r8),parameter::CA=1006.0
 real(r8),parameter::CV=1860.0
 real(r8),parameter::CR=1900.0
 real(r8),parameter::VONKRM=0.4
 real(r8),parameter::VDIFF=0.0000212
 real(r8),parameter::P0=101300.0
 real(r8),parameter::TKL=0.57
 real(r8),parameter::TKI=2.2
 real(r8),parameter::TKA=0.025
 real(r8),parameter::TKR=0.05
 real(r8),parameter::TKSA=8.8
 real(r8),parameter::TKSI=2.92
 real(r8),parameter::TKCL=2.92
 real(r8),parameter::TKOM=0.25

 !used in SOILTK
 real(r8),parameter::GAOM=0.5
 real(r8),parameter::GASA=0.144
 real(r8),parameter::GASI=0.144
 real(r8),parameter::GAC=0.125

 !
 real(r8),parameter::RESMA=-53.72
 real(r8),parameter::RESMB=1.32
 real(r8),parameter::RESMC=0.0
 real(r8),parameter::RESTKA=0.007
 real(r8),parameter::RESTKB=4.0

 ! COMMON BLOCK LWRCOF
 real(r8),parameter::STEFAN=5.6697E-08
 real(r8),parameter::EMATM1=0.261
 real(r8),parameter::EMATM2=0.000777
 real(r8),parameter::EMITC=0.95
 real(r8),parameter::EMITR=0.95
 real(r8),parameter::EMITSP=0.9
 real(r8),parameter::EMITS=0.95

end module constvar_mod
