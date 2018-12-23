module swrcoe_mod
 use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
 implicit none

 real(r8),parameter::SOLCON=1360
 real(r8)::DIFATM
 real(r8)::DIFRES
 real(r8)::SNOCOF  !SNOCOF IS MINIMUM DEPTH OF SNOW FOR COMPLETE GROUND COVER
 real(r8)::SNOEXP

 data DIFATM/0.76/ DIFRES/0.667/ SNOCOF/0.03/  SNOEXP/1.0/
end module swrcoe_mod
