module radcan_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use dims_mod
  implicit none
  real(r8),dimension(NX,NY,NCMAX-1),target::TDIRCC2d,TDIFFC2d
  real(r8),dimension(NX,NY,NPMAX+1,NCMAX-1),target::DIRKL2d,DIFKL2d
  real(r8),dimension(NX,NY,NPMAX),target::FBDU2d,FDDU2d
  real(r8),dimension(NX,NY,NPMAX,NCMAX-1),target::tlclwr2d
end module radcan_mod

