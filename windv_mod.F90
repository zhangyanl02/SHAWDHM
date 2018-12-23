module windv_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use dims_mod
  implicit none
  real(r8),dimension(NX,NY),target::ZH2d,ZM2d,ZERO2d,USTAR2d,STABLE2d,CONRH2d,CONRV2d,ZMSUB2d,ZHSUB2d,ZERSUB2d
  real(r8),dimension(NX,NY,NCMAX),target::WINDC2d
  real(r8),dimension(NX,NY,NRMAX),target::WINDR2d
end module windv_mod


