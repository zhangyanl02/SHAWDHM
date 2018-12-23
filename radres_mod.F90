module radres_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use dims_mod
  implicit none
  real(r8),dimension(NX,NY,NRMAX),target::TDIREC2d,TDIFFU2d
end module radres_mod
