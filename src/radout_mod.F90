module radout_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use dims_mod
  implicit none
  real(r8),dimension(NX,NY,NCMAX+NRMAX),target::dir2d,down2d,up2d
  real(r8),dimension(NX,NY,2),target::dwnlwr2d,uplwr2d,dwnswr2d,upswr2d
  real(r8),dimension(NX,NY),target::albsrf2d
  real(r8),dimension(NX,NY),target::SW_on_soil2d,SW_on_snow2d
end module radout_mod
