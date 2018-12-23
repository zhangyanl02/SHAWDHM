module forcing_mod
  use shr_kind_mod, only: r8 => shr_kind_r8
  use dims_mod
  implicit none
  SAVE
  real(r8),dimension(NX,NY,24)::SUNHOR2d,TMPDAY2d,WINDAY2d,HUMDAY2d,PRECIP2d,SNODEN2d,SOITMP2d,VLCDAY2d
  real(r8),dimension(NSMAX),24)::SOILXT2d
end module forcing_mod
