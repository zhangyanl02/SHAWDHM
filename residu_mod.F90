module residu_mod
 use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8
 use dims_mod
 implicit none
 real(r8),dimension(NX,NY,NRMAX),target::EVAP2d,EVAPK2d
end module residu_mod
