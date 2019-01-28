module writeit_mod
 use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8
 use dims_mod
 implicit none
 real(r8),dimension(NX,NY),target::hflux12d,rh2d,hnc2d,xlenc2d,contk2d
 integer(i4),dimension(NX,NY),target::nnc2d
 real(r8),dimension(NX,NY,NPMAX,NCMAX-1),target::tleaf2d
end module writeit_mod
