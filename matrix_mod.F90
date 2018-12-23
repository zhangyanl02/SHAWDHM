module matrix_mod
 use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8
 use dims_mod
 implicit none
 public
 real(r8),dimension(NX,NY,maxlayer),target::A12d,B12d,C12d,D12d,A22d,B22d,C22d,D22d
end module matrix_mod
