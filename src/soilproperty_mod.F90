module soilproperty_mod
 use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4
 use dims_mod
 public
 real(r8),dimension(NX,NY,NSMAX),target::B2d,ENTRY2d,SAT2d,thfc2d,thr2d,alpha2d,n2d,l2d
 real(r8),dimension(NX,NY,NSMAX),target::satk2d,RHOB2d,vapcof2d,vapexp2d,sand2d,silt2d,clay2d,om2d
 real(r8),dimension(NX,NY,NSALTMAX,NSMAX),target::nsalt2d,SALTKQ2d
end module soilproperty_mod
