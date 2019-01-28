module clayrs_mod
    use shr_kind_mod, only: r8 => shr_kind_r8
    use dims_mod
    implicit none
    real(r8),dimension(NX,NY,NPMAX,NCMAX-1),target::DRYCAN2d,CANLAI2d,RLEAF2d
    real(r8),dimension(NX,NY,NPMAX),target::TOTLAI2d,TOTROT2d
    integer(i4),dimension(NX,NY,NPMAX),target::IEVAP2d
    real(r8),dimension(NX,NY,NPMAX,NSMAX),target::RROOT2d,ROOTDN2d
    real(r8),dimension(NX,NY,NCMAX),target::ZZC2d

    real(r8),dimension(NX,NY,NPMAX),target::PXYLEM2d
    real(r8),dimension(NX,NY,NPMAX,NCMAX-1),target::RHCAN2d,ETCAN2d
    integer(i4),dimension(NX,NY,NPMAX),target::INIT2d
end module clayrs_mod
