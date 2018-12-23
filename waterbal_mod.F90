module waterbal_mod
    use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8,r4=> shr_kind_r4
    use dims_mod
    implicit none

    real(r8),dimension(NX,NY),target::RAIN2d,DPCAN2d,DCAN2d,DSNOW2d,DRES2d,DSOIL2d,&
    TRUNOF2d,POND22d,TPERC2d,TETSUM2d,TEVAP2d,CUMVAP2d,SWE2d
    DATA RAIN2d /TOTGRID*0.0/
    DATA DPCAN2d /TOTGRID*0.0/
    DATA DCAN2d /TOTGRID*0.0/
    DATA DSNOW2d /TOTGRID*0.0/
    DATA DRES2d /TOTGRID*0.0/
    DATA DSOIL2d /TOTGRID*0.0/
    DATA TRUNOF2d /TOTGRID*0.0/
    DATA POND22d /TOTGRID*0.0/
    DATA TPERC2d /TOTGRID*0.0/
    DATA TETSUM2d /TOTGRID*0.0/
    DATA TEVAP2d /TOTGRID*0.0/
    DATA CUMVAP2d /TOTGRID*0.0/
end module waterbal_mod
