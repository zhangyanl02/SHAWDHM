module controlpara_mod
  use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4
  use dims_mod
  implicit none
    public

!   namelist jobname
    character(200)::hydro_para_dir,output_dir,forcing_dir,lai_dir,result2_dir,result1_dir,simulation_dir,shadows_dir

!   timepara
    integer(i4)::mtstep,nhrpdt,jstart,hrstar,yrstar,jend,hrend,yrend
    real(r8)::hrnoon

!   control
    real(r8)::toler,height,pondmx,canma,canmb,wcandt1,snotmp
    integer(i4)::inph2o,mwatrxt,mpltgro,mzcinp,ivlcbc,itmpbc,restart,shadeef,nsalt
    integer(i4)::nplant

    integer(i4)::debug
    real(r8)::wt,wdt,dtime
    real(r8),dimension(nx,ny),target::dt2d,wt2d,wdt2d


    integer(i4),dimension(6)::level
    integer(i4),dimension(15)::lvlout
    real(r8)::wcmax
    integer(i4)::inital,deflate
    integer(i4)::maxstep
    integer(i4)::hydro_module
    integer(i4)::timezone
    integer(i4)::maxiter
    integer(i4)::daily_output
    integer(i4)::point_simulation
    integer(i4)::ntimes
    real(r8)::GroundTempGradients


    integer(i4)::T_soil_out,VLC_soil_out,VIC_soil_out,T_Residue_out,TSP_out,ZSP_out,DGL_out,EVAP_out
    integer(i4)::ET_out,Runoff_out,surfacerunoff_out,SWE_out,THFLUX_out,infill_out,snowmelt_out

  !inital: intial time
  !mtstep: meteo data time step 0:hourly,1:daily
  !mzcinp: 0:model will generate spacing of nodes within the canopy ||||not eq 0->user input of root distribution and node spacing within the canopy (assumed to be constant for entire simulation
  !mpltgro:  2:set flag for user input of node spacing   | else: model will generate spacing of nodes within the canopy
  !inph2o: not equal 1:input soil moisture is water content ||| 1:input soil moisture is matric potential
  !mwatrxt: 1->read soil sink data
  !ivlcbc: gt 0-> unit gradient is assumed for water flux at bottom boundary;||| 0->input water content specified for water flux at bottom boundar
  !itmpbc: le 0 ->lower temperature boundary specified from input file|||1->estimate temp at lower boundary by method of hirota et al 2002|||2->zero heat flux at bottom boundary;

end module controlpara_mod
