module dims_mod
  use shr_kind_mod, only: r8 => shr_kind_r8,i4 => shr_kind_i4,i8 => shr_kind_i8
  implicit none
  public
  integer(i4), parameter ::  NX = 176
  integer(i4), parameter ::  NY = 144
  integer(i4), parameter ::  NSMAX = 30
  integer(i4), parameter ::  NCMAX = 11
  integer(i4), parameter ::  NRMAX = 10
  integer(i4), parameter ::  NSPMAX = 100
  integer(i4), parameter ::  NSALTMAX = 10
  integer(i4), parameter ::  NPMAX = 8
  integer(i4), parameter ::  TOTGRID = NX*NY
  integer(i4), parameter ::  nvtyps = 1
  integer(i4), parameter ::  maxlayer=200
  real(r8)   , parameter ::  DDX = 2000.0, DDY=2000.0
end module dims_mod
