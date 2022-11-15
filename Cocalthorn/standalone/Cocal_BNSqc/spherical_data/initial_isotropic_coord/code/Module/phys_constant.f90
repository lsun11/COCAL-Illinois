module phys_constant
  implicit none
  Integer, parameter :: long = SELECTED_REAL_KIND(15,307)
  real(long), parameter :: pi = 3.14159265358979d+0
  real(long), parameter :: g = 6.67428d-08
  real(long), parameter :: c = 2.99792458d+10
  real(long), parameter :: msol   = 1.98892d+33
  real(long), parameter :: solmas = 1.98892d+33
  integer, parameter :: nnrg = 500, nntg = 128, nnpg = 256, nnlg = 30
  integer, parameter :: nnrf = 200, nntf = 128, nnpf = 256, nnlf = 30
  integer, parameter :: nnpeos = 10
end module phys_constant
