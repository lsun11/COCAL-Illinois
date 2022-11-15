module grid_parameter_cartesian
  use phys_constant, only : long
  implicit none
  integer, parameter :: nnx = 2401, nny = 2401, nnz = 2401
!  integer, parameter :: nnx = 801, nny = 801, nnz = 801
  integer :: nx, ny, nz, nstar   ! Cartesian coordinate grid points
  integer :: nxhalf, nx_mid      ! nx_mid is a half point for odd nx
  integer :: nx0, ny0, nz0       ! For fish eye
  real(long) :: facco, haba
  character(4) :: chpr
contains
subroutine read_parameter_cartesian
  implicit none
  open(1,file='rnspar_cartesian.dat',status='old')
  read(1,'(4i5)') nx, ny, nz, nstar
  read(1,'(3x,a7)') chpr
  read(1,'(3i5)') nx0, ny0, nz0
  read(1,'(1p,2e11.3)') facco, haba
  close(1)
  nxhalf = nx/2
  nx_mid = nx/2 + 1
end subroutine read_parameter_cartesian
end module grid_parameter_cartesian
