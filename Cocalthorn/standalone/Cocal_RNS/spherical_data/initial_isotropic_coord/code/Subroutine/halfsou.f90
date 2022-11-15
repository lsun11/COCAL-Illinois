subroutine halfsou(soug)
!
  use phys_constant, only : nnrg
  use grid_parameter_1D, only : nrg
  use coordinate_grav_r_1D, only : rg, hrg
  implicit none
!
  real(8), external      :: fn_lagint
  real(8), intent(inout) :: soug(0:nnrg)
  real(8) :: x4(4), f4(4)
  real(8) :: hhrr, sougb(0:nnrg)
  integer :: irg, irg0
!
  sougb(0:nrg) = soug(0:nrg)
!
  do irg = 1, nrg
    hhrr = hrg(irg)
    irg0 = min0(max0(irg-2,0),nrg-3) - 1
    x4(1:4) = rg(irg0+1:irg0+4)
    f4(1:4) = sougb(irg0+1:irg0+4)
    soug(irg) = fn_lagint(x4,f4,hhrr)
  end do
!
  soug(0) = 0.0d0
!
end subroutine halfsou
