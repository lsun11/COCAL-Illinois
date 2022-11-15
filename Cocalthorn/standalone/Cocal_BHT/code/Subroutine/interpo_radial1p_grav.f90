subroutine interpo_radial1p_grav(grv,val,rv,it,ip)
  use phys_constant, only : long
  use grid_parameter, only : nrg
  use coordinate_grav_r, only : rg
  implicit none
  real(long), external :: lagint_4th
  real(long), pointer :: grv(:,:,:)
  real(long), intent(out) :: val
  real(long), intent(in)  :: rv
  integer, intent(in)     :: it, ip
  real(long) :: x(4), f(4)
  integer :: irg, ir0
!
  do irg = 0, nrg-1
    if (rv.le.rg(irg)) then 
      ir0 = min0(max0(0,irg-2),nrg-3)
      exit
    end if
  end do
  x(1:4) = rg(ir0:ir0+3)
  f(1:4) = grv(ir0:ir0+3,it,ip)
  val = lagint_4th(x,f,rv)
!
end subroutine interpo_radial1p_grav
