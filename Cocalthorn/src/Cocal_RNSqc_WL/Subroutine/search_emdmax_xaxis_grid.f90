subroutine search_emdmax_xaxis_grid(iremax)
  use phys_constant,  only  : long
  use grid_parameter, only  : nrf, ntfeq, npfxzp
  use def_matter, only  :   emd
  use make_array_1d
  implicit none
  real(long) :: emdmax
  integer    :: iremax
  real(long), pointer :: fnc_x(:)
!
! search maximum emd on positive x-axis
!
  call alloc_array1d(fnc_x,0,nrf)
  fnc_x(0:nrf)=emd(0:nrf,ntfeq,npfxzp)
  iremax = maxloc(fnc_x,DIM=1) - 1
  deallocate(fnc_x)
!
!!  emdmax = emd(0,ntfxy,0)
!!  iremax = 0
!!  do ir = 0, nrf
!!    if (emd(ir,ntfxy,0).gt.emdmax) then
!!      emdmax = emd(ir,ntfxy,0)
!!      iremax = ir
!!    end if
!!  end do
!
end subroutine search_emdmax_xaxis_grid
