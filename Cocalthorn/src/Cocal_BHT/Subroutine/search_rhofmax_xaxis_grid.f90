subroutine search_rhofmax_xaxis_grid(iremax)
  use phys_constant,  only  : long
  use grid_parameter, only  : nrf, ntfxy
  use def_matter, only  :  rhof 
  implicit none
  real(long) :: rhofmax
  integer    :: ir, iremax
!
! search maximum emd on positive x-axis
  rhofmax = rhof(0,ntfxy,0)
  iremax = 0
  do ir = 0, nrf
    if (rhof(ir,ntfxy,0).gt.rhofmax) then
      rhofmax = rhof(ir,ntfxy,0)
      iremax = ir
    end if
  end do
!
end subroutine search_rhofmax_xaxis_grid
