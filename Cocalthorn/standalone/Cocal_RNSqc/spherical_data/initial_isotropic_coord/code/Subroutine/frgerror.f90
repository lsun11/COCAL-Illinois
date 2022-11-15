subroutine frgerror(backg,rnewg,epsmaxg,irger)
!
! --- Compute non-variables for inteporation.
!
  use phys_constant, only : nnrg
  use grid_parameter_1D, only : nrgin
  use grid_parameter_1D, only : nrg
  implicit none
!
  real(8), intent(in)    :: rnewg(0:nnrg),backg(0:nnrg)
  real(8), intent(out)   :: epsmaxg
  integer, intent(inout) :: irger
  real(8) :: devi, edet, edetb, error
  integer :: irg
!
!
! --- Set improved values for quantities on GR-coordinate and 
! --- convergence check.  
!
  epsmaxg = 0.0d0
  do irg = 0, nrg
    edet  = rnewg(irg)
    edetb = backg(irg)
    devi  = dabs(rnewg(irg)) + dabs(backg(irg))
!      
    if (irg <= nrgin.and.devi >= 1.0d-8) then
      error = dabs(2.d0*(edet - edetb))/devi
      if(error  >  epsmaxg) then
        epsmaxg = error
        irger = irg
      end if
    end if
  end do
!
end subroutine frgerror
