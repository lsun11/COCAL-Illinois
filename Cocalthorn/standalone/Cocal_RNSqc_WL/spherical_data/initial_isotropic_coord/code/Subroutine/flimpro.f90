subroutine flimpro(back,rnew,fffac,epsmax,irerr,iterr,iperr,isw)
!
! --- Compute non-variables for inteporation.
!
  use phys_constant, only : nnrg
  use grid_parameter_1D, only : nrf
  implicit none
!
  real(8), intent(inout) :: rnew(0:nnrg), back(0:nnrg), fffac, epsmax
  integer, intent(inout) :: irerr
  integer, intent(in)    :: isw
  real(8) :: devi, edet, edetb, error
  integer :: ir, iterr, iperr
!
! --- Set improved values for quantities on fluid-coordinate and 
! --- convergence check.  
!
  epsmax = 0.0d0
  do ir = 0, nrf - isw
    rnew(ir) = fffac*rnew(ir) + (1.d0-fffac)*back(ir)
    edet  = rnew(ir)
    edetb = back(ir)
    devi  = dabs(rnew(ir)) + dabs(back(ir))
    if (dabs(edet+edetb) >= 1.0d-8) then
      error = dabs(2.d0*(edet - edetb))/devi
      if(error  >  epsmax) then
        epsmax = error
        irerr = ir
      end if
    end if
  end do
!
end subroutine flimpro
