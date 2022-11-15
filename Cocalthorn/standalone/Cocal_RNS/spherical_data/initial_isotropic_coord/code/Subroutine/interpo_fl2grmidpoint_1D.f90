subroutine interpo_fl2grmidpoint_1D(soug)
  use phys_constant, only : nnrg
  use grid_parameter_1D, only : nrg, nrf
  use coordinate_grav_r_1D, only : rg, hrg
  implicit none
  real(8), external :: fn_lagint, fn_lagint_2nd
  real(8), intent(inout) :: soug(0:nnrg)
  real(8) :: x(4), f(4), x2(2), f2(2), flv(0:nnrg)
  real(8) :: rrgg, small = 1.0d-14
  integer :: irf, irg, ir0, irf0, ir
!
  flv(0:nrg) = soug(0:nrg)
  soug(0:nrg)= 0.0d0
  do irg = 1, nrg
    if (hrg(irg).gt.rg(nrf)) exit
    do irf = 0, nrf
      if (hrg(irg).le.rg(irf)) then 
        irf0= irf-1
        ir0 = min0(max0(0,irf-2),nrf-3)
        exit
      end if
    end do
    if (irf0.eq.nrf-1) then
      x2(1:2) = rg(irf0:irf0+1)
      f2(1:2) = flv(irf0:irf0+1)
      rrgg = hrg(irg)
      soug(irg) = fn_lagint_2nd(x2,f2,rrgg)
    else 
      x(1:4) = rg(ir0:ir0+3)
      f(1:4) = flv(ir0:ir0+3)
      rrgg = hrg(irg)
      soug(irg) = fn_lagint(x,f,rrgg)
    end if
  end do
!!
!    open(26,file='frgaux2.dat',status='unknown')
!    do ir = 0, nrg
!      write(26,*) rg(ir), flv(ir)
!    end do
!      write(26,*) ' '
!    do ir = 1, nrg
!      write(26,*) hrg(ir), soug(ir)
!    end do
!    close(26)
!stop
end subroutine interpo_fl2grmidpoint_1D
