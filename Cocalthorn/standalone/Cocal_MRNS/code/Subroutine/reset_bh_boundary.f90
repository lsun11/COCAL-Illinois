subroutine reset_bh_boundary(char_bc)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : psi, bvxd, alph, bvyd
  implicit none
  character(len=1), intent(in) :: char_bc
! Reset boundary condition for BH
!  psi(0,0:ntg,0:npg) = 0.0d0
!  psi(nrg,0:ntg,0:npg) = 1.0d0
  if (char_bc.eq.'d') then
    psi(0,0:ntg,0:npg) = bvxd(0,0:ntg,0:npg)
    alph(0,0:ntg,0:npg) = bvyd(0,0:ntg,0:npg)
  endif
  psi(nrg,0:ntg,0:npg) = bvxd(nrg,0:ntg,0:npg)  ! this is psi=1+m/2.0r1+m/2.0r2
  alph(nrg,0:ntg,0:npg) = bvyd(nrg,0:ntg,0:npg) ! this is alph=(1-m/2.0r1-m/2.0r2)/psi
!
end subroutine reset_bh_boundary
