module radial_perm_fn_grav_1D
  use phys_constant, only : nnrg
  implicit none
  real(8) :: hfsn(0:nnrg, 0:nnrg)	!def_grfsn
contains
subroutine calc_hfsn
!
  use grid_parameter_1D, only : nrg
  use coordinate_grav_r_1D, only : hrg, rg
  implicit none
  integer :: irr, ir
!
! -------------------------------------
! --- Computation of functoin hfsn. ---
! -------------------------------------
!
!    hfsn(irr,ir) --> for r' < r
!                     for r < r'
!
  hfsn(0:nrg,0:nrg) = 0.0d0
!
  do ir  = 1, nrg
    do irr = 1, nrg
      if (hrg(irr) < rg(ir)) then
        hfsn(irr,ir) = 1.0d0/rg(ir)
      else
        hfsn(irr,ir) = 1.0d0/hrg(irr)
      end if
    end do
  end do
!
  hfsn(1:nrg,0) = 1.0d0/hrg(1:nrg)
!
end subroutine calc_hfsn
end module radial_perm_fn_grav_1D
