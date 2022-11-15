subroutine gauge_shift
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use def_matter,  only : rs
  use def_metric, only : psi, bvxu, bvyu, bvzu
  use interface_poisson_solver
  use make_array_3d
  use interface_interpo_linear_type0
  use interface_grgrad_4th_gridpoint
  use interface_grgrad_midpoint_type0
  implicit none
!
  real(long), pointer :: sou(:,:,:), gauge(:,:,:)
  real(long) :: psigc, bvxgc, bvygc, bvzgc, 
  &             divbv, dxpsi, dypsi, dzpsi, bvudpsi, psi6inv
  real(long) :: dgadx, dgady, dgadz
  integer :: irg, itg, ipg
!
  call alloc_array3d(sou,0,nrg,0,ntg,0,npg)
  call alloc_array3d(gauge,0,nrg,0,ntg,0,npg)
!
! --- Impose extended Dirac gauge condition.
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        call interpo_linear_type0(bvxgc,bvxu,irg,itg,ipg)
        call interpo_linear_type0(bvygc,bvyu,irg,itg,ipg)
        call interpo_linear_type0(bvzgc,bvzu,irg,itg,ipg)
        call grgrad_midpoint_type0(psi,dxpsi,dypsi,dzpsi,irg,itg,ipg)
        divbv   = cdivbv(irg,itg,ipg)
        bvudpsi = bvxgc*dxpsi + bvygc*dypsi + bvzgc*dzpsi
        sou(irg,itg,ipg)=psigc**6*(-divbv-6.0d0/psigc*bvudpsi)
!
      end do
    end do
  end do
!
  call poisson_solver(sou,gauge)
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        psi6inv = 1.0d0/psi(irg,itg,ipg)**6
        call grgrad_4th_gridpoint(gauge,dgadx,dgady,dgadz,irg,itg,ipg)
!
        bvxu(irg,itg,ipg) = bvxu(irg,itg,ipg) + psi6inv*dgadx
        bvyu(irg,itg,ipg) = bvyu(irg,itg,ipg) + psi6inv*dgady
        bvzu(irg,itg,ipg) = bvzu(irg,itg,ipg) + psi6inv*dgadz
!
      end do
    end do
  end do
!
! -- For get_Module_Subrotine_names.sh script
! call gauge_shift.f90_kerr_schild
!
  deallocate(sou)
  deallocate(gauge)
end subroutine gauge_shift
