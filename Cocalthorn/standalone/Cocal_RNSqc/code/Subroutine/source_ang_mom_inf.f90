subroutine source_ang_mom_inf(sous,irg)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : tfkij
  use def_metric_excurve_grid, only : tfkij_grid
  use coordinate_grav_r, only : rg
  use def_vector_irg, only : hvec_irg_cm_phig, hvec_irg_cbh_xg
  implicit none
  real(long), pointer :: sous(:,:)
  real(long) :: ni, vphi_cm, Aij, work(2,2)
  integer    :: irg, itg, ipg, ia, ib
!
  call calc_vector_irg(2,irg)
!
!  write (6,*) 'In source_ang_mom_inf    irg=', irg

  do ipg = 1, npg
    do itg = 1, ntg
      sous(itg,ipg)=0.0d0
      do ib = 1, 3
        do ia = 1, 3
!         tfkij is calculated at the midpoints. Here theta,phi must be midpoints
!         but irg=mass_ir is an endpoint.
!          Aij = 0.5d0*(tfkij(irg,itg,ipg,ia,ib) + tfkij(irg+1,itg,ipg,ia,ib))
          work(1:2,1:2) = tfkij_grid(irg, itg-1:itg,ipg-1:ipg, ia,ib)
          call interpo_linear1p_type0_2Dsurf(Aij,work)

          ni = hvec_irg_cbh_xg(itg,ipg,ib)/rg(irg)
          vphi_cm = hvec_irg_cm_phig(itg,ipg,ia)
!
          sous(itg,ipg) = sous(itg,ipg) + Aij*vphi_cm*ni
        end do
      end do
    end do
  end do
!
end subroutine source_ang_mom_inf

