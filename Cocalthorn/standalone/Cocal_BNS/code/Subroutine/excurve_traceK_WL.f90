subroutine excurve_traceK_WL
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg, hrg
  use def_metric, only : psi, alph, bvxu, bvyu, bvzu, trk
  use def_metric_excurve_grid, only : trk_grid
  use def_shift_derivatives,      only : cdivbv
  use def_shift_derivatives_grid, only : cdivbv_grid
  use interface_grgrad_midpoint_type0
  use interface_grgrad_gridpoint_type0
  use interface_interpo_linear_type0
  implicit none
  real(8) :: psigc, alphgc, bvxgc, bvygc, bvzgc, alphinv, alpsinv, &
  &          divbv, dxpsi, dypsi, dzpsi, bvudpsi
  integer :: ipg, irg, itg
!
! --- Compute trace of extringic curvature.  
!
!itg = ntg/2-2; ipg = 2
!open(16,file='test_traceK',status='unknown')
!  do irg = 1, nrg
!    write(16,'(1p,12e20.12)') rg(irg), hrg(irg), &
!    &   trk(irg,itg,ipg), trk_grid(irg,itg,ipg)
!  end do
!close(16)
!
! --- Mid points
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        call interpo_linear_type0(alphgc,alph,irg,itg,ipg)
        call interpo_linear_type0(bvxgc,bvxu,irg,itg,ipg)
        call interpo_linear_type0(bvygc,bvyu,irg,itg,ipg)
        call interpo_linear_type0(bvzgc,bvzu,irg,itg,ipg)
        call grgrad_midpoint_type0(psi,dxpsi,dypsi,dzpsi,irg,itg,ipg)
        alphinv = 1.0d0/alphgc
        alpsinv = 1.0d0/(alphgc*psigc)
        divbv   = cdivbv(irg,itg,ipg)
        bvudpsi = bvxgc*dxpsi + bvygc*dypsi + bvzgc*dzpsi
!
        trk(irg,itg,ipg) = alphinv*divbv + 6.0d0*alpsinv*bvudpsi
!
      end do
    end do
  end do
!
! --- Grid points
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        psigc  =  psi(irg,itg,ipg)
        alphgc = alph(irg,itg,ipg)
        bvxgc  = bvxu(irg,itg,ipg)
        bvygc  = bvyu(irg,itg,ipg)
        bvzgc  = bvzu(irg,itg,ipg)
        call grgrad_gridpoint_type0(psi,dxpsi,dypsi,dzpsi,irg,itg,ipg)
        alphinv = 1.0d0/alphgc
        alpsinv = 1.0d0/(alphgc*psigc)
        divbv   = cdivbv_grid(irg,itg,ipg)
        bvudpsi = bvxgc*dxpsi + bvygc*dypsi + bvzgc*dzpsi
!
        trk_grid(irg,itg,ipg) = alphinv*divbv + 6.0d0*alpsinv*bvudpsi
!
      end do
    end do
  end do
!
!itg = ntg/2-2; ipg = 2
!open(16,file='test_traceK2',status='unknown')
!  do irg = 1, nrg
!    write(16,'(1p,12e20.12)') rg(irg), hrg(irg), &
!    &   trk(irg,itg,ipg), trk_grid(irg,itg,ipg)
!  end do
!close(16)
!irg = 2; itg = ntg/2-2; ipg = 2
!write(6,'(1p,3e16.8)') trk(irg,itg,ipg), trk_grid(irg,itg,ipg)
!stop
end subroutine excurve_traceK_WL
