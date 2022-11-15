subroutine gauge_shift
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use def_matter,  only : rs
  use def_metric, only : psi, alph, bvxu, bvyu, bvzu
  use def_shift_derivatives, only : cdivbv
  use def_vector_x, only : hvec_xg
  use interface_poisson_solver_1bh_homosol
  use make_array_2d
  use make_array_3d
  use interface_grgrad_4th_gridpoint
  use interface_grgrad_midpoint_type0
  use interface_interpo_linear_type0
  implicit none
  real(long), pointer :: sou(:,:,:), gauge(:,:,:)
  real(long), pointer :: sou_bhsurf(:,:), dsou_bhsurf(:,:)
  real(long), pointer :: sou_outsurf(:,:), dsou_outsurf(:,:)
  real(long) :: psigc, alphgc, bvxgc, bvygc, bvzgc, alphinv, alpsinv, &
  &             divbv, dxpsi, dypsi, dzpsi, bvudpsi, x, y, z, traceK, &
  &             psi6inv
  real(long) :: dgadx, dgady, dgadz
  integer :: irg, itg, ipg
!
  call alloc_array3d(sou,0,nrg,0,ntg,0,npg)
  call alloc_array3d(gauge,0,nrg,0,ntg,0,npg)
  call alloc_array2d(sou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(sou_outsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_outsurf,0,ntg,0,npg)
!
! --- Impose extended Dirac gauge condition.
!
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
        x = hvec_xg(irg,itg,ipg,1)
        y = hvec_xg(irg,itg,ipg,2)
        z = hvec_xg(irg,itg,ipg,3)
        call kerr_schild_traceK(x,y,z,traceK)
        divbv   = cdivbv(irg,itg,ipg)
        bvudpsi = bvxgc*dxpsi + bvygc*dypsi + bvzgc*dzpsi
        sou(irg,itg,ipg)=psigc**6*(alphgc*traceK-divbv-6.0d0/psigc*bvudpsi)
!
!!if (itg.eq.ntg/2-2.and.ipg.eq.2) write(6,*) sou(irg,itg,ipg), &
!if (itg.eq.1.and.ipg.eq.2) write(6,*) x, y, z
!if (itg.eq.1.and.ipg.eq.2) write(6,*) sou(irg,itg,ipg), &
!&   alphgc*traceK, -divbv -6.0d0/psigc*bvudpsi
!if (itg.eq.1.and.ipg.eq.2) write(6,*) alphgc, traceK
!
      end do
    end do
  end do
!
   sou_bhsurf( 0:ntg,0:npg) = 0.0d0
  dsou_bhsurf( 0:ntg,0:npg) = 0.0d0
   sou_outsurf(0:ntg,0:npg) = 0.0d0
  dsou_outsurf(0:ntg,0:npg) = 0.0d0
  call poisson_solver_1bh_homosol('nd',sou, &
  &                         sou_bhsurf,dsou_bhsurf, & 
  &                        sou_outsurf,dsou_outsurf,gauge)
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
  deallocate(sou)
  deallocate(gauge)
  deallocate(sou_bhsurf)
  deallocate(dsou_bhsurf)
  deallocate(sou_outsurf)
  deallocate(dsou_outsurf)
end subroutine gauge_shift
