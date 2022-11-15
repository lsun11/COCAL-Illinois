subroutine source_EMenergy_axisym_WL(sou_MtorB,sou_MpolB,sou_MeleE)
  use phys_constant,  only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_metric,     only : psi
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_faraday_tensor, only : fijd_grid, fiju_grid, fidfiu_grid
  use def_SEM_tensor_EMF, only : rhoH_EMF, jmd_EMF, smijd_EMF, trsm_EMF
  use interface_interpo_linear_type0_meridian
  use make_array_3d
  implicit none
  integer :: irg, itg, ipg
  real(long), pointer :: sou_MtorB(:,:,:), sou_MpolB(:,:,:), sou_MeleE(:,:,:)
  real(long), pointer :: sou_MtorB_grid(:,:,:), sou_MpolB_grid(:,:,:), &
  &                      sou_MeleE_grid(:,:,:)
  real(long) :: pi4inv, pi8inv, val
  real(long) :: psigc, psigc4inv, psigc8inv
  real(long) :: fidfiugc, fijdgc(3,3), fijugc(3,3)
!
! --- Compute components of faraday tensor
! --- whose values are assigned on the mid points. 
!
  call alloc_array3d(sou_MtorB_grid,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sou_MpolB_grid,0,nrg,0,ntg,0,npg)
  call alloc_array3d(sou_MeleE_grid,0,nrg,0,ntg,0,npg)
!
  pi4inv = 1.0d0/(4.0d0*pi)
  pi8inv = 1.0d0/(8.0d0*pi)
!
  fijdgc(1:3,1:3) = 0.0d0
  fijugc(1:3,1:3) = 0.0d0
  ipg = 0
  do itg = 0, ntg
    do irg = 0, nrg
!
      psigc  = psi(irg,itg,ipg)
      psigc4inv = 1.0d0/psigc**4
      psigc8inv = 1.0d0/psigc**8
!
      fidfiugc = fidfiu_grid(irg,itg,ipg)
      fijdgc(1,2) = fijd_grid(irg,itg,ipg,1) ; fijdgc(2,1) = - fijdgc(1,2)
      fijdgc(1,3) = fijd_grid(irg,itg,ipg,2) ; fijdgc(3,1) = - fijdgc(1,3)
      fijdgc(2,3) = fijd_grid(irg,itg,ipg,3) ; fijdgc(3,2) = - fijdgc(2,3)
      fijugc(1,2) = fiju_grid(irg,itg,ipg,1) ; fijugc(2,1) = - fijugc(1,2)
      fijugc(1,3) = fiju_grid(irg,itg,ipg,2) ; fijugc(3,1) = - fijugc(1,3)
      fijugc(2,3) = fiju_grid(irg,itg,ipg,3) ; fijugc(3,2) = - fijugc(2,3)
!
      sou_MtorB_grid(irg,itg,ipg) = pi8inv*psigc8inv* fijdgc(1,3)*fijugc(1,3)
      sou_MpolB_grid(irg,itg,ipg) = pi8inv*psigc8inv*(fijdgc(1,2)*fijugc(1,2) &
      &                                             + fijdgc(3,2)*fijugc(3,2))
      sou_MeleE_grid(irg,itg,ipg) = pi8inv*psigc4inv* fidfiugc
!
    end do
  end do
!
  ipg = 0
  do itg = 1, ntg
    do irg = 1, nrg
      call interpo_linear_type0_meridian(val,sou_MtorB_grid,irg,itg,ipg)
      sou_MtorB(irg,itg,ipg) = val
      call interpo_linear_type0_meridian(val,sou_MpolB_grid,irg,itg,ipg)
      sou_MpolB(irg,itg,ipg) = val
      call interpo_linear_type0_meridian(val,sou_MeleE_grid,irg,itg,ipg)
      sou_MeleE(irg,itg,ipg) = val
    end do 
  end do 
!
  do ipg = 1, npg
    sou_MtorB(0:nrg,0:ntg,ipg) = sou_MtorB(0:nrg,0:ntg,0)
    sou_MpolB(0:nrg,0:ntg,ipg) = sou_MpolB(0:nrg,0:ntg,0)
    sou_MeleE(0:nrg,0:ntg,ipg) = sou_MeleE(0:nrg,0:ntg,0)
  end do 
!
  deallocate(sou_MtorB_grid)
  deallocate(sou_MpolB_grid)
  deallocate(sou_MeleE_grid)
!
end subroutine source_EMenergy_axisym_WL
