subroutine source_virial_gravity_CF(sou_Wgra)
  use phys_constant,  only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : psi, alph, tfkijkij, trk, &
  &                      bvxu, bvyu, bvzu, bvxd, bvyd, bvzd
!
  use interface_grgrad_midpoint
  use interface_interpo_linear_type0
  use make_array_3d
  use make_array_4d
  implicit none
  integer :: ia, ib, irg, itg, ipg
  real(long) :: psifc, psifc6, alpfc, pi4inv
  real(long) :: alpgc, psigc, bvxgc, bvygc, bvzgc, &
  &             psigc2, psigc6, alp2inv, psal2inv, aijaij, trkgc, &
  &             dapsi(3), daalph(3), dpsi2, dalph2, bvdal
  real(long), pointer :: sou_Wgra(:,:,:)
  real(long), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
  real(long), pointer :: grada(:,:,:,:), gradp(:,:,:,:)
!
! --- Compute source terms for virial relations.
! --- source for gravitational energies is defined on mid points.
!
  call alloc_array3d(dfdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdz,1,nrg,1,ntg,1,npg)
  call alloc_array4d(grada,1,nrg,1,ntg,1,npg,1,3)
  call alloc_array4d(gradp,1,nrg,1,ntg,1,npg,1,3)
!
  call grgrad_midpoint(psi,dfdx,dfdy,dfdz)
  gradp(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  gradp(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  gradp(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
  call grgrad_midpoint(alph,dfdx,dfdy,dfdz)
  grada(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  grada(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  grada(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
  pi4inv = 1.0d0/(4.0d0*pi)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        call interpo_linear_type0(alpgc,alph,irg,itg,ipg)
        call interpo_linear_type0(bvxgc,bvxu,irg,itg,ipg)
        call interpo_linear_type0(bvygc,bvyu,irg,itg,ipg)
        call interpo_linear_type0(bvzgc,bvzu,irg,itg,ipg)
!
        psigc2 = psigc**2
        psigc6 = psigc**6
        alp2inv = 1.0d0/alpgc**2
        psal2inv= (psigc/alpgc)**2
        aijaij = tfkijkij(irg,itg,ipg)
        trkgc  = trk(irg,itg,ipg)
        dapsi(1:3) = gradp(irg,itg,ipg,1:3)
        daalph(1:3) = grada(irg,itg,ipg,1:3)
!
        dpsi2  = 0.0d0
        dalph2 = 0.0d0
        dpsi2  =  dapsi(1)**2 +  dapsi(2)**2 +  dapsi(3)**2
        dalph2 = daalph(1)**2 + daalph(2)**2 + daalph(3)**2
        bvdal = bvxgc*daalph(1) +  bvygc*daalph(2) +  bvzgc*daalph(3)
!
        sou_Wgra(irg,itg,ipg) = pi4inv*(2.0d0*dpsi2 - psal2inv*dalph2 &
        &                     + 0.75d0*psigc6*(aijaij - 2.0d0*trkgc/3.0d0) &
        &                     + psigc6*alp2inv*trkgc*bvdal)
!
      end do
    end do
  end do
!
  deallocate(dfdx)
  deallocate(dfdy)
  deallocate(dfdz)
  deallocate(grada)
  deallocate(gradp)
!
end subroutine source_virial_gravity_CF
