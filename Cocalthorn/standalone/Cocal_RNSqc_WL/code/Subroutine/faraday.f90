subroutine faraday
  use phys_constant,  only : long
  use make_array_3d
  use grid_parameter, only : nrg, ntg, npg
  use def_metric,  only : alph, bvxu, bvyu, bvzu
  use def_emfield,  only : alva, vaxd, vayd, vazd
  use def_faraday_tensor,  only : fxd, fyd, fzd, fijd
  use def_emfield_derivatives, only : Lie_bAxd, Lie_bAyd, Lie_bAzd
  use interface_grgrad_midpoint
  use interface_interpo_linear_type0
  implicit none
  integer :: irg, itg, ipg
  real(long) :: alpgc, ainvh, bvxugc, bvyugc, bvzugc
  real(long) :: vaxdgc, vaydgc, vazdgc
  real(long) :: dbvxdxgc, dbvxdygc, dbvxdzgc, &
  &             dbvydxgc, dbvydygc, dbvydzgc, & 
  &             dbvzdxgc, dbvzdygc, dbvzdzgc 
  real(long) :: dvaxdxgc, dvaxdygc, dvaxdzgc, & 
  &             dvaydxgc, dvaydygc, dvaydzgc, &
  &             dvazdxgc, dvazdygc, dvazdzgc 
  real(long) :: dalvadxgc, dalvadygc, dalvadzgc
  real(long) :: lie_bAx, lie_bAy, lie_bAz
  real(long), pointer :: dbvxdx(:,:,:), dbvxdy(:,:,:), dbvxdz(:,:,:)
  real(long), pointer :: dbvydx(:,:,:), dbvydy(:,:,:), dbvydz(:,:,:)
  real(long), pointer :: dbvzdx(:,:,:), dbvzdy(:,:,:), dbvzdz(:,:,:)
  real(long), pointer :: dvaxdx(:,:,:), dvaxdy(:,:,:), dvaxdz(:,:,:)
  real(long), pointer :: dvaydx(:,:,:), dvaydy(:,:,:), dvaydz(:,:,:)
  real(long), pointer :: dvazdx(:,:,:), dvazdy(:,:,:), dvazdz(:,:,:)
  real(long), pointer :: dalvadx(:,:,:), dalvady(:,:,:), dalvadz(:,:,:)
!
! --- Compute components of faraday tensor
! --- whose values are assigned on the mid points. 
!
  call alloc_array3d(dbvxdx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dbvxdy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dbvxdz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dbvydx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dbvydy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dbvydz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dbvzdx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dbvzdy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dbvzdz, 0, nrg, 0, ntg, 0, npg)
  call grgrad_midpoint(bvxu,dbvxdx,dbvxdy,dbvxdz)
  call grgrad_midpoint(bvyu,dbvydx,dbvydy,dbvydz)
  call grgrad_midpoint(bvzu,dbvzdx,dbvzdy,dbvzdz)
!
  call alloc_array3d(dvaxdx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dvaxdy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dvaxdz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dvaydx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dvaydy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dvaydz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dvazdx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dvazdy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dvazdz, 0, nrg, 0, ntg, 0, npg)
  call grgrad_midpoint(vaxd,dvaxdx,dvaxdy,dvaxdz)
  call grgrad_midpoint(vayd,dvaydx,dvaydy,dvaydz)
  call grgrad_midpoint(vazd,dvazdx,dvazdy,dvazdz)
!
  call alloc_array3d(dalvadx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dalvady, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(dalvadz, 0, nrg, 0, ntg, 0, npg)
  call grgrad_midpoint(alva,dalvadx,dalvady,dalvadz)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(alpgc,alph,irg,itg,ipg)
        ainvh = 1.0d0/alpgc    ! 1/\alpha
        call interpo_linear_type0(bvxugc,bvxu,irg,itg,ipg)
        call interpo_linear_type0(bvyugc,bvyu,irg,itg,ipg)
        call interpo_linear_type0(bvzugc,bvzu,irg,itg,ipg)
        call interpo_linear_type0(vaxdgc,vaxd,irg,itg,ipg)
        call interpo_linear_type0(vaydgc,vayd,irg,itg,ipg)
        call interpo_linear_type0(vazdgc,vazd,irg,itg,ipg)
        dbvxdxgc = dbvxdx(irg,itg,ipg) ; dvaxdxgc = dvaxdx(irg,itg,ipg)
        dbvxdygc = dbvxdy(irg,itg,ipg) ; dvaxdygc = dvaxdy(irg,itg,ipg)
        dbvxdzgc = dbvxdz(irg,itg,ipg) ; dvaxdzgc = dvaxdz(irg,itg,ipg)
        dbvydxgc = dbvydx(irg,itg,ipg) ; dvaydxgc = dvaydx(irg,itg,ipg)
        dbvydygc = dbvydy(irg,itg,ipg) ; dvaydygc = dvaydy(irg,itg,ipg)
        dbvydzgc = dbvydz(irg,itg,ipg) ; dvaydzgc = dvaydz(irg,itg,ipg)
        dbvzdxgc = dbvzdx(irg,itg,ipg) ; dvazdxgc = dvazdx(irg,itg,ipg)
        dbvzdygc = dbvzdy(irg,itg,ipg) ; dvazdygc = dvazdy(irg,itg,ipg)
        dbvzdzgc = dbvzdz(irg,itg,ipg) ; dvazdzgc = dvazdz(irg,itg,ipg)
        dalvadxgc = dalvadx(irg,itg,ipg)
        dalvadygc = dalvady(irg,itg,ipg)
        dalvadzgc = dalvadz(irg,itg,ipg)
!
        lie_bAx = bvxugc*dvaxdxgc + bvyugc*dvaxdygc + bvzugc*dvaxdzgc &
        &       + vaxdgc*dbvxdxgc + vaydgc*dbvydxgc + vazdgc*dbvzdxgc
        lie_bAy = bvxugc*dvaydxgc + bvyugc*dvaydygc + bvzugc*dvaydzgc &
        &       + vaxdgc*dbvxdygc + vaydgc*dbvydygc + vazdgc*dbvzdygc
        lie_bAz = bvxugc*dvazdxgc + bvyugc*dvazdygc + bvzugc*dvazdzgc &
        &       + vaxdgc*dbvxdzgc + vaydgc*dbvydzgc + vazdgc*dbvzdzgc
!
        Lie_bAxd(irg,itg,ipg) = lie_bAx
        Lie_bAyd(irg,itg,ipg) = lie_bAy
        Lie_bAzd(irg,itg,ipg) = lie_bAz
!
        fxd(irg,itg,ipg) = ainvh*(lie_bAx - dalvadxgc)
        fyd(irg,itg,ipg) = ainvh*(lie_bAy - dalvadygc)
        fzd(irg,itg,ipg) = ainvh*(lie_bAz - dalvadzgc)
        fijd(irg,itg,ipg,1) = dvaydxgc - dvaxdygc
        fijd(irg,itg,ipg,2) = dvazdxgc - dvaxdzgc
        fijd(irg,itg,ipg,3) = dvazdygc - dvaydzgc
!
      end do
    end do
  end do
!
  deallocate(dalvadx)
  deallocate(dalvady)
  deallocate(dalvadz)
  deallocate(dvaxdx)
  deallocate(dvaxdy)
  deallocate(dvaxdz)
  deallocate(dvaydx)
  deallocate(dvaydy)
  deallocate(dvaydz)
  deallocate(dvazdx)
  deallocate(dvazdy)
  deallocate(dvazdz)
!
  deallocate(dbvxdx)
  deallocate(dbvxdy)
  deallocate(dbvxdz)
  deallocate(dbvydx)
  deallocate(dbvydy)
  deallocate(dbvydz)
  deallocate(dbvzdx)
  deallocate(dbvzdy)
  deallocate(dbvzdz)
!
end subroutine faraday
