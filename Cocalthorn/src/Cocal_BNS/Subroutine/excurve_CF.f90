subroutine excurve_CF(cobj)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric,  only : alph, bvxd, bvyd, bvzd, tfkij, tfkijkij, trk
  use interface_grgrad_midpoint_r3rd_nsbh
  use interface_interpo_linear_type0
  use make_array_3d
  implicit none
  integer :: info
  integer :: ia, ib
  integer :: irg, itg, ipg
  real(long) :: fa23, diver
  real(long) :: alpgc, ainvh, cdivbv
  real(long), pointer :: dbvxdx(:,:,:), dbvxdy(:,:,:), dbvxdz(:,:,:)
  real(long), pointer :: dbvydx(:,:,:), dbvydy(:,:,:), dbvydz(:,:,:)
  real(long), pointer :: dbvzdx(:,:,:), dbvzdy(:,:,:), dbvzdz(:,:,:)
  character(len=2), intent(in) :: cobj
!
! --- Compute components of extringic curvature 
! --- whose values are assigned on the mid points. 
!
!test
  info = 0
!test
!
  fa23 = 2.0d0/3.0d0
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

  call grgrad_midpoint_r3rd_nsbh(bvxd,dbvxdx,dbvxdy,dbvxdz,cobj)
  call grgrad_midpoint_r3rd_nsbh(bvyd,dbvydx,dbvydy,dbvydz,cobj)
  call grgrad_midpoint_r3rd_nsbh(bvzd,dbvzdx,dbvzdy,dbvzdz,cobj)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(alpgc,alph,irg,itg,ipg)
        ainvh = 0.5d0/alpgc    ! 1/2\alpha
        cdivbv = dbvxdx(irg,itg,ipg) + dbvydy(irg,itg,ipg) &
               + dbvzdz(irg,itg,ipg)
        diver = fa23*cdivbv
!
        tfkij(irg,itg,ipg,1,1) = ainvh*(2.0d0*dbvxdx(irg,itg,ipg) &
        - diver)
        tfkij(irg,itg,ipg,2,2) = ainvh*(2.0d0*dbvydy(irg,itg,ipg) & 
        - diver)
        tfkij(irg,itg,ipg,3,3) = ainvh*(2.0d0*dbvzdz(irg,itg,ipg) & 
        - diver)
        tfkij(irg,itg,ipg,1,2) = ainvh*(dbvydx(irg,itg,ipg) & 
        + dbvxdy(irg,itg,ipg))
        tfkij(irg,itg,ipg,1,3) = ainvh*(dbvzdx(irg,itg,ipg) & 
        + dbvxdz(irg,itg,ipg))
        tfkij(irg,itg,ipg,2,3) = ainvh*(dbvzdy(irg,itg,ipg) & 
        + dbvydz(irg,itg,ipg))
        tfkij(irg,itg,ipg,2,1) = tfkij(irg,itg,ipg,1,2)
        tfkij(irg,itg,ipg,3,1) = tfkij(irg,itg,ipg,1,3)
        tfkij(irg,itg,ipg,3,2) = tfkij(irg,itg,ipg,2,3)
! 
        tfkijkij(irg,itg,ipg) = 0.0d0
        trk(irg,itg,ipg) = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            tfkijkij(irg,itg,ipg) = tfkijkij(irg,itg,ipg) & 
            + tfkij(irg,itg,ipg,ia,ib)*tfkij(irg,itg,ipg,ia,ib)
          end do
        end do
!
        if (tfkijkij(irg,itg,ipg) /= 0.0d0) info = 1
!
      end do
    end do
  end do
!
  if (info /= 1) write(6,*) ' ### Warning K_ij = 0 *** '
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
!!irg = nrg-3; itg = ntg/2; ipg = 2
!!write(6,*) 'shift sample'
!!write(6,*) bvxd(irg,itg,ipg), bvyd(irg,itg,ipg), bvzd(irg,itg,ipg)
!!write(6,*) 'excurve sample'
!!write(6,*) tfkij(irg,itg,ipg,1,1),tfkij(irg,itg,ipg,1,2),tfkij(irg,itg,ipg,1,3)
!!write(6,*) tfkij(irg,itg,ipg,2,1),tfkij(irg,itg,ipg,2,2),tfkij(irg,itg,ipg,2,3)
!!write(6,*) tfkij(irg,itg,ipg,3,1),tfkij(irg,itg,ipg,3,2),tfkij(irg,itg,ipg,3,3)
!!write(6,*) tfkijkij(irg,itg,ipg)
!
end subroutine excurve_CF
