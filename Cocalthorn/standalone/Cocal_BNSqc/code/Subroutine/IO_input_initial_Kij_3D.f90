subroutine IO_input_initial_Kij_3D
  use phys_constant, only : long, nnrg, nntg, nnpg
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_excurve_grid, only : tfkij_grid 
  use grid_temporary
  use interface_interpo_3D_initial_4th
  use interface_interpo_3D_initial_4th_surface
  use make_array_2d
  use make_array_3d
  implicit none
  real(long), pointer :: kxx_tmp(:,:,:),  kxy_tmp(:,:,:),  kxz_tmp(:,:,:), &
  &                      kyy_tmp(:,:,:),  kyz_tmp(:,:,:),  kzz_tmp(:,:,:)
  integer :: ir, it, ip, nrgmax, ntgmax, npgmax
!
!
  open(13,file='rnsgra_Kij_3D.las',status='old')
  read(13,'(5i5)') nrgtmp, ntgtmp, npgtmp
  nrgmax = max0(nrg,nrgtmp)
  ntgmax = max0(ntg,ntgtmp)
  npgmax = max0(npg,npgtmp)
  call alloc_array3d(kxx_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(kxy_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(kxz_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(kyy_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(kyz_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(kzz_tmp,0,nrgmax,0,ntgmax,0,npgmax)

  do ip = 0, npgtmp
    do it = 0, ntgtmp
      do ir = 0, nrgtmp
        read(13,'(1p,6e23.15)') kxx_tmp(ir,it,ip), &
        &                       kxy_tmp(ir,it,ip), &
        &                       kxz_tmp(ir,it,ip), &
        &                       kyy_tmp(ir,it,ip), &
        &                       kyz_tmp(ir,it,ip), &
        &                       kzz_tmp(ir,it,ip)
      end do
    end do
  end do
  close(13)

!
! --- Interpolation
!
  if (nrg.ne.nrgtmp.or.ntg.ne.ntgtmp.or.npg.ne.npgtmp) then
    write(6,*) "Interpolating to new grid..."
    call interpo_3D_initial_4th(kxx_tmp,'grco')
    call interpo_3D_initial_4th(kxy_tmp,'grco')
    call interpo_3D_initial_4th(kxz_tmp,'grco')
    call interpo_3D_initial_4th(kyy_tmp,'grco')
    call interpo_3D_initial_4th(kyz_tmp,'grco')
    call interpo_3D_initial_4th(kzz_tmp,'grco')
  end if
!
!
!     NEEDS WORK    ************************************8
!
!
!

  tfkij_grid(0:nrg,0:ntg,0:npg,1,1) = kxx_tmp(0:nrg,0:ntg,0:npg)
  tfkij_grid(0:nrg,0:ntg,0:npg,1,2) = kxy_tmp(0:nrg,0:ntg,0:npg)
  tfkij_grid(0:nrg,0:ntg,0:npg,1,3) = kxz_tmp(0:nrg,0:ntg,0:npg)
  tfkij_grid(0:nrg,0:ntg,0:npg,2,2) = kyy_tmp(0:nrg,0:ntg,0:npg)
  tfkij_grid(0:nrg,0:ntg,0:npg,2,3) = kyz_tmp(0:nrg,0:ntg,0:npg)
  tfkij_grid(0:nrg,0:ntg,0:npg,3,3) = kzz_tmp(0:nrg,0:ntg,0:npg)
!
  deallocate(kxx_tmp)
  deallocate(kxy_tmp)
  deallocate(kxz_tmp)
  deallocate(kyy_tmp)
  deallocate(kyz_tmp)
  deallocate(kzz_tmp)
!
end subroutine IO_input_initial_Kij_3D
