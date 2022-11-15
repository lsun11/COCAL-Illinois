subroutine IO_input_initial_3D_WL_BHT
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
!  use def_metric_dihiju, only : dihixu_grid, dihiyu_grid, dihizu_grid
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          gaugex, gaugey, gaugez
  use def_CTT_decomposition, only : wvxd, wvyd, wvzd, sigt
!  use def_metric_excurve_grid, only : trk_grid
  use grid_temporary, only : nrgtmp, ntgtmp, npgtmp, rgtmp, thgtmp, phigtmp
  use interface_interpo_3D_initial_4th
  use make_array_3d
  implicit none
  real(long), pointer :: gx_tmp(:,:,:),   gy_tmp(:,:,:),   gz_tmp(:,:,:)
  real(long), pointer :: sigt_tmp(:,:,:), wvxd_tmp(:,:,:), wvyd_tmp(:,:,:), wvzd_tmp(:,:,:)
  real(long), pointer :: hxxd_tmp(:,:,:), hxyd_tmp(:,:,:), hxzd_tmp(:,:,:),&
  &                      hyyd_tmp(:,:,:), hyzd_tmp(:,:,:), hzzd_tmp(:,:,:)
  integer :: irg, itg, ipg, nrgmax, ntgmax, npgmax
!
  nrgmax = max0(nrg,nrgtmp)
  ntgmax = max0(ntg,ntgtmp)
  npgmax = max0(npg,npgtmp)
  call alloc_array3d(hxxd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(hxyd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(hxzd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(hyyd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(hyzd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(hzzd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
!
!  call alloc_array3d(trk_tmp,  0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(gx_tmp,   0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(gy_tmp,   0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(gz_tmp,   0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(wvxd_tmp, 0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(wvyd_tmp, 0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(wvzd_tmp, 0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(sigt_tmp, 0,nrgmax,0,ntgmax,0,npgmax)
!
! --- Metric potentials.
  open(13,file='rnsgra_hij_3D.ini',status='unknown')
  read(13,'(5i5)') nrgtmp, ntgtmp, npgtmp
  do ipg = 0, npgtmp
    do itg = 0, ntgtmp
      do irg = 0, nrgtmp
        read(13,'(1p,6e23.15)') hxxd_tmp(irg,itg,ipg), &
        &                       hxyd_tmp(irg,itg,ipg), &
        &                       hxzd_tmp(irg,itg,ipg), &
        &                       hyyd_tmp(irg,itg,ipg), &
        &                       hyzd_tmp(irg,itg,ipg), &
        &                       hzzd_tmp(irg,itg,ipg)
      end do
    end do
  end do
  close(13)
!
! --- Gauge potentials.
  open(14,file='rnsgra_gauge_3D.ini',status='unknown')
  read(14,'(5i5)') nrgtmp, ntgtmp, npgtmp
  do ipg = 0, npgtmp
    do itg = 0, ntgtmp
      do irg = 0, nrgtmp
        read(14,'(1p,10e23.15)')  gx_tmp(irg,itg,ipg), &
        &                         gy_tmp(irg,itg,ipg), &
        &                         gz_tmp(irg,itg,ipg), &
        &                       wvxd_tmp(irg,itg,ipg), &
        &                       wvyd_tmp(irg,itg,ipg), &
        &                       wvzd_tmp(irg,itg,ipg), &
        &                       sigt_tmp(irg,itg,ipg)
      end do
    end do
  end do
  close(14)
!
  if (nrg.ne.nrgtmp.or.ntg.ne.ntgtmp.or.npg.ne.npgtmp) then
    write(6,*) "Interpolating to new grid..."
    call interpo_3D_initial_4th(hxxd_tmp,'grco')
    call interpo_3D_initial_4th(hxyd_tmp,'grco')
    call interpo_3D_initial_4th(hxzd_tmp,'grco')
    call interpo_3D_initial_4th(hyyd_tmp,'grco')
    call interpo_3D_initial_4th(hyzd_tmp,'grco')
    call interpo_3D_initial_4th(hzzd_tmp,'grco')
!    call interpo_3D_initial_4th(trk_tmp, 'grco')
    call interpo_3D_initial_4th(gx_tmp,  'grco')
    call interpo_3D_initial_4th(gy_tmp,  'grco')
    call interpo_3D_initial_4th(gz_tmp,  'grco')
    call interpo_3D_initial_4th(wvxd_tmp,'grco')
    call interpo_3D_initial_4th(wvyd_tmp,'grco')
    call interpo_3D_initial_4th(wvzd_tmp,'grco')
    call interpo_3D_initial_4th(sigt_tmp,'grco')
  end if
!
  hxxd(0:nrg,0:ntg,0:npg) = hxxd_tmp(0:nrg,0:ntg,0:npg)
  hxyd(0:nrg,0:ntg,0:npg) = hxyd_tmp(0:nrg,0:ntg,0:npg)
  hxzd(0:nrg,0:ntg,0:npg) = hxzd_tmp(0:nrg,0:ntg,0:npg)
  hyyd(0:nrg,0:ntg,0:npg) = hyyd_tmp(0:nrg,0:ntg,0:npg)
  hyzd(0:nrg,0:ntg,0:npg) = hyzd_tmp(0:nrg,0:ntg,0:npg)
  hzzd(0:nrg,0:ntg,0:npg) = hzzd_tmp(0:nrg,0:ntg,0:npg)
!
!  trk_grid(0:nrg,0:ntg,0:npg)    = trk_tmp(0:nrg,0:ntg,0:npg)
  gaugex(0:nrg,0:ntg,0:npg) = gx_tmp(0:nrg,0:ntg,0:npg)
  gaugey(0:nrg,0:ntg,0:npg) = gy_tmp(0:nrg,0:ntg,0:npg)
  gaugez(0:nrg,0:ntg,0:npg) = gz_tmp(0:nrg,0:ntg,0:npg)
!
  wvxd(0:nrg,0:ntg,0:npg) = wvxd_tmp(0:nrg,0:ntg,0:npg)
  wvyd(0:nrg,0:ntg,0:npg) = wvyd_tmp(0:nrg,0:ntg,0:npg)
  wvzd(0:nrg,0:ntg,0:npg) = wvzd_tmp(0:nrg,0:ntg,0:npg)
  sigt(0:nrg,0:ntg,0:npg) = sigt_tmp(0:nrg,0:ntg,0:npg)
!
  deallocate(hxxd_tmp)
  deallocate(hxyd_tmp)
  deallocate(hxzd_tmp)
  deallocate(hyyd_tmp)
  deallocate(hyzd_tmp)
  deallocate(hzzd_tmp)
!  deallocate(trk_tmp)
  deallocate(gx_tmp)
  deallocate(gy_tmp)
  deallocate(gz_tmp)
  deallocate(wvxd_tmp)
  deallocate(wvyd_tmp)
  deallocate(wvzd_tmp)
  deallocate(sigt_tmp)
!
end subroutine IO_input_initial_3D_WL_BHT
