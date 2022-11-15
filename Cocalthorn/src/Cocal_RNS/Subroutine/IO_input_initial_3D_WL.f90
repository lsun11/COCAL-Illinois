subroutine IO_input_initial_3D_WL
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          gaugex, gaugey, gaugez
  use grid_temporary, only : nrgtmp, ntgtmp, npgtmp, rgtmp, thgtmp, phigtmp
  use interface_interpo_3D_initial_4th
  use make_array_3d
  implicit none
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
! --- Metric potentials.
  open(13,file='rnsgra_hij_3D.ini',status='unknown')
  read(13,'(5i5)') nrgtmp, ntgtmp, npgtmp
  do ipg = 0, npgtmp
    do itg = 0, ntgtmp
      do irg = 0, nrgtmp
        read(13,'(1p,6e20.12)') hxxd_tmp(irg,itg,ipg), &
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
  if (nrg.ne.nrgtmp.or.ntg.ne.ntgtmp.or.npg.ne.npgtmp) then
    call interpo_3D_initial_4th(hxxd_tmp,'grco')
    call interpo_3D_initial_4th(hxyd_tmp,'grco')
    call interpo_3D_initial_4th(hxzd_tmp,'grco')
    call interpo_3D_initial_4th(hyyd_tmp,'grco')
    call interpo_3D_initial_4th(hyzd_tmp,'grco')
    call interpo_3D_initial_4th(hzzd_tmp,'grco')
  end if
!
  hxxd(0:nrg,0:ntg,0:npg) = hxxd_tmp(0:nrg,0:ntg,0:npg)
  hxyd(0:nrg,0:ntg,0:npg) = hxyd_tmp(0:nrg,0:ntg,0:npg)
  hxzd(0:nrg,0:ntg,0:npg) = hxzd_tmp(0:nrg,0:ntg,0:npg)
  hyyd(0:nrg,0:ntg,0:npg) = hyyd_tmp(0:nrg,0:ntg,0:npg)
  hyzd(0:nrg,0:ntg,0:npg) = hyzd_tmp(0:nrg,0:ntg,0:npg)
  hzzd(0:nrg,0:ntg,0:npg) = hzzd_tmp(0:nrg,0:ntg,0:npg)
!
  gaugex(0:nrg,0:ntg,0:npg) = 0.0d0
  gaugey(0:nrg,0:ntg,0:npg) = 0.0d0
  gaugez(0:nrg,0:ntg,0:npg) = 0.0d0
!
  deallocate(hxxd_tmp)
  deallocate(hxyd_tmp)
  deallocate(hxzd_tmp)
  deallocate(hyyd_tmp)
  deallocate(hyzd_tmp)
  deallocate(hzzd_tmp)
!
end subroutine IO_input_initial_3D_WL
