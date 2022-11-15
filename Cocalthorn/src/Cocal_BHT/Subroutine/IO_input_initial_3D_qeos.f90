subroutine IO_input_initial_3D_qeos
  use phys_constant, only : long, nnrg, nntg, nnpg
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd
  use def_matter, only : rhof, rs, omef
  use def_matter_parameter, only : ome, ber, radi
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig 
  use grid_temporary
  use interface_interpo_3D_initial_4th
  use interface_interpo_3D_initial_4th_surface
  use make_array_2d
  use make_array_3d
  implicit none
  real(long), pointer :: rho_tmp(:,:,:),  omef_tmp(:,:,:), rs_tmp(:,:)
  real(long), pointer :: psi_tmp(:,:,:),  alph_tmp(:,:,:), &
  &                      bvxd_tmp(:,:,:), bvyd_tmp(:,:,:), bvzd_tmp(:,:,:)
  integer :: ir, it, ip, nrgmax, ntgmax, npgmax, nrfmax, ntfmax, npfmax
!
! --- Coordinate grids.
  open(14,file='rnsgrids_3D.ini',status='old')
  read(14,'(5i5)') nrgtmp, ntgtmp, npgtmp
  do ir = 0, nrgtmp
    read(14,'(1p,6e20.12)') rgtmp(ir)
  end do
  do it = 0, ntgtmp
    read(14,'(1p,6e20.12)') thgtmp(it)
  end do
  do ip = 0, npgtmp
    read(14,'(1p,6e20.12)') phigtmp(ip)
  end do
  close(14)
!
  nrgmax = max0(nrg,nrgtmp)
  ntgmax = max0(ntg,ntgtmp)
  npgmax = max0(npg,npgtmp)
  call alloc_array3d( psi_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(alph_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(bvxd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(bvyd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(bvzd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
!
! --- Matter
  open(12,file='rnsflu_3D.ini',status='old')
  read(12,'(5i5)') nrftmp, ntftmp, npftmp
!
  nrfmax = max0(nrf,nrftmp)
  ntfmax = max0(ntf,ntftmp)
  npfmax = max0(npf,npftmp)
  call alloc_array3d( rho_tmp,0,nrfmax,0,ntfmax,0,npfmax)
  call alloc_array2d(  rs_tmp,         0,ntfmax,0,npfmax)
  call alloc_array3d(omef_tmp,0,nrfmax,0,ntfmax,0,npfmax)
!
  do ip = 0, npftmp
    do it = 0, ntftmp
      do ir = 0, nrftmp
        read(12,'(1p,6e20.12)') rho_tmp(ir,it,ip), rs_tmp(it,ip), &
        &                      omef_tmp(ir,it,ip)
      end do
    end do
  end do
  read(12,'(1p,6e20.12)') ome, ber, radi
  close(12)
!
! --- Metric potentials.
  open(13,file='rnsgra_3D.ini',status='old')
  read(13,'(5i5)') nrgtmp, ntgtmp, npgtmp
  do ip = 0, npgtmp
    do it = 0, ntgtmp
      do ir = 0, nrgtmp
        read(13,'(1p,6e20.12)')  psi_tmp(ir,it,ip), &
        &                       alph_tmp(ir,it,ip), &
        &                       bvxd_tmp(ir,it,ip), &
        &                       bvyd_tmp(ir,it,ip), &
        &                       bvzd_tmp(ir,it,ip)
      end do
    end do
  end do
  read(13,'(1p,6e20.12)') ome, ber, radi
  close(13)
!
! --- Interpolation
!
  if (nrf.ne.nrftmp.or.ntf.ne.ntftmp.or.npf.ne.npftmp) then
    call interpo_3D_initial_4th( rho_tmp,'sfco')
    call interpo_3D_initial_4th(omef_tmp,'sfco')
  end if
  if (ntf.ne.ntftmp.or.npf.ne.npftmp) then
    call interpo_3D_initial_4th_surface(rs_tmp)
  end if
  if (nrg.ne.nrgtmp.or.ntg.ne.ntgtmp.or.npg.ne.npgtmp) then
    call interpo_3D_initial_4th( psi_tmp,'grco')
    call interpo_3D_initial_4th(alph_tmp,'grco')
    call interpo_3D_initial_4th(bvxd_tmp,'grco')
    call interpo_3D_initial_4th(bvyd_tmp,'grco')
    call interpo_3D_initial_4th(bvzd_tmp,'grco')
  end if
!
  rhof(0:nrf,0:ntf,0:npf) =  rho_tmp(0:nrf,0:ntf,0:npf)
    rs(0:ntf,0:npf)       =   rs_tmp(0:ntf,0:npf)
  omef(0:nrf,0:ntf,0:npf) = omef_tmp(0:nrf,0:ntf,0:npf)
   psi(0:nrg,0:ntg,0:npg) =  psi_tmp(0:nrg,0:ntg,0:npg)
  alph(0:nrg,0:ntg,0:npg) = alph_tmp(0:nrg,0:ntg,0:npg)
  bvxd(0:nrg,0:ntg,0:npg) = bvxd_tmp(0:nrg,0:ntg,0:npg)
  bvyd(0:nrg,0:ntg,0:npg) = bvyd_tmp(0:nrg,0:ntg,0:npg)
  bvzd(0:nrg,0:ntg,0:npg) = bvzd_tmp(0:nrg,0:ntg,0:npg)
!
  deallocate( rho_tmp)
  deallocate(  rs_tmp)
  deallocate(omef_tmp)
  deallocate( psi_tmp)
  deallocate(alph_tmp)
  deallocate(bvxd_tmp)
  deallocate(bvyd_tmp)
  deallocate(bvzd_tmp)
!
end subroutine IO_input_initial_3D_qeos
