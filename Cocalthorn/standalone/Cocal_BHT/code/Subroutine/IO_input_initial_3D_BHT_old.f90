subroutine IO_input_initial_3D_BHT
  use phys_constant, only : long, nnrg, nntg, nnpg
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : alph, psi, bvxd, bvyd, bvzd, alps
  use def_matter, only : emdg, omeg
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
  real(long), pointer :: emdg_tmp(:,:,:),  omeg_tmp(:,:,:)
  real(long), pointer :: psi_tmp(:,:,:),  alph_tmp(:,:,:), &
  &                      bvxd_tmp(:,:,:), bvyd_tmp(:,:,:), bvzd_tmp(:,:,:)
  integer :: ir, it, ip, nrgmax, ntgmax, npgmax
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
  read(12,'(5i5)') nrgtmp, ntgtmp, npgtmp
!
  call alloc_array3d(emdg_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(omeg_tmp,0,nrgmax,0,ntgmax,0,npgmax)
!
  do ip = 0, npgtmp
    do it = 0, ntgtmp
      do ir = 0, nrgtmp
        read(12,'(1p,6e20.12)') emdg_tmp(ir,it,ip), omeg_tmp(ir,it,ip)
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
  if (nrg.ne.nrgtmp.or.ntg.ne.ntgtmp.or.npg.ne.npgtmp) then
    call interpo_3D_initial_4th( psi_tmp,'grco')
    call interpo_3D_initial_4th(alph_tmp,'grco')
    call interpo_3D_initial_4th(bvxd_tmp,'grco')
    call interpo_3D_initial_4th(bvyd_tmp,'grco')
    call interpo_3D_initial_4th(bvzd_tmp,'grco')
    call interpo_3D_initial_4th(emdg_tmp,'grco')
    call interpo_3D_initial_4th(omeg_tmp,'grco')
  end if
!
  emdg(0:nrg,0:ntg,0:npg) = emdg_tmp(0:nrg,0:ntg,0:npg)
  omeg(0:nrg,0:ntg,0:npg) = omeg_tmp(0:nrg,0:ntg,0:npg)
   psi(0:nrg,0:ntg,0:npg) =  psi_tmp(0:nrg,0:ntg,0:npg)
  alph(0:nrg,0:ntg,0:npg) = alph_tmp(0:nrg,0:ntg,0:npg)
  bvxd(0:nrg,0:ntg,0:npg) = bvxd_tmp(0:nrg,0:ntg,0:npg)
  bvyd(0:nrg,0:ntg,0:npg) = bvyd_tmp(0:nrg,0:ntg,0:npg)
  bvzd(0:nrg,0:ntg,0:npg) = bvzd_tmp(0:nrg,0:ntg,0:npg)
!
  alps(0:nrg,0:ntg,0:npg) = alph(0:nrg,0:ntg,0:npg)*psi(0:nrg,0:ntg,0:npg)

  deallocate(emdg_tmp)
  deallocate(omeg_tmp)
  deallocate( psi_tmp)
  deallocate(alph_tmp)
  deallocate(bvxd_tmp)
  deallocate(bvyd_tmp)
  deallocate(bvzd_tmp)
!
end subroutine IO_input_initial_3D_BHT
