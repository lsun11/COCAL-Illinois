subroutine IO_input_initial_3D_EMF
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_emfield, only  : va, vaxd, vayd, vazd
  use def_faraday_tensor, only  : fxd_grid, fyd_grid, fzd_grid, fijd_grid
  use grid_temporary, only : nrgtmp, ntgtmp, npgtmp, rgtmp, thgtmp, phigtmp
  use interface_interpo_3D_initial_4th
  use make_array_3d
  implicit none
  real(long), pointer :: va_tmp(:,:,:)
  real(long), pointer :: vaxd_tmp(:,:,:), vayd_tmp(:,:,:), vazd_tmp(:,:,:)
  real(long), pointer ::  fxd_tmp(:,:,:),  fyd_tmp(:,:,:),  fzd_tmp(:,:,:)
  real(long), pointer :: fij1_tmp(:,:,:), fij2_tmp(:,:,:), fij3_tmp(:,:,:)
  integer :: irg, itg, ipg, nrgmax, ntgmax, npgmax
!
  nrgmax = max0(nrg,nrgtmp)
  ntgmax = max0(ntg,ntgtmp)
  npgmax = max0(npg,npgtmp)
  call alloc_array3d(  va_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(vaxd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(vayd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(vazd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d( fxd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d( fyd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d( fzd_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(fij1_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(fij2_tmp,0,nrgmax,0,ntgmax,0,npgmax)
  call alloc_array3d(fij3_tmp,0,nrgmax,0,ntgmax,0,npgmax)
!                                                                              
! --- EM 1-form and Faraday tensor
  open(13,file='rnsEMF_3D.ini',status='unknown')
  read(13,'(5i5)') nrgtmp, ntgtmp, npgtmp
  do ipg = 0, npgtmp
    do itg = 0, ntgtmp
      do irg = 0, nrgtmp
        read(13,'(1p,6e20.12)')   va_tmp(irg,itg,ipg), &
        &                       vaxd_tmp(irg,itg,ipg), &
        &                       vayd_tmp(irg,itg,ipg), &
        &                       vazd_tmp(irg,itg,ipg)
      end do
    end do
  end do
  close(13)
!
  open(13,file='rnsEMF_faraday_3D.ini',status='unknown')
  read(13,'(5i5)') nrgtmp, ntgtmp, npgtmp
  do ipg = 0, npgtmp
    do itg = 0, ntgtmp
      do irg = 0, nrgtmp
        read(13,'(1p,6e20.12)')  fxd_tmp(irg,itg,ipg), &
        &                        fyd_tmp(irg,itg,ipg), &
        &                        fzd_tmp(irg,itg,ipg), &
        &                       fij1_tmp(irg,itg,ipg), &
        &                       fij2_tmp(irg,itg,ipg), &
        &                       fij3_tmp(irg,itg,ipg)
      end do
    end do
  end do
  close(13)
!
  if (nrg.ne.nrgtmp.or.ntg.ne.ntgtmp.or.npg.ne.npgtmp) then
    call interpo_3D_initial_4th(  va_tmp,'grco')
    call interpo_3D_initial_4th(vaxd_tmp,'grco')
    call interpo_3D_initial_4th(vayd_tmp,'grco')
    call interpo_3D_initial_4th(vazd_tmp,'grco')
    call interpo_3D_initial_4th( fxd_tmp,'grco')
    call interpo_3D_initial_4th( fyd_tmp,'grco')
    call interpo_3D_initial_4th( fzd_tmp,'grco')
    call interpo_3D_initial_4th(fij1_tmp,'grco')
    call interpo_3D_initial_4th(fij2_tmp,'grco')
    call interpo_3D_initial_4th(fij3_tmp,'grco')
  end if
!
    va(0:nrg,0:ntg,0:npg)        =   va_tmp(0:nrg,0:ntg,0:npg)
  vaxd(0:nrg,0:ntg,0:npg)        = vaxd_tmp(0:nrg,0:ntg,0:npg)
  vayd(0:nrg,0:ntg,0:npg)        = vayd_tmp(0:nrg,0:ntg,0:npg)
  vazd(0:nrg,0:ntg,0:npg)        = vazd_tmp(0:nrg,0:ntg,0:npg)
   fxd_grid(0:nrg,0:ntg,0:npg)   =  fxd_tmp(0:nrg,0:ntg,0:npg)
   fyd_grid(0:nrg,0:ntg,0:npg)   =  fyd_tmp(0:nrg,0:ntg,0:npg)
   fzd_grid(0:nrg,0:ntg,0:npg)   =  fzd_tmp(0:nrg,0:ntg,0:npg)
  fijd_grid(0:nrg,0:ntg,0:npg,1) = fij1_tmp(0:nrg,0:ntg,0:npg)
  fijd_grid(0:nrg,0:ntg,0:npg,2) = fij2_tmp(0:nrg,0:ntg,0:npg)
  fijd_grid(0:nrg,0:ntg,0:npg,3) = fij3_tmp(0:nrg,0:ntg,0:npg)
!
  deallocate(va_tmp)
  deallocate(vaxd_tmp)
  deallocate(vayd_tmp)
  deallocate(vazd_tmp)
  deallocate(fxd_tmp)
  deallocate(fyd_tmp)
  deallocate(fzd_tmp)
  deallocate(fij1_tmp)
  deallocate(fij2_tmp)
  deallocate(fij3_tmp)
!
end subroutine IO_input_initial_3D_EMF
