subroutine IO_input_initial_3D_4ve
  use phys_constant,  only : long
  use grid_parameter, only : nrf, ntf, npf
  use def_matter, only  : utf, uxf, uyf, uzf
  use grid_temporary, only : nrftmp, ntftmp, npftmp, rgtmp, thgtmp, phigtmp
  use interface_interpo_3D_initial_4th
  use make_array_3d
  implicit none
  real(long), pointer  :: utf_tmp(:,:,:), &
  &                       uxf_tmp(:,:,:), uyf_tmp(:,:,:), uzf_tmp(:,:,:)
  integer :: ir, it, ip, nrfmax, ntfmax, npfmax
!
  nrfmax = max0(nrf,nrftmp)
  ntfmax = max0(ntf,ntftmp)
  npfmax = max0(npf,npftmp)
  call alloc_array3d(utf_tmp,0,nrfmax,0,ntfmax,0,npfmax)
  call alloc_array3d(uxf_tmp,0,nrfmax,0,ntfmax,0,npfmax)
  call alloc_array3d(uyf_tmp,0,nrfmax,0,ntfmax,0,npfmax)
  call alloc_array3d(uzf_tmp,0,nrfmax,0,ntfmax,0,npfmax)
!
! --- 4 velocity on fluid coordinate
  open(13,file='rns4ve_3D.ini',status='old')
  read(13,'(5i5)') nrftmp, ntftmp, npftmp
  do ip = 0, npftmp
    do it = 0, ntftmp
      do ir = 0, nrftmp
        read(13,'(1p,6e20.12)')  utf_tmp(ir,it,ip), &
        &                        uxf_tmp(ir,it,ip), &
        &                        uyf_tmp(ir,it,ip), &
        &                        uzf_tmp(ir,it,ip)
      end do
    end do
  end do
  close(13)
!
  if (nrf.ne.nrftmp.or.ntf.ne.ntftmp.or.npf.ne.npftmp) then
    call interpo_3D_initial_4th(utf_tmp,'sfco')
    call interpo_3D_initial_4th(uxf_tmp,'sfco')
    call interpo_3D_initial_4th(uyf_tmp,'sfco')
    call interpo_3D_initial_4th(uzf_tmp,'sfco')
  end if
!
  utf(0:nrf,0:ntf,0:npf) = utf_tmp(0:nrf,0:ntf,0:npf)
  uxf(0:nrf,0:ntf,0:npf) = uxf_tmp(0:nrf,0:ntf,0:npf)
  uyf(0:nrf,0:ntf,0:npf) = uyf_tmp(0:nrf,0:ntf,0:npf)
  uzf(0:nrf,0:ntf,0:npf) = uzf_tmp(0:nrf,0:ntf,0:npf)
!
  deallocate(utf_tmp)
  deallocate(uxf_tmp)
  deallocate(uyf_tmp)
  deallocate(uzf_tmp)
!
end subroutine IO_input_initial_3D_4ve
