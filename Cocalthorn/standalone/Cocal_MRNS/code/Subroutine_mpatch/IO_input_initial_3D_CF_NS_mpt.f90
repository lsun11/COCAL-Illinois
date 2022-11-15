subroutine IO_input_initial_3D_CF_NS_mpt(impt)
  use phys_constant, only : long, nnrg, nntg, nnpg
  use def_metric
  use def_matter, only : emd, rs
  use def_matter_parameter, only : ome, ber, radi
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig 
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp
  real(long) :: rgtmp(0:nnrg), thgtmp(0:nntg), phigtmp(0:nnpg)
  integer,intent(in)  :: impt
  character(len=1) :: np(3) = (/'1', '2', '3'/)

!
! --- Matter
  if (impt==1 .or. impt==2)  then
    open(12,file='bnsflu_3D_mpt'//np(impt)//'.ini',status='old')
    read(12,'(5i5)') nrtmp, nttmp, nptmp
    do ip = 0, nptmp
      do it = 0, nttmp
        do ir = 0, nrtmp
          read(12,'(1p,6e20.12)') emd(ir,it,ip)
        end do
      end do
    end do
    read(12,'(1p,6e20.12)') ome, ber, radi
    close(12)
!
    open(15,file='bnssur_3D_mpt'//np(impt)//'.ini',status='old')
    read(15,'(5i5)') nttmp, nptmp
    do ip = 0, nptmp
      do it = 0, nttmp
        read(15,'(1p,6e20.12)')   rs(it,ip)
      end do
    end do
    close(15)
  end if
!
! --- Metric potentials.
  open(13,file='bnsgra_3D_mpt'//np(impt)//'.ini',status='old')
  read(13,'(5i5)') nrtmp, nttmp, nptmp
  do ip = 0, nptmp
    do it = 0, nttmp
      do ir = 0, nrtmp
        read(13,'(1p,6e20.12)')  psi(ir,it,ip), &
    &                           alph(ir,it,ip), &
    &                           bvxd(ir,it,ip), &
    &                           bvyd(ir,it,ip), &
    &                           bvzd(ir,it,ip)
      end do
    end do
  end do
  close(13)
  alps(0:nrtmp,0:nttmp,0:nptmp) = alph(0:nrtmp,0:nttmp,0:nptmp)*psi(0:nrtmp,0:nttmp,0:nptmp)
!
! --- Coordinate grids.
  open(14,file='bnsgrids_3D_mpt'//np(impt)//'.ini',status='old')
  read(14,'(5i5)') nrtmp, nttmp, nptmp
  do ir = 0, nrtmp
    read(14,'(1p,6e20.12)') rgtmp(ir)
  end do
  do it = 0, nttmp
    read(14,'(1p,6e20.12)') thgtmp(it)
  end do
  do ip = 0, nptmp
    read(14,'(1p,6e20.12)') phigtmp(ip)
  end do
  close(14)

end subroutine IO_input_initial_3D_CF_NS_mpt
