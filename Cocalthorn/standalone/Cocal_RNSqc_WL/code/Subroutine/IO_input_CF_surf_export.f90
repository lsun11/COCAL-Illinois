subroutine IO_input_CF_surf_export(filenm,rs)
  use phys_constant, only : long, nnrg, nntg, nnpg
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp
  real(8), pointer ::  rs(:,:)
  character(len=*) :: filenm
!
! --- Star surface
  open(15,file=trim(filenm),status='old')
  read(15,'(5i5)')   nttmp,  nptmp
  do ip = 0, nptmp
    do it = 0, nttmp
      read(15,'(1p,6e20.12)')   rs(it,ip)
    end do
  end do
  close(15)
!
end subroutine IO_input_CF_surf_export
