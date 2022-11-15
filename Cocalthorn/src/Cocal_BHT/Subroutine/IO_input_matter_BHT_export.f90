subroutine IO_input_matter_BHT_export(filenm, emdg, omeg, ome, ber, radi)
  use phys_constant, only : long
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp
  real(8), pointer :: emdg(:,:,:), omeg(:,:,:)
  real(8) :: ome, ber, radi, emdc
  character(len=*) :: filenm
!
  write(6,*) "Reading emdg, omeg..."   
! --- Matter
  open(12,file=trim(filenm),status='old')
  read(12,'(5i5)')  nrtmp, nttmp, nptmp
  do ip = 0, nptmp
    do it = 0, nttmp
      do ir = 0, nrtmp
        read(12,'(1p,6e23.15)') emdg(ir,it,ip),  omeg(ir,it,ip)
      end do
    end do
  end do
  read(12,'(1p,6e23.15)') ome, ber, radi
!  read(12,'(1p,6e23.15)') emdc
  close(12)
!
end subroutine IO_input_matter_BHT_export
