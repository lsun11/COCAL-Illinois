subroutine IO_output_3D_general(filename,fg,mg,pot)
  use phys_constant, only : long
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp, ir0,it0,ip0
  real(long), pointer  :: pot(:,:,:)
  character(len=1), intent(in) :: fg, mg
  character(len=30), intent(in) :: filename
!
  if(fg == 'f') then
    nrtmp=nrf;  nttmp=ntf;  nptmp=npf
  else
    nrtmp=nrg;  nttmp=ntg;  nptmp=npg
  end if

  if (mg=='m') then
    ir0=1;  it0=1;  ip0=1
  else
    ir0=0;  it0=0;  ip0=0
  end if

  open(12,file=filename, status='unknown')
  write(12,'(5i5)')  nrtmp, nttmp, nptmp
  do ip = ip0, nptmp
    do it = it0, nttmp
      do ir = ir0, nrtmp
        write(12,'(1p,e20.12)') pot(ir,it,ip)
      end do
    end do
  end do
  close(12)
!
end subroutine IO_output_3D_general
