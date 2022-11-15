subroutine IO_output_1D_general(filename,fg,mg,pot,jr,jt,jp)
  use phys_constant, only : long
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use coordinate_grav_r, only : rg, hrg
  use coordinate_grav_theta, only : thg, hthg
  use coordinate_grav_phi, only : phig, hphig
  use make_array_1d
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp, ii, ir0,it0,ip0
  real(long), pointer  :: pot(:,:,:)
  real(long), pointer  :: ra(:), tha(:), pha(:)
  integer, intent(in) :: jr, jt, jp
  character(len=1), intent(in) :: fg,mg
  character(len=30), intent(in) :: filename
!
  call alloc_array1d(ra,0,nrg)
  call alloc_array1d(tha,0,ntg)
  call alloc_array1d(pha,0,npg)


  if(fg == 'f') then
    nrtmp=nrf;  nttmp=ntf;  nptmp=npf
  else
    nrtmp=nrg;  nttmp=ntg;  nptmp=npg
  end if

  if (mg=='m') then
    ir0=1;  it0=1;  ip0=1
    ra(1:nrg)  = hrg(1:nrg)
    tha(1:ntg) = hthg(1:ntg)
    pha(1:npg) = hphig(1:npg)
  else
    ir0=0;  it0=0;  ip0=0
    ra(0:nrg)  = rg(0:nrg)
    tha(0:ntg) = thg(0:ntg)
    pha(0:npg) = phig(0:npg)
  end if

  open(12,file=filename, status='unknown')
  if (jr==-1)  then
    do ii=ir0,nrtmp
      write(12,'(1p,2e20.12)') ra(ii), pot(ii,jt,jp)
    end do
  end if
  if (jt==-1)  then
    do ii=it0,nttmp
      write(12,'(1p,2e20.12)') tha(ii), pot(jr,ii,jp)
    end do
  end if
  if (jp==-1)  then
    do ii=ip0,nptmp
      write(12,'(1p,2e20.12)') pha(ii), pot(jr,jt,ii)
    end do
  end if
  close(12)
!
  deallocate(ra)
  deallocate(tha)
  deallocate(pha)
!
end subroutine IO_output_1D_general
