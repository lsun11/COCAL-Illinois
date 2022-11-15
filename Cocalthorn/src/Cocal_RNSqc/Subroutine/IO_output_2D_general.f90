subroutine IO_output_2D_general(filename,fg,mg,pot,pl)
  use phys_constant, only : long
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  use trigonometry_grav_theta, only : sinthg, costhg, hsinthg, hcosthg
  use trigonometry_grav_phi,   only : sinphig, cosphig, hsinphig, hcosphig
  use coordinate_grav_r, only : rg, hrg
  use make_array_1d
  implicit none
  integer :: ir, it, ip, nrtmp, nttmp, nptmp, ii, ir0,it0,ip0
  real(long), pointer  :: pot(:,:,:)
  real(long), pointer  :: ra(:), sintha(:), costha(:), sinpha(:), cospha(:)
  real(long) :: xx, yy, zz
  character(len=1), intent(in) :: fg,mg
  character(len=2), intent(in) :: pl
  character(len=30), intent(in) :: filename
!
  call alloc_array1d(ra,0,nrg)
  call alloc_array1d(sintha,0,ntg)
  call alloc_array1d(costha,0,ntg)
  call alloc_array1d(sinpha,0,npg)
  call alloc_array1d(cospha,0,npg)

  if(fg == 'f') then
    nrtmp=nrf;  nttmp=ntf;  nptmp=npf
  else
    nrtmp=nrg;  nttmp=ntg;  nptmp=npg
  end if

  if (mg=='m') then
    ir0=1;  it0=1;  ip0=1
    ra(1:nrg)     = hrg(1:nrg)
    sintha(1:ntg) = hsinthg(1:ntg)
    costha(1:ntg) = hcosthg(1:ntg)
    sinpha(1:npg) = hsinphig(1:npg)
    cospha(1:npg) = hcosphig(1:npg)
  else
    ir0=0;  it0=0;  ip0=0
    ra(0:nrg)     = rg(0:nrg)
    sintha(0:ntg) = sinthg(0:ntg)
    costha(0:ntg) = costhg(0:ntg)
    sinpha(0:npg) = sinphig(0:npg)
    cospha(0:npg) = cosphig(0:npg)
  end if

  open(12,file=filename, status='unknown')
  if (pl=='xy' .or. pl=='yx')  then
    it = nttmp/2
    do ip = ip0, nptmp
      do ir = ir0, nrtmp
        xx = ra(ir)*sintha(it)*cospha(ip)
        yy = ra(ir)*sintha(it)*sinpha(ip)
        write(12,'(1p,3e23.15)') xx, yy, pot(ir,it,ip)
      end do
    end do
  end if
!
  if (pl=='yz' .or. pl=='zy')  then
    ip = nptmp/4
    do it = it0, nttmp
      do ir = ir0, nrtmp
        yy = ra(ir)*sintha(it)*sinpha(ip)
        zz = ra(ir)*costha(it)
        write(12,'(1p,3e23.15)') yy, zz, pot(ir,it,ip)
      end do
    end do
!
    ip = 3*nptmp/4
    do it = it0, nttmp
      do ir = ir0, nrtmp
        yy = ra(ir)*sintha(it)*sinpha(ip)
        zz = ra(ir)*costha(it)
        write(12,'(1p,3e23.15)') yy, zz, pot(ir,it,ip)
      end do
    end do
  end if
!
  if (pl=='zx' .or. pl=='xz')  then
    ip = 0
    do it = it0, nttmp
      do ir = ir0, nrtmp
        xx = ra(ir)*sintha(it)*cospha(ip)
        zz = ra(ir)*costha(it)
        write(12,'(1p,3e23.15)') xx, zz, pot(ir,it,ip)
      end do
    end do
!
    ip = nptmp/2
    do it = it0, nttmp
      do ir = ir0, nrtmp
        xx = ra(ir)*sintha(it)*cospha(ip)
        zz = ra(ir)*costha(it)
        write(12,'(1p,3e23.15)') xx, zz, pot(ir,it,ip)
      end do
    end do
  end if
  close(12)
!
  deallocate(ra)
  deallocate(sintha)
  deallocate(costha)
  deallocate(sinpha)
  deallocate(cospha)
!
end subroutine IO_output_2D_general
