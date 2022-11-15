subroutine IO_input_grav_export_Kij(filenm,kxx,kxy,kxz,kyy,kyz,kzz)
  use phys_constant, only : long
  implicit none
  integer :: irg, itg, ipg, nrtmp, nttmp, nptmp
  real(8), pointer :: kxx(:,:,:), kxy(:,:,:), kxz(:,:,:), &
      &               kyy(:,:,:), kyz(:,:,:), kzz(:,:,:)
  character(len=*) :: filenm
!
  write(6,*) "Reading Kij..."
! --- Metric potentials.
  open(13,file=trim(filenm),status='old')
  read(13,'(5i5)')  nrtmp, nttmp, nptmp
  do ipg = 0, nptmp
    do itg = 0, nttmp
      do irg = 0, nrtmp
        read(13,'(1p,6e23.15)')  kxx(irg,itg,ipg), &
        &                        kxy(irg,itg,ipg), &
        &                        kxz(irg,itg,ipg), &
        &                        kyy(irg,itg,ipg), &
        &                        kyz(irg,itg,ipg), &
        &                        kzz(irg,itg,ipg)
      end do
    end do
  end do
  close(13)
!
end subroutine IO_input_grav_export_Kij
