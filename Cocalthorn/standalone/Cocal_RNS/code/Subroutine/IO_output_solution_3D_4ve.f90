subroutine IO_output_solution_3D_4ve
  use phys_constant, only : long
  use def_matter, only  : utf, uxf, uyf, uzf
  use grid_parameter, only  : nrf, ntf, npf
  implicit none
  integer :: ir, it, ip
!
! --- 4 velocity on fluid coordinate
  open(13,file='rns4ve_3D.las',status='unknown')
  write(13,'(5i5)') nrf, ntf, npf
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        write(13,'(1p,6e20.12)')  utf(ir,it,ip), &
     &                            uxf(ir,it,ip), &
     &                            uyf(ir,it,ip), &
     &                            uzf(ir,it,ip)
      end do
    end do
  end do
  close(13)
!
end subroutine IO_output_solution_3D_4ve
