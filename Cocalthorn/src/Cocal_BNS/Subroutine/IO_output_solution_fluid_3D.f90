subroutine IO_output_solution_fluid_3D
  use def_matter, only    : emd, rs, rhof, utf, uxf, uyf, uzf
  use def_matter_parameter, only : ome, ber, radi
  use grid_parameter, only : nrf, ntf, npf
  implicit none
  integer :: ir, it, ip
!
! --- Matter
  open(12,file='rnsflu_all_3D.las',status='unknown')
  write(12,'(5i5)') nrf, ntf, npf
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        write(12,'(1p,7e20.12)') emd(ir,it,ip), rs(it,ip), rhof(ir,it,ip), &
        &      utf(ir,it,ip), uxf(ir,it,ip), uyf(ir,it,ip), uzf(ir,it,ip)
      end do
    end do
  end do
  write(12,'(1p,6e20.12)') ome, ber, radi
  close(12)
!
end subroutine IO_output_solution_fluid_3D
