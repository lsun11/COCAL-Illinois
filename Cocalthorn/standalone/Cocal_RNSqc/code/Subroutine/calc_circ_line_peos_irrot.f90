subroutine calc_circ_line_peos_irrot(impt)
  use phys_constant, only  : long, pi
  use grid_parameter, only : ntfeq, npfxzp, npfyzp, nrf, ntf, npf
  use coordinate_grav_r, only : hrg, rg
  use def_matter_parameter, only : radi
  use make_array_1d
  use make_array_2d
  use coordinate_grav_phi, only : phig
  use def_quantities, only : circ_line_xy, circ_line_yz, circ_line_zx, circ_shift_xy
  use interface_source_circ_line_peos_irrot
  use interface_source_circ_line_peos_shift
  use interface_line_int_fluid
  implicit none
  real(long) :: line_int, circ_line_ip
  real(long),pointer :: souf_circ_line_ph_peos(:,:), souf_circ_line_th_peos(:,:)
  integer ::  irg, irs, ii, ia, ib, impt, iphi
  character(len=1) :: np(3) = (/'1', '2', '3'/)
!
  write(6,'(a90)') "---------------------------------------------------------------------------------------------"
  do ia = nrf,nrf 
    ib = ntfeq    ! equatorial plane
    call alloc_array2d(souf_circ_line_ph_peos, 1, npf, 1, 3)
    call source_circ_line_peos_irrot("ph", ia, ib, souf_circ_line_ph_peos)
    call line_int_fluid("ph", ia, ib, souf_circ_line_ph_peos, line_int)
    circ_line_xy = line_int*radi    ! along positive z-axis
!
    souf_circ_line_ph_peos = 0.0d0
    call source_circ_line_peos_shift("ph", ia, ib, souf_circ_line_ph_peos)
    call line_int_fluid("ph", ia, ib, souf_circ_line_ph_peos, line_int)
    circ_shift_xy = line_int*radi    ! along positive z-axis
    deallocate(souf_circ_line_ph_peos)
    write(6,'(a8,i3,a42,1p,e23.15)') 'For irf=', ia, '  Circulation (line integral) on xy plane=', circ_line_xy
    write(6,'(a8,i3,a42,1p,e23.15)') 'For irf=', ia, '  Shift circulation on xy plane          =', circ_shift_xy
  end do 

!  write(6,'(a90)') "---------------------------------------------------------------------------------------------"
  do ia = nrf,nrf 
    ib = npfxzp   ! phi = 0.0
    call alloc_array2d(souf_circ_line_th_peos, 1, 2*ntf, 1, 3)
    call source_circ_line_peos_irrot("th", ia, ib, souf_circ_line_th_peos)
    call line_int_fluid("th", ia, ib, souf_circ_line_th_peos, line_int)
    circ_line_zx = line_int*radi    ! along positive y-axis
    deallocate(souf_circ_line_th_peos)
    write(6,'(a8,i3,a42,1p,e23.15)') 'For irf=', ia, '  Circulation (line integral) on zx plane=', circ_line_zx
  end do 

!  write(6,'(a90)') "---------------------------------------------------------------------------------------------"
  do ia = nrf,nrf 
    ib = npfyzp   ! phi = pi/2
    call alloc_array2d(souf_circ_line_th_peos, 1, 2*ntf, 1, 3)
    call source_circ_line_peos_irrot("th", ia, ib, souf_circ_line_th_peos)
    call line_int_fluid("th", ia, ib, souf_circ_line_th_peos, line_int)
    circ_line_yz = -line_int*radi   ! along positive x-axis
    deallocate(souf_circ_line_th_peos)
    write(6,'(a8,i3,a42,1p,e23.15)') 'For irf=', ia, '  Circulation (line integral) on yz plane=', circ_line_yz
  end do 

! Output all meridional circulations
  open(12,file='circ_line_mpt'//np(impt)//'.dat', status='unknown')
  call alloc_array2d(souf_circ_line_th_peos, 1, 2*ntf, 1, 3)
  do iphi=0,npf/2
    do ia = nrf,nrf
      souf_circ_line_th_peos = 0.0d0
      ib = iphi    !npfyzp   ! phi = pi/2
      call source_circ_line_peos_irrot("th", ia, ib, souf_circ_line_th_peos)
      call line_int_fluid("th", ia, ib, souf_circ_line_th_peos, line_int)
      circ_line_ip = -line_int*radi   ! along positive x-axis
!      write(6,'(a8,i3,a42,1p,e23.15)') 'For irf=', ia, '  Circulation (line
!      integral) on yz plane=', circ_line_ip
      write(12,'(1p,6e20.12)')  phig(iphi), circ_line_ip
    end do
  end do
  deallocate(souf_circ_line_th_peos)
  close(12)

end subroutine calc_circ_line_peos_irrot
