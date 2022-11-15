subroutine calc_circ_surf_peos_irrot(impt)
  use phys_constant, only  : long, pi
  use grid_parameter, only : ntfeq, npfxzp, npfyzp, nrf, ntf, npf
  use coordinate_grav_r, only : hrg, rg
  use def_matter_parameter, only : radi
  use make_array_1d
  use make_array_2d
  use coordinate_grav_phi, only : phig
  use def_quantities, only : circ_surf_xy, circ_surf_yz, circ_surf_zx
  use interface_source_circ_surf_peos_irrot
  use interface_plane_surf_int_fluid
  implicit none
  real(long) :: surf_int, circ_surf_ip
  real(long),pointer :: souf_circ_surf_ph_peos(:,:), souf_circ_surf_th_peos(:,:)
  integer ::  irg, irs, ii, ia, ib, impt, iphi
  character(len=1) :: np(3) = (/'1', '2', '3'/)
!
  write(6,'(a90)') "---------------------------------------------------------------------------------------------"
  do ia = nrf,nrf  
    ib = ntfeq    ! equatorial plane
    call alloc_array2d(souf_circ_surf_ph_peos, 1, ia, 1, npf)
    call source_circ_surf_peos_irrot("ph", ia, ib, souf_circ_surf_ph_peos)
    call plane_surf_int_fluid("ph", ia, ib, souf_circ_surf_ph_peos, surf_int)
    circ_surf_xy = surf_int*radi
    deallocate(souf_circ_surf_ph_peos)
    write(6,'(a8,i3,a42,1p,e23.15)') 'For irf=', ia, '  Circulation (surf integral) on xy plane=',circ_surf_xy
  end do

!  write(6,'(a90)') "---------------------------------------------------------------------------------------------"
  do ia = nrf,nrf  
    ib = npfxzp   ! phi = 0.0
    call alloc_array2d(souf_circ_surf_th_peos, 1, ia, 1, 2*ntf)
    call source_circ_surf_peos_irrot("th", ia, ib, souf_circ_surf_th_peos)
    call plane_surf_int_fluid("th", ia, ib, souf_circ_surf_th_peos, surf_int)
    circ_surf_zx = surf_int*radi
    deallocate(souf_circ_surf_th_peos)
    write(6,'(a8,i3,a42,1p,e23.15)') 'For irf=', ia, '  Circulation (surf integral) on zx plane=',circ_surf_zx
  end do

!  write(6,'(a90)') "---------------------------------------------------------------------------------------------"
  do ia = nrf,nrf  
    ib = npfyzp   ! phi = pi/2
    call alloc_array2d(souf_circ_surf_th_peos, 1, ia, 1, 2*ntf)
    call source_circ_surf_peos_irrot("th", ia, ib, souf_circ_surf_th_peos)
    call plane_surf_int_fluid("th", ia, ib, souf_circ_surf_th_peos, surf_int)
    circ_surf_yz = -surf_int*radi
    deallocate(souf_circ_surf_th_peos)
    write(6,'(a8,i3,a42,1p,e23.15)') 'For irf=', ia, '  Circulation (surf integral) on yz plane=',circ_surf_yz
  end do

! Output all meridional circulations
  open(12,file='circ_surf_mpt'//np(impt)//'.dat', status='unknown')
  do iphi=0,npf/2
    do ia = nrf,nrf
      ib = iphi    !npfyzp   ! phi = pi/2
      call alloc_array2d(souf_circ_surf_th_peos, 1, ia, 1, 2*ntf)
      call source_circ_surf_peos_irrot("th", ia, ib, souf_circ_surf_th_peos)
      call plane_surf_int_fluid("th", ia, ib, souf_circ_surf_th_peos, surf_int)
      circ_surf_ip = -surf_int*radi
      write(12,'(1p,6e20.12)')  phig(iphi), circ_surf_ip
      deallocate(souf_circ_surf_th_peos)
    end do
  end do
  close(12)

end subroutine calc_circ_surf_peos_irrot
