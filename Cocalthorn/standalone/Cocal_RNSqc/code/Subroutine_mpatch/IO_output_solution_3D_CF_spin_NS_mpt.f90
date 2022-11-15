subroutine IO_output_solution_3D_CF_spin_NS_mpt(impt,csol)
  use phys_constant, only : long
  use grid_parameter, only : r_surf
  use def_metric, only  : alph, psi, bvxd, bvyd, bvzd
  use def_matter, only  : emd, rs, vep
  use def_matter_parameter
  use def_velocity_rot
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig 
  use grid_parameter, only  :   nrg, ntg, npg, nrf, ntf, npf
  implicit none
  integer :: ir, it, ip
  integer,intent(in)  :: impt
  character(len=4), intent(in) :: csol
  character(len=1) :: np(3) = (/'1', '2', '3'/)
!
! --- Matter
  if (impt==1 .or. impt==2) then
    open(11,file='r_surf_mpt'//np(impt)//'.dat',status='unknown')
    write(11,'(1p,1e23.15)') r_surf
    close(11)

    open(12,file='bnsflu_3D_mpt'//np(impt)//csol ,status='unknown')
    write(12,'(5i5)') nrf, ntf, npf
    do ip = 0, npf
      do it = 0, ntf
        do ir = 0, nrf
          write(12,'(1p,6e20.12)') emd(ir,it,ip), vep(ir,it,ip), wxspf(ir,it,ip), &
                &   wyspf(ir,it,ip), wzspf(ir,it,ip)
        end do
      end do
    end do
    write(12,'(1p,6e20.12)') ome, ber, radi
    write(12,'(1p,6e20.12)') confpow, omespx, omespy, omespz
    close(12)
!
    open(15,file='bnssur_3D_mpt'//np(impt)//csol ,status='unknown')
    write(15,'(5i5)') ntf, npf
    do ip = 0, npf
      do it = 0, ntf
        write(15,'(1p,6e20.12)')   rs(it,ip)
      end do
    end do
    close(15)
  end if
!
! --- Metric potentials.
  open(13,file='bnsgra_3D_mpt'//np(impt)//csol ,status='unknown')
  write(13,'(5i5)') nrg, ntg, npg
  do ip = 0, npg
    do it = 0, ntg
      do ir = 0, nrg
        write(13,'(1p,6e20.12)')  psi(ir,it,ip), &
    &                            alph(ir,it,ip), &
    &                            bvxd(ir,it,ip), &
    &                            bvyd(ir,it,ip), &
    &                            bvzd(ir,it,ip)
      end do
    end do
  end do
  close(13)
!
! --- Coordinate grids.
  open(14,file='bnsgrids_3D_mpt'//np(impt)//csol ,status='unknown')
  write(14,'(5i5)') nrg, ntg, npg
  do ir = 0, nrg
    write(14,'(1p,6e20.12)') rg(ir)
  end do
  do it = 0, ntg
    write(14,'(1p,6e20.12)') thg(it)
  end do
  do ip = 0, npg
    write(14,'(1p,6e20.12)') phig(ip)
  end do
  close(14)
!
end subroutine IO_output_solution_3D_CF_spin_NS_mpt
