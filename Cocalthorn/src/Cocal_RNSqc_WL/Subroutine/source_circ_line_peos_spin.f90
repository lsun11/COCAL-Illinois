subroutine source_circ_line_peos_spin(cline, ia, ib, souf)
  use phys_constant, only :  long
  use grid_parameter
  use def_matter, only : vepxf, vepyf, vepzf
  use def_matter_parameter, only : ome, ber, confpow
  use def_metric, only : psi
  use def_metric_on_SFC_CF, only : psif
  use def_velocity_rot
!  use coordinate_grav_r, only : rg
!  use def_vector_phi, only : vec_phif
  use interface_interpo_gr2fl
  use make_array_2d
  implicit none
  character(len=2), intent(in) :: cline
  integer,          intent(in) :: ia, ib
  real(long), pointer :: souf(:,:)
  real(long), pointer :: fph(:,:), fth(:,:)
  real(long) :: omefc, jomef_intfc, ui(3), vphif(3)
  real(long) :: hh, ut, pre, rho, ene, qq, psif4
  integer    :: irf, itf, ipf, irs, j, it
!
  call interpo_gr2fl(psi , psif )

  if (cline=="ph") then         ! Circulation on any parallel, radius
    call alloc_array2d(fph, 0, npf, 1, 3)
    irf = ia
    itf = ib
    souf = 0.0d0;   fph=0.0d0
    do ipf = 0, npf
      psif4 = psif(irf,itf,ipf)**(4+confpow)

!     hu_i = partial_i(Phi)
      fph(ipf,1) = vepxf(irf,itf,ipf) + psif4*wxspf(irf,itf,ipf)
      fph(ipf,2) = vepyf(irf,itf,ipf) + psif4*wyspf(irf,itf,ipf)
      fph(ipf,3) = vepzf(irf,itf,ipf) + psif4*wzspf(irf,itf,ipf)
    end do
    do ipf = 1, npf
      souf(ipf,1) = 0.5d0*( fph(ipf,1) + fph(ipf-1,1) )
      souf(ipf,2) = 0.5d0*( fph(ipf,2) + fph(ipf-1,2) )
      souf(ipf,3) = 0.5d0*( fph(ipf,3) + fph(ipf-1,3) )
    end do 
    deallocate(fph)
  end if

  if (cline=="th") then         ! Circulation on any meridian, radius
    call alloc_array2d(fth, 0, 2*ntf, 1, 3)
    irf = ia
    ipf = ib
    souf = 0.0d0;   fth=0.0d0
    do it = 0, 2*ntf
      itf = it
      if (it > ntf) then
        itf = ntf - (it - ntf)  ! gridpoint 
        ipf = ib + npf/2        ! gridpoint
      end if
!      write(6,'(a3,4i5)') "it=", it,irf,itf,ipf
      psif4 = psif(irf,itf,ipf)**(4+confpow)

      fth(it,1) = vepxf(irf,itf,ipf) + psif4*wxspf(irf,itf,ipf)
      fth(it,2) = vepyf(irf,itf,ipf) + psif4*wyspf(irf,itf,ipf)
      fth(it,3) = vepzf(irf,itf,ipf) + psif4*wzspf(irf,itf,ipf)
    end do 
    do it = 1, 2*ntf
      souf(it,1) = 0.5d0*( fth(it,1) + fth(it-1,1) )
      souf(it,2) = 0.5d0*( fth(it,2) + fth(it-1,2) )
      souf(it,3) = 0.5d0*( fth(it,3) + fth(it-1,3) )
    end do 
    deallocate(fth)
  end if


end subroutine source_circ_line_peos_spin
