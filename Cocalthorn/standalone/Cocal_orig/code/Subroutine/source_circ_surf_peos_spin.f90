subroutine source_circ_surf_peos_spin(cline, ia, ib, souf)
  use phys_constant, only :  long
  use grid_parameter
  use def_matter, only : vepxf, vepyf, vepzf
  use def_matter_parameter, only : ome, ber, confpow
  use def_metric, only : psi
  use def_metric_on_SFC_CF, only : psif
  use def_velocity_rot
!  use coordinate_grav_r, only : rg
!  use def_vector_phi, only : vec_phif
  use trigonometry_grav_phi, only : sinphig, cosphig
  use interface_interpo_gr2fl
  use interface_flgrad_midpoint_type0_parallel
  use interface_flgrad_midpoint_type0_meridian
  use make_array_2d
  implicit none
  character(len=2), intent(in) :: cline
  integer,          intent(in) :: ia, ib
  real(long), pointer :: souf(:,:)
  real(long), pointer :: fx(:,:), fy(:,:), fz(:,:)
  real(long), pointer :: gx(:,:), gy(:,:), gz(:,:)
  real(long) :: dfxdx, dfxdy, dfxdz
  real(long) :: dfydx, dfydy, dfydz
  real(long) :: dfzdx, dfzdy, dfzdz
  real(long) :: omefc, jomef_intfc, ui(3), vphif(3)
  real(long) :: hh, ut, pre, rho, ene, qq, psif4
  integer    :: irf, itf, ipf, irs, j, ir, it
!
  call interpo_gr2fl(psi , psif )

  if (cline=="ph") then   ! Circulation on any parallel plane, radius
    call alloc_array2d(fx, 0, ia, 0, npf)
    call alloc_array2d(fy, 0, ia, 0, npf)
!    call alloc_array2d(fz, 0, ia, 0, npf)
    itf = ib
    souf = 0.0d0;  fx=0.0d0;  fy=0.0d0
    do irf = 0, ia
      do ipf = 0, npf
        psif4 = psif(irf,itf,ipf)**(4+confpow)

!       hu_i = partial_i(Phi)
        fx(irf,ipf) = vepxf(irf,itf,ipf) + psif4*wxspf(irf,itf,ipf)
        fy(irf,ipf) = vepyf(irf,itf,ipf) + psif4*wyspf(irf,itf,ipf)
!        fz(irf,ipf) = hh*ui(3)
      end do
    end do
    do irf = 1, ia
      do ipf = 1, npf
        call flgrad_midpoint_type0_parallel(fx,dfxdx,dfxdy,dfxdz,irf,itf,ipf)
        call flgrad_midpoint_type0_parallel(fy,dfydx,dfydy,dfydz,irf,itf,ipf)
        souf(irf,ipf) = dfydx - dfxdy
      end do
    end do 
    deallocate(fx); deallocate(fy); !deallocate(fz);
  end if

  if (cline=="th") then   ! Circulation on any meridian plane, radius
    call alloc_array2d(fx, 0, ia, 0, ntf)
    call alloc_array2d(fy, 0, ia, 0, ntf)
    call alloc_array2d(fz, 0, ia, 0, ntf)
    call alloc_array2d(gx, 0, ia, 0, ntf)
    call alloc_array2d(gy, 0, ia, 0, ntf)
    call alloc_array2d(gz, 0, ia, 0, ntf)
    souf = 0.0d0;  
    fx=0.0d0;  fy=0.0d0; fz=0.0d0
    gx=0.0d0;  gy=0.0d0; gz=0.0d0

    do irf = 0, ia
      ipf = ib
      do itf = 0, ntf      !  Half of the disk 
        psif4 = psif(irf,itf,ipf)**(4+confpow)

        fx(irf,itf) = vepxf(irf,itf,ipf) + psif4*wxspf(irf,itf,ipf)
        fy(irf,itf) = vepyf(irf,itf,ipf) + psif4*wyspf(irf,itf,ipf)
        fz(irf,itf) = vepzf(irf,itf,ipf) + psif4*wzspf(irf,itf,ipf)
      end do

      ipf = ib + npf/2   !  The other half of the disk

      do it = ntf, 2*ntf
        itf = ntf - (it - ntf)  ! gridpoint: ntf -> 0
        psif4 = psif(irf,itf,ipf)**(4+confpow)

        gx(irf,itf) = vepxf(irf,itf,ipf) + psif4*wxspf(irf,itf,ipf)
        gy(irf,itf) = vepyf(irf,itf,ipf) + psif4*wyspf(irf,itf,ipf)
        gz(irf,itf) = vepzf(irf,itf,ipf) + psif4*wzspf(irf,itf,ipf)
      end do
    end do

!   souf(1:ia, 1:2*ntf)
    do ir = 1, ia
      irf = ir
      ipf = ib
      do it = 1, ntf        
        itf = it
        call flgrad_midpoint_type0_meridian(fx,dfxdx,dfxdy,dfxdz,irf,itf,ipf)
        call flgrad_midpoint_type0_meridian(fy,dfydx,dfydy,dfydz,irf,itf,ipf)
        call flgrad_midpoint_type0_meridian(fz,dfzdx,dfzdy,dfzdz,irf,itf,ipf)

!       Normal vector to phi=const. plane is ( -sin(ph), cos(ph), 0)
        souf(ir,it) = -sinphig(ib) * (dfzdy - dfydz) + &
                    &  cosphig(ib) * (dfxdz - dfzdx) 
      end do

      ipf = ib + npf/2
      do it = ntf+1, 2*ntf  
        itf = ntf - (it - ntf) + 1     ! midpoint  ntf -> 1
        call flgrad_midpoint_type0_meridian(gx,dfxdx,dfxdy,dfxdz,irf,itf,ipf)
        call flgrad_midpoint_type0_meridian(gy,dfydx,dfydy,dfydz,irf,itf,ipf)
        call flgrad_midpoint_type0_meridian(gz,dfzdx,dfzdy,dfzdz,irf,itf,ipf)

        souf(ir,it) = -sinphig(ib) * (dfzdy - dfydz) + &
                    &  cosphig(ib) * (dfxdz - dfzdx)
      end do
    end do
    deallocate(fx); deallocate(fy); deallocate(fz);
    deallocate(gx); deallocate(gy); deallocate(gz);
  end if

end subroutine source_circ_surf_peos_spin
