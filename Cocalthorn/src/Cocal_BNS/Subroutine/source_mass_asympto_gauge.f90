subroutine source_mass_asympto_gauge(sousf1,irg)
  use phys_constant, only  : long
  use grid_parameter, only : ntg, npg
  use def_formulation, only : chgra, chope
  use def_kerr_schild, only : Fxu_grid_ks, Fyu_grid_ks, Fzu_grid_ks
  use trigonometry_grav_theta, only : hsinthg, hcosthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use make_array_2d
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: fx(:,:), fy(:,:), fz(:,:), sousf1(:,:)
  integer    :: irg
  integer    :: itg, ipg
  real(long) :: valx, valy, valz, xx, yy, zz
!
  call alloc_array2d(fx, 0,ntg, 0,npg)
  call alloc_array2d(fy, 0,ntg, 0,npg)
  call alloc_array2d(fz, 0,ntg, 0,npg)
!
  if (chgra .eq. 'k') then
    do ipg = 0, npg
      do itg = 0, ntg
        fx(itg,ipg) = Fxu_grid_ks(irg,itg,ipg)
        fy(itg,ipg) = Fyu_grid_ks(irg,itg,ipg)
        fz(itg,ipg) = Fzu_grid_ks(irg,itg,ipg)
      end do
    end do
  else    !  this part can change to anything
    fx = 0.0d0
    fy = 0.0d0
    fz = 0.0d0
  end if
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(valx,fx,itg,ipg)
      call interpo_linear_type0_2Dsurf(valy,fy,itg,ipg)
      call interpo_linear_type0_2Dsurf(valz,fz,itg,ipg)

      xx = hsinthg(itg)*hcosphig(ipg)
      yy = hsinthg(itg)*hsinphig(ipg)
      zz = hcosthg(itg)

      sousf1(itg,ipg) = xx*valx + yy*valy + zz*valz
    end do
  end do
!
  deallocate(fx)
  deallocate(fy)
  deallocate(fz)
end subroutine source_mass_asympto_gauge
