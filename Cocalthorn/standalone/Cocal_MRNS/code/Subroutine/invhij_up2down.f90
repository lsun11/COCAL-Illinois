subroutine invhij_up2down
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  implicit none
  real(8) :: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz, &
  &          hod1, hod2, hod3, detgm, detgmi, &
  &          gmxxd, gmxyd, gmxzd, gmyxd, gmyyd, gmyzd, &
  &          gmzxd, gmzyd, gmzzd
  integer :: ipg, itg, irg
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        hxx = hxxu(irg,itg,ipg)
        hxy = hxyu(irg,itg,ipg)
        hxz = hxzu(irg,itg,ipg)
        hyx = hxy
        hyy = hyyu(irg,itg,ipg)
        hyz = hyzu(irg,itg,ipg)
        hzx = hxz
        hzy = hyz
        hzz = hzzu(irg,itg,ipg)
!
        hod1 = hxx + hyy + hzz
        hod2 = hxx*hyy + hxx*hzz + hyy*hzz &
           &     - hxy*hyx - hxz*hzx - hyz*hzy
        hod3 = hxx*hyy*hzz + hxy*hyz*hzx + hxz*hyx*hzy &
           &     - hxx*hyz*hzy - hxy*hyx*hzz - hxz*hyy*hzx
        detgm  = 1.0d0 + hod1 + hod2 + hod3
        detgmi = 1.0d0/detgm
!
        hod1  = + hyy + hzz
        hod2  = + hyy*hzz - hyz*hzy
        gmxxd = (1.0d0 + hod1 + hod2)*detgmi
        hod1  = - hxy
        hod2  = + hxz*hzy - hxy*hzz
        gmxyd = (hod1 + hod2)*detgmi
        hod1  = - hxz
        hod2  = + hxy*hyz - hxz*hyy
        gmxzd = (hod1 + hod2)*detgmi
        hod1  = + hxx + hzz
        hod2  = + hxx*hzz - hxz*hzx
        gmyyd = (1.0d0 + hod1 + hod2)*detgmi
        hod1  = - hyz
        hod2  = + hxz*hyx - hxx*hyz
        gmyzd = (hod1 + hod2)*detgmi
        hod1  = + hxx + hyy
        hod2  = + hxx*hyy - hxy*hyx
        gmzzd = (1.0d0 + hod1 + hod2)*detgmi
        gmyxd = gmxyd
        gmzxd = gmxzd
        gmzyd = gmyzd
!
        hxxd(irg,itg,ipg) = gmxxd - 1.0d0
        hxyd(irg,itg,ipg) = gmxyd
        hxzd(irg,itg,ipg) = gmxzd
        hyyd(irg,itg,ipg) = gmyyd - 1.0d0
        hyzd(irg,itg,ipg) = gmyzd
        hzzd(irg,itg,ipg) = gmzzd - 1.0d0
!
      end do
    end do
  end do
!
end subroutine invhij_up2down
