subroutine adjusthij
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : psi
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  implicit none
  real(8) :: detgm, psi4, psinew, psiold, twelve, &
  &          gmxx, gmxy, gmxz, gmyx, gmyy, gmyz, gmzx, gmzy, gmzz, &
  &          hod1, hod2, hod3, &
  &          hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
  integer :: ipg, itg, irg
!
!write(6,*) psi(10,ntg/2,0), hxxd(10,ntg/2,0), hyyd(10,ntg/2,0)
  twelve = 1.0d0/12.0d0
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        hxx = hxxd(irg,itg,ipg)
        hxy = hxyd(irg,itg,ipg)
        hxz = hxzd(irg,itg,ipg)
        hyx = hxy
        hyy = hyyd(irg,itg,ipg)
        hyz = hyzd(irg,itg,ipg)
        hzx = hxz
        hzy = hyz
        hzz = hzzd(irg,itg,ipg)
!
        gmxx = 1.0d0 + hxx
        gmxy = hxy
        gmxz = hxz
        gmyx = hxy
        gmyy = 1.0d0 + hyy
        gmyz = hyz
        gmzx = hxz
        gmzy = hyz
        gmzz = 1.0d0 + hzz
!
        hod1 = hxx + hyy + hzz
        hod2 = hxx*hyy + hxx*hzz + hyy*hzz &
        &    - hxy*hyx - hxz*hzx - hyz*hzy
        hod3 = hxx*hyy*hzz + hxy*hyz*hzx + hxz*hyx*hzy &
        &    - hxx*hyz*hzy - hxy*hyx*hzz - hxz*hyy*hzx
        detgm = 1.0d0 + hod1 + hod2 + hod3
!
        psiold = psi(irg,itg,ipg)
        psinew = psiold*detgm**twelve
        psi4 = (psiold/psinew)**4
!
!if(irg.eq.10.and.itg.eq.ntg/2.and.ipg.eq.0) &
!write(6,*) psiold, psinew, detgm
        psi(irg,itg,ipg) = psinew
        hxxd(irg,itg,ipg) = gmxx*psi4 - 1.0d0
        hxyd(irg,itg,ipg) = gmxy*psi4
        hxzd(irg,itg,ipg) = gmxz*psi4
        hyyd(irg,itg,ipg) = gmyy*psi4 - 1.0d0
        hyzd(irg,itg,ipg) = gmyz*psi4
        hzzd(irg,itg,ipg) = gmzz*psi4 - 1.0d0
!
      end do
    end do
  end do
!
end subroutine adjusthij
