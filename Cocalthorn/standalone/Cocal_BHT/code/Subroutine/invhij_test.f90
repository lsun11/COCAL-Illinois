subroutine invhij
  use grid_parameter, only : nrg, ntg, npg
  use def_metric
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu, trh
  use interface_IO_output_1D_general
  implicit none
  real(8) :: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz, &
  &          hod1, hod2, hod3, detgm, detgmi, &
  &          gmxxu, gmxyu, gmxzu, gmyxu, gmyyu, gmyzu, &
  &          gmzxu, gmzyu, gmzzu, ps4oal2, bbxx, bbxy, bbxz, bbyy, & 
  &          bbyz, bbzz
  integer :: ipg, itg, irg, iter_count
  character(30) :: char1, char2, char3, char4, char5
!
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
        gmxxu = (1.0d0 + hod1 + hod2)*detgmi
        hod1  = - hxy
        hod2  = + hxz*hzy - hxy*hzz
        gmxyu = (hod1 + hod2)*detgmi
        hod1  = - hxz
        hod2  = + hxy*hyz - hxz*hyy
        gmxzu = (hod1 + hod2)*detgmi
        hod1  = + hxx + hzz
        hod2  = + hxx*hzz - hxz*hzx
        gmyyu = (1.0d0 + hod1 + hod2)*detgmi
        hod1  = - hyz
        hod2  = + hxz*hyx - hxx*hyz
        gmyzu = (hod1 + hod2)*detgmi
        hod1  = + hxx + hyy
        hod2  = + hxx*hyy - hxy*hyx
        gmzzu = (1.0d0 + hod1 + hod2)*detgmi
        gmyxu = gmxyu
        gmzxu = gmxzu
        gmzyu = gmyzu
!
!        ps4oal2 = psi(irg,itg,ipg)**4/alph(irg,itg,ipg)**2
!        bbxx = ps4oal2*bvxu(irg,itg,ipg)*bvxu(irg,itg,ipg)
!        bbxy = ps4oal2*bvxu(irg,itg,ipg)*bvyu(irg,itg,ipg)
!        bbxz = ps4oal2*bvxu(irg,itg,ipg)*bvzu(irg,itg,ipg)
!        bbyy = ps4oal2*bvyu(irg,itg,ipg)*bvyu(irg,itg,ipg)
!        bbyz = ps4oal2*bvyu(irg,itg,ipg)*bvzu(irg,itg,ipg)
!        bbzz = ps4oal2*bvzu(irg,itg,ipg)*bvzu(irg,itg,ipg)

        hxxu(irg,itg,ipg) = gmxxu - 1.0d0 !- bbxx
        hxyu(irg,itg,ipg) = gmxyu         !- bbxy
        hxzu(irg,itg,ipg) = gmxzu         !- bbxz
        hyyu(irg,itg,ipg) = gmyyu - 1.0d0 !- bbyy
        hyzu(irg,itg,ipg) = gmyzu         !- bbyz
        hzzu(irg,itg,ipg) = gmzzu - 1.0d0 !- bbzz
!
        trh(irg,itg,ipg)  = hxx + hyy + hzz
      end do
    end do
  end do
!
!  iter_count=1
!  write(char1, '(i5)') iter_count
!  char2 = adjustl(char1)
!  char3 = 'hxxu' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hxxu, -1, ntg/2, 0)
!  char3 = 'hxyu' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hxyu, -1, ntg/2, 0)
!  char3 = 'hxzu' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hxzu, -1, ntg/2, 0)
!  char3 = 'hyyu' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hyyu, -1, ntg/2, 0)
!  char3 = 'hyzu' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hyzu, -1, ntg/2, 0)
!  char3 = 'hzzu' // trim(char2) // '.txt'
!  call IO_output_1D_general(char3, 'g', 'g', hzzu, -1, ntg/2, 0)
!
end subroutine invhij
