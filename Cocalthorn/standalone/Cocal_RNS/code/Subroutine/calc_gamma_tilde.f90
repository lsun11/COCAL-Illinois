subroutine calc_gamma_tilde(iter_count)
  use phys_constant, only :  long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : psi
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use interface_IO_output_1D_general
  use make_array_3d
  implicit none
  real(long) :: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz, &
     &          hod1, hod2, hod3, detgm, detgmi, &
     &          gmxxu, gmxyu, gmxzu, gmyxu, gmyyu, gmyzu, &
     &          gmzxu, gmzyu, gmzzu
  integer :: ipg, itg, irg, iter_count
  real(long), pointer :: gamt(:,:,:) 
  character(30) :: char1, char2, char3, char4, char5

  call alloc_array3d(gamt,0,nrg,0,ntg,0,npg)
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

        gamt(irg,itg,ipg)  = 1.0d0 + hod1 + hod2 + hod3

      end do
    end do
  end do
!
  write(char1, '(i5)') iter_count
  char2 = adjustl(char1)
  char3 = 'gamt_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', gamt, -1, ntg/2,0)
  char3 = 'gamt_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', gamt, -1, ntg/2,npg/4)
  char3 = 'gamt_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', gamt, -1, 0,0)
!
  deallocate(gamt)
!
end subroutine calc_gamma_tilde
