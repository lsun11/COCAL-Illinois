subroutine calc_dihij_gridpoint_bhex(iter_count)
  use phys_constant, only :  long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_metric_dihiju,  only : dihixu, dihiyu, dihizu
  use interface_grgrad_4th_gridpoint_bhex
  use interface_IO_output_1D_general
  use make_array_3d
  implicit none
  real(long), pointer :: Fxu(:,:,:), Fyu(:,:,:), Fzu(:,:,:)
  real(long) :: dhxxux, dhxxuy, dhxxuz, dhxyux, dhxyuy, dhxyuz, &
  &             dhxzux, dhxzuy, dhxzuz, dhyyux, dhyyuy, dhyyuz, &
  &             dhyzux, dhyzuy, dhyzuz, dhzzux, dhzzuy, dhzzuz, &
  &             hxidiv, hyidiv, hzidiv
  integer :: ipg, irg, itg, iter_count
  character(30) :: char1, char2, char3, char4, char5
!
  call alloc_array3d(Fxu,0,nrg,0,ntg,0,npg)
  call alloc_array3d(Fyu,0,nrg,0,ntg,0,npg)
  call alloc_array3d(Fzu,0,nrg,0,ntg,0,npg)

  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        call grgrad_4th_gridpoint_bhex(hxxu,dhxxux,dhxxuy,dhxxuz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hxyu,dhxyux,dhxyuy,dhxyuz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hxzu,dhxzux,dhxzuy,dhxzuz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hyyu,dhyyux,dhyyuy,dhyyuz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hyzu,dhyzux,dhyzuy,dhyzuz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hzzu,dhzzux,dhzzuy,dhzzuz,irg,itg,ipg)
        Fxu(irg,itg,ipg) = dhxxux + dhxyuy + dhxzuz
        Fyu(irg,itg,ipg) = dhxyux + dhyyuy + dhyzuz
        Fzu(irg,itg,ipg) = dhxzux + dhyzuy + dhzzuz    
      end do
    end do
  end do
!
  write(char1, '(i5)') iter_count
  char2 = adjustl(char1)

  char3 = 'Fxu_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', Fxu, -1, ntg/2,0)
  char3 = 'Fxu_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', Fxu, -1, ntg/2,npg/4)
  char3 = 'Fxu_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', Fxu, -1, 0,0)

  char3 = 'Fyu_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', Fyu, -1, ntg/2,0)
  char3 = 'Fyu_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', Fyu, -1, ntg/2,npg/4)
  char3 = 'Fyu_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', Fyu, -1, 0,0)

  char3 = 'Fzu_xa' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', Fzu, -1, ntg/2,0)
  char3 = 'Fzu_ya' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', Fzu, -1, ntg/2,npg/4)
  char3 = 'Fzu_za' // trim(char2) // '.txt'
  call IO_output_1D_general(char3, 'g', 'g', Fzu, -1, 0,0)
!
  deallocate(Fxu)
  deallocate(Fyu)
  deallocate(Fzu)
!
end subroutine calc_dihij_gridpoint_bhex
