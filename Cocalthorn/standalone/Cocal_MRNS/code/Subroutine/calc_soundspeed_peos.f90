subroutine calc_soundspeed_peos
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, ntfeq, ntfpolp, npfxzp, npfyzp
  use def_matter, only : emd
  use make_array_3d
  use interface_IO_output_1D_general
  implicit none
  real(long), pointer  :: c_s(:,:,:)
  real(long) :: q, ss
  integer :: irf, itf, ipf
  character(30) :: char1 
!
  call alloc_array3d(c_s, 0, nrf, 0, ntf, 0, npf)

  open(30,file='rns_ss_gt_1.dat',status='unknown')
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        q = emd(irf,itf,ipf)
        call peos_sound_speed(q,ss)
        c_s(irf,itf,ipf) = ss
        if (ss > 1)  write(30,'(3i5,1p,e23.15)') irf, itf, ipf, ss
      end do
    end do
  end do
  close(30)
!
  char1 = 'rns_sound_speed_pxaxis.dat'
  call IO_output_1D_general(char1,'f','g',c_s,-1, ntfeq, npfxzp)
!
  deallocate(c_s)
!
end subroutine calc_soundspeed_peos
