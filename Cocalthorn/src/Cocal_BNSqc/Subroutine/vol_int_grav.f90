subroutine vol_int_grav(soug,vol)
  use phys_constant, only  :   long
  use grid_parameter, only  :   nrg, ntg, npg
  use weight_midpoint_grav, only : hwrtpg
  implicit none
  integer     ::   ir,it,ip
  real(long)  ::   vol, hsoug
  real(long), pointer  ::   soug(:,:,:)
  vol = 0.0d0
  do ip = 1, npg
    do it = 1, ntg
      do ir = 1, nrg
        vol = vol + soug(ir,it,ip)*hwrtpg(ir,it,ip)
      end do
    end do
  end do
end subroutine vol_int_grav
