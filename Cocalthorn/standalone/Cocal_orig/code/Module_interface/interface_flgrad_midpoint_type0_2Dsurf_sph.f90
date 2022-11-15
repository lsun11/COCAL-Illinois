module interface_flgrad_midpoint_type0_2Dsurf_sph
  implicit none
  interface 
    subroutine flgrad_midpoint_type0_2Dsurf_sph(fnc,dfdth,dfdphi,itf,ipf)
      real(8), pointer :: fnc(:,:)
      real(8) :: dfdth,dfdphi
      integer :: itf, ipf
    end subroutine flgrad_midpoint_type0_2Dsurf_sph
  end interface
end module interface_flgrad_midpoint_type0_2Dsurf_sph
