module interface_calc_surface_normal_tangent_midpoint
  implicit none
  interface 
    subroutine calc_surface_normal_tangent_midpoint(rs,nv,tv,itf,ipf)
      real(8), pointer :: rs(:,:), nv(:), tv(:)
      integer, intent(in)  :: itf,ipf
    end subroutine calc_surface_normal_tangent_midpoint
  end interface
end module interface_calc_surface_normal_tangent_midpoint
