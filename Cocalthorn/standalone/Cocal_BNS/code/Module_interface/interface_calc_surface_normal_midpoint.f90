module interface_calc_surface_normal_midpoint
  implicit none
  interface 
    subroutine calc_surface_normal_midpoint(rs,rnx,rny,rnz,it,ip)
      real(8), pointer :: rs(:,:)
      real(8), intent(out) :: rnx, rny, rnz
      integer, intent(in)  :: it,ip
    end subroutine calc_surface_normal_midpoint
  end interface
end module interface_calc_surface_normal_midpoint
