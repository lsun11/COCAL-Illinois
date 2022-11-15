module interface_source_quad_pole_qeos
  implicit none
  interface 
    subroutine source_quad_pole_qeos(souf5d, iquad)
      real(8), pointer     :: souf5d(:,:,:,:,:)
      integer :: iquad
    end subroutine source_quad_pole_qeos
  end interface
end module interface_source_quad_pole_qeos
