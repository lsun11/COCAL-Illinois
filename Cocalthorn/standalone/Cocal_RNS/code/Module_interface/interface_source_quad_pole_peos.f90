module interface_source_quad_pole_peos
  implicit none
  interface 
    subroutine source_quad_pole_peos(souf5d, iquad)
      real(8), pointer     :: souf5d(:,:,:,:,:)
      integer :: iquad
    end subroutine source_quad_pole_peos
  end interface
end module interface_source_quad_pole_peos
