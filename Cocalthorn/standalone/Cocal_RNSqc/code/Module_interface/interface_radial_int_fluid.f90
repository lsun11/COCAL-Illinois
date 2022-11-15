module interface_radial_int_fluid
  implicit none
  interface 
    subroutine radial_int_fluid(sou,radius,it,ip)
      real(8), pointer     :: sou(:)
      real(8), intent(out) :: radius
      integer, intent(in)  :: it, ip
    end subroutine radial_int_fluid
  end interface
end module interface_radial_int_fluid
