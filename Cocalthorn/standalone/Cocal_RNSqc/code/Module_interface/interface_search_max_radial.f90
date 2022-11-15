module interface_search_max_radial
  implicit none
  interface 
    subroutine search_max_radial(iremax,fnc,fncmax,xsol)
      real(8), pointer :: fnc(:)
      integer :: iremax
      real(8) :: fncmax, xsol
    end subroutine search_max_radial
  end interface
end module interface_search_max_radial
