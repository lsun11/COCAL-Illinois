module interface_grdr_gridpoint_type0
  implicit none
  interface 
    subroutine grdr_gridpoint_type0(fnc,deriv,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: deriv
      integer :: irg, itg, ipg
    end subroutine grdr_gridpoint_type0
  end interface
end module interface_grdr_gridpoint_type0
