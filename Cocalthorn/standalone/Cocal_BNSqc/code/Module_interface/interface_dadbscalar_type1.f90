module interface_dadbscalar_type1
  implicit none
  interface 
    subroutine dadbscalar_type1(fnc,d2fdxdx,d2fdxdy,d2fdxdz,d2fdydx,d2fdydy,d2fdydz,d2fdzdx,d2fdzdy,d2fdzdz,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: d2fdxdx,d2fdxdy,d2fdxdz,d2fdydx,d2fdydy,d2fdydz,d2fdzdx,d2fdzdy,d2fdzdz
      integer :: irg, itg, ipg
    end subroutine dadbscalar_type1
  end interface
end module interface_dadbscalar_type1
