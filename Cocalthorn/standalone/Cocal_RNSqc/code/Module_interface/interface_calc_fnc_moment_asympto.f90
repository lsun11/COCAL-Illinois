module interface_calc_fnc_moment_asympto
  implicit none
  interface 
    subroutine calc_fnc_moment_asympto(fnc,fnc_moment)
      real(8), pointer :: fnc(:,:,:)
      real(8) :: fnc_moment(3)
    end subroutine calc_fnc_moment_asympto
  end interface
end module interface_calc_fnc_moment_asympto
