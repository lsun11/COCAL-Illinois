module interface_interpolation_fillup_binary_parity
  implicit none
  interface
    subroutine interpolation_fillup_binary_parity(fnc,par)
      real(8), pointer    :: fnc(:,:,:)
      real(8), intent(in) :: par
    end subroutine interpolation_fillup_binary_parity
  end interface
end module interface_interpolation_fillup_binary_parity
