module interface_interpolation_fillup_binary_parity_mpt
  implicit none
  interface
    subroutine interpolation_fillup_binary_parity_mpt(fnc,fnc_ex,parchar)
      real(8), pointer    :: fnc(:,:,:), fnc_ex(:,:,:)
      character(len=2) :: parchar
    end subroutine interpolation_fillup_binary_parity_mpt
  end interface
end module interface_interpolation_fillup_binary_parity_mpt
