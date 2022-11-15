module interface_interpolation_fillup_binary_mpt
  implicit none
  interface
    subroutine interpolation_fillup_binary_mpt(fnc,fnc_ex)
      real(8), pointer :: fnc(:,:,:), fnc_ex(:,:,:)
    end subroutine interpolation_fillup_binary_mpt
  end interface
end module interface_interpolation_fillup_binary_mpt
