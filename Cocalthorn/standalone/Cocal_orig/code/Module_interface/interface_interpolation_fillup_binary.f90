module interface_interpolation_fillup_binary
  implicit none
  interface
    subroutine interpolation_fillup_binary(fnc)
      real(8), pointer :: fnc(:,:,:)
    end subroutine interpolation_fillup_binary
  end interface
end module interface_interpolation_fillup_binary
