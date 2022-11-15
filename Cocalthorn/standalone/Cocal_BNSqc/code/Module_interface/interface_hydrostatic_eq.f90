module interface_hydrostatic_eq
  implicit none
  interface 
    subroutine hydrostatic_eq(emd)
      real(8), pointer :: emd(:,:,:)
    end subroutine hydrostatic_eq
  end interface
end module interface_hydrostatic_eq
