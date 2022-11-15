module interface_hydrostatic_eq_WL_peos
  implicit none
  interface 
    subroutine hydrostatic_eq_WL_peos(emd)
      real(8), pointer :: emd(:,:,:)
    end subroutine hydrostatic_eq_WL_peos
  end interface
end module interface_hydrostatic_eq_WL_peos
