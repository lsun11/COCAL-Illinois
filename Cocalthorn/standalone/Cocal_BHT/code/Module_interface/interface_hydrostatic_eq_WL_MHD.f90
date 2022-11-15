module interface_hydrostatic_eq_WL_MHD
  implicit none
  interface 
    subroutine hydrostatic_eq_WL_MHD(emd,utf,uxf,uyf,uzf)
      real(8), pointer :: emd(:,:,:)
      real(8), pointer :: utf(:,:,:), uxf(:,:,:), uyf(:,:,:), uzf(:,:,:)
    end subroutine hydrostatic_eq_WL_MHD
  end interface
end module interface_hydrostatic_eq_WL_MHD
