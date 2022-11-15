module interface_calc_surface
  implicit none
  interface 
    subroutine calc_surface(rsnew,emd)
      real(8), pointer :: emd(:,:,:)
      real(8), pointer :: rsnew(:,:)
    end subroutine calc_surface
  end interface
end module interface_calc_surface
