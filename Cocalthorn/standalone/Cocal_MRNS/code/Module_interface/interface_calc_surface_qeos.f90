module interface_calc_surface_qeos
  implicit none
  interface 
    subroutine calc_surface_qeos(rsnew,rhof)
      real(8), pointer :: rhof(:,:,:)
      real(8), pointer :: rsnew(:,:)
    end subroutine calc_surface_qeos
  end interface
end module interface_calc_surface_qeos
