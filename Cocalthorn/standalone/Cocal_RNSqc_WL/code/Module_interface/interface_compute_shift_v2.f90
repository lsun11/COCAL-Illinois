module interface_compute_shift_v2
  implicit none
  interface 
    subroutine compute_shift_v2(Bfun,potx,poty,potz)
      real(8), pointer :: Bfun(:,:,:), potx(:,:,:), poty(:,:,:), potz(:,:,:) 
    end subroutine compute_shift_v2
  end interface
end module interface_compute_shift_v2
