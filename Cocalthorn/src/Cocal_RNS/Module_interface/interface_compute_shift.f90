module interface_compute_shift
  implicit none
  interface 
    subroutine compute_shift(potx,poty,potz,gvec,bfnc)      
      real(8), pointer :: bfnc(:,:,:) 
      real(8), pointer :: gvec(:,:,:,:)
      real(8), pointer :: potx(:,:,:), poty(:,:,:), potz(:,:,:)
    end subroutine compute_shift
  end interface
end module interface_compute_shift
