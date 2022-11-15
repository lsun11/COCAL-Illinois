module interface_compute_dBfun
  implicit none
  interface 
    subroutine compute_dBfun(Bfun,dBfundx,dBfundy,dBfundz)
      real(8), pointer :: Bfun(:,:,:), dBfundx(:,:,:), dBfundy(:,:,:), dBfundz(:,:,:) 
    end subroutine compute_dBfun
  end interface
end module interface_compute_dBfun
