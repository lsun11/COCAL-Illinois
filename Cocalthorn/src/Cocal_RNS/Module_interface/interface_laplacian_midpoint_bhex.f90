module interface_laplacian_midpoint_bhex
  implicit none
  interface 
    subroutine laplacian_midpoint_bhex(fnc,lapfnc)
      real(8), pointer :: fnc(:,:,:), lapfnc(:,:,:)
    end subroutine laplacian_midpoint_bhex
  end interface
end module interface_laplacian_midpoint_bhex
