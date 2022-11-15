module interface_copy_to_hgfn_and_gfnsf
  implicit none
  interface 
    subroutine copy_to_hgfn_and_gfnsf(hgfn_bc,gfnsf_bc)
      real(8), pointer :: hgfn_bc(:,:,:), gfnsf_bc(:,:,:)
    end subroutine copy_to_hgfn_and_gfnsf
  end interface
end module interface_copy_to_hgfn_and_gfnsf
