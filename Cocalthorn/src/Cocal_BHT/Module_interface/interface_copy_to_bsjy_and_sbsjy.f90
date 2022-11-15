module interface_copy_to_bsjy_and_sbsjy
  implicit none
  interface 
    subroutine copy_to_bsjy_and_sbsjy(bsjy_bc,sbsjy_bc,sbsjyp_bc)
      real(8), pointer :: bsjy_bc(:,:,:,:)
      real(8), pointer :: sbsjy_bc(:,:,:), sbsjyp_bc(:,:,:)
    end subroutine copy_to_bsjy_and_sbsjy
  end interface
end module interface_copy_to_bsjy_and_sbsjy
