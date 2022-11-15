module interface_update_grfield_inAH
  implicit none
  interface 
    subroutine update_grfield_inAH(pot,grfield,convf)
      real(8), pointer :: pot(:,:,:)
      real(8), pointer :: grfield(:,:,:)
      real(8), intent(in) :: convf
    end subroutine update_grfield_inAH
  end interface
end module interface_update_grfield_inAH
