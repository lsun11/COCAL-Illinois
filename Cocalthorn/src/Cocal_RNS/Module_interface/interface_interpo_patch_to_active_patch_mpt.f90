module interface_interpo_patch_to_active_patch_mpt
  implicit none
  interface
    subroutine interpo_patch_to_active_patch_mpt(fnc,cfn,ra,tha,phia)
      real(8), pointer     :: fnc(:,:,:)
      real(8), intent(out) :: cfn
      real(8), intent(in)  :: ra, tha, phia
    end subroutine interpo_patch_to_active_patch_mpt
  end interface
end module interface_interpo_patch_to_active_patch_mpt
