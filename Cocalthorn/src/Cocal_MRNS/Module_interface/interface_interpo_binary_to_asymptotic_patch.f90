module interface_binary_to_asymptotic_patch
  implicit none
  interface
    subroutine interpo_binary_to_asymptotic_patch(fnc,cfn,ira,ita,ipa)
      real(8), pointer     :: fnc(:,:,:)
      real(8), intent(out) :: cfn
      integer, intent(in)  :: ira, ita, ipa
    end subroutine interpo_binary_to_asymptotic_patch
  end interface
end module interface_interpo_binary_to_asymptotic_patch
