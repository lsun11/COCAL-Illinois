module interface_source_ang_mom_smarr
  implicit none
  interface 
    subroutine source_ang_mom_smarr(sous)
      real(8), pointer  :: sous(:,:)
    end subroutine source_ang_mom_smarr
  end interface
end module interface_source_ang_mom_smarr
