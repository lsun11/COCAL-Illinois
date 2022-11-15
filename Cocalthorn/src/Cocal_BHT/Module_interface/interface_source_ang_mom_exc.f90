module interface_source_ang_mom_exc
  implicit none
  interface 
    subroutine source_ang_mom_exc(sous)
      real(8), pointer  :: sous(:,:)
    end subroutine source_ang_mom_exc
  end interface
end module interface_source_ang_mom_exc
