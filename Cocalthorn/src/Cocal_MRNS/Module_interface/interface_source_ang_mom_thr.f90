module interface_source_ang_mom_thr
  implicit none
  interface 
    subroutine source_ang_mom_thr(sous)
      real(8), pointer  :: sous(:,:)
    end subroutine source_ang_mom_thr
  end interface
end module interface_source_ang_mom_thr
