module interface_source_ang_mom_inf
  implicit none
  interface 
    subroutine source_ang_mom_inf(sous, irg)
      real(8), pointer  :: sous(:,:)
      integer :: irg
    end subroutine source_ang_mom_inf
  end interface
end module interface_source_ang_mom_inf
