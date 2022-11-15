module interface_source_ang_mom_qeos
  implicit none
  interface 
    subroutine source_ang_mom_qeos(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_ang_mom_qeos
  end interface
end module interface_source_ang_mom_qeos
