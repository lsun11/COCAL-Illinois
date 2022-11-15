module interface_source_ang_mom_peos
  implicit none
  interface 
    subroutine source_ang_mom_peos(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_ang_mom_peos
  end interface
end module interface_source_ang_mom_peos
