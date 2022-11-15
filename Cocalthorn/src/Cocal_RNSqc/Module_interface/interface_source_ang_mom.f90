module interface_source_ang_mom
  implicit none
  interface 
    subroutine source_ang_mom(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_ang_mom
  end interface
end module interface_source_ang_mom
