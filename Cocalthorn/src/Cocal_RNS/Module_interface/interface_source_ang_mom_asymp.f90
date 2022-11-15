module interface_source_ang_mom_asymp
  implicit none
  interface 
    subroutine source_ang_mom_asymp(sousf,irg)
      real(8), pointer     :: sousf(:,:)
      integer, intent(in)  :: irg
    end subroutine source_ang_mom_asymp
  end interface
end module interface_source_ang_mom_asymp
