module interface_error_monitor_surface
  implicit none
  interface 
    subroutine error_monitor_surface(pot,pot_bak,char,itmo,ipmo)
      real(8), pointer     :: pot(:,:), pot_bak(:,:)
      character(len=5)     :: char
      integer              :: itmo, ipmo
    end subroutine error_monitor_surface
  end interface
end module interface_error_monitor_surface
