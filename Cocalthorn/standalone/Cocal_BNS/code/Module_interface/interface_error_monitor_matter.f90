module interface_error_monitor_matter
  implicit none
  interface 
    subroutine error_monitor_matter(pot,pot_bak,char,irmo,itmo,ipmo)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      character(len=5)     :: char
      integer              :: irmo, itmo, ipmo
    end subroutine error_monitor_matter
  end interface
end module interface_error_monitor_matter
