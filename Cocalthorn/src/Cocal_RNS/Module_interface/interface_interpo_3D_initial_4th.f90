module interface_interpo_3D_initial_4th
  implicit none
  interface 
    subroutine interpo_3D_initial_4th(fnc,char_co)
      real(8), pointer :: fnc(:,:,:)
      character(len=4) :: char_co 
    end subroutine interpo_3D_initial_4th
  end interface
end module interface_interpo_3D_initial_4th
