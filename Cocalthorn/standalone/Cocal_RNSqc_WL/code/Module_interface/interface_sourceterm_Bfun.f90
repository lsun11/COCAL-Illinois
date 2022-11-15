module interface_sourceterm_Bfun
  implicit none
  interface 
    subroutine sourceterm_Bfun(sou,Bfun,potx,poty,potz)
      real(8), pointer :: sou(:,:,:), Bfun(:,:,:), potx(:,:,:), poty(:,:,:), potz(:,:,:)
    end subroutine sourceterm_Bfun
  end interface
end module interface_sourceterm_Bfun
