module interface_source_mp_minus_madm
  implicit none
  interface 
    subroutine source_mp_minus_madm(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_mp_minus_madm
  end interface
end module interface_source_mp_minus_madm
