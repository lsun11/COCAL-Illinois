module interface_source_mp_minus_madm_peos
  implicit none
  interface 
    subroutine source_mp_minus_madm_peos(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_mp_minus_madm_peos
  end interface
end module interface_source_mp_minus_madm_peos
