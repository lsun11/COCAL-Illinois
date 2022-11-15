module interface_source_EMenergy_axisym_WL
  implicit none
  interface 
    subroutine source_EMenergy_axisym_WL(sou_MtorB,sou_MpolB,sou_MeleE)
      real(8), pointer :: sou_MtorB(:,:,:), sou_MpolB(:,:,:), sou_MeleE(:,:,:)
    end subroutine source_EMenergy_axisym_WL
  end interface
end module interface_source_EMenergy_axisym_WL
