module interface_index_vec_down2up
  implicit none
  interface 
    subroutine index_vec_down2up(vecxu,vecyu,veczu,vecxd,vecyd,veczd)
      real(8), pointer :: vecxu(:,:,:), vecyu(:,:,:), veczu(:,:,:), &
      &                   vecxd(:,:,:), vecyd(:,:,:), veczd(:,:,:)
    end subroutine index_vec_down2up
  end interface
end module interface_index_vec_down2up
