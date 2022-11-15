module interface_index_vec_up2down
  implicit none
  interface 
    subroutine index_vec_up2down(vecxd,vecyd,veczd,vecxu,vecyu,veczu)
      real(8), pointer :: vecxu(:,:,:), vecyu(:,:,:), veczu(:,:,:), &
      &                   vecxd(:,:,:), vecyd(:,:,:), veczd(:,:,:)
    end subroutine index_vec_up2down
  end interface
end module interface_index_vec_up2down
