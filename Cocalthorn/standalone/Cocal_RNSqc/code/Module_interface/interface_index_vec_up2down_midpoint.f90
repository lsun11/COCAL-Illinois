module interface_index_vec_up2down_midpoint
  implicit none
  interface 
    subroutine index_vec_up2down_midpoint(vecxd,vecyd,veczd,vecxu,vecyu,veczu)
      real(8), pointer :: vecxu(:,:,:), vecyu(:,:,:), veczu(:,:,:), &
      &                   vecxd(:,:,:), vecyd(:,:,:), veczd(:,:,:)
    end subroutine index_vec_up2down_midpoint
  end interface
end module interface_index_vec_up2down_midpoint
