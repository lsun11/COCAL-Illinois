module interface_correct_matter_source_midpoint
  implicit none
  interface 
    subroutine correct_matter_source_midpoint(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine correct_matter_source_midpoint
  end interface
end module interface_correct_matter_source_midpoint
