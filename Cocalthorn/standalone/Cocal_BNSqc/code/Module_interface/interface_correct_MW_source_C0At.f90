module interface_correct_MW_source_C0At
  implicit none
  interface 
    subroutine correct_MW_source_C0At(sou,intp)
      real(8), pointer :: sou(:,:,:)
      integer :: intp
    end subroutine correct_MW_source_C0At
  end interface
end module interface_correct_MW_source_C0At
