module interface_grgrad_vec_covariant
  implicit none
  interface 
    subroutine grgrad_vec_covariant(index,fnx,fny,fnz,pdfnx,pdfny,pdfnz, &
               &                                      cdfnx,cdfny,cdfnz)
      real(8), pointer :: fnx(:,:,:), fny(:,:,:), fnz(:,:,:)
      real(8), pointer :: pdfnx(:,:,:,:), pdfny(:,:,:,:), pdfnz(:,:,:,:)
      real(8), pointer :: cdfnx(:,:,:,:), cdfny(:,:,:,:), cdfnz(:,:,:,:)
      character(len=1) :: index
    end subroutine grgrad_vec_covariant
  end interface
end module interface_grgrad_vec_covariant
