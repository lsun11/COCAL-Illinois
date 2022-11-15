module interface_grgrad_vec_partial
  implicit none
  interface 
    subroutine grgrad_vec_partial(fnx,fny,fnz,pdfnx,pdfny,pdfnz)
      real(8), pointer :: fnx(:,:,:), fny(:,:,:), fnz(:,:,:)
      real(8), pointer :: pdfnx(:,:,:,:), pdfny(:,:,:,:), pdfnz(:,:,:,:)
    end subroutine grgrad_vec_partial
  end interface
end module interface_grgrad_vec_partial
