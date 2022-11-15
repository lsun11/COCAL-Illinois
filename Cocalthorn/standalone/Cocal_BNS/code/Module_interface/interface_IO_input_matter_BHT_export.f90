module interface_IO_input_matter_BHT_export
  implicit none
  interface 
    subroutine IO_input_matter_BHT_export(filenm,emdg,omeg,ome,ber,radi)
      real(8), pointer :: emdg(:,:,:), omeg(:,:,:)
      real(8) ::  ome,ber,radi
      character(len=*) :: filenm
    end subroutine IO_input_matter_BHT_export
  end interface
end module interface_IO_input_matter_BHT_export
