module interface_length_conversion
  implicit none
  interface 
    subroutine length_conversion(i,inlen,outlen)
      real(8) :: inlen, outlen
      integer :: i
    end subroutine length_conversion
  end interface
end module interface_length_conversion
