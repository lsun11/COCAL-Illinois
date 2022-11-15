module interface_minv
  implicit none
  interface 
    subroutine minv(aa,bb,nn,nnz)
      integer :: nn, nnz
      real(8) :: aa(nnz,nnz),bb(nnz)
    end subroutine minv
  end interface
end module interface_minv
