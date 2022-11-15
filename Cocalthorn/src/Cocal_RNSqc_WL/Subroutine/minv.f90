subroutine minv(aa,bb,nn,nnz)
  use phys_constant, only  :   long
  implicit none
  integer      ::    nn, nnz
  real(long)   ::    aa(nnz,nnz),bb(nnz)
  real(long)   ::    fmax, aab, sig
  integer      ::    nnp1, nnm1
  integer      ::    i, j, k, kp1, ifmax
!
! solution of sm*aa=ss
  nnp1 = nn + 1
  aa(1:nn,nnp1) = bb(1:nn)
!
  nnm1 = nn - 1  
  do k = 1, nnm1
    kp1 = k + 1
!
    fmax = 1.e-10
    ifmax = 1
!
    do i = k, nn
      aab = abs(aa(i,k))
      if(aab > fmax) then
         fmax = aab
         ifmax = i
      end if
    end do
    if(ifmax /= k) then
      sig = 1.0d0
      if(aa(k,k)*aa(ifmax,k) < 0.) sig = -1.0d0
      aa(k,1:nnp1) = aa(k,1:nnp1) + sig*aa(ifmax,1:nnp1)
    end if
!
    aa(k,kp1:nnp1) = aa(k,kp1:nnp1)/aa(k,k)
! what the previous part did is to make aa(k,k) = 1.      
    do i = kp1, nnp1
      aa(kp1:nn, i) = aa(kp1:nn, i) - aa(kp1:nn, k)*aa(k, i)
    end do
! now make aa(kp1:nn,k) = 0. aa(n,n) is still unchanged.      
!
  end do
!
  aa(nn,nnp1) = aa(nn,nnp1)/aa(nn,nn)
  do i = 1, nnm1
    k   = nn - i
    kp1 = k + 1
    do j = kp1, nn
      aa(k,nnp1) = aa(k,nnp1) - aa(k,j)*aa(j,nnp1)
    end do
  end do
! Finish the inverse remember all diagonal of a is 1.      
!
  bb (1:nn) = aa(1:nn,nnp1)
!
end subroutine minv
