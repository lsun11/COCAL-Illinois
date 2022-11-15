function dfdx_nth(no,x,f,v)
  use phys_constant, only : long
  implicit none
  integer, intent(in) :: no
  real(long), intent(in) :: v
  real(long), intent(in) :: x(no+1),f(no+1)
  real(long) :: dfdx_nth
  real(long) :: deno, nume, prodxv
  real(long) :: dx(0:no,0:no), xv(0:no), wex(0:no)
  integer :: ir(0:no), iset(no), comb(no-1)
  integer :: i,j,k
!
  do i=0,no
    ir(i) = i+1
  end do

  do i=0,no
    do j= i+1,no
      dx(i,j) = x(ir(i)) - x(ir(j))
      dx(j,i) =-dx(i,j)
    end do
  end do

  do i=0,no
    xv(i) = v - x(ir(i))
  end do
!
  do i=0,no
    deno = 1.0d0;  k=1
    do j=0,no
      if (i.ne.j) then
        deno = deno*dx(i,j)
        iset(k)=j;  k=k+1
      end if
    end do
    nume=0.0d0
    comb(1)=1
    call combinations(1,2,1)
!   iterate(1,n-r+1,1) : combinations of r over n elements
!   n-r+1=no - (no-1) + 1 = 2
!
!   write(6,*)  "nume=", nume
    wex(i) = nume/deno
  end do

  dfdx_nth =0.0d0
  do i=0,no
    dfdx_nth = dfdx_nth + wex(i)*f(ir(i))
  end do

contains 

RECURSIVE subroutine combinations(is,ie,jj) 
  integer, intent(in) :: is,ie
  integer :: ii,jj,m
  do ii=is,ie
    comb(jj)=ii
    if(jj.lt.(no-1))    call combinations(ii+1,ie+1,jj+1)
    if(jj.eq.(no-1)) then
!      WRITE(6,'(10i3)')  (iset(comb(m)), m=1,no-1)
      prodxv=1.0d0
      do m=1,no-1
        prodxv = prodxv*xv(iset(comb(m)))
      end do
      nume = nume + prodxv
!      write(6,'(1p,e23.15)') nume
    end if
  end do
END subroutine combinations
end function dfdx_nth

