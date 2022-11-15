!  Associated Legendre function and factorials
!______________________________________________
module legendre_fn_grav
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  plmg(:,:,:), hplmg(:,:,:)
  real(long), pointer  ::  facnmg(:,:), epsig(:) 
contains
subroutine allocate_legendre
  use phys_constant,  only : long
  use grid_parameter, only : ntg, nlg     ! nlg : maximum value of multipole
  use make_array_3d
  use make_array_2d
  use make_array_1d
  implicit none
!
! --- allocate arrays for legendre expansion
!
  call alloc_array3d(plmg, 0, nlg, 0, nlg, 0, ntg)
  call alloc_array3d(hplmg, 0, nlg, 0, nlg, 1, ntg)
  call alloc_array2d(facnmg, 0, nlg, 0, nlg)
  call alloc_array1d(epsig, 0, nlg)
!
end subroutine allocate_legendre
SUBROUTINE legendre
  use phys_constant,  only : nnlg, long
  use grid_parameter, only : ntg, nlg     ! nlg : maximum value of multipole
  use coordinate_grav_theta, only : thg, hthg
  use make_array_3d
  use make_array_2d
  use make_array_1d
  implicit none
  integer              ::  it, mm, nn
  real(long)           ::  fmm, fnn, fnmfm
  real(long)           ::  pnag(0:nnlg,0:nnlg), theta
!
! --- computation of associated Legendre polynomials.
!
  do it = 0, ntg
    theta = thg(it)
    call legendre_theta(pnag,theta)
    plmg(0:nlg,0:nlg,it) = pnag(0:nlg,0:nlg)
  end do
  do it = 1, ntg
    theta = hthg(it)
    call legendre_theta(pnag,theta)
    hplmg(0:nlg,0:nlg,it) = pnag(0:nlg,0:nlg)
  end do
! 
! --- Computation of factors, facnmg & epsig. ---
!
  do nn = 0, nlg
    do mm = 0, nlg
      facnmg(nn,mm) = 0.0d+0
    end do
  end do
!
  facnmg(0,0) = 1.0d0
  do nn = 1, nlg
    fnn  = real(nn)
    facnmg(nn,0) = 1.0d0
    do mm = 1, nn
      fmm  = real(mm)
      fnmfm= fnn-fmm + 1.0d0 
      facnmg(nn,mm) = facnmg(nn,mm-1)/(fnn+fmm)/fnmfm
    end do
  end do
!
  do mm = 0, nlg 
    epsig(mm) = 2.0d+0 
    if (mm == 0) epsig(mm) = 1.0d+0 
  end do
!   
end subroutine legendre
SUBROUTINE legendre_theta(pnag,theta)
  use phys_constant,  only : nnlg, long
  use grid_parameter, only : nlg
  implicit none
  integer     ::  it, mm, nn, kk
  real(long)  ::  fmm, fkk, cc, ss, q1, q2
  real(long)  ::  pnag(0:nnlg,0:nnlg), theta, fac(0:nnlg)
!
  cc = cos(theta)
  ss = sin(theta)
! 
  pnag(0:nlg,0:nlg)  = 0.0d0
!      
  fac(0) = 1.0d0
  do mm = 1, nlg
    fmm = real(mm)
    fac(mm) = (2.0d0*fmm-1.0d0) * fac(mm-1)
  end do
! 
  pnag(0,0) = 1.d0
  do mm = 1, nlg
    pnag(mm,mm) = fac(mm) * (-ss)**mm 
  end do
!
  do mm = 0, nlg-1
    fmm = real(mm)
    pnag(mm+1,mm) = (2.0d0*fmm + 1.0d0)* cc*  pnag(mm,mm)
  end do
! 
  do mm = 0, nlg-2
    fmm = real(mm)
    do kk = 2, nlg-mm
      fkk = real(kk)
      q1 = ( 2.0d0 * fmm + 2.0d0 * fkk - 1.0d0 ) / fkk
      q2 = ( 2.0d0 * fmm + fkk - 1.0d0 ) / fkk
      pnag(mm+kk,mm) = q1 * cc * pnag(mm+kk-1,mm)- q2 * pnag(mm+kk-2,mm)
    end do
  end do
end subroutine legendre_theta
end module legendre_fn_grav
