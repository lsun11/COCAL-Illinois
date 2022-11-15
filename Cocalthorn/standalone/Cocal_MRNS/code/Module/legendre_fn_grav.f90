!  Associated Legendre function and factorials
!______________________________________________
module legendre_fn_grav
  use phys_constant, only : long
  implicit none
  real(long), pointer ::    plmg(:,:,:),    hplmg(:,:,:)
  real(long), pointer ::  dtplmg(:,:,:),  hdtplmg(:,:,:)
  real(long), pointer ::   yplmg(:,:,:),   hyplmg(:,:,:)
  real(long), pointer :: dtyplmg(:,:,:), hdtyplmg(:,:,:)
  real(long), pointer :: facnmg(:,:), epsig(:) 
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
  call alloc_array3d(dtplmg, 0, nlg, 0, nlg, 0, ntg)
  call alloc_array3d(hdtplmg, 0, nlg, 0, nlg, 1, ntg)
  call alloc_array3d(yplmg, 0, nlg, 0, nlg, 0, ntg)
  call alloc_array3d(hyplmg, 0, nlg, 0, nlg, 1, ntg)
  call alloc_array3d(dtyplmg, 0, nlg, 0, nlg, 0, ntg)
  call alloc_array3d(hdtyplmg, 0, nlg, 0, nlg, 1, ntg)
  call alloc_array2d(facnmg, 0, nlg, 0, nlg)
  call alloc_array1d(epsig, 0, nlg)
!
end subroutine allocate_legendre

SUBROUTINE legendre
  use phys_constant,  only : nnlg, long, pi
  use grid_parameter, only : ntg, nlg     ! nlg : maximum value of multipole
  use coordinate_grav_theta, only : thg, hthg
  use make_array_3d
  use make_array_2d
  use make_array_1d
  implicit none
  integer              ::  it, mm, nn
  real(long)           ::  fmm, fnn, fnmfm, yplm_fac, theta
  real(long)           ::  pnag(0:nnlg,0:nnlg), dtpnag(0:nnlg,0:nnlg)
!
! --- computation of associated Legendre polynomials.
!
!  call alloc_array3d(plmg, 0, nlg, 0, nlg, 0, ntg)
!  call alloc_array3d(hplmg, 0, nlg, 0, nlg, 1, ntg)
!  call alloc_array3d(dtplmg, 0, nlg, 0, nlg, 0, ntg)
!  call alloc_array3d(hdtplmg, 0, nlg, 0, nlg, 1, ntg)
!  call alloc_array3d(yplmg, 0, nlg, 0, nlg, 0, ntg)
!  call alloc_array3d(hyplmg, 0, nlg, 0, nlg, 1, ntg)
!  call alloc_array3d(dtyplmg, 0, nlg, 0, nlg, 0, ntg)
!  call alloc_array3d(hdtyplmg, 0, nlg, 0, nlg, 1, ntg)
  do it = 0, ntg
    theta = thg(it)
    call legendre_theta(pnag,dtpnag,theta)
    plmg(0:nlg,0:nlg,it) = pnag(0:nlg,0:nlg)
    dtplmg(0:nlg,0:nlg,it) = dtpnag(0:nlg,0:nlg)
  end do
  do it = 1, ntg
    theta = hthg(it)
    call legendre_theta(pnag,dtpnag,theta)
    hplmg(0:nlg,0:nlg,it) = pnag(0:nlg,0:nlg)
    hdtplmg(0:nlg,0:nlg,it) = dtpnag(0:nlg,0:nlg)
  end do
! 
! --- Computation of factors, facnmg & epsig. ---
!
!  call alloc_array2d(facnmg, 0, nlg, 0, nlg)
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
!  call alloc_array1d(epsig, 0, nlg)
  do mm = 0, nlg 
    epsig(mm) = 2.0d+0 
    if (mm == 0) epsig(mm) = 1.0d+0 
  end do
!   
  do nn = 0, nlg
    do mm = 0, nn
      do it = 0, ntg
        yplm_fac = (epsig(mm)*(2.0d0*dble(nn)+1.0d0)/(4.0d0*pi)* &
        &          facnmg(nn,mm))**0.5d0
           yplmg(nn,mm,it) = yplm_fac*   plmg(nn,mm,it)
         dtyplmg(nn,mm,it) = yplm_fac* dtplmg(nn,mm,it)
        if (it.eq.0) cycle
          hyplmg(nn,mm,it) = yplm_fac*  hplmg(nn,mm,it)
        hdtyplmg(nn,mm,it) = yplm_fac*hdtplmg(nn,mm,it)
      end do
    end do
  end do
!
end subroutine legendre
SUBROUTINE legendre_theta(pnag,dtpnag,theta)
  use phys_constant,  only : nnlg, long
  use grid_parameter, only : nlg
  implicit none
  integer     ::  it, mm, nn, kk
  real(long)  ::  fmm, fkk, cc, ss, q1, q2
  real(long)  ::  pnag(0:nnlg,0:nnlg), dtpnag(0:nnlg,0:nnlg)
  real(long)  ::  theta, fac(0:nnlg)
!
  cc = cos(theta)
  ss = sin(theta)
! 
  pnag(0:nlg,0:nlg)  = 0.0d0
  dtpnag(0:nlg,0:nlg)  = 0.0d0
!      
  fac(0) = 1.0d0
  do mm = 1, nlg
    fmm = real(mm)
    fac(mm) = (2.0d0*fmm-1.0d0) * fac(mm-1)
  end do
! 
  pnag(0,0) = 1.d0
  dtpnag(0,0) = 0.d0
  do mm = 1, nlg
    pnag(mm,mm) = fac(mm) * (-ss)**mm 
  end do
!
  mm = 1
  dtpnag(mm,mm) = fac(mm) * dble(mm)*(-cc)
  do mm = 2, nlg
    dtpnag(mm,mm) = fac(mm) * dble(mm)*(-ss)**(mm-1)*(-cc)
  end do
!
  do mm = 0, nlg-1
    fmm = real(mm)
      pnag(mm+1,mm) = (2.0d0*fmm + 1.0d0)* cc*  pnag(mm,mm)
    dtpnag(mm+1,mm) = (2.0d0*fmm + 1.0d0)*(cc*dtpnag(mm,mm) &
    &                                     -ss*  pnag(mm,mm))
  end do
! 
  do mm = 0, nlg-2
    fmm = real(mm)
    do kk = 2, nlg-mm
      fkk = real(kk)
      q1 = ( 2.0d0 * fmm + 2.0d0 * fkk - 1.0d0 ) / fkk
      q2 = ( 2.0d0 * fmm + fkk - 1.0d0 ) / fkk
        pnag(mm+kk,mm) = q1* cc*  pnag(mm+kk-1,mm) - q2*  pnag(mm+kk-2,mm)
      dtpnag(mm+kk,mm) = q1*(cc*dtpnag(mm+kk-1,mm) - ss * pnag(mm+kk-1,mm)) &
     &                                             - q2*dtpnag(mm+kk-2,mm)
    end do
  end do
end subroutine legendre_theta
end module legendre_fn_grav
