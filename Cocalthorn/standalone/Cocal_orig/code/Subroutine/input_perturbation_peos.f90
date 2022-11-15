subroutine input_perturbation_peos
  use phys_constant, only : long, pi
  use coordinate_grav_r, only : rg
  use grid_parameter, only : nrf, ntf, npf
  use trigonometry_grav_theta, only : sinthg , costhg
  use trigonometry_grav_phi,   only : sinphig, cosphig, cosmpg
  use def_matter, only : rs, emd
  use def_matter_parameter, only : ome, radi
  use make_array_3d
  implicit none
  integer    :: irf, itf, ipf, iper
  real(long) :: xx, yy, zz, ww, rr, phia, ampl, aa
  real(long) :: qq, hh, pre, rho, ene, abin, abct
  real(long), pointer :: rhopf(:,:,:)
!
! Assume rho2 = rho*aa and  q = pre/rho with
!        pre  = kappa*rho**gamma. Then  
!        pre2 = kappa*rho2**gamma = pre*aa**gamma which implies 
!        q2   = q*aa**(gamma-1)
!
  call alloc_array3d(rhopf, 0, nrf, 0, ntf, 0, npf)

! Amplitude and type of perturbation
  ampl=0.2d0;    iper=4
!  ampl=0.15d0;    iper=1

  do irf = 0, nrf
    do itf = 0, ntf
      do ipf = 0, npf
        xx = rg(irf)*rs(itf,ipf)*sinthg(itf)*cosphig(ipf)
        yy = rg(irf)*rs(itf,ipf)*sinthg(itf)*sinphig(ipf)
        zz = rg(irf)*rs(itf,ipf)*costhg(itf)
        rr = xx*xx + yy*yy
        ww = xx*xx - yy*yy
        if (rr .ne. 0)  then
          phia = atan2(yy,xx)
          if (yy<0)  phia = 2.0d0*pi + phia
        end if 

        qq = emd(irf,itf,ipf)
!        call peos_q2hprhoab(qq, hh, pre, rho, ene, abin, abct)
        call peos_q2hprho(qq, hh, pre, rho, ene)
 
        if (iper==1)  aa = 1.0d0 - ampl + ampl*cosmpg(2,ipf)
        if (iper==2)  aa = 1.0d0 - ampl + ampl*cos(ome*phia/radi)
        if (iper==4)  aa = 1.0d0 + ampl*ww

        rhopf(irf,itf,ipf) = rho*aa
!        emd(irf,itf,ipf) = qq*aa**(abin-1.0d0)
        emd(irf,itf,ipf) = qq*aa
      end do
    end do
  end do
!
  irf=ntf/2;    ipf=npf/2 
  do irf = nrf, 0, -1
    write(6,'(a13,1p,3e25.15)') "x, rho, emd =", rg(irf)*rs(itf,ipf), rhopf(irf,itf,ipf), emd(irf,itf,ipf)
  end do
  irf=ntf/2;    ipf=0  
  do irf = 0, nrf
    write(6,'(a13,1p,3e25.15)') "x, rho, emd =", rg(irf)*rs(itf,ipf), rhopf(irf,itf,ipf), emd(irf,itf,ipf)
  end do
  irf=ntf/2;    ipf=3*npf/4 
  do irf = nrf, 0, -1
    write(6,'(a13,1p,3e25.15)') "y, rho, emd =", rg(irf)*rs(itf,ipf), rhopf(irf,itf,ipf), emd(irf,itf,ipf)
  end do
  irf=ntf/2;    ipf=npf/4  
  do irf = 0, nrf
    write(6,'(a13,1p,3e25.15)') "y, rho, emd =", rg(irf)*rs(itf,ipf), rhopf(irf,itf,ipf), emd(irf,itf,ipf)
  end do
!
  deallocate(rhopf)
end subroutine input_perturbation_peos
