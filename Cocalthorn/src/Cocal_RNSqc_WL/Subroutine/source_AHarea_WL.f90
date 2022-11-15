subroutine source_AHarea_WL(sousg)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_horizon, only : ahz
  use def_bh_parameter
!  use def_kerr_schild, only : kerr_a, reh_ks
  use def_metric, only  : psi
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use coordinate_grav_r, only : rg  
  use trigonometry_grav_theta, only : hsinthg, hcosthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use coordinate_grav_theta,   only : dthg
  use interface_interpo_linear_surface_type0
  use interface_interpo_linear_type0_2Dsurf
  use interface_grdtp_2Dsurf_midpoint_type0
  use make_array_2d
  implicit none
  real(long), pointer  :: sousg(:,:), jacob(:,:), gamij(:,:), ps(:,:)
  real(long), pointer  :: hxx(:,:), hxy(:,:), hxz(:,:), hyy(:,:), hyz(:,:), hzz(:,:)
  real(long) :: psim, psim4, gmxxd, gmxyd, gmxzd, gmyyd, gmyzd, gmzzd
  real(long) :: gam_thth, gam_phph, gam_thph, x1, x2, f1, f2, r2pa2 
  real(long) :: rc, drcdth, drcdph, a2, r2, dcosangdth, dcosangdph
  real(long) :: sinph12, cosph12, cosang
  integer :: irg, itg, ipg, ia,ib, ieh
!
! Computes the area element dS of a two dimensional surface, ahz(itg,ipg), in the WL formalism.
! dS = sqrt(g_thth * g_phph - g_thph**2) dth dph
!
!
  call alloc_array2d(jacob,1,3,1,3)
  call alloc_array2d(gamij,1,3,1,3)
  call alloc_array2d(ps, 0,ntg,0,npg)
  call alloc_array2d(hxx,0,ntg,0,npg)
  call alloc_array2d(hxy,0,ntg,0,npg)
  call alloc_array2d(hxz,0,ntg,0,npg)
  call alloc_array2d(hyy,0,ntg,0,npg)
  call alloc_array2d(hyz,0,ntg,0,npg)
  call alloc_array2d(hzz,0,ntg,0,npg)
!

  do ipg = 0, npg
    do itg = 0, ntg

      if (rg(0)>ahz(itg,ipg)) then
!        write(6,*) " *** Warning *** : AH inside excised region: (itg,ipg)=", itg,ipg, &
!           &       "     r_ah=", ahz(itg,ipg)
        ieh=0
      else
        do irg=0,nrg
          if (rg(irg)<=ahz(itg,ipg) .and. ahz(itg,ipg)<rg(irg+1)) then
            ieh=irg
            exit
          end if
        end do
        if (ipg==0 .and. itg==0) write(6,*)  "(itg,ipg)=(0,0) : ahz(itg,ipg), ieh, rg(irg), rg(irg+1)", &
                         &                   ahz(itg,ipg), ieh, rg(ieh),rg(ieh+1)
      end if

      x1=rg(ieh);           x2=rg(ieh+1);      
      f1= psi(ieh,itg,ipg); f2= psi(ieh+1,itg,ipg);     ps(itg,ipg)=f2-(x2-ahz(itg,ipg))*(f2-f1)/(x2-x1);
      f1=hxxd(ieh,itg,ipg); f2=hxxd(ieh+1,itg,ipg);    hxx(itg,ipg)=f2-(x2-ahz(itg,ipg))*(f2-f1)/(x2-x1);
      f1=hxyd(ieh,itg,ipg); f2=hxyd(ieh+1,itg,ipg);    hxy(itg,ipg)=f2-(x2-ahz(itg,ipg))*(f2-f1)/(x2-x1);
      f1=hxzd(ieh,itg,ipg); f2=hxzd(ieh+1,itg,ipg);    hxz(itg,ipg)=f2-(x2-ahz(itg,ipg))*(f2-f1)/(x2-x1);
      f1=hyyd(ieh,itg,ipg); f2=hyyd(ieh+1,itg,ipg);    hyy(itg,ipg)=f2-(x2-ahz(itg,ipg))*(f2-f1)/(x2-x1);
      f1=hyzd(ieh,itg,ipg); f2=hyzd(ieh+1,itg,ipg);    hyz(itg,ipg)=f2-(x2-ahz(itg,ipg))*(f2-f1)/(x2-x1);
      f1=hzzd(ieh,itg,ipg); f2=hzzd(ieh+1,itg,ipg);    hzz(itg,ipg)=f2-(x2-ahz(itg,ipg))*(f2-f1)/(x2-x1);
    end do
  end do
!  write(6,*) hxx(10,10)
!  write(6,*) hxxd(irg,10,10), hxxd(irg+1,10,10)

  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(rc,   ahz,itg,ipg)
      call interpo_linear_type0_2Dsurf(psim,  ps,itg,ipg)
      call interpo_linear_type0_2Dsurf(gmxxd,hxx,itg,ipg)
      call interpo_linear_type0_2Dsurf(gmxyd,hxy,itg,ipg)
      call interpo_linear_type0_2Dsurf(gmxzd,hxz,itg,ipg)
      call interpo_linear_type0_2Dsurf(gmyyd,hyy,itg,ipg)
      call interpo_linear_type0_2Dsurf(gmyzd,hyz,itg,ipg)
      call interpo_linear_type0_2Dsurf(gmzzd,hzz,itg,ipg)

      psim4 = psim*psim*psim*psim

      gamij(1,1) = (1.0d0+gmxxd)*psim4
      gamij(1,2) = (      gmxyd)*psim4
      gamij(1,3) = (      gmxzd)*psim4
      gamij(2,2) = (1.0d0+gmyyd)*psim4
      gamij(2,3) = (      gmyzd)*psim4
      gamij(3,3) = (1.0d0+gmzzd)*psim4

      gamij(2,1) = gamij(1,2)
      gamij(3,1) = gamij(1,3)
      gamij(3,2) = gamij(2,3)


      call grdtp_2Dsurf_midpoint_type0(ahz,drcdth,drcdph,itg,ipg)
 
      jacob(1,1) = drcdth*hsinthg(itg)*hcosphig(ipg) + rc*hcosthg(itg)*hcosphig(ipg)
      jacob(1,2) = drcdph*hsinthg(itg)*hcosphig(ipg) - rc*hsinthg(itg)*hsinphig(ipg)
!
      jacob(2,1) = drcdth*hsinthg(itg)*hsinphig(ipg) + rc*hcosthg(itg)*hsinphig(ipg)
      jacob(2,2) = drcdph*hsinthg(itg)*hsinphig(ipg) + rc*hsinthg(itg)*hcosphig(ipg)
!
      jacob(3,1) = drcdth*hcosthg(itg) - rc*hsinthg(itg)
      jacob(3,2) = drcdph*hcosthg(itg)
     
      gam_thth = 0.0d0
      gam_phph = 0.0d0
      gam_thph = 0.0d0
      do ia=1,3
        do ib=1,3
          gam_thth = gam_thth + jacob(ia,1)*gamij(ia,ib)*jacob(ib,1)
          gam_phph = gam_phph + jacob(ia,2)*gamij(ia,ib)*jacob(ib,2)
          gam_thph = gam_thph + jacob(ia,1)*gamij(ia,ib)*jacob(ib,2)
        end do
      end do

!     NOTE: The source here contains the typical weight r**2 sin(theta)
      sousg(itg,ipg) = dsqrt(gam_thth*gam_phph - gam_thph*gam_thph)
    end do
  end do
!

  deallocate(jacob)
  deallocate(gamij)
  deallocate(ps)
  deallocate(hxx)
  deallocate(hxy)
  deallocate(hxz)
  deallocate(hyy)
  deallocate(hyz)
  deallocate(hzz)
!
end subroutine source_AHarea_WL
