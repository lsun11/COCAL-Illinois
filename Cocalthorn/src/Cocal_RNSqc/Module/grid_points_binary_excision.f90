module grid_points_binary_excision
  use phys_constant, only : long
  implicit none
  integer, pointer :: irg_exin(:,:),  itg_exin(:,:),  ipg_exin(:,:)
  integer, pointer :: irg_exout(:,:), itg_exout(:,:), ipg_exout(:,:)
  real(long), pointer :: rg_exin(:,:),  thg_exin(:,:),  phig_exin(:,:)
  real(long), pointer :: rg_exout(:,:), thg_exout(:,:), phig_exout(:,:)
  integer, pointer :: ihrg_exin(:,:),  ihtg_exin(:,:),  ihpg_exin(:,:)
  integer, pointer :: ihrg_exout(:,:), ihtg_exout(:,:), ihpg_exout(:,:)
  real(long), pointer :: hrg_exin(:,:),  hthg_exin(:,:),  hphig_exin(:,:)
  real(long), pointer :: hrg_exout(:,:), hthg_exout(:,:), hphig_exout(:,:)
!
  real(long), pointer :: rb(:,:,:), thb(:,:,:), phib(:,:,:)
  real(long), pointer :: hrb(:,:,:), hthb(:,:,:), hphib(:,:,:)
  integer :: ntg_exin_min, ntg_exout_max, npg_exin_max, npg_exout_min
!
contains
!
subroutine allocate_grid_points_binary_excision
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use make_array_2d
  use make_array_3d
  use make_int_array_2d
  implicit none
!
  call alloc_int_array2d(irg_exin, 0,ntg,0,npg)
  call alloc_int_array2d(irg_exout,0,ntg,0,npg)
  call alloc_int_array2d(itg_exin, 0,nrg,0,npg)
  call alloc_int_array2d(itg_exout,0,nrg,0,npg)
  call alloc_int_array2d(ipg_exin, 0,nrg,0,ntg)
  call alloc_int_array2d(ipg_exout,0,nrg,0,ntg)
  call alloc_array2d(rg_exin, 0,ntg,0,npg)
  call alloc_array2d(rg_exout,0,ntg,0,npg)
  call alloc_array2d(thg_exin, 0,nrg,0,npg)
  call alloc_array2d(thg_exout,0,nrg,0,npg)
  call alloc_array2d(phig_exin, 0,nrg,0,ntg)
  call alloc_array2d(phig_exout,0,nrg,0,ntg)
!
  call alloc_int_array2d(ihrg_exin, 0,ntg,0,npg)
  call alloc_int_array2d(ihrg_exout,0,ntg,0,npg)
  call alloc_int_array2d(ihtg_exin, 0,nrg,0,npg)
  call alloc_int_array2d(ihtg_exout,0,nrg,0,npg)
  call alloc_int_array2d(ihpg_exin, 0,nrg,0,ntg)
  call alloc_int_array2d(ihpg_exout,0,nrg,0,ntg)
  call alloc_array2d(hrg_exin, 0,ntg,0,npg)
  call alloc_array2d(hrg_exout,0,ntg,0,npg)
  call alloc_array2d(hthg_exin, 0,nrg,0,npg)
  call alloc_array2d(hthg_exout,0,nrg,0,npg)
  call alloc_array2d(hphig_exin, 0,nrg,0,ntg)
  call alloc_array2d(hphig_exout,0,nrg,0,ntg)
!
  call alloc_array3d(rb,0,nrg,0,ntg,0,npg)
  call alloc_array3d(thb,0,nrg,0,ntg,0,npg)
  call alloc_array3d(phib,0,nrg,0,ntg,0,npg)
  call alloc_array3d(hrb,1,nrg,1,ntg,1,npg)
  call alloc_array3d(hthb,1,nrg,1,ntg,1,npg)
  call alloc_array3d(hphib,1,nrg,1,ntg,1,npg)
!
  irg_exin( 0:ntg,0:npg) = 0
  irg_exout(0:ntg,0:npg) = 0
  itg_exin( 0:nrg,0:npg) = 0
  itg_exout(0:nrg,0:npg) = 0
  ipg_exin( 0:nrg,0:ntg) = 0
  ipg_exout(0:nrg,0:ntg) = 0
  rg_exin( 0:ntg,0:npg) = 0.0d0
  rg_exout(0:ntg,0:npg) = 0.0d0
  thg_exin( 0:nrg,0:npg) = 0.0d0
  thg_exout(0:nrg,0:npg) = 0.0d0
  phig_exin( 0:nrg,0:ntg) = 0.0d0
  phig_exout(0:nrg,0:ntg) = 0.0d0
!
  ihrg_exin( 0:ntg,0:npg) = 0
  ihrg_exout(0:ntg,0:npg) = 0
  ihtg_exin( 0:nrg,0:npg) = 0
  ihtg_exout(0:nrg,0:npg) = 0
  ihpg_exin( 0:nrg,0:ntg) = 0
  ihpg_exout(0:nrg,0:ntg) = 0
  hrg_exin( 0:ntg,0:npg) = 0.0d0
  hrg_exout(0:ntg,0:npg) = 0.0d0
  hthg_exin( 0:nrg,0:npg) = 0.0d0
  hthg_exout(0:nrg,0:npg) = 0.0d0
  hphig_exin( 0:nrg,0:ntg) = 0.0d0
  hphig_exout(0:nrg,0:ntg) = 0.0d0
!
  rb(0:nrg,0:ntg,0:npg) = 0.0d0
  thb(0:nrg,0:ntg,0:npg) = 0.0d0
  phib(0:nrg,0:ntg,0:npg) = 0.0d0
  hrb(1:nrg,1:ntg,1:npg) = 0.0d0
  hthb(1:nrg,1:ntg,1:npg) = 0.0d0
  hphib(1:nrg,1:ntg,1:npg) = 0.0d0
!
end subroutine allocate_grid_points_binary_excision
!
! --- Determin coordinates and grid number at the excision radius.
!
subroutine calc_grid_points_binary_excision
  use phys_constant, only : long, pi
  use coordinate_grav_r, only : rg, hrg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig, dphig
  use grid_parameter, only : nrg, ntg, npg, ntgxy, npgxzp, npgxzm
  use def_binary_parameter, only : sepa
  use grid_parameter_binary_excision, only : ex_radius
  use trigonometry_grav_theta, only :  sinthg, costhg, hsinthg, hcosthg
  use trigonometry_grav_phi, only :  sinphig, cosphig, hsinphig, hcosphig
  implicit none
!
  real(long) :: small = 1.0d-10
  real(long) :: rrgg, sth, cth, sphi, cphi
  real(long) :: point_line_dis2, root_det
  real(long) :: x0, y0, z0, rg_exc_dis2, xxyy2
  integer    :: irg, itg, ipg
!
!-- radial direction.
!
  do ipg = 0, npg
    do itg = 0, ntg
!
      irg_exin(itg,ipg) = 0
      irg_exout(itg,ipg) = 0
      rg_exin(itg,ipg) = 0.0d0
      rg_exout(itg,ipg) = 0.0d0
!
      sth = sinthg(itg)
      cphi = cosphig(ipg)
      point_line_dis2 = sepa**2*(1.0 - sth**2*cphi**2)
!
      if (point_line_dis2.lt.ex_radius**2) then
!
        root_det = sqrt(ex_radius**2 - point_line_dis2)
        rg_exin(itg,ipg)  = sepa*sth*cphi - root_det
        rg_exout(itg,ipg) = sepa*sth*cphi + root_det
!
        do irg = 0, nrg - 1
          if (rg_exin(itg,ipg).ge.rg(irg)-small.and. &
          &   rg_exin(itg,ipg).lt.rg(irg+1)-small) irg_exin(itg,ipg)=irg 
          if (rg_exout(itg,ipg).gt.rg(irg)+small.and. &
          &   rg_exout(itg,ipg).le.rg(irg+1)+small) irg_exout(itg,ipg)=irg+1
        end do
      end if
!
      if (itg.eq.0.or.ipg.eq.0) cycle    ! for half points
!
      ihrg_exin(itg,ipg) = 1
      ihrg_exout(itg,ipg) = 1
      hrg_exin(itg,ipg) = 0.0d0
      hrg_exout(itg,ipg) = 0.0d0
!
      sth = hsinthg(itg)
      cphi = hcosphig(ipg)
      point_line_dis2 = sepa**2*(1.0 - sth**2*cphi**2)
!
      if (point_line_dis2.lt.ex_radius**2) then
!
        root_det = sqrt(ex_radius**2 - point_line_dis2)
        hrg_exin(itg,ipg)  = sepa*sth*cphi - root_det
        hrg_exout(itg,ipg) = sepa*sth*cphi + root_det
!
        do irg = 0, nrg - 1  ! compare with rg, not with hrg
        if (hrg_exin(itg,ipg).ge.rg(irg).and. &
       &    hrg_exin(itg,ipg).lt.rg(irg+1)) ihrg_exin(itg,ipg)=irg+1
        if (hrg_exout(itg,ipg).gt.rg(irg).and. &
       &    hrg_exout(itg,ipg).le.rg(irg+1))ihrg_exout(itg,ipg)=irg+1
        end do
      end if
    end do
  end do
!
!-- theta direction
!
  do ipg = 0, npg
    do irg = 0, nrg
      itg_exin(irg,ipg) = 0
      itg_exout(irg,ipg) = 0
      thg_exin(irg,ipg) = 0.0d0
      thg_exout(irg,ipg) = 0.0d0
!
      if (irg.eq.0) cycle
      if (ipg.ge.npgxzp.and.ipg.le.npgxzm) cycle
!
      rrgg = rg(irg)
      cphi = cosphig(ipg)
!
      sth = (rrgg**2 + sepa**2 - ex_radius**2)/(2.0d0*sepa*rrgg*cphi)
      if (sth.le.1.0d0) then
        thg_exin(irg,ipg) = dasin(sth)
        thg_exout(irg,ipg) = pi - dasin(sth)
!
        do itg = 0, ntg-1
          if (thg_exin(irg,ipg).ge.thg(itg).and. &
         &    thg_exin(irg,ipg).lt.thg(itg+1)) itg_exin(irg,ipg) = itg
        end do
        itg_exout(irg,ipg) = ntg - itg_exin(irg,ipg)
      end if
    end do
  end do
!
!! for half points
  do ipg = 1, npg
    do irg = 1, nrg
      ihtg_exin(irg,ipg) = 1
      ihtg_exout(irg,ipg) = 1
      hthg_exin(irg,ipg) = 0.0d0
      hthg_exout(irg,ipg) = 0.0d0
!
      rrgg = hrg(irg)
      cphi = hcosphig(ipg)
!
      sth = (rrgg**2 + sepa**2 - ex_radius**2)/(2.0d0*sepa*rrgg*cphi)
      if (sth.le.1.0d0) then
        hthg_exin(irg,ipg) = dasin(sth)
        hthg_exout(irg,ipg) = pi - dasin(sth)
!
        do itg = 0, ntg-1
          if (hthg_exin(irg,ipg).ge.thg(itg).and. &
        &     hthg_exin(irg,ipg).lt.thg(itg+1)) ihtg_exin(irg,ipg) = itg+1
          if (hthg_exout(irg,ipg).ge.thg(itg).and. &
        &     hthg_exout(irg,ipg).lt.thg(itg+1)) ihtg_exout(irg,ipg)=itg+1
        end do
      end if
    end do
  end do
!
!-- phi direction
!
  do itg = 0, ntg
    do irg = 0, nrg
      ipg_exin(irg,itg) = 0
      ipg_exout(irg,itg) = npg
      phig_exin(irg,itg) = 0.0d0
      phig_exout(irg,itg) = 2.0d0*pi
!
      if (irg.eq.0.or.itg.eq.0.or.itg.eq.ntg) cycle 
!
      rrgg = rg(irg)
      sth = sinthg(itg)
!
      cphi = (rrgg**2 + sepa**2 - ex_radius**2)/(2.0d0*sepa*rrgg*sth)
      if (cphi.lt.1.0d0) then
        phig_exin(irg,itg) = dacos(cphi)
        phig_exout(irg,itg) = 2.0d0*pi - dacos(cphi)
        do ipg = 0, npg-1
          if (phig_exin(irg,itg).gt.phig(ipg).and. &
         &    phig_exin(irg,itg).le.phig(ipg+1)) ipg_exin(irg,itg)=ipg+1
          if (phig_exout(irg,itg).ge.phig(ipg).and. &
         &    phig_exout(irg,itg).lt.phig(ipg+1)) ipg_exout(irg,itg)=ipg
        end do
      end if
    end do
  end do
!
!! for half points
  do itg = 1, ntg
    do irg = 1, nrg
      ihpg_exin(irg,itg) = 1
      ihpg_exout(irg,itg) = npg
      hphig_exin(irg,itg) = 0.0d0
      hphig_exout(irg,itg) = 0.0d0
!see weight_midpoint_binary_excision.f90
!no good      hphig_exin(irg,itg) = 0.5d0*dphig
!no good      hphig_exout(irg,itg) = 2.0d0*pi - 0.5d0*dphig
!
      rrgg = hrg(irg)
      sth = hsinthg(itg)
!
      cphi = (rrgg**2 + sepa**2 - ex_radius**2)/(2.0d0*sepa*rrgg*sth)
      if (cphi.lt.1.0d0) then
        hphig_exin(irg,itg) = dacos(cphi)
        hphig_exout(irg,itg) = 2.0d0*pi- dacos(cphi)
        do ipg = 0, npg-1
          if (hphig_exin(irg,itg).gt.phig(ipg).and. &
         &    hphig_exin(irg,itg).le.phig(ipg+1)) ihpg_exin(irg,itg)=ipg+1
          if (hphig_exout(irg,itg).ge.phig(ipg).and. &
         &    hphig_exout(irg,itg).lt.phig(ipg+1)) ihpg_exout(irg,itg)=ipg+1
        end do
      end if
    end do
  end do
!
! --  coordinates of the central grid seen from the center of excision.
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        rb(irg,itg,ipg) = 0.0d0
        thb(irg,itg,ipg) = 0.0d0
        phib(irg,itg,ipg) = 0.0d0
!
        rrgg = rg(irg)
        sth = sinthg(itg)
        cth = costhg(itg)
        sphi = sinphig(ipg)
        cphi = cosphig(ipg)
        x0 = rrgg*sth*cphi
        y0 = rrgg*sth*sphi
        z0 = rrgg*cth
!
!        rg_exc_dis2 = rrgg**2 + sepa**2 - 2.0d0*rrgg*sepa*sth*cphi
        rg_exc_dis2 = (x0 - sepa)**2 + y0**2 + z0**2
        xxyy2       = (x0 - sepa)**2 + y0**2
!
        if (rg_exc_dis2.gt.0.0d0) &
      & rb(irg,itg,ipg)  = sqrt(rg_exc_dis2)
        thb(irg,itg,ipg) = atan2(sqrt(xxyy2),z0)
        phib(irg,itg,ipg) = dmod(2.0d0*pi+datan2(y0,x0-sepa),2.0d0*pi)
!
      end do
    end do
  end do
!
! -- coordinates of the central half grid seen from the center of excision.
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        hrb(irg,itg,ipg) = 0.0d0
        hthb(irg,itg,ipg) = 0.0d0
        hphib(irg,itg,ipg) = 0.0d0
!
        rrgg = hrg(irg)
        sth = hsinthg(itg)
        cth = hcosthg(itg)
        sphi = hsinphig(ipg)
        cphi = hcosphig(ipg)
        x0 = rrgg*sth*cphi
        y0 = rrgg*sth*sphi
        z0 = rrgg*cth
!
        rg_exc_dis2 = (x0 - sepa)**2 + y0**2 + z0**2
        xxyy2       = (x0 - sepa)**2 + y0**2
!
        if (rg_exc_dis2.gt.0.0d0) &
      & hrb(irg,itg,ipg)  = sqrt(rg_exc_dis2)
        hthb(irg,itg,ipg) = atan2(sqrt(xxyy2),z0)
        hphib(irg,itg,ipg) = dmod(2.0d0*pi+datan2(y0,x0-sepa),2.0d0*pi)
!
      end do
    end do
  end do
!
  ntg_exin_min = ntg
  ntg_exout_max = 0
  do irg = 0, nrg
    if (itg_exin(irg,npgxzp).gt.0) &
    & ntg_exin_min = min0(itg_exin(irg,npgxzp),ntg_exin_min)
    if (itg_exout(irg,npgxzp).gt.0) &
    & ntg_exout_max = max0(itg_exout(irg,npgxzp),ntg_exout_max)
  end do
  npg_exin_max = 0
  npg_exout_min = npg
  do irg = 0, nrg
    if (ipg_exin(irg,ntgxy).gt.0) &
    & npg_exin_max = max0(ipg_exin(irg,ntgxy),npg_exin_max)
    if (ipg_exout(irg,ntgxy).lt.npg) &
    & npg_exout_min = min0(ipg_exout(irg,ntgxy),npg_exout_min)
  end do
!
  write(6,*) 'ntg_exin_min = ' , ntg_exin_min, &
  &    ',     ntg_exout_max = ', ntg_exout_max
  write(6,*) 'npg_exin_max = ' , npg_exin_max, &
  &    ',     npg_exout_min = ', npg_exout_min
!
end subroutine calc_grid_points_binary_excision
end module grid_points_binary_excision
