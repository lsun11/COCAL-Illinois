! r_coordinate for the field
!______________________________________________
module coordinate_grav_r_1D
  use phys_constant, only : nnrg, long
  use grid_parameter_1D, only : nrg, nrgin, rgin, rgmid, rgout, nrf, nrg_1, r_surf
  implicit none
  real(long)  :: rg(0:nnrg), rginv(0:nnrg) ! radial grid points.  inv means 1/r
  real(long)  :: hrg(nnrg), hrginv(nnrg)   ! mid points.
  real(long)  :: drg(nnrg), drginv(nnrg)   ! intervals.
contains
subroutine grid_r_1D
! Grid points of r-coordinate
! --- Set up of meshes  ---  ( r(1) = dr ( not the center ) )
!
  implicit none
  real(long)  ::  drdr, drdrinv
  real(long) :: ratrr  ! ratio between subsequent step outside rgmid (or rather nrgin). \delta_{j+1} = k \delta_{j}
  real(long) :: rvdom, ratrrb, alge, dalge, error, rdet, rdetf
  real(long) :: rfdom, reso_min, drmin, drmininv, ratrf, ratrfb
  integer    :: nrout, ir, itry
!
  drdr =  1.0d0/dble(nrf)
  drdrinv = 1.0e0/drdr
  rvdom = rgout - rgmid
  nrout = nrg - nrgin
!
    drg(1:nrgin)    = drdr
    drginv(1:nrgin) = drdrinv
! 
  if (rgin.ge.1.0d-13.and.rgin.lt.1.0d0) then
    rfdom = 1.0d0 - rgin
    ratrf = 0.1d0
    do
      ratrfb = ratrf
      alge = -(ratrf**nrf - 1.0d0 - rfdom*drdrinv*(ratrf-1.0e0))
      dalge = dble(nrf)*ratrf**(nrf-1) - rfdom*drdrinv
      ratrf = ratrf + alge/dalge
      error = 2.d0*(ratrf - ratrfb)/(ratrf + ratrfb)
      if (abs(error) <= 1.0d-14) exit
    end do
    if (ratrf.gt.1.0d0) stop 'ratrf not good'
    drg(nrf) = drdr
    drginv(nrf) = 1.0d0/drg(nrf)
    do ir = nrf-1, 1, -1
      drg(ir) = drg(ir+1)*ratrf
      drginv(ir) = 1.0d0/drg(ir)
    end do
    do ir = nrf+1, nrgin
      drg(ir) = drg(nrf)
      drginv(ir) = 1.0d0/drg(ir)
    end do
  end if
!
do itry = 0, 3
  if (itry.eq.0) ratrr = 2.0e-01   ! good for finer grid
  if (itry.eq.1) ratrr = 2.0e0   ! good for finer grid
  if (itry.eq.2) ratrr = 10.0e0  ! good for coarser grid
  if (itry.eq.3) ratrr = 50.0e0  ! good for coarsest grid
!
  DO
    ratrrb = ratrr
    alge = -(ratrr**(nrout+1)-ratrr-rvdom*drdrinv*(ratrr-1.0e0))
    dalge = dble(nrout+1)*ratrr**nrout - 1.0e0 - rvdom*drdrinv
    ratrr = ratrr + alge/dalge
    error = 2.d0*(ratrr - ratrrb)/(ratrr + ratrrb)
!
    IF (abs(error)<=1.d-14) EXIT
  END DO
! r-coordinate. interval dr
  DO ir = nrgin+1, nrg
    drg(ir) = ratrr*drg(ir-1)
    drginv(ir) = 1.0e0/drg(ir)
  END DO
! r-coordinate. grid ponit and mid point.
  rg(0) = rgin
  if (rg(0).le.1.0d-14) then
    rginv(0) = 0.0e0
  else
    rginv(0) = 1.0d0/rg(0)
  end if
  DO ir = 1, nrg
    rg(ir) = rg(ir-1) + drg(ir)
    hrg(ir) = (rg(ir) + rg(ir-1))*0.5e0
    rginv(ir) = 1.0e0/rg(ir)
    hrginv(ir) = 1.0e0/hrg(ir)
!      write(6,*)ir,rg(ir)
  END DO
!Œ‹‰Ê“¯‚¶‚É‚·‚é‚½‚ßŸŽè‚É‚Â‚¯‚½‚µ
  write(6,"(' GR ', 4es23.15)") error, ratrr
  write(6,"('    ', 4es23.15)") rgout, rg(nrg)
!‚Â‚¯‚½‚µI‚í‚è
!      WRITE(6,'(1p,2e12.4)') error, ratrr
!      WRITE(6,'(1p,2e12.4)') rgout, rg(nrg)
  rdet  = (rg(nrg) - rgout)/rgout
  rdetf = 0.0d0
  if (rgin.lt.1.0d0) rdetf =  rg(nrf) - 1.0d0
  IF(dabs(rdet)< 1.e-10) exit
  IF(dabs(rdet)>=1.e-10)then
    WRITE(6,*) ' bad coordinate itry =', itry
    WRITE(6,*) ' bad coordinate GR ', rdet, ratrr
    WRITE(6,*) ' bad coordinate FLUID ', rdetf, ratrf
    if (itry.eq.3) STOP
  END IF
end do
end subroutine grid_r_1D

subroutine grid_r_1D_surf   ! this is grid_r_bns

  implicit none
  real(long)  ::  drdr, drdrinv
  real(long) :: ratrr  ! ratio between subsequent step outside rgmid (or rather nrgin). \delta_{j+1} = k \delta_{j}
  real(long) :: rvdom, ratrrb, alge, dalge, error, rdet, rdetf, rdet1
  real(long) :: rfdom, reso_min, drmin, drmininv, ratrf, ratrfb
  integer    :: nrout, ir, itry, nrf_bh
  real(long) :: dr1, dr1inv, tolerant_val = 5.0d-14
  real(long) :: dr2, dr2inv
  integer    :: nr_count, count, nr3
  character(len=3) :: char_bh

!
! Region I r<=rgmid  for mpt3.
  drdr =  1.0d0/dble(nrf)
  drdrinv = 1.0d0/drdr
!
  drg(1:nrgin+1)    = drdr
  drginv(1:nrgin+1) = drdrinv
!
! the following if block only for mpt1,2:   Sets up the grid for r<3rgmid
! for mpt3 we have only the above lines i.e.  r<rgmid
!
  if (rgin==0.0d0) then
!
    drdr = r_surf/dble(nrf)
    drdrinv = 1.0d0/drdr

!   Region S
    drg(1:nrf)    = drdr
    drginv(1:nrf) = drdrinv

!   Region II
    dr2 =  1.0e0/dble(nrg_1)            ! nrg_1 plays the role of nrf in BBH
    drg(nrg_1+1:nrgin)    = dr2
    drginv(nrg_1+1:nrgin) = 1.0e0/dr2

    if(r_surf<1.0d0) then
      if(nrf==nrg_1) then
        write(6,*) "nrf==nrg_1 but r_surf<1.0...exiting"
        stop
      end if
      do itry = 0, 1
        if (itry.eq.0) ratrr = 2.0e-01 ! good for finer grid
        if (itry.eq.1) ratrr = 2.0e0   ! good for finer grid
!
        dr1 = drg(nrf)                 ! drg(nrf+1)= k drg(nrf)
        dr1inv = 1.0d0/dr1
        rvdom = 1.0d0 - r_surf
        nrout = nrg_1 - nrf

        DO
          ratrrb = ratrr
          alge = -(ratrr**(nrout+1)-ratrr-rvdom*dr1inv*(ratrr-1.0e0))
          dalge = dble(nrout+1)*ratrr**nrout - 1.0e0 - rvdom*dr1inv
          ratrr = ratrr + alge/dalge
          error = 2.d0*(ratrr - ratrrb)/(ratrr + ratrrb)
!
          IF (abs(error)<=1.d-14) EXIT
        END DO
!
!       r-coordinate. interval dr
        do ir = nrf+1, nrg_1
          drg(ir) = ratrr*drg(ir-1)
          drginv(ir) = 1.0e0/drg(ir)
        end do

        rg(0) = rgin
        DO ir = 1, nrgin
          rg(ir) = rg(ir-1) + drg(ir)
!          write(6,*) ir, rg(ir), drg(ir)
        end do
!
        rdet  = (rg(nrgin) - rgmid)/rgmid
        rdetf = 0.0d0
        if (rgin.lt.1.0d0) rdet1 =  rg(nrg_1) - 1.0d0
        if (rgin.lt.1.0d0) rdetf =  rg(nrf) - r_surf
        IF(dabs(rdet)< 1.e-10) then
          exit
        end if
        IF(dabs(rdet)>=1.e-10)then
          WRITE(6,*) ' bad coordinate itry =', itry
          WRITE(6,*) ' bad coordinate GR : rdet1, rdet, ratrr ', rdet1, rdet, ratrr
          WRITE(6,*) ' bad coordinate FLUID : rdetf ', rdetf
          if (itry.eq.1) STOP
        END IF
      end do
    else
      if(nrf.ne.nrg_1)  then
        write(6,*) "nrf is not equal to nrg_1, but r_surf=1...exiting"
        stop
      end if
    end if                   !  r_surf < 1.0
!
!   region III.
    drdr =  drg(nrgin)       !  In bhex first interval in region III: drg(nrgin+1)= k drg(nrgin)
                             !  If we want to change drdr we must also modify
                             !  LINE (*)
    drdrinv = 1.0d0/drdr
    rvdom = 2.0d0*rgmid
    nr3 = nrg_1               !  plays the role of nrf in BBH
    ratrr = 2.0e0
!
    DO
      ratrrb = ratrr
      alge = -(ratrr**(nr3+1)-ratrr-rvdom*drdrinv*(ratrr-1.0e0))
      dalge = dble(nr3+1)*ratrr**nr3 - 1.0e0 - rvdom*drdrinv
      ratrr = ratrr + alge/dalge
      error = 2.d0*(ratrr - ratrrb)/(ratrr + ratrrb)
!
      IF (abs(error)<=1.d-14) EXIT
    END DO

!   r-coordinate. interval dr:   drg(nrgin+1)= k drdr
    drg(nrgin+1)    = ratrr*drdr
    drginv(nrgin+1) = 1.0e0/drg(nrgin+1)

    DO ir = nrgin+2, nrgin+nr3
      drg(ir) = ratrr*drg(ir-1)               !   LINE (*)
      drginv(ir) = 1.0e0/drg(ir)
    END DO
!
    rg(0) = rgin
    DO ir = 1, nrgin+nr3
      rg(ir) = rg(ir-1) + drg(ir)
    end do
  end if                         !  rgin==0
!
!  Now for r > 3rgmid if mpt1,2
!  or      r > rgmid  if mpt3
!
do itry = 0, 3
  if (itry.eq.0) ratrr = 2.0e-01 ! good for finer grid
  if (itry.eq.1) ratrr = 2.0e0   ! good for finer grid
  if (itry.eq.2) ratrr = 10.0e0  ! good for coarser grid
  if (itry.eq.3) ratrr = 50.0e0  ! good for coarsest grid
!
  drdr =  1.0d0/dble(nrf)
  drdrinv = 1.0d0/drdr
  rvdom = rgout - rgmid
  nrout = nrg - nrgin
  if (rgin==0.0d0) then
    drdr =  drg(nrgin+nr3)
    drdrinv = 1.0d0/drdr
    rvdom = rgout - 3.0d0*rgmid
    nrout = nrg - (nrgin + nr3)
  end if
!
  DO
    ratrrb = ratrr
    alge = -(ratrr**(nrout+1)-ratrr-rvdom*drdrinv*(ratrr-1.0e0))
    dalge = dble(nrout+1)*ratrr**nrout - 1.0e0 - rvdom*drdrinv
    ratrr = ratrr + alge/dalge
    error = 2.d0*(ratrr - ratrrb)/(ratrr + ratrrb)
!
    IF (abs(error)<=1.d-14) EXIT
  END DO

! r-coordinate. interval dr
  nr_count = nrgin+1
  if (rgin==0.0d0) nr_count = nrgin+nr3+1
  DO ir = nr_count, nrg
    drg(ir) = ratrr*drg(ir-1)
    drginv(ir) = 1.0e0/drg(ir)
  END DO
! r-coordinate. grid ponit and mid point.
  rg(0) = rgin
  if (rg(0).le.1.0d-14) then
    rginv(0) = 0.0e0
  else
    rginv(0) = 1.0d0/rg(0)
  end if
  DO ir = 1, nrg
    rg(ir) = rg(ir-1) + drg(ir)
    hrg(ir) = (rg(ir) + rg(ir-1))*0.5e0
    rginv(ir) = 1.0e0/rg(ir)
    hrginv(ir) = 1.0e0/hrg(ir)
!      write(6,*)ir,rg(ir)
  END DO
!      WRITE(6,'(1p,2e12.4)') error, ratrr
!      WRITE(6,'(1p,2e12.4)') rgout, rg(nrg)
  rdet  = (rg(nrg) - rgout)/rgout
  rdetf = 0.0d0
  if (rgin.lt.1.0d0) rdet1 =  rg(nrg_1) - 1.0d0
  if (rgin.lt.1.0d0) rdetf =  rg(nrf) - r_surf
  IF(dabs(rdet)< 1.e-10) then
    write(6,*) '1D good r grid.'
    WRITE(6,*) 'good coordinate itry =', itry
    WRITE(6,*) 'good coordinate GR : rdet1, rdet, ratrr ', rdet1, rdet, ratrr
    WRITE(6,*) 'good coordinate FLUID : rdetf ', rdetf
    exit
  end if
  IF(dabs(rdet)>=1.e-10)then
    WRITE(6,*) ' bad coordinate itry =', itry
    WRITE(6,*) ' bad coordinate GR : rdet1, rdet, ratrr', rdet1, rdet, ratrr
    WRITE(6,*) ' bad coordinate FLUID : rdetf ', rdetf
    if (itry.eq.3) STOP
  END IF
end do
!
end subroutine grid_r_1D_surf
!
!
subroutine grid_r_1D_surf_const

  implicit none
  real(long)  ::  drdr, drdrinv
  real(long) :: ratrr  ! ratio between subsequent step outside rgmid (or rather nrgin). \delta_{j+1} = k \delta_{j}
  real(long) :: rvdom, ratrrb, alge, dalge, error, rdet, rdetf, rdet1
  real(long) :: rfdom, reso_min, drmin, drmininv, ratrf, ratrfb
  integer    :: nrout, ir, itry, nrf_bh, nrbns
  real(long) :: dr1, dr1inv, tolerant_val = 5.0d-14
  real(long) :: dr2, dr2inv, rbns
  integer    :: nr_count, count, nr3
  character(len=3) :: char_bh

!
! Region I r<=rgmid  for mpt3.
  drdr =  1.0d0/dble(nrf)
  drdrinv = 1.0d0/drdr
!
  drg(1:nrgin+1)    = drdr
  drginv(1:nrgin+1) = drdrinv
!
! the following if block only for mpt1,2:   Sets up the grid for r<3rgmid
! for mpt3 we have only the above lines i.e.  r<rgmid
!
  if (rgin==0.0d0) then
    nrbns = 7*nrf
    rbns  = 7.0d0
    rvdom = rgout - rbns
    nrout = nrg - nrbns
!
    if (nrout<0)  then
      write(6,*) "It must be nrg > 7*nrf....exiting"
      stop
    end if
!   r-coordinate. interval dr
    DO ir = 1, nrbns
      drg(ir) = drdr
      drginv(ir) = 1.0e0/drg(ir)
    END DO
!
    rg(0) = rgin
    DO ir = 1, nrbns
      rg(ir) = rg(ir-1) + drg(ir)
!      write(6,*) ir, rg(ir), drg(ir)
    end do
!
    rdet  = (rg(nrbns) - rbns)/rbns
    IF(dabs(rdet)< 1.e-10) then
      WRITE(6,*) 'good coordinate GR : rdet, ratrr ', rdet
    end if
    IF(dabs(rdet)>=1.e-10)then
      WRITE(6,*) ' bad coordinate GR : rdet, ratrr ', rdet
      stop
    END IF
  end if
!
!  Now for r > 3rgmid if mpt1,2
!  or      r > rgmid  if mpt3
!
do itry = 0, 3
  if (itry.eq.0) ratrr = 2.0e-01 ! good for finer grid
  if (itry.eq.1) ratrr = 2.0e0   ! good for finer grid
  if (itry.eq.2) ratrr = 10.0e0  ! good for coarser grid
  if (itry.eq.3) ratrr = 50.0e0  ! good for coarsest grid
!
  drdr =  1.0d0/dble(nrf)
  drdrinv = 1.0d0/drdr
  rvdom = rgout - rgmid
  nrout = nrg - nrgin
  if (rgin==0.0d0) then
    drdr =  drg(nrbns)
    drdrinv = 1.0d0/drdr
    rvdom = rgout - rbns
    nrout = nrg - nrbns
  end if
!
  DO
    ratrrb = ratrr
    alge = -(ratrr**(nrout+1)-ratrr-rvdom*drdrinv*(ratrr-1.0e0))
    dalge = dble(nrout+1)*ratrr**nrout - 1.0e0 - rvdom*drdrinv
    ratrr = ratrr + alge/dalge
    error = 2.d0*(ratrr - ratrrb)/(ratrr + ratrrb)
!
    IF (abs(error)<=1.d-14) EXIT
  END DO

  write(6,*) 'region II: ratrr=', ratrr
! r-coordinate. interval dr
  nr_count = nrgin+1
  if (rgin==0.0d0) nr_count = nrbns+1
  DO ir = nr_count, nrg
    drg(ir) = ratrr*drg(ir-1)
    drginv(ir) = 1.0e0/drg(ir)
  END DO
! r-coordinate. grid ponit and mid point.
  rg(0) = rgin
  if (rg(0).le.1.0d-14) then
    rginv(0) = 0.0e0
  else
    rginv(0) = 1.0d0/rg(0)
  end if
  DO ir = 1, nrg
    rg(ir) = rg(ir-1) + drg(ir)
    hrg(ir) = (rg(ir) + rg(ir-1))*0.5e0
    rginv(ir) = 1.0e0/rg(ir)
    hrginv(ir) = 1.0e0/hrg(ir)
!      write(6,*)ir,rg(ir)
  END DO
!
!
  rdet  = (rg(nrg) - rgout)/rgout
  rdetf = 0.0d0
  if (rgin.lt.1.0d0) rdetf =  rg(nrf) - 1.0d0
  IF(dabs(rdet)< 1.e-10) then
    WRITE(6,*) 'good coordinate itry =', itry
    WRITE(6,*) 'good coordinate GR : rdet, ratrr ', rdet, ratrr
    WRITE(6,*) 'good coordinate FLUID : rdetf, ratrf', rdetf, ratrf
    exit
  end if
  IF(dabs(rdet)>=1.e-10)then
    WRITE(6,*) ' bad coordinate itry =', itry
    WRITE(6,*) ' bad coordinate GR : rdet, ratrr', rdet, ratrr
    WRITE(6,*) ' bad coordinate FLUID : rdetf, ratrf', rdetf, ratrf
    if (itry.eq.3) STOP
  END IF
end do
!
end subroutine grid_r_1D_surf_const

end module coordinate_grav_r_1D
