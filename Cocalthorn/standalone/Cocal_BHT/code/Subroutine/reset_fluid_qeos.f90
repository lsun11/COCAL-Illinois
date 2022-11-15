subroutine reset_fluid_qeos
  use def_matter, only : rhof, rs
  use def_matter_parameter, only : rhoc_qs, rhos_qs
  use grid_parameter, only : nrf, ntf, npf, &
  &                          ntfeq, ntfxy, npfyzp, npfxzp, npfxzm, &
  &                          ratio, NS_shape, EQ_point
  implicit none
  integer :: it, ip, iremax
! For parameter equations
  rs(ntfeq,npfxzp) = 1.0d0; rs(ntfeq,npf) = 1.0d0
! EQ_point = XZ : Solve parameter eqs at the equator and pole
  if (EQ_point.eq.'XZ') rs(0,0:npf) = ratio   !! rs(ntf,0:npf) = ratio !! 
! EQ_point = XY : Solve parameter eqs at two equatorial points
  if (EQ_point.eq.'XY') rs(ntfeq,npfyzp) = ratio
! Reset central value of Emden fn.
  call search_rhofmax_xaxis_grid(iremax)
  if (iremax.eq.0) then
    rhof(iremax,0:ntf,0:npf) = rhoc_qs
  else
    rhof(iremax,ntfeq,0:npf) = rhoc_qs
  end if
! Reset surface values of Emden fn.
!  emd(nrf,0:ntf,0:npf) = emdc*1.0d-05
  rhof(nrf,0:ntf,0:npf) = rhos_qs
! Reset ip = 0 eq ip = npf
  rhof(0:nrf,0:ntf,npf) = rhof(0:nrf,0:ntf,0)
  rs(0:ntf,npf) = rs(0:ntf,0)
!
! Impose symmetry of the shape
! NS_shape = JB : tri-axial configuration
! NS_shape = ML : axi-symmetric configuration
!
  if (NS_shape.eq.'JB') then
! Impose tri-axial and equatorial symmetry
  do it = 0, ntfeq
    do ip = 0, npfyzp
      rs(it,npfxzm-ip) = rs(it,ip)
      rs(it,npfxzm+ip) = rs(it,ip)
      rs(it,npf   -ip) = rs(it,ip)
      rs(ntf-it,       ip) = rs(it,ip)
      rs(ntf-it,npfxzm-ip) = rs(it,ip)
      rs(ntf-it,npfxzm+ip) = rs(it,ip)
      rs(ntf-it,npf   -ip) = rs(it,ip)
      rhof(0:nrf,it,npfxzm-ip) = rhof(0:nrf,it,ip)
      rhof(0:nrf,it,npfxzm+ip) = rhof(0:nrf,it,ip)
      rhof(0:nrf,it,npf   -ip) = rhof(0:nrf,it,ip)
      rhof(0:nrf,ntf-it,       ip) = rhof(0:nrf,it,ip)
      rhof(0:nrf,ntf-it,npfxzm-ip) = rhof(0:nrf,it,ip)
      rhof(0:nrf,ntf-it,npfxzm+ip) = rhof(0:nrf,it,ip)
      rhof(0:nrf,ntf-it,npf   -ip) = rhof(0:nrf,it,ip)
    end do
  end do
  else if (NS_shape.eq.'ML') then
! Impose axi-symmetry and equatorial symmetry
!test  do it = 0, ntf
  do it = 0, ntfeq
    do ip = 1, npf
      rs(it,ip) = rs(it,0)
      rs(ntf-it,ip) = rs(it,0)
      rhof(0:nrf,it,ip) = rhof(0:nrf,it,0)
      rhof(0:nrf,ntf-it,ip) = rhof(0:nrf,it,0)
    end do
  end do
  end if
!
! ### surface is forced to be smaller than 1.0
! ### can be chosen by select_reset_fluid_radius_le_1.sh
!le1  do it = 0, ntf
!le1    do ip = 0, npf
!le1      if (rs(it,ip).gt.1.0d0) then
!le1        rs(it,ip) = 1.0d0
!le1      end if
!le1    end do
!le1  end do
!!
end subroutine reset_fluid_qeos
