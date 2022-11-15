subroutine reset_fluid_BNS
  use def_matter, only : emd, rs
  use def_matter_parameter, only : emdc
  use grid_parameter, only : nrf, ntf, npf, &
  &                          ntfeq, ntfxy, npfyzp, npfxzp, npfxzm, &
  &                          ratio, NS_shape, EQ_point, r_surf
  implicit none
  integer :: it, ip
! For parameter equations
  rs(ntfeq,npfxzp) = 1.0d0;  rs(ntfeq,npf) = 1.0d0

! EQ_point = XZ : Solve parameter eqs at the equator and pole
!  if (EQ_point.eq.'XZ') rs(0,0:npf) = ratio   !! rs(ntf,0:npf) = ratio !! 
! EQ_point = XY : Solve parameter eqs at two equatorial points
!  if (EQ_point.eq.'XY') rs(ntfeq,npfyzp) = ratio

! Reset central value of Emden fn.
  emd(0,0:ntf,0:npf) = emdc

! Reset surface values of Emden fn.
  emd(nrf,0:ntf,0:npf) = emdc*1.0d-14

! Reset ip = 0 eq ip = npf
  emd(0:nrf,0:ntf,npf) = emd(0:nrf,0:ntf,0)
  rs(0:ntf,npf) = rs(0:ntf,0)  
!
! ### surface is forced to be smaller than 1.0
  do it = 0, ntf
    do ip = 0, npf
      if (rs(it,ip).gt.1.0d0) then
        rs(it,ip) = 1.0d0
      end if
    end do
  end do
!
end subroutine reset_fluid_BNS
