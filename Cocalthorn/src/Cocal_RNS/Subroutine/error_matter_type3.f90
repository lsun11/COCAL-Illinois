subroutine error_matter_type3(pot,pot_bak,error,flag,cpot)
  use phys_constant,  only : long
  use def_matter, only : rs
  use coordinate_grav_r, only : rg
  use grid_parameter, only : nrf, ntf, npf, eps
  implicit none
  real(long), pointer :: pot(:,:,:), pot_bak(:,:,:)
  real(long), intent(out) :: error
  integer,    intent(out) :: flag
  character(len=4), intent(in) :: cpot
  real(long) :: error_pot = 0.0d0, small = 1.0d-14
  integer    :: irf, itf, ipf, ire, ite, ipe
!
  error = 0.0d0
  flag = 0
  do irf = 0, nrf-1
    do itf = 0, ntf
      do ipf = 0, npf
        error_pot = 2.0d0*dabs(pot(irf,itf,ipf) -     pot_bak(irf,itf,ipf)) &
      &                 /(dabs(pot(irf,itf,ipf))+dabs(pot_bak(irf,itf,ipf)) &
      &                 + small)
        if (error_pot > eps) flag = 1
!        error = dmax1(error,error_pot)
        if(error_pot > error) then
          error = error_pot
          ire = irf
          ite = itf
          ipe = ipf
        end if
      end do
    end do
  end do
!  write(6,*) ''
  write(6,'(1a19,3i5,a10,1p,e14.6,a2,a4,a1,1p,e20.12,a8,1p,e14.6)') '-- max error grid =',ire,ite,ipe, &
    &            "  rf coor=", rg(ire)*rs(ite,ipe), "  ", cpot, "=", pot(ire,ite,ipe), "  Error=", error

end subroutine error_matter_type3
